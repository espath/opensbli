#!/usr/bin/env python
# Import all the functions from opensbli
from opensbli import *
from sympy import sin, log, cos, pi, tanh
from opensbli.utilities.helperfunctions import substitute_simulation_parameters, print_iteration_ops
from stats import *

# STEP 0 Create the equations required for the numerical solution
simulation_parameters = {
'Re'    :   '400.0',
'gama'  :   '1.4',
'Minf'  :   '0.4',
'Pr'    :   '0.7',
'dt'    :   '0.0001',
'niter' :   '250000',
'block0np0' :   '360',
'block0np1' :   '355',
'block0np2' :   '300',
'Delta0block0'  :   '6.0/block0np0',
'Delta1block0'  :   '2.0/(block0np1-1)',
'Delta2block0'  :   '3.0/block0np2',
"c0"    :   '-1',
"c1"    :   '0',
"c2"    :   '0',
"lx0"   :   "6.0",
"lx2"   :   "3.0",
"stretch"   :   "1.7",
"Twall" :   "1.0",
"theta" :   "0.1",
'lambda0_TVD'       : 'dt/Delta0block0',
'lambda1_TVD'       : 'dt/Delta1block0',
'lambda2_TVD'       : 'dt/Delta2block0',
"kappa_TVD"        : '1.25',
'shock_factor'      :   '1',
'inv_rfact0_block0'     :   '1.0/Delta0block0',
'inv_rfact1_block0'     :   '1.0/Delta1block0',
'inv_rfact2_block0'     :   '1.0/Delta2block0',
'TENO_CT'       :   '1e-6',
'aCF'          : '100.0', # Hyperbolic tangent forcing parameter
}

# Isothermal-Isothermal or Isothermal-Adiabatic walls
isothermal_adiabatic = False
# Turn on/off DRP explicit filtering
explicit_filters = True
# Select which scheme to use
#### Direct application of WENO/TENO ####
weno = False
teno = False
#### Central + WENO/TVD filter step ####
TVD = False # if False -> Central + WENO filter is used
# Statistics?
stats = True
# Problem dimension
ndim = 3
# # Constants that are used
constants = ["Re", "Pr", "gama", "Minf", "c_j", "aCF"]
# symbol for the coordinate system in the equations
coordinate_symbol = "x"
conservative = True
einstein_eq = EinsteinEquation()
# Choose which set of equations to use
if weno or teno:
    if weno:
        sc1 = "**{\'scheme\':\'Weno\'}"
    else:
        sc1 = "**{\'scheme\':\'Teno\'}"
    # Define the compresible Navier-Stokes equations in Einstein notation.
    mass = "Eq(Der(rho,t), - Conservative(rhou_j,x_j,%s))" % sc1
    momentum = "Eq(Der(rhou_i,t) , -Conservative(rhou_i*u_j + KD(_i,_j)*p,x_j , %s) + Der(tau_i_j,x_j) )" % sc1
    energy = "Eq(Der(rhoE,t), - Conservative((p+rhoE)*u_j,x_j, %s) - Der(q_j,x_j) + Der(u_i*tau_i_j ,x_j) )" % sc1
    stress_tensor = "Eq(tau_i_j, (mu/Re)*(Der(u_i,x_j)+ Der(u_j,x_i) - (2/3)* KD(_i,_j)* Der(u_k,x_k)))"
    heat_flux = "Eq(q_j, (-mu/((gama-1)*Minf*Minf*Pr*Re))*Der(T,x_j))"
    # Substitutions
    substitutions = [stress_tensor, heat_flux]
    base_eqns = [mass, momentum, energy]
    # Expand the base equations
    for i, base in enumerate(base_eqns):
        base_eqns[i] = einstein_eq.expand(base, ndim, coordinate_symbol, substitutions, constants)
    base_eqns = flatten(base_eqns)
    mass, momentum, energy = base_eqns[0], base_eqns[1:ndim+1], base_eqns[-1]
else:
    NS = NS_Split('KGP', ndim, constants, coordinate_symbol=coordinate_symbol, conservative=conservative, viscosity='dynamic', energy_formulation='enthalpy', debug=False)
    mass, momentum, energy = NS.mass, NS.momentum, NS.energy

# Add counter-flow forcing terms
CF_forcing_momentum = "Eq(rhou_i, -KD(_i,_j)*c_j*phi)"
CF_forcing_energy = "Eq(rhoE, - Dot(c_j*phi, u_j))"
CF_mom = einstein_eq.expand(CF_forcing_momentum, ndim, coordinate_symbol, [], constants)
CF_energy = einstein_eq.expand(CF_forcing_energy, ndim, coordinate_symbol, [], constants)
# Add the forcing terms to momentum and energy equations
for i in range(ndim):
    momentum[i] = OpenSBLIEq(momentum[i].lhs, momentum[i].rhs + CF_mom[i].rhs)
energy = OpenSBLIEq(energy.lhs, energy.rhs + CF_energy.rhs)


# Expand the simulation equations, for this create a simulation equations class
simulation_eq = SimulationEquations()
simulation_eq.add_equations(mass)
simulation_eq.add_equations(momentum)
simulation_eq.add_equations(energy)
# Constituent relations used in the system
velocity = "Eq(u_i, rhou_i/rho)"
pressure = "Eq(p, (gama-1)*(rhoE - rho*(1/2)*(KD(_i,_j)*u_i*u_j)))"
temperature = "Eq(T, p*gama*Minf*Minf/(rho))"
viscosity = "Eq(mu, (T**0.7))"
enthalpy = "Eq(H, (rhoE + p) / rho)"
speed_of_sound = "Eq(a, (gama*p/rho)**0.5)"
Force = "Eq(phi, tanh(aCF*DataObject('x1')))"
internal_energy = "Eq(e, p / (rho*(gama-1)))"
# Expand the constituent relations and them to the constituent relations class
constituent = ConstituentRelations()  # Instantiate constituent relations object
for input_eqn in [velocity, pressure, temperature, enthalpy, viscosity, speed_of_sound, Force, internal_energy]:
    eqns = einstein_eq.expand(input_eqn, ndim, coordinate_symbol, [], constants)
    constituent.add_equations(eqns)

# Create a simulation block
block = SimulationBlock(ndim, block_number=0, conservative=conservative)

# Transform equations for non uniform grid distribution
metriceq = MetricsEquation()
metriceq.generate_transformations(ndim, coordinate_symbol, [(False, False), (True, False), (False, False)], 2)
simulation_eq.apply_metrics(metriceq)

## Select the numerical schemes to use
schemes = {}
if weno or teno:
    scheme_halo_depth = [-4, 4]
else:
    scheme_halo_depth = [-2, 2]
# Scheme selection
if weno:
    Avg = SimpleAverage([0, 1])
    LF = LFWeno(order=5, formulation='Z', averaging=Avg, flux_type='LLF')
    # LF = HLLCWeno(order=5, formulation='Z', averaging=Avg, flux_type='HLLC-LM')
    # Add to schemes
    schemes[LF.name] = LF
elif teno:
    Avg = RoeAverage([0, 1])
    LF = LFTeno(order=6, averaging=Avg, flux_type='LLF')
    schemes[LF.name] = LF

# Central scheme 
fns = 'u0 u1 u2 T'
cent = StoreSome(4, fns)
schemes[cent.name] = cent
# Time-stepping
rk = RungeKuttaLS(4)
schemes[rk.name] = rk
block.set_discretisation_schemes(schemes)

# Define the variables used for creating boundary conditions and the initialisation
# dx and dy of the grid
dx, dy, dz = block.deltas
# Indices for the grid location
i, j, k = block.grid_indexes
# Some constants used
gama, Minf, Re = symbols('gama Minf Re', **{'cls': ConstantObject})
q_vector = flatten(simulation_eq.time_advance_arrays)

##### Boundary condition selection
boundaries = []
# Periodic boundaries in x0 direction
direction = 0
boundaries += [PeriodicBC(direction, side=0, halos=scheme_halo_depth)]
boundaries += [PeriodicBC(direction, side=1, halos=scheme_halo_depth)]
# Isothermal wall in x1 direction
Twall = ConstantObject("Twall")
wall_energy = [Eq(q_vector[-1], Twall*q_vector[0] / (gama * Minf**2.0 * (gama - S.One)))]
lower_wall_eq = wall_energy[:]
direction = 1
boundaries += [IsothermalWallBC(direction, 0, lower_wall_eq)]
# boundaries += [IsothermalWall_ZeroPressureGradBC(direction, 0, None)]
# Side 1 (top) boundary
if isothermal_adiabatic:
    boundaries += [AdiabaticWall_CarpenterBC(direction, 1)]
else:
    boundaries += [IsothermalWallBC(direction, 1, lower_wall_eq)]
    # boundaries += [IsothermalWall_ZeroPressureGradBC(direction, 1, None)]
# Periodic boundaries in x2 direction
direction = 2 # spanwise
boundaries += [PeriodicBC(direction, side=0, halos=scheme_halo_depth)]
boundaries += [PeriodicBC(direction, side=1, halos=scheme_halo_depth)]
# set the boundaries for the block
block.set_block_boundaries(boundaries)

# Create the grid and intial conditions
# Arrays to store x and y coordinates, i.e (x0 and x1)
x, y, z = symbols('x0:%d' % ndim, **{'cls': DataObject})
grid_equations = []
# Equations for generating the grid, simple equispacing grid, later change to stretched
#grid_equations += [Eq(x, i * dx), Eq(y, -S.One + j * dy), Eq(z, k * dz)]
ny = ConstantObject('block0np1')
stretch = ConstantObject('stretch')
Ly = 2
stretched_eqn =  0.5*Ly*(1.0-((tanh(stretch*(1.0-2.0*(j/(ny-1.0)))))/(tanh(stretch))))-1.0
grid_equations += [Eq(x, i * dx), Eq(y, stretched_eqn), Eq(z, k * dz)]
initial_equations = []

lx, ly, lz = symbols('lx0:%s' % ndim, **{'cls': ConstantObject})
sx, sy, sz = symbols('sx0:%s' % ndim, **{'cls': GridVariable})
cx, cy, cz = symbols('cx0:%s' % ndim, **{'cls': GridVariable})
ubar, amp, vonkar, b = symbols('ubar amp vonkar b', **{'cls': GridVariable})
Re = ConstantObject("Re")
# Set the constants for initialisation
initial_equations += [Eq(b, 5.5), Eq(vonkar, 2.5)]
# Equation for umean
initial_equations += [Eq(ubar, Piecewise(((1 - Abs(y)) * (Re), (1 - Abs(y)) * Re < 10.0),
                         (vonkar * log((1 - Abs(y)) * (Re)) + b, True)))]
# Amplitude of the disturbances
initial_equations += [Eq(amp, 0.1 * (vonkar * log((Re)) + b))]
# sin disturbances in each direction, 16, 8
initial_equations += [Eq(sx, sin(4.0 * pi * x / lx)), Eq(sy, sin(pi * y)),
                         Eq(sz, sin(2.0 * pi * z / lz))]

subs_dict = {}
# cosine disturbances in each direction
initial_equations += [Eq(cx, cos(4.0 * pi * x / lx)), Eq(cy, S.One + \
                         cos(pi * y)), Eq(cz, cos(2.0 * pi * z / lz))]
# Substitutions required to differentiate the equations symbolically
for eqn in initial_equations:
    subs_dict[eqn.lhs] = eqn.rhs
d, u, v, w, p = symbols('d u0 u1 u2 p', **{'cls': GridVariable})
# Perturbations
initial_equations += [Eq(u, ubar + amp * (lx / 2.0) * cx * sy * sz)]
u_eqn = initial_equations[-1].rhs.subs(subs_dict)
initial_equations += [Eq(v, - amp * sx * cy * sz)]
v_eqn = initial_equations[-1].rhs.subs(subs_dict)
initial_equations += [Eq(w, - amp * (lz / 2.0) * sx * sy * cz)]
w_eqn = initial_equations[-1].rhs.subs(subs_dict)

# Check the initial condition
from sympy import diff
output = Eq(symbols('divergence'), diff(u_eqn, x) + diff(v_eqn, y)  + diff(w_eqn, z))
print("Continuity check on the initial condition:")
pprint("du/dx + dv/dy + dw/dz")
pprint(output)
initial_equations += [Eq(p, S.One / (gama * Minf**2.0))]
initial_equations += [Eq(d, S.One)]
# Now evaluate the Conservative vector
initial_equations += [Eq(q_vector[0], d)]
initial_equations += [Eq(q_vector[1], d*u), Eq(q_vector[2], d*v), Eq(q_vector[3], d*w)]
initial_equations += [Eq(q_vector[4], p / (gama - S.One) + 0.5 * d * (u**2 + v**2 + w**2))]

# Instantiate a grid based initialisation classes
initial = GridBasedInitialisation()
initial.add_equations(grid_equations + initial_equations)

# Set the equation classes for the block (list)
block.set_equations([constituent, simulation_eq, initial, metriceq])

# WENO/TVD filter if not using direct application of WENO/TENO
if not weno and not teno:
    if TVD:
        WF = TVDFilter(block, airfoil=False)
        block.set_equations(WF.equation_classes)
    else:
        WF = WENOFilter(block, order=7, formulation='Z', flux_type='LLF', airfoil=False, metrics=metriceq, optimize=True)
        block.set_equations(WF.equation_classes)
# Additional filtering if needed
if explicit_filters:
    # Dispersion relation preserving filters
    DRP = ExplicitFilter(block, [0,1,2], width=9, filter_type='DRP', optimized=False, sigma=0.1, airfoil=False, multi_block=None)
    block.set_equations(DRP.equation_classes)

# Add statistics calculations
if stats:
    # Create the statistics equations, this shows another way of writing the equations
    q_vector = flatten(simulation_eq.time_advance_arrays)
    # Create equation classes for statistics, see stats.py in the current folder.
    # Later an automatic way of creating statistics equations will be provided
    stat_equation_classes = favre_averaged_stats(ndim, q_vector)
    statistics = [DataObject('u2u2mean'), DataObject('rhou2mean'), DataObject('rhou2u1mean'), DataObject('u0u0mean'), 
        DataObject('rhou1u0mean'), DataObject('E_mean'), DataObject('u1u0mean'), DataObject('u1u1mean'), DataObject('rhou2u0mean'),
        DataObject('rhou0mean'), DataObject('rhou1mean'), DataObject('pp_mean'), DataObject('rhou2u2mean'), DataObject('u2mean'), 
        DataObject('M_mean'), DataObject('u2u0mean'), DataObject('p_mean'), DataObject('a_mean'), DataObject('T_mean'), 
        DataObject('rhou0u0mean'), DataObject('rhomean'), DataObject('mu_mean'), DataObject('u2u1mean'), DataObject('TT_mean'), 
        DataObject('rhou1u1mean'), DataObject('u0mean'), DataObject('u1mean'), DataObject('D11')] 
    kwargs = {'iotype': "Write", 'name': "stats_output.h5"}
    stat_arrays = statistics
    stats_hdf5 = iohdf5(arrays=stat_arrays, **kwargs)
else:
    stat_equation_classes = []
block.set_equations(stat_equation_classes)

# STEP 4 add io for the block
kwargs = {'iotype': "Write"}
output_arrays = simulation_eq.time_advance_arrays + [x, y, z, DataObject('D11')]
output_hdf5 = iohdf5(arrays=output_arrays, save_every=5000, **kwargs)
block.setio([output_hdf5])
try:
    block.setio([stats_hdf5])
except:
    pass

# Add 2D I/O slicing option
grid_slice_hdf5_side = iohdf5_slices(blocknumber=0, **{'iotype': "Init"})
# x-y side view coordinates
coords = [([DataObject('x0'), DataObject('x1')], 2, 'block0np2/2')]
grid_slice_hdf5_side.add_slices(coords)
# Q vector slices written out in time
slices_hdf5_side = iohdf5_slices(save_every=500, blocknumber=0, **{'iotype': "Write"})
# x-y side view
if not weno and not teno and not TVD:
    slice_arrays = [DataObject('%s' % i) for i in ['rho', 'rhou0', 'rhou1', 'rhou2', 'rhoE']] + [DataObject('WENO_filter')]
else:
    slice_arrays = [DataObject('%s' % i) for i in ['rho', 'rhou0', 'rhou1', 'rhou2', 'rhoE']]
slices = [(slice_arrays, 2, 'block0np2/2')]
slices_hdf5_side.add_slices(slices)
grid_slice_hdf5_surfaces = iohdf5_slices(blocknumber=1, **{'iotype': "Init"})
# Surface coordinates
coords = [([DataObject('x0'), DataObject('x2')], 1, 20)]
coords += [([DataObject('x0'), DataObject('x2')], 1, 'block0np1 - 20')]
grid_slice_hdf5_surfaces.add_slices(coords)
# Q vector slices written out in time
slices_hdf5_surfaces = iohdf5_slices(save_every=500, blocknumber=0, **{'iotype': "Write"})
# Wall normal planes 
slices = [(slice_arrays, 1, 20)]
slices += [(slice_arrays, 1, 'block0np1 - 20')]
slices_hdf5_surfaces.add_slices(slices)
block.setio([grid_slice_hdf5_side, slices_hdf5_side, slices_hdf5_surfaces])

# Perform the symbolic discretisation of the equations
block.discretise()

# Apply a periodic BC for WENO filter: # Note -> has to be done after block.discretise()
if not teno and not weno and not TVD:
    WF.update_periodic_boundary(block, halos=[-4,4])

#Add some full [-5,5] halo swaps over the periodic directions only when the filter is called
if explicit_filters:
    def create_exchange_calls_codes(block, dsets):
        kernels = []
        arrays = [block.location_dataset(a) for a in flatten(dsets)]
        for direction in [0,2]:
            for side in [0,1]:
                BC = PeriodicBC(direction, side, halos=[-5, 5], full_depth=True)
                kernels += [BC.apply(arrays, block)]
        return kernels

    # Make some full swaps for interfaces before filtering
    filter_swaps = create_exchange_calls_codes(block, ['rho', 'rhou0', 'rhou1', 'rhou2', 'rhoE'])
    for no, eq in enumerate(block.list_of_equation_classes):
        if isinstance(eq, UserDefinedEquations):
            if eq.full_swap:
                eq.Kernels += filter_swaps

# Simulation monitors
# Density at entrance plane, lower wall, midspan
arrays = ['rho']
probe_locations = [(0, 0, '(block0np2-1)/2')]
# u velocity at entrance plane, 1 point off the wall, midspan
arrays += ['u0']
probe_locations += [(0, 1, '(block0np2-1)/2')]
# Temperature at entrance plane, 1st point off lower wall, midspan
arrays += ['T']
probe_locations += [(0, 1, '(block0np2-1)/2')]
SM = SimulationMonitor(arrays, probe_locations, block, print_frequency=100)
# create an algorithm from the numerical solution
alg = TraditionalAlgorithmRK(block)#, simulation_monitor=SM)

# Set the simulation data type: if not set "Double" is default
SimulationDataType.set_datatype(Double)
mixed_precision = False
if mixed_precision:
    # Define custom precision options
    mixed_precision_config = {
    'q_vector' : ([], FloatC),
    'RK_arrays' : ([], FloatC),
    'casting' : False # Explicit casting strategy
    }
    # Write the OPSC compatible code for the numerical solution
    OPSC(alg, OPS_V2=True, mixed_precision_config=mixed_precision_config)
else:
    OPSC(alg, OPS_V2=True)
# Add the simulation constants to the OPS C code
substitute_simulation_parameters(simulation_parameters.keys(), simulation_parameters.values())
print_iteration_ops(NaN_check='rho')
