#!/usr/bin/env python
from opensbli import *
from sympy import sin, exp, pi, tan, cos
import copy
from opensbli.multiblock.algorithm import TraditionalAlgorithmRKMB
from opensbli.postprocess.airfoil import *
from sympy.functions.elementary.piecewise import Piecewise, ExprCondPair
import os
# Disable the gmpy library for this case to avoid deepcopy issues
os.environ['MPMATH_NOGMPY'] = '1'

# Simulation parameters
simulation_parameters = {
# Physical parameters
'gama'      :   '1.4',
'Minf'      :   '0.72',
'Pr'        :   '0.71',
'Re'        :   '5.0e5',
'dt'        :   '5.0e-5',
'niter'     :   '1000000',
'sigma_filt'        :   '0.01',
'SuthT'     :   '110.4',
'RefT'      :   '273.15',
'stat_frequency'        :   '100',
'shock_factor'      :   '1.0',
'Twall'     :   '1.0',
'Lz'        :   '0.05',
'span_factor'       :   'Lz/0.05',
# Block 0
'block0np0'     :   '701',
'block0np1'     :   '681',
'block0np2'     :   '50',
'Delta0block0'      :   '5.0/(block0np0 - 1.0)',
'Delta1block0'      :   '22.5/(block0np1 - 1.0)',
'Delta2block0'      :   '0.05/block0np2',
'inv_rfact0_block0'     :   '1.0/Delta0block0',
'inv_rfact1_block0'     :   '1.0/Delta1block0',
'inv_rfact2_block0'     :   '1.0/Delta2block0',
# Block 1
'block1np0'     :   '2249',
'block1np1'     :   '681',
'block1np2'     :   '50',
'Delta0block1'      :   '2.0461756979465546/(block1np0 - 1.0)',
'Delta1block1'      :   '22.5/(block1np1 - 1.0)',
'Delta2block1'      :   '0.05/block1np2',
'inv_rfact0_block1'     :   '1.0/Delta0block1',
'inv_rfact1_block1'     :   '1.0/Delta1block1',
'inv_rfact2_block1'     :   '1.0/Delta2block1',
# Block 2
'block2np0'     :   '701',
'block2np1'     :   '681',
'block2np2'     :   '50',
'Delta0block2'      :   '5.0/(block2np0 - 1.0)',
'Delta1block2'      :   '22.5/(block2np1 - 1.0)',
'Delta2block2'      :   '0.05/block2np2',
'inv_rfact0_block2'     :   '1.0/Delta0block2',
'inv_rfact1_block2'     :   '1.0/Delta1block2',
'inv_rfact2_block2'     :   '1.0/Delta2block2',
# Add forcing modes
'tripA'     :   '0.075', # trip amplitude
'tripSigma'     :   '0.00833',
'xts'       :   '0.1',
'xtp'       :   '0.1',
'omega_0'       :   '26', # temporal frequency
'omega_1'       :   '88',
'omega_2'       :   '200',
'k_0'       :   '3.0', # wavenumbers
'k_1'       :   '4.0',
'k_2'       :   '4.0',
'phi_0'     :   '0.0', # phase shift
'phi_1'     :   'M_PI',
'phi_2'     :   '-M_PI/2',
}

# Define the problem
import itertools
def create_exchange_calls_codes(multiblock_descriptor, dsets):
    kernels = []
    for block in multiblock_descriptor.blocks:
        arrays = [block.location_dataset(a) for a in flatten(dsets)]
        kernels += block.apply_interface_bc(arrays, multiblock_descriptor, full_halo_swap=True)
    return kernels

ndim = 3
nblocks = 3
# Set non-conservative LHS to reduce array storage
conservative = True
multi_block = MultiBlock(ndim, nblocks, conservative=conservative)
SimulationDataType.set_datatype(Double)

# # Constants that are used
constants = ["Re", "Pr", "gama", "Minf", "RefT", "SuthT"]
# Define coordinate direction symbol (x) this will be x_i, x_j, x_k
coordinate_symbol = "x"
metriceq = MetricsEquation()
metriceq.generate_transformations(ndim, coordinate_symbol, [(True, True), (True, True), (False, False)], 2, latex_debug=False)
#Create an optional substitutions dictionary, this will be used to modify the equations when parsed
optional_subs_dict = metriceq.metric_subs
Einstein_expansion = EinsteinEquation()
Einstein_expansion.optional_subs_dict = optional_subs_dict
metric_vel = "Eq(U_i, D_i_j*u_j)"
eqns = Einstein_expansion.expand(metric_vel, ndim, coordinate_symbol, [], constants)
for eq in eqns:
    Einstein_expansion.optional_subs_dict[eq.lhs] = eq.rhs

NS = NS_Split('KGP', ndim, constants, coordinate_symbol=coordinate_symbol, conservative=conservative, viscosity='dynamic', energy_formulation='enthalpy', debug=False)

mass, momentum, energy = NS.mass, NS.momentum, NS.energy
# Expand the simulation equations, for this create a simulation equations class
simulation_eq = SimulationEquations()
simulation_eq.add_equations(mass)
simulation_eq.add_equations(momentum)
simulation_eq.add_equations(energy)

# Formulas for the variables used in the equations
constituent = ConstituentRelations()
if conservative:
    velocity = "Eq(u_i, rhou_i/rho)"
    eqns = Einstein_expansion.expand(velocity, ndim, coordinate_symbol, [], constants)
    constituent.add_equations(eqns)
    pressure = "Eq(p, (gama-1)*(rhoE - (1/2)*rho*(KD(_i,_j)*u_i*u_j)))"
    enthalpy = "Eq(H, (rhoE + p) / rho)"
else:
    pressure = "Eq(p, (gama-1)*(Et - (1/2)*(KD(_i,_j)*u_i*u_j)))"
    enthalpy = "Eq(H, Et + p / rho)"

temperature = "Eq(T, p*gama*Minf*Minf/(rho))"
viscosity = "Eq(mu, (T**(1.5)*(1.0+SuthT/RefT)/(T+SuthT/RefT)))"

eqns = Einstein_expansion.expand(pressure, ndim, coordinate_symbol, [], constants)
constituent.add_equations(eqns)
eqns = Einstein_expansion.expand(temperature, ndim, coordinate_symbol, [], constants)
constituent.add_equations(eqns)
eqns = Einstein_expansion.expand(viscosity, ndim, coordinate_symbol, [], constants)
constituent.add_equations(eqns)
eqns = Einstein_expansion.expand(enthalpy, ndim, coordinate_symbol, [], constants)
constituent.add_equations(eqns)

# Transform the equations into curvilinear form
simulation_eq.apply_metrics(metriceq)
# Specify the numerical schemes
schemes = {}
rk = RungeKuttaLS(4)
schemes[rk.name] = rk
# cent = Central(4)
cent = StoreSome(4, 'u0 u1 u2 T', merged=True)
schemes[cent.name] = cent
multi_block.set_discretisation_schemes(schemes)

# Initial conditions
## Need to change for swept cases
d, u0, u1, u2, p = symbols("d, u0:3, p", **{'cls':GridVariable})
gama, Minf = symbols("gama, Minf", **{'cls':ConstantObject})
initial_equations = []
initial_equations += [Eq(d, 1.0)]
initial_equations += [Eq(u0, 1.0)]
initial_equations += [Eq(u1, 0.0)]
initial_equations += [Eq(u2, 0.0)]
initial_equations += [Eq(p, 1.0/(gama*Minf**2.0))]

# Set the q vector values for initial condition and farfield boundaries
q_vector = flatten(simulation_eq.time_advance_arrays)
if conservative:
    initial_equations += [Eq(q_vector[0], d)]
    initial_equations += [Eq(q_vector[1], d*u0)]
    initial_equations += [Eq(q_vector[2], d*u1)]
    initial_equations += [Eq(q_vector[3], d*u2)]
    initial_equations += [Eq(q_vector[4], p/(gama-1.0) + 0.5* d *(u0**2+u1**2+ u2**2))]
else:
    initial_equations += [Eq(q_vector[0], d)]
    initial_equations += [Eq(q_vector[1], u0)]
    initial_equations += [Eq(q_vector[2], u1)]
    initial_equations += [Eq(q_vector[3], u2)]
    initial_equations += [Eq(q_vector[4], p/(d*(gama-1.0)) + 0.5*(u0**2+u1**2+ u2**2))]

temp_x2 = [Eq(GridVariable('temp'), DataObject('x2'))]
initial = GridBasedInitialisation()
initial.add_equations(copy.deepcopy(initial_equations) + temp_x2)
multi_block.set_equations([initial])

# block 0 boundary conditions
mb_bcs = {0:None, 1:None, 2:None}
# Boundary conditions for block 0
# The boundary conditions are [InterfaceBC, outflow] in x0 direction and [SharedInterfaceBC, Inflow]  in x1 direction 
# Matching boundaries are located at are [1,0,0] and [2, 1, 0]
block0_bc = []
direction = 0
side = 0
block0_bc.append(InterfaceBC(direction=0, side=0,  halos=[-2,2], name="block0_to_block1", match=(1, 0, 0, True)))
block0_bc.append(ExtrapolationBC(direction=0, side=1, order=0))
block0_bc.append(SharedInterfaceBC(direction=1, side=0, halos=[-4,4], name="block0_to_block2", match=(2, 1, 0, True)))
block0_bc.append(DirichletBC(direction=1, side=1, equations=initial_equations))
block0_bc.append(PeriodicBC(direction=2, side=0, halos=[-2,2], corners=False))
block0_bc.append(PeriodicBC(direction=2, side=1, halos=[-2,2], corners=False))
mb_bcs[0] = block0_bc

# Boundary conditions for block 1
#The boundary conditions are [InterfaceBC, InterfaceBC] in x0 direction and [wall, Inflow]  in x1 direction 
# Matching boundaries are located at are [0,0,0] and [2, 0, 0]
block1_bc = []
block1_bc.append(InterfaceBC(direction=0, side=0, halos=[-2,2], name="block1_to_block0", match=(0, 0, 0, True)))
block1_bc.append(InterfaceBC(direction=0, side=1, halos=[-2,2], name="block1_to_block2", match=(2, 0, 0, False)))

## Wall condition
# Isothermal wall in x1 direction
gama, Minf, Twall = symbols('gama Minf Twall', **{'cls': ConstantObject})
# Boundary-layer tripping
tripped = True
direction, side = 1, 0
if tripped:
    from sympy import pi
    Amp, sigma, xts, xtp, Lz, span_factor = symbols('tripA tripSigma xts xtp Lz span_factor', **{'cls':ConstantObject})
    # Time dependence
    # current_iter = multi_block.get_block(nblocks-1).get_temporal_schemes[0].iteration_number # Current iteration number
    current_iter = Globalvariable("iter", integer=True)
    dt, omega0, omega1, omega2 = symbols('dt omega_0 omega_1 omega_2', **{'cls': ConstantObject})
    t = dt*current_iter 
    # Spatial Modes
    k0, k1, k2 = symbols('k_0 k_1 k_2', **{'cls': ConstantObject})
    phi0, phi1, phi2 = symbols('phi_0 phi_1 phi_2', **{'cls': ConstantObject})
    # Coordinate arrays
    x0, z0 = DataObject('x0'), DataObject('x2')
    conditional_expressions = []
    # Suction side trip
    SS_trip = Amp*(exp(-(x0 - xts)**2  / (2*sigma**2))*(sin(2*pi*k0 * z0 * span_factor/Lz)*sin(omega0*t + phi0) + sin(2*pi*k1 * z0 * span_factor/Lz)*sin(omega1*t + phi1) + sin(2*pi*k2 * z0 * span_factor/Lz)*sin(omega2*t + phi2)))
    # Pressure side trip
    PS_trip = Amp*(exp(-(x0 - xtp)**2  / (2*sigma**2))*(sin(2*pi*k0 * z0 * span_factor/Lz)*sin(omega0*t + phi0) + sin(2*pi*k1 * z0 * span_factor/Lz)*sin(omega1*t + phi1) + sin(2*pi*k2 * z0 * span_factor/Lz)*sin(omega2*t + phi2)))
    # Index in x direction, assuming anti-clockwise grid configuration here
    idx = multi_block.get_block(1).grid_indexes[0] # x index on block 1 (airfoil block)
    expr_condition_pairs = Piecewise((SS_trip, idx > ConstantObject('block1np0')/2), (PS_trip, idx < ConstantObject('block1np0')/2),  (0, True))
    v = OpenSBLIEq(GridVariable('v'), expr_condition_pairs)
    rhov_wall = OpenSBLIEq(DataObject('rhou1'), DataObject('rho')*GridVariable('v'))
    # Adding the non-zero V component to the calculation of rhoE at the wall
    rhoE_wall = OpenSBLIEq(DataObject('rhoE'), DataObject('rho')*Twall/(gama*(gama-1.0)*Minf**2.0) + 0.5*DataObject('rho')*GridVariable('v')**2)
    # Equations to set rhov and rhoE on the wall
    wall_eqns = [rhov_wall, rhoE_wall]
    wall_energy = [rhoE_wall]
    block1_bc.append(ForcingStripBC(direction, 0, v, wall_eqns, corners=False, multi_block=True))
else:
    # Energy on the wall is set
    if conservative:
        wall_energy = [Eq(q_vector[-1], q_vector[0]*Twall / (gama * Minf**2.0 * (gama - S.One)))]
    else:
        wall_energy = [Eq(q_vector[-1], Twall / (gama * Minf**2.0 * (gama - S.One)))]
    direction = 1
    lower_wall_eq = wall_energy[:]

    block1_bc.append(IsothermalWallBC(direction=1, side=0, corners=False, equations=lower_wall_eq, multi_block=True))


## Farfield
block1_bc.append(DirichletBC(direction=1, side=1, equations=initial_equations))
block1_bc.append(PeriodicBC(direction=2, side=0, halos=[-2,2], corners=False))
block1_bc.append(PeriodicBC(direction=2, side=1, halos=[-2,2], corners=False))
mb_bcs[1] = block1_bc

# Boundary conditions for block 2
# The boundary conditions are [InterfaceBC, outflow] in x0 direction and  SharedInterfaceBC, Inflow]  in x1 direction 
# Matching boundaries are located at are [1,0,1] and [0, 1, 0]
block2_bc = []
block2_bc.append(InterfaceBC(direction=0, side=0,  halos=[-2,2], name="block2_to_block1", match=(1, 0, 1, False)))
block2_bc.append(ExtrapolationBC(direction=0, side=1, order=0))
block2_bc.append(SharedInterfaceBC(direction=1, side=0,  halos=[-4,4], name="block2_to_block0", match=(0, 1, 0, True)))
block2_bc.append(DirichletBC(direction=1, side=1, equations=initial_equations))
block2_bc.append(PeriodicBC(direction=2, side=0, halos=[-2,2], corners=False))
block2_bc.append(PeriodicBC(direction=2, side=1, halos=[-2,2], corners=False))
mb_bcs[2] = block2_bc
# Set the multi block boundary conditions
multi_block.set_block_boundaries(mb_bcs)

# Set the equations on the blocks
multi_block.set_equations([simulation_eq, constituent, metriceq])

# Add statsistics gathering
stats = False
if stats:
    # Create the statistics equations, this shows another way of writing the equations
    from airfoil_stats import favre_averaged_stats
    q_vector = flatten(simulation_eq.time_advance_arrays)
    stat_equation_classes, stats_arrays = favre_averaged_stats(ndim, q_vector, conservative=conservative)
else:
    stat_equation_classes, stats_arrays = [], []

multi_block.set_equations(stat_equation_classes)

# Add filters to each block
filters = {0:[], 1:[], 2:[]}
shock_filters = []
for no, block in enumerate(multi_block.blocks):
    if no == 1: # Main aerofoil block, C-mesh. Don't filter near the aerofoil
        WF = WENOFilter(block, order=5, metrics=metriceq, airfoil=True, optimize=True, flux_type='LLF')
        shock_filters.append(WF)
        filters[no] += [WF.equation_classes]

# Add DRP filters for freestream
for no, block in enumerate(multi_block.blocks):
    filters[no] += [ExplicitFilter(block, [0,1,2], width=9, filter_type='DRP', optimized=False, sigma=0.3333333, airfoil=True, multi_block=multi_block).equation_classes]

# Add a binomial filter on the outlet boundary to kill reflections
for no, block in enumerate(multi_block.blocks):
    i, j, k = block.grid_indexes[0], block.grid_indexes[1], block.grid_indexes[2]
    if no == 0:
        grid_condition = i >= 794
        filters[no] += [BinomialFilter(block, order=6, grid_condition=grid_condition).equation_classes]
        grid_condition = j >= 475
        filters[no] += [BinomialFilter(block, order=6, grid_condition=grid_condition).equation_classes]
    elif no == 1:
        grid_condition = j >= 475
        filters[no] += [BinomialFilter(block, order=6, grid_condition=grid_condition).equation_classes]
    elif no == 2:
        grid_condition = i >= 794
        filters[no] += [BinomialFilter(block, order=6, grid_condition=grid_condition).equation_classes]
        grid_condition = j >= 475
        filters[no] += [BinomialFilter(block, order=6, grid_condition=grid_condition).equation_classes]
    
multi_block.set_filters(filters)

# HDF5 input/output
x,y,z = symbols("x0, x1, x2", **{'cls':DataObject})
kwargs = {'iotype': "Write"}
q_hdf5 = iohdf5(save_every=50000, **kwargs)
q_hdf5.add_arrays(simulation_eq.time_advance_arrays)
# Read in the grid file
kwargs = {'iotype': "Read"}
grid_hdf5 = iohdf5(**kwargs)
grid_hdf5.add_arrays([x, y, z])
# Stats HDF5 and write metrics to the grid file
metrics_hdf5 = iohdf5(arrays=metriceq.grid_der_wks, **{'position': "init", 'iotype': 'Write', 'name': "metrics.h5"})
# Set the I/O on the block
multi_block.setio([q_hdf5, grid_hdf5, metrics_hdf5])

# Block 0 slicing
grid_slice_hdf5_0 = iohdf5_slices(blocknumber=0, **{'iotype': "Init"})
# x-y view
coords = [([DataObject('x0'), DataObject('x1')], 2, 'block0np2/2')]
grid_slice_hdf5_0.add_slices(coords)
# Q vector slices written out in time
slices_hdf5_0 = iohdf5_slices(save_every=1000, blocknumber=0, **{'iotype': "Write"})
# x-y view
slices = [(q_vector, 2, 'block0np2/2')]
slices_hdf5_0.add_slices(slices)
# Block 2 slicing
grid_slice_hdf5_2 = iohdf5_slices(blocknumber=2, **{'iotype': "Init"})
# x-y view
coords = [([DataObject('x0'), DataObject('x1')], 2, 'block2np2/2')]
grid_slice_hdf5_2.add_slices(coords)
# Q vector slices written out in time
slices_hdf5_2 = iohdf5_slices(save_every=1000, blocknumber=2, **{'iotype': "Write"})
# Surface above the cylinder, 5 points off the wall
slices = [(q_vector, 2, 'block2np2/2')]
slices_hdf5_2.add_slices(slices)

# Block 1 slicing (airfoil block)
grid_slice_hdf5_1 = iohdf5_slices(blocknumber=1, **{'iotype': "Init"})
# x-y view
coords = [([DataObject('x0'), DataObject('x1')], 2, 'block1np2/2')]
# Airfoil surface coordinates
coords += [([DataObject('x0'), DataObject('x2')], 1, 0)]
coords += [([DataObject('x0'), DataObject('x2')], 1, 1)]
coords += [([DataObject('x0'), DataObject('x2')], 1, 2)]
coords += [([DataObject('x0'), DataObject('x2')], 1, 3)]
coords += [([DataObject('x0'), DataObject('x2')], 1, 4)]
coords += [([DataObject('x0'), DataObject('x2')], 1, 5)]

grid_slice_hdf5_1.add_slices(coords)
# Q vector slices written out in time
slices_hdf5_1 = iohdf5_slices(save_every=1000, blocknumber=1, **{'iotype': "Write"})
# x-y view
slices = [(q_vector, 2, 'block1np2/2')]
# Airfoil 
slices += [(q_vector, 1, 0)]
slices += [(q_vector, 1, 1)]
slices += [(q_vector, 1, 2)]
slices += [(q_vector, 1, 3)]
slices += [(q_vector, 1, 4)]
slices += [(q_vector, 1, 5)]
slices_hdf5_1.add_slices(slices)

# Set both HDF5 slicing objects
multi_block.setio([grid_slice_hdf5_0, slices_hdf5_0, grid_slice_hdf5_1, slices_hdf5_1, grid_slice_hdf5_2, slices_hdf5_2])

# Perform the discretization
multi_block.discretise()

# Add a periodic boundary condition call for WENO filters
for i, block in enumerate(multi_block.blocks):
    if i == 1:
        shock_filters[0].update_periodic_boundary(block, halos=[-4,4])

# Add the wake treatment kernelss
wake_ker = generate_wake_kernel(q_vector, multi_block, wall_energy[0])
# Sponge zones for outer boundaries
# Outlet
outlet_sponge_block0 = generate_outlet_sponge(q_vector, multi_block.get_block(0), Lx=5.0, npoints=10)
outlet_sponge_block2 = generate_outlet_sponge(q_vector, multi_block.get_block(2), Lx=5.0, npoints=10)
# Farfield
farfield_sponge_block0 = generate_farfield_sponge(q_vector, multi_block.get_block(0), Ly=22.5, npoints=10)
farfield_sponge_block1 = generate_farfield_sponge(q_vector, multi_block.get_block(1), Ly=22.5, npoints=10)
farfield_sponge_block2 = generate_farfield_sponge(q_vector, multi_block.get_block(2), Ly=22.5, npoints=10)

# Add wake exchanges and kernels to block2 boundary conditions
b = multi_block.get_block(2)
for no, eq in enumerate(b.list_of_equation_classes):
    # Add the sponge kernels after updating the residuals for all blocks, before time advancement
    if isinstance(eq, SimulationEquations):
        eq.Kernels += [outlet_sponge_block0, outlet_sponge_block2, farfield_sponge_block0, farfield_sponge_block1, farfield_sponge_block2]
        eq.boundary_kernels += wake_ker

# Make some full swaps for interfaces before filtering
filter_swaps = create_exchange_calls_codes(multi_block, simulation_eq.time_advance_arrays)
for block in multi_block.blocks:
    for no, eq in enumerate(block.list_of_equation_classes):
        if isinstance(eq, UserDefinedEquations):
            if eq.full_swap:
                eq.Kernels += filter_swaps

# Add some full [-5,5] halo swaps over the periodic directions only when the filter is called
def create_periodic_BCs(multi_block, dsets):
    kernels = []
    for block in multi_block.blocks:
        arrays = [block.location_dataset(a) for a in flatten(dsets)]
        for direction in [2]:
            for side in [0,1]:
                BC = PeriodicBC(direction, side, halos=[-4,4], corners=False)
                kernels += [BC.apply(arrays, block)]
    return kernels

# Periodic boundary condition for DRP filters
if conservative:
    dsets = ['rho', 'rhou0', 'rhou1', 'rhou2', 'rhoE']
else:
    dsets = ['rho', 'u0', 'u1', 'u2', 'Et']
DRP_periodic = create_periodic_BCs(multi_block, dsets)
for block in multi_block.blocks:
    for no, eq in enumerate(block.list_of_equation_classes):
        if isinstance(eq, UserDefinedEquations):
            if eq.order == 0:
                eq.Kernels += DRP_periodic

# Create the OPS C code
alg = TraditionalAlgorithmRKMB(multi_block)
OPSC(alg, OPS_diagnostics=1)
# Add the simulation constants to the OPS C code
substitute_simulation_parameters(simulation_parameters.keys(), simulation_parameters.values())
print_iteration_ops(NaN_check='rho', every=100, nblocks=nblocks)
