#!/usr/bin/env python
from opensbli import *
import copy
from opensbli.utilities.katzer_init import Initialise_Katzer
from opensbli.utilities.helperfunctions import substitute_simulation_parameters, output_hdf5, print_iteration_ops
from opensbli.utilities.oblique_shock import ShockConditions
from sympy import tan, pi, tanh, sinh, cosh, exp, cos, sin
import time_averaging
from forcing import forcing_setup

simulation_parameters = {
'gama'      :   '1.4',
'Minf'      :   '1.5',
'Pr'        :   '0.72',
'Re'        :   '750.0',
'Twall'     :   '1.3809973268575328',
'dt'        :   '0.01',
'niter'     :   '300000',
'Lx0'       :   '375.0',
'Lx1'       :   '140.0',
'Lx2'       :   '27.32',
'block0np0'     :   '2050',
'block0np1'     :   '325',
'block0np2'     :   '200',
'Delta0block0'      :   'Lx0/(block0np0-1)',
'Delta1block0'      :   'Lx1/(block0np1-1)',
'Delta2block0'      :   'Lx2/(block0np2)',
'SuthT'     :   '110.4',
'RefT'      :   '202.17',
'eps'       :   '1e-30',
'bx'        :   '1.5',
'by'        :   '5.0',
'TENO_CT'   :   '1e-7',
'epsilon'       :   '1.0e-16',
}

# Select forcing method examples: 1. Single most unstable mode as acoustic freestream perturbation or 2. Many modes forced at the wall
single_mode = True

single_mode_parameters = {
'bta'       :   '0.23',
'omega'     :   '0.1011',
'xF'        :   '20.0',
'yF'        :   '1.0',
'A'     :   '0.1',
}

multi_mode_parameters = {
'xa' : '10.0', # start and end of tripping strip in x
'xb' : '30.0',
'omega' : '0.12313',
'A' : '0.25',
}

# Instatiate equation classes
eq = EinsteinEquation()
# Define coordinate direction symbol (x) this will be x_i, x_j, x_k
coordinate_symbol = "x"
constants = ["Re", "Pr", "gama", "Minf", "SuthT", "RefT"]
# Define the problem
ndim = 3
# Pure TENO, pure WENO
teno = True
weno = False
# TVD or WENO filter method
TVD  = False
stress_tensor = "Eq(tau_i_j, (mu/Re)*(Der(u_i,x_j)+ Der(u_j,x_i) - (2/3)* KD(_i,_j)* Der(u_k,x_k)))"
heat_flux = "Eq(q_j, (-mu/((gama-1)*Minf*Minf*Pr*Re))*Der(T,x_j))"
# Substitutions
substitutions = [stress_tensor, heat_flux]
if teno or weno:
    if teno:
        sc1 = "**{\'scheme\':\'Teno\'}"
    else:
        sc1 = "**{\'scheme\':\'Weno\'}"
    # Define the compresible Navier-Stokes equations in Einstein notation.
    mass = "Eq(Der(rho,t), - Conservative(rhou_j,x_j,%s))" % sc1
    momentum = "Eq(Der(rhou_i,t) , -Conservative(rhou_i*u_j + KD(_i,_j)*p,x_j , %s) + Der(tau_i_j,x_j) )" % sc1
    energy = "Eq(Der(rhoE,t), - Conservative((p+rhoE)*u_j,x_j, %s) - Der(q_j,x_j) + Der(u_i*tau_i_j ,x_j) )" % sc1
    base_eqns = [mass, momentum, energy]
    # Expand the base equations
    for i, base in enumerate(base_eqns):
        base_eqns[i] = eq.expand(base, ndim, coordinate_symbol, substitutions, constants)
else:
    NS = NS_Split('KGP', ndim, constants, coordinate_symbol=coordinate_symbol, conservative=True, viscosity='dynamic', energy_formulation='enthalpy', debug=False)
    mass, momentum, energy = NS.mass, NS.momentum, NS.energy
    base_eqns = [mass, momentum, energy]

    
# Formulas for the variables used in the equations
velocity = "Eq(u_i, rhou_i/rho)"
pressure = "Eq(p, (gama-1)*(rhoE - rho*(1/2)*(KD(_i,_j)*u_i*u_j)))"
speed_of_sound = "Eq(a, (gama*p/rho)**0.5)"
temperature = "Eq(T, p*gama*Minf*Minf/(rho))"
viscosity = "Eq(mu, (T**(1.5)*(1.0+SuthT/RefT)/(T+SuthT/RefT)))"
constituent_eqns = [velocity, pressure, speed_of_sound, temperature, viscosity]
# Expand the constituent relations
for i, CR in enumerate(constituent_eqns):
    constituent_eqns[i] = eq.expand(CR, ndim, coordinate_symbol, substitutions, constants)
if not weno and not teno:
    ## append constituent eqns for split_forms
    internal_energy = "Eq(e, p / (rho*(gama-1)))"
    enthalpy = "Eq(H, (rhoE + p) / rho)"
    eqn1 = eq.expand(internal_energy, ndim, coordinate_symbol, [], constants)
    eqn2 = eq.expand(enthalpy, ndim, coordinate_symbol, [], constants)
    constituent_eqns.append(eqn1)
    constituent_eqns.append(eqn2)
    halos = [-2, 2]
else:
    halos = [-4, 4]

# Create a simulation block
block = SimulationBlock(ndim, block_number=0)
# Local dictionary for the block
local_dict = {"block": block, "GridVariable": GridVariable, "DataObject": DataObject}

# Metric transformation
metriceq = MetricsEquation()
metriceq.generate_transformations(ndim, coordinate_symbol, [(True, False), (True, False), (False, False)], 2)

schemes = {}
if teno:
    teno_order = 6
    Avg = RoeAverage([0, 1])
    LF = LFTeno(teno_order, averaging=Avg, flux_type='LLF')
    schemes[LF.name] = LF
elif weno:
    weno_order = 7
    Avg = RoeAverage([0, 1])
    LF = LFWeno(weno_order, averaging=Avg, formulation='Z', flux_type='LLF')
    #LF = HLLCWeno(order=5, formulation='Z', averaging=Avg, flux_type='HLLC-LM')
    schemes[LF.name] = LF

# Central scheme
fns = 'u0 u1 u2 T'
cent = StoreSome(4, fns)
schemes[cent.name] = cent
rk = RungeKuttaLS(3, formulation='SSP')
schemes[rk.name] = rk
block.set_discretisation_schemes(schemes)

# Create SimulationEquations and Constituent relations, add the expanded equations
simulation_eq = SimulationEquations()
constituent = ConstituentRelations()

if single_mode: # one mode, acoustic source
    # Create the body force after specifying a time stepping scheme
    A, B, xF, yF, bta, omega, dt = symbols('A B xF yF bta omega dt', **{'cls': ConstantObject})
    forcing_const = ["A", "B", "xF", "yF", "bta", "omega", "dt"]

    current_iter = block.get_temporal_schemes[0].iteration_number
    x0, x1, x2 = symbols('x0 x1 x2', **{'cls': DataObject})
    # damping = (1 - exp(-x1 / 0.5))
    body_force = Eq(DataObject('BF'), A*exp(-((x0-xF)**2+(x1-yF)**2))*cos(bta*x2)*sin(omega*dt*current_iter))

    # Add forcing as an acoustic term to the continuity equation
    for i, eq in enumerate(base_eqns):
        if i == 0:
            base_eqns[i] = Eq(eq.lhs, eq.rhs + body_force.lhs)
    # Add the body forcing term to the constituent relations
    constituent.add_equations(body_force)


# Add the equations to solve
for eqn in base_eqns:
    simulation_eq.add_equations(eqn)
for eqn in constituent_eqns:
    constituent.add_equations(eqn)
# Apply metric transformation to the simulation equations
simulation_eq.apply_metrics(metriceq)

# Define the boundary conditions
boundaries = [[0, 0] for t in range(ndim)]
# Left pressure extrapolation at x= 0, inlet conditions
direction, side = 0, 0
if teno or weno:
    boundaries[direction][side] = InletPressureExtrapolateBC(direction, side)
else:
    boundaries[direction][side] = InletTransferBC(direction, side)
# Right extrapolation at outlet
direction, side = 0, 1
boundaries[direction][side] = ExtrapolationBC(direction, side, order=0)

# Bottom wall boundary
wall_const = ["Minf", "Twall"]
for con in wall_const:
    local_dict[con] = ConstantObject(con)
direction, side = 1, 0
if single_mode:
    rhoE_wall = parse_expr("Eq(DataObject(rhoE), DataObject(rho)*Twall/(gama*(gama-1.0)*Minf**2.0))", local_dict=local_dict)
    wall_eqns = [rhoE_wall]
    boundaries[direction][side] = IsothermalWallBC(direction, side, wall_eqns)
else: # wall forcing
    evaluations = forcing_setup(20, direction, block, 0)
    rhov_wall = parse_expr("Eq(DataObject(rhou1), DataObject(rho)*GridVariable(u_norm))", local_dict=local_dict)
    rhoE_wall = parse_expr("Eq(DataObject(rhoE), DataObject(rho)*Twall/(gama*(gama-1.0)*Minf**2.0) + 0.5*DataObject(rho)*GridVariable(u_norm)**2)", local_dict=local_dict)
    # Equations to set rhov and rhoE on the wall
    wall_eqns = [rhov_wall, rhoE_wall]
    # Bottom forced wall
    boundaries[direction][side] = ForcingStripBC(direction, side, evaluations, wall_eqns)

# Top dirichlet shock generator condition
direction, side = 1, 1
boundaries[direction][side] = ExtrapolationBC(direction, side, order=0)
# Periodic direction 2
direction = 2
for side in [0,1]:
    boundaries[direction][side] = PeriodicBC(direction, side, halos=halos)
block.set_block_boundaries(boundaries)

# Generate initial condition, similarity solution and grid coordinates
Re, xMach, Tinf = 750.0, 1.5, 202.17
local_dict['Lx1'] = ConstantObject('Lx1')
local_dict['by'] = ConstantObject('by')
# Ensure the grid size passed to the initialisation routine matches the grid sizes used in the simulation parameters
grid_const = ["Lx0", "bx", "block0np0", "Lx1", "by", "Lx2"]
for con in grid_const:
    local_dict[con] = ConstantObject(con)
    CTD.add_constant(ConstantObject(con))
gridx0 = parse_expr("Eq(DataObject(x0), Lx0*(1.0 - sinh(bx*block.deltas[0]*(block0np0 - 1.0 - block.grid_indexes[0])/Lx0)/sinh(bx)))", local_dict=local_dict)
gridx1 = parse_expr("Eq(DataObject(x1), Lx1*sinh(by*block.deltas[1]*block.grid_indexes[1]/Lx1)/sinh(by))", local_dict=local_dict)
gridx2 = parse_expr("Eq(DataObject(x2), block.deltas[2]*block.grid_indexes[2])", local_dict=local_dict)
coordinate_evaluation = [gridx0, gridx1, gridx2]
polynomial_directions = [(False, DataObject('x0')), (True, DataObject('x1')), (False, DataObject('x2'))]
n_poly_coefficients = 50

initial = Initialise_Katzer(polynomial_directions, n_poly_coefficients, Re, xMach, Tinf, coordinate_evaluations=coordinate_evaluation)

# Set I/O options
q_vector = simulation_eq.time_advance_arrays
kwargs = {'iotype': "Write"}
h5 = iohdf5(arrays=q_vector, save_every=5000, **kwargs)
h5.add_arrays([DataObject('x0'), DataObject('x1'), DataObject('x2'), DataObject('D11')])
block.setio(h5)

## 2D HDF5 slicing output:
# Add 2D I/O slicing option
# Block 0 slicing
grid_slices = iohdf5_slices(blocknumber=0, **{'iotype': "Init"})
# x-y side view coordinates
coords = [([DataObject('x0'), DataObject('x1')], 2, 'block0np2/2')]
grid_slices.add_slices(coords)
# Q vector slices written out in time
slices_hdf5_side = iohdf5_slices(save_every=1000, blocknumber=0, **{'iotype': "Write"})
# x-y side view
slices = [(q_vector, 2, 'block0np2/2')]
slices_hdf5_side.add_slices(slices)
# Surface coordinates
coords = [([DataObject('x0'), DataObject('x2')], 1, 20)]
grid_slices.add_slices(coords)
# Q vector slices written out in time
slices_hdf5_surfaces = iohdf5_slices(save_every=1000, blocknumber=0, **{'iotype': "Write"})
# Wall normal planes 
slices = [(q_vector, 1, 20)]
slices_hdf5_surfaces.add_slices(slices)

# Set both HDF5 slicing objects
block.setio([grid_slices, slices_hdf5_side, slices_hdf5_surfaces])

# Settings to turn on statistics gathering 
stats = True
if stats:
    stats_class = time_averaging.get_stats_classes()
    stats_arrays = time_averaging.get_arrays()
    kwargs = {'iotype': "Write", "name": "stats_output.h5"}
    h5_stats = iohdf5(**kwargs)
    h5_stats.add_arrays(stats_arrays)
    block.setio(h5_stats)

# Set equations on the block and discretise
block.set_equations([simulation_eq, constituent, initial, metriceq] + stats_class)
block.discretise()
# Create an algorithm and write the OPS C code
alg = TraditionalAlgorithmRK(block)
SimulationDataType.set_datatype(Double)
OPSC(alg)
# Add the simulation constants to the OPS C code
if single_mode:
    simulation_parameters.update(single_mode_parameters)
else:
    simulation_parameters.update(multi_mode_parameters)

substitute_simulation_parameters(simulation_parameters.keys(), simulation_parameters.values())
print_iteration_ops(NaN_check='rho')
