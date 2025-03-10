#!/usr/bin/env python

# Import all the functions from opensbli
from opensbli import *
from opensbli.utilities.helperfunctions import substitute_simulation_parameters

# Input parameters for the simulation
simulation_parameters = {
'Re'	: '90.0',
'gama'	: '1.4',
'Minf'	: '0.01',
'Pr'	: '0.72',
'dt'	: '0.0001',
'niter'	: '100000',
'block0np0'	: '16',
'block0np1'	: '81',
'Delta0block0'	: '2.0*M_PI/block0np0',
'Delta1block0'	: '2.0/(block0np1-1)',
"c0"	: '-1',
"c1"	: '0',
"SuthT"	: "110.4",
"RefT"	: "273.0",
"Twall"	: "1.0",
}

# Problem dimension
ndim = 2
# # Constants that are used
constants = ["Re", "Pr", "gama", "Minf", "c_j", "RefT", "SuthT"]
# # symbol for the coordinate system in the equations
coordinate_symbol = "x"
# symbol for the coordinate system in the equations
conservative = True
NS = NS_Split('KGP', ndim, constants, coordinate_symbol=coordinate_symbol, conservative=conservative, viscosity='dynamic')

mass, momentum, energy = NS.mass, NS.momentum, NS.energy
# Add channel forcing term
for i, eqn in enumerate(momentum):
	momentum[i] = OpenSBLIEq(eqn.lhs, eqn.rhs - ConstantObject("c%d" % i))
energy = OpenSBLIEq(energy.lhs, energy.rhs - (ConstantObject("c0")*DataObject('u0')))
# Expand the simulation equations, for this create a simulation equations class
simulation_eq = SimulationEquations()
simulation_eq.add_equations(mass)
simulation_eq.add_equations(momentum)
simulation_eq.add_equations(energy)
# Constituent relations used in the system
velocity = "Eq(u_i, rhou_i/rho)"
pressure = "Eq(p, (gama-1)*(rhoE - rho*(1/2)*(KD(_i,_j)*u_i*u_j)))"
temperature = "Eq(T, p*gama*Minf*Minf/(rho))"
# exact solution for mu=const
viscosity = "Eq(mu, 1.0)"
#viscosity = "Eq(mu, (T**(1.5)*(1.0+SuthT/RefT)/(T+SuthT/RefT)))"

# Instantiate EinsteinEquation class for expanding the Einstein indices in the equations
einstein_eq = EinsteinEquation()
# Expand the constituent relations and them to the constituent relations class
constituent = ConstituentRelations()  # Instantiate constituent relations object
# Expand momentum add the expanded equations to the constituent relations
eqns = einstein_eq.expand(velocity, ndim, coordinate_symbol, [], constants)
constituent.add_equations(eqns)
# Expand pressure add the expanded equations to the constituent relations
eqns = einstein_eq.expand(pressure, ndim, coordinate_symbol, [], constants)
constituent.add_equations(eqns)
# Expand temperature add the expanded equations to the constituent relations
eqns = einstein_eq.expand(temperature, ndim, coordinate_symbol, [], constants)
constituent.add_equations(eqns)
# Expand temperature add the expanded equations to the constituent relations
eqns = einstein_eq.expand(viscosity, ndim, coordinate_symbol, [], constants)
constituent.add_equations(eqns)

# Create a simulation block
block = SimulationBlock(ndim, block_number=0, conservative=conservative)
# Define the variables used for creating boundary conditions and the initialisation
# dx and dy of the grid
dx, dy = block.deltas
# Indices for the grid location
i, j = block.grid_indexes
# Some constants used
gama, Re, Pr, Minf = symbols('gama Re Pr Minf ', **{'cls': ConstantObject})
""" Conservative vector is the time advancement arrays of the simulation equations.
the order follows the order in which they are added to the simulation equations
class, i.e. arrays of density, momentum (components), energy in the present case"""
q_vector = flatten(simulation_eq.time_advance_arrays)

# STEP 1
# Set the boundary conditions on the block
boundaries = []
# For laminar channel flow case the boundaries are periodic in x and walls in y
# Periodic boundaries in x0 direction
direction = 0
boundaries += [PeriodicBC(direction, side=0)]
boundaries += [PeriodicBC(direction, side=1)]

# Isothermal wall in x1 direction
local_dict = {"block": block, "GridVariable": GridVariable, "DataObject": DataObject}

# Energy on the wall is set
wall_const = ["Minf", "Twall"]
for con in wall_const:
    local_dict[con] = ConstantObject(con)
# Isothermal wall condition
rhoE_wall = [parse_expr("Eq(DataObject(rhoE), DataObject(rho)*Twall/(gama*(gama-1.0)*Minf**2.0))", local_dict=local_dict)]
# Side 0 (bottom wall) boundary
boundaries += [IsothermalWallBC(direction=1, side=0, equations=rhoE_wall)]

# Side 1 (top) boundary
# Select one adiabatic wall one isothermal wall. Otherwise two isothermal walls
mixed_wall_condition = False

if mixed_wall_condition:
	boundaries += [AdiabaticWallBC(direction=1, side=1)]
else:
	boundaries += [IsothermalWallBC(direction=1, side=1, equations=rhoE_wall)]

# set the boundaries for the block
block.set_block_boundaries(boundaries)

# The equation classes used for the block, these are
# simulation equations, constituent relations, grid and initial conditions

# Create the grid and intial conditions
# Arrays to store x and y coordinates, i.e (x0 and x1)
x, y = symbols('x0:%d' % ndim, **{'cls': DataObject})
# Equations for generating the grid, simple equispacing grid
grid_equations = [Eq(x, i * dx), Eq(y, -1.0 + j * dy)]

# Initialisation equations
initial_equations = []
# local varibales for temperature and pressure
density, pressure, ulam = symbols('r p u0', **{'cls': GridVariable})
# Equations for pressure and temperature
initial_equations += [Eq(ulam,0.5*Re*(1.0-y**2))]
initial_equations += [Eq(pressure, 1.0 / (gama * Minf**2.0))]
initial_equations += [Eq(density, 1.0/(1.0 + (gama-1.0)*Pr*Re**2*Minf**2*(1.0-y**4)/12.0))]

# Initialise the conservative vector
initial_equations += [Eq(q_vector[0], density)]
initial_equations += [Eq(q_vector[1], density*ulam)]
initial_equations += [Eq(q_vector[2], 0.0)]
initial_equations += [Eq(q_vector[3], pressure/(gama - 1.0)+0.5*density*ulam**2)]

# Instantiate a grid based initialisation classes
initial = GridBasedInitialisation()
initial.add_equations(grid_equations + initial_equations)

# STEP 2
# Set the equation classes for the block (list)
block.set_equations([constituent, simulation_eq, initial])


# STEP 3
# Create the dictionary of schemes
schemes = {}
# Central scheme for spatial discretisation and add to the schemes dictionary
fns = 'u0 u1 T'
# cent = Central(4)
cent = StoreSome(4, fns)
schemes[cent.name] = cent
# RungeKutta scheme for temporal discretisation and add to the schemes dictionary
rk = RungeKuttaLS(3)
schemes[rk.name] = rk
# Set the discretisation schemes to be used (a python dictionary)
block.set_discretisation_schemes(schemes)

# STEP 4 add io for the block
kwargs = {'iotype': "Write"}
output_arrays = simulation_eq.time_advance_arrays + [x, y]
output_hdf5 = iohdf5(arrays=output_arrays, **kwargs)
block.setio([output_hdf5])

# STEP 6
# Perform the symbolic discretisation of the equations
block.discretise()

# STEP 7
# create an algorithm from the numerical solution
alg = TraditionalAlgorithmRK(block)

# STEP 8
# set the simulation data type: if not set "Double" is default
SimulationDataType.set_datatype(Double)

# STEP 9
# Write the OPSC compatible code for the numerical solution
OPSC(alg)

# STEP 10
# Add the simulation constants to the OPS C code
substitute_simulation_parameters(simulation_parameters.keys(), simulation_parameters.values())
print_iteration_ops(NaN_check='rho')
