#!/usr/bin/env python

# Import all the functions from opensbli
from opensbli import *
from sympy import sin, cos, sinh, tanh, exp, pi
#import copy
from opensbli.utilities.helperfunctions import substitute_simulation_parameters

simulation_parameters = {
'Re'        :   '200.0',   
'gama'      :   '1.4',   
'Minf'      :   '0.4',   
'Pr'        :   '0.72',   
'Sc'        :   '1.0',   
'dt'        :   '0.02',   
'niter'     :   '1000',   
'block0np0'     :   '81',   
'block0np1'     :   '121',   
'Delta0block0'      :   '10.0/(block0np0-1)',   
'Delta1block0'      :   '1.0/(block0np1-1)',   
'Ly'        :   '20.0',   
'stretch'       :   '2.0',
}

# Number of dimensions of the system to be solved
ndim = 2
stats = False
# Define the compresible Navier-Stokes equations in Einstein notation, by default the scheme is Central no need to
# Specify the schemes
mass = "Eq(Der(rho,t), - Skew(rho*u_j,x_j))"
momentum = "Eq(Der(rhou_i,t) , - Skew(rhou_i*u_j, x_j) - Der(p,x_i)  + Der(tau_i_j,x_j))"
energy = "Eq(Der(rhoE,t), - Skew(rhoE*u_j,x_j) - Conservative(p*u_j,x_j) - Der(q_j,x_j) + Der(u_i*tau_i_j ,x_j))"
scalar = "Eq(Der(rhof,t), - Skew(rhof*u_j,x_j) + Der(1.0/(Re*Sc)*Der(f,x_j),x_j))"

# Substitutions used in the equations
stress_tensor = "Eq(tau_i_j, (1.0/Re)*(Der(u_i,x_j)+ Der(u_j,x_i)- (2/3)* KD(_i,_j)* Der(u_k,x_k)))"
heat_flux = "Eq(q_j, -(1.0/((gama-1)*Minf*Minf*Pr*Re))*Der(T,x_j))"

substitutions = [stress_tensor, heat_flux]

# Constants that are used
constants = ["Re", "Pr", "Sc", "gama", "Minf", "mu"]

# symbol for the coordinate system in the equations
coordinate_symbol = "x"

# Constituent relations used in the system
velocity = "Eq(u_i, rhou_i/rho)"
pressure = "Eq(p, (gama-1)*(rhoE - rho*(1/2)*(KD(_i,_j)*u_i*u_j)))"
temperature = "Eq(T, p*gama*Minf*Minf/(rho))"
mixturefraction = "Eq(f, rhof/rho)"

# Instantiate EinsteinEquation class for expanding the Einstein indices in the equations
einstein_eq = EinsteinEquation()

# Expand the simulation equations, for this create a simulation equations class
simulation_eq = SimulationEquations()

# Expand mass and add the expanded equations to the simulation equations
eqns = einstein_eq.expand(mass, ndim, coordinate_symbol, substitutions, constants)
simulation_eq.add_equations(eqns)

# Expand momentum add the expanded equations to the simulation equations
eqns = einstein_eq.expand(momentum, ndim, coordinate_symbol, substitutions, constants)
simulation_eq.add_equations(eqns)

# Expand energy equation add the expanded equations to the simulation equations
eqns = einstein_eq.expand(energy, ndim, coordinate_symbol, substitutions, constants)
simulation_eq.add_equations(eqns)

# Expand energy equation add the expanded equations to the simulation equations
eqns = einstein_eq.expand(scalar, ndim, coordinate_symbol, substitutions, constants)
simulation_eq.add_equations(eqns)

# Expand the constituent relations and them to the constituent relations class
constituent = ConstituentRelations()  # Instantiate constituent relations object

# Expand momentum add the expanded equations to the constituent relations
eqns = einstein_eq.expand(velocity, ndim, coordinate_symbol, substitutions, constants)
constituent.add_equations(eqns)

# Expand pressure add the expanded equations to the constituent relations
eqns = einstein_eq.expand(pressure, ndim, coordinate_symbol, substitutions, constants)
constituent.add_equations(eqns)

# Expand temperature add the expanded equations to the constituent relations
eqns = einstein_eq.expand(temperature, ndim, coordinate_symbol, substitutions, constants)
constituent.add_equations(eqns)

# Expand temperature add the expanded equations to the constituent relations
eqns = einstein_eq.expand(mixturefraction, ndim, coordinate_symbol, substitutions, constants)
constituent.add_equations(eqns)

# Write the expanded equations to a Latex file with a given name and titile
latex = LatexWriter()
latex.open('equations.tex', "Einstein Expansion of the simulation equations")
latex.write_string("Simulation equations\n")
for index, eq in enumerate(flatten(simulation_eq.equations)):
    latex.write_expression(eq)

latex.write_string("Constituent relations\n")
for index, eq in enumerate(flatten(constituent.equations)):
    latex.write_expression(eq)

latex.close()

# Create a simulation block
block = SimulationBlock(ndim, block_number=0)

# Local dictionary for parsing the expressions
local_dict = {"block": block, "GridVariable": GridVariable, "DataObject": DataObject}

# initial conditions
dx, dy = block.deltas
i, j = block.grid_indexes
gama, Minf = symbols('gama Minf', **{'cls': ConstantObject})
q_vector=flatten(simulation_eq.time_advance_arrays)

x, y = symbols('x0:%d' % ndim, **{'cls': DataObject})
grid_equations= []
nx = ConstantObject('block0np0')
ny = ConstantObject('block0np1')
stretch = ConstantObject('stretch')
Lx=nx*dx
Ly=ConstantObject('Ly')
stretch_eqn=0.5*Ly*sinh(stretch*(j-(ny-1)/2)/((ny-1)/2))/sinh(stretch)
grid_equations += [Eq(x, i*dx), Eq(y,stretch_eqn)]

# initial conditions: in contrast to YSD this is currently missing Sutherland variable viscosity and the u component of the forcing 
# also note that Re is here based on freestream velocity (-1 to 1 here,compared to -0.5 to 0.5 in YSD). Also no subharmonic here yet.
# To compare with YSD output at t/2 ie t=40 there is 20 here, forcing amplitudes x2 here to be equivalent

initial_equations = []
rho, uref, vpert, p, T ,fref = symbols('rho, u, v, p, T, f', **{'cls': GridVariable})

initial_equations += [Eq(uref,tanh(2.0*y))]
initial_equations += [Eq(fref,0.5*(1.0+tanh(2.0*y)))]
initial_equations += [Eq(vpert,0.1*cos(2.0*pi*x/Lx)*exp(-y**2.0/10.0))]
initial_equations += [Eq(T,1.0+Minf**2*(gama-1.0)/2.0*(1.0-uref**2))]
initial_equations += [Eq(p,1.0/(gama*Minf**2.0))]
initial_equations += [Eq(rho,gama*Minf**2*p/T)]

initial_equations += [Eq(q_vector[0],rho)]
initial_equations += [Eq(q_vector[1],rho*uref)]
initial_equations += [Eq(q_vector[2],rho*vpert)]
initial_equations += [Eq(q_vector[3],p/(gama-1.0)+0.5*rho*(uref**2+vpert**2))]
initial_equations += [Eq(q_vector[4],rho*fref)]

# initial conditions
initial = GridBasedInitialisation()
initial.add_equations(grid_equations + initial_equations)

# metrics
metriceq =  MetricsEquation()
#metriceq.generate_transformations(ndim, coordinate_symbol, [(True, True), (True, True)], 2)
metriceq.generate_transformations(ndim, coordinate_symbol, [(False, False), (True, False)], 2)
simulation_eq.apply_metrics(metriceq)

# Create a schemes dictionary to be used for discretisation
schemes = {}
fns = 'u0 u1 T'
# cent = StoreSome(4, fns)
cent = Central(4)
schemes[cent.name] = cent
rk = RungeKuttaLS(3)
schemes[rk.name] = rk

# Create boundaries, one for each side per dimension, so in total 6 BC's for 3D'
boundaries = []
direction=0
boundaries += [PeriodicBC(direction, side=0)]
boundaries += [PeriodicBC(direction, side=1)]
direction=1
boundaries += [SymmetryBC(direction, 0)]
boundaries += [SymmetryBC(direction, 1)]
block.set_block_boundaries(boundaries)


# set the IO class to write out arrays
kwargs = {'iotype': "Write"}
h5 = iohdf5(save_every=200, **kwargs)
h5.add_arrays(simulation_eq.time_advance_arrays + [x, y])
block.setio(copy.deepcopy(h5))

# set the equations to be solved on the block
block.set_equations([copy.deepcopy(constituent), copy.deepcopy(simulation_eq), initial, metriceq])
# set the discretisation schemes
block.set_discretisation_schemes(schemes)

# Discretise the equations on the block
block.discretise()

# create an algorithm from the discretised computations
alg = TraditionalAlgorithmRK(block)

# set the simulation data type, for more information on the datatypes see opensbli.core.datatypes
SimulationDataType.set_datatype(Double)

# Write the code for the algorithm
OPSC(alg)
# Add the simulation constants to the OPS C code
substitute_simulation_parameters(simulation_parameters.keys(), simulation_parameters.values())
print_iteration_ops(NaN_check='rho')
