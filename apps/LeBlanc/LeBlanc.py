#!/usr/bin/env python
# Import all the functions from opensbli
from opensbli import *
import copy
from opensbli.utilities.helperfunctions import substitute_simulation_parameters

# Input parameters for the simulation
simulation_parameters = {
    "gama"                 : "5.0/3.0",
    "dt"                   : "0.0001",
    "niter"                : "ceil(6.0/0.0001)",
    "block0np0"            : "900",
    "Delta0block0"         : "9.0/(block0np0-1)",
    "eps"                  : "1e-15",
    "TENO_CT"              : "1e-5",
    "inv_rfact0_block0"    : "'1.0/Delta0block0"
}

# Direct application of shock-capturing scheme, otherwise central scheme with filter-step example
teno = False
weno = True
ndim = 1
# Define all the constants in the equations
constants = ["gama"]
conservative = True
# Define coordinate direction symbol (x) this will be x_i, x_j, x_k
coordinate_symbol = "x"
# Substitutions
substitutions = []
eq = EinsteinEquation()

if teno or weno:
    if teno:
        sc1 = "**{\'scheme\':\'Teno\'}"
    else:
        sc1 = "**{\'scheme\':\'Weno\'}"
    # Define the compresible Navier-Stokes equations in Einstein notation.
    a = "Conservative(rhou_j,x_j,%s)" % sc1
    mass = "Eq(Der(rho,t), - %s)" % (a)
    a = "Conservative(rhou_i*u_j + KD(_i,_j)*p,x_j , %s)" % sc1
    momentum = "Eq(Der(rhou_i,t) , -%s  )" % (a)
    a = "Conservative((p+rhoE)*u_j,x_j, %s)" % sc1
    energy = "Eq(Der(rhoE,t), - %s  )" % (a)

    mass = eq.expand(mass, ndim, coordinate_symbol, substitutions, constants)
    momentum = eq.expand(momentum, ndim, coordinate_symbol, substitutions, constants)
    energy = eq.expand(energy, ndim, coordinate_symbol, substitutions, constants)
else:
    NS = NS_Split('KGP', ndim, constants, coordinate_symbol=coordinate_symbol, conservative=conservative, viscosity='inviscid', energy_formulation='enthalpy', debug=False)
    mass, momentum, energy = NS.mass, NS.momentum, NS.energy

# Expand the simulation equations, for this create a simulation equations class
simulation_eq = SimulationEquations()
simulation_eq.add_equations(mass)
simulation_eq.add_equations(momentum)
simulation_eq.add_equations(energy)

# Constituent relations
if conservative:
    pressure = "Eq(p, (gama-1)*(rhoE - (1/2)*rho*(KD(_i,_j)*u_i*u_j)))"
    velocity = "Eq(u_i, rhou_i/rho)"
    enthalpy = "Eq(H, (rhoE + p) / rho)"
else:
    pressure = "Eq(p, rho*(gama-1)*(Et - (1/2)*(KD(_i,_j)*u_i*u_j)))"
    enthalpy = "Eq(H, Et + p / rho)"

speed_of_sound = "Eq(a, (gama*p/rho)**0.5)"


constituent = ConstituentRelations()
eqns = eq.expand(velocity, ndim, coordinate_symbol, substitutions, constants)
constituent.add_equations(eqns)
eqns = eq.expand(pressure, ndim, coordinate_symbol, substitutions, constants)
constituent.add_equations(eqns)
eqns = eq.expand(speed_of_sound, ndim, coordinate_symbol, substitutions, constants)
constituent.add_equations(eqns)
eqns = eq.expand(enthalpy, ndim, coordinate_symbol, substitutions, constants)
constituent.add_equations(eqns)

block = SimulationBlock(ndim, block_number=0)

# Initial conditions
initial = GridBasedInitialisation()
# x = "GridVariable(x0)"
x0 = "Eq(GridVariable(x0), block.deltas[0]*block.grid_indexes[0])"
x0_dset = "Eq(DataObject(x0), GridVariable(x0))"

d = "Eq(GridVariable(d), Piecewise((1, x0 < 3),(1e-3, True)))"
u0 = "Eq(GridVariable(u0), Piecewise((0, x0 < 3),(0, True)))"
p = "Eq(GridVariable(p), Piecewise(((2/3)*1e-1, x0 < 3), ((2/3)*1e-10, True)))"

rho = "Eq(DataObject(rho), d)"
rhou0 = "Eq(DataObject(rhou0), d*u0)"
rhoE = "Eq(DataObject(rhoE), p/(gama-1.0) + 0.5* d *(u0**2.0))"

# Temp arrays for testing
kappa = "Eq(DataObject(kappa), 0)"
q0 = "Eq(DataObject(q0), 0)"
q1 = "Eq(DataObject(q1), 0)"
q2 = "Eq(DataObject(q2), 0)"
eqns = [x0, x0_dset, u0, p, d, rho, rhou0, rhoE, kappa, q0, q1, q2]

local_dict = {"block": block, "GridVariable": GridVariable, "DataObject": DataObject}
initial_equations = [parse_expr(eq, local_dict=local_dict) for eq in eqns]
initial = GridBasedInitialisation()
initial.add_equations(initial_equations)


# LeBlanc boundary condition values left side
arrays = flatten(simulation_eq.time_advance_arrays)
subs_dict = {Symbol('x0'): 0}
boundary_eqns = [x0, u0, p, d, rho, rhou0, rhoE]
boundary_eqns = [parse_expr(eq, local_dict=local_dict) for eq in boundary_eqns]
left_eqns = [eq.subs(subs_dict) for eq in boundary_eqns]

subs_dict = {Symbol('x0'): 9.0}

right_eqns = [eq.subs(subs_dict) for eq in boundary_eqns]

boundaries = []
# Create boundaries, one for each side per dimension
for direction in range(ndim):
    boundaries += [DirichletBC(direction, 0, left_eqns)]
    boundaries += [DirichletBC(direction, 1, right_eqns)]

schemes = {}
# Spatial scheme
if teno or weno:
    Avg = RoeAverage([0, 1])
    # Avg = SimpleAverage([0, 1])
    if teno:
        # LF = LFTeno(order=6, averaging=Avg, flux_type='LLF', flux_split=False)
        LF = HLLCTeno(order=6, averaging=Avg, flux_type='HLLC-LM')
    else:
        LF = LFWeno(order=3, formulation='JS', averaging=Avg, flux_type='LLF', flux_split=True)
        # LF = HLLCWeno(order=5, formulation='JS', averaging=Avg, flux_type='HLLC-LM')
    # Add to schemes
    schemes[LF.name] = LF
else:
    fns = 'u0'
    cent = StoreSome(4, fns)
    # cent = Central(4)
    schemes[cent.name] = cent
# Time-stepping
rk = RungeKuttaLS(3, formulation='SSP')
schemes[rk.name] = rk

block.set_block_boundaries(boundaries)
kwargs = {'iotype': "Write"}
h5 = iohdf5(**kwargs)
h5.add_arrays(simulation_eq.time_advance_arrays)
h5.add_arrays([DataObject('x0')])

if not teno and not weno:
    # WENO filter for shock-capturing
    WF = WENOFilter(block, order=3, dissipation_sensor='Ducros', flux_type='LLF', airfoil=False)
    block.set_equations(WF.equation_classes)
    h5.add_arrays([ DataObject('kappa')])

block.setio(copy.deepcopy(h5))
block.set_equations([copy.deepcopy(constituent), copy.deepcopy(simulation_eq), initial])
block.set_discretisation_schemes(schemes)

block.discretise()

alg = TraditionalAlgorithmRK(block)
SimulationDataType.set_datatype(Double)
OPSC(alg)
# Add the simulation constants to the OPS C code
substitute_simulation_parameters(simulation_parameters.keys(), simulation_parameters.values())
print_iteration_ops(NaN_check='rho')
