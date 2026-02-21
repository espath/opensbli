#!/usr/bin/env python
# Import all the functions from opensbli
from opensbli import *
import copy
from opensbli.utilities.helperfunctions import substitute_simulation_parameters

# Direct application of shock-capturing scheme, otherwise central scheme with filter-step example
weno = False
teno = False
ndim = 1
# Define all the constants in the equations
constants = ["gama", "Minf"]
# Define coordinate direction symbol (x) this will be x_i, x_j, x_k
coordinate_symbol = "x"
# Substitutions
conservative = True
substitutions = []
eq = EinsteinEquation()

if weno or teno:
    if weno:
        sc1 = "**{'scheme':'Weno'}"
    else:
        sc1 = "**{'scheme':'Teno'}"

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
    NS = NS_Split(
        "Feiereisen",
        ndim,
        constants,
        coordinate_symbol=coordinate_symbol,
        conservative=conservative,
        viscosity="inviscid",
        energy_formulation="enthalpy",
        debug=False,
    )
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
x0 = "Eq(DataObject(x0), block.deltas[0]*block.grid_indexes[0])"
d = "Eq(GridVariable(d), Piecewise((1.0, DataObject(x0) < 0.5), (0.125, True)))"
p = "Eq(GridVariable(p), Piecewise((1.0, DataObject(x0) < 0.5), (0.10, True)))"
u0 = "Eq(GridVariable(u0), Piecewise((0.0, DataObject(x0) < 0.5),(0, True)))"
rho = "Eq(DataObject(rho), d)"
rhou0 = "Eq(DataObject(rhou0), d*u0)"
rhoE = "Eq(DataObject(rhoE), p/(gama-1.0) + 0.5* d *(u0**2.0))"
eqns = [x0, u0, p, d, rho, rhou0, rhoE]
local_dict = {"block": block, "GridVariable": GridVariable, "DataObject": DataObject}
initial_equations = [parse_expr(eq, local_dict=local_dict) for eq in eqns]
initial = GridBasedInitialisation()
initial.add_equations(initial_equations)

# Left boundary condition
d, u, p = symbols("d u0 p", **{"cls": GridVariable})
left_eqns = [OpenSBLIEq(d, 1.0), OpenSBLIEq(u, 0.0), OpenSBLIEq(p, 1.0)] + [
    parse_expr(eq, local_dict=local_dict) for eq in [rho, rhou0, rhoE]
]
# Right boundary condition
right_eqns = [OpenSBLIEq(d, 0.125), OpenSBLIEq(u, 0.0), OpenSBLIEq(p, 0.1)] + [
    parse_expr(eq, local_dict=local_dict) for eq in [rho, rhou0, rhoE]
]

boundaries = []
# Create boundaries, one for each side per dimension
for direction in range(ndim):
    boundaries += [DirichletBC(direction, 0, left_eqns)]
    boundaries += [DirichletBC(direction, 1, right_eqns)]

schemes = {}
# Spatial scheme
if weno:
    Avg = SimpleAverage([0, 1])
    # LF = LFWeno(order=7, formulation='Z', averaging=Avg, flux_type='LLF')
    LF = HLLCWeno(order=5, formulation="Z", averaging=Avg, flux_type="HLLC-LM")
    # Add to schemes
    schemes[LF.name] = LF
elif teno:
    Avg = SimpleAverage([0, 1])
    LF = LFTeno(order=5, averaging=Avg, flux_type="GLF")
    # LF = HLLCTeno(order=6, averaging=Avg, flux_type='HLLC-LM')
    # Add to schemes
    schemes[LF.name] = LF
else:
    fns = "u0"
    cent = StoreSome(4, fns)
    schemes[cent.name] = cent

# Time-stepping
rk = RungeKuttaLS(3, formulation="SSP")
schemes[rk.name] = rk

block.set_block_boundaries(boundaries)
kwargs = {"iotype": "Write"}
h5 = iohdf5(**kwargs)
h5.add_arrays(simulation_eq.time_advance_arrays)
h5.add_arrays([DataObject("x0")])

# if not weno:
# h5.add_arrays([DataObject('kappa'), DataObject('q0'), DataObject('q1'), DataObject('q2')])
block.setio(copy.deepcopy(h5))

TVD = False
if not weno and not teno:
    if TVD:
        TVD_filter = TVDFilter(block, airfoil=False)
        block.set_equations(TVD_filter.equation_classes)
    else:
        WF = WENOFilter(
            block,
            order=5,
            formulation="Z",
            flux_type="LLF",
            airfoil=False,
            optimize=False,
        )
        block.set_equations(WF.equation_classes)


block.set_equations([copy.deepcopy(constituent), copy.deepcopy(simulation_eq), initial])
block.set_discretisation_schemes(schemes)
block.discretise()

alg = TraditionalAlgorithmRK(block)
SimulationDataType.set_datatype(Double)
OPSC(alg)
constants = [
    "gama",
    "dt",
    "niter",
    "block0np0",
    "Delta0block0",
    "eps",
    "inv_rfact0_block0",
    "lambda0_TVD",
    "TENO_CT",
]
values = [
    "1.4",
    "0.0002",
    # dt=2e-4, niter=1000 -> final time t=0.2
    "1000",
    "200",
    "1.0/(block0np0-1)",
    "1.0e-16",
    "1.0/Delta0block0",
    "dt/Delta0block0",
    "1.0e-5",
]
substitute_simulation_parameters(constants, values)
print_iteration_ops(NaN_check="rho")
