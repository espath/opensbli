#!/usr/bin/env python
# Import all the functions from opensbli
from opensbli import *
import copy
import os
from pathlib import Path
import h5py
from opensbli.utilities.helperfunctions import substitute_simulation_parameters
from sympy import pi, sin, cos, Abs, sqrt

grid_path = Path(__file__).with_name("data.h5")
with h5py.File(grid_path, "r") as _f:
    _shape = _f["opensbliblock00"]["x0_B0"].shape
if len(_shape) != 2:
    raise RuntimeError(f"Expected 2D grid in {grid_path}, found shape {_shape}")

# OpenSBLI direction-0 aligns with second dataset axis in this case.
mesh_np0_default = str(_shape[1])
mesh_np1_default = str(_shape[0])

simulation_parameters = {
'Re'        :   os.getenv('RE_VALUE', '300.0'),
'gama'      :   os.getenv('GAMMA_VALUE', '1.4'),
'Minf'      :   os.getenv('MINF_VALUE', '1.5'),
'Pr'        :   os.getenv('PR_VALUE', '0.71'),
'dt'        :   os.getenv('DT_VALUE', '0.0001'),
'niter'     :   os.getenv('NITER_VALUE', '5000000'),
'block0np0'     :   os.getenv('NP0_VALUE', mesh_np0_default),
'block0np1'     :   os.getenv('NP1_VALUE', mesh_np1_default),
'Delta0block0'      :   'M_PI/(block0np0-1)',
'Delta1block0'      :   '242.2/(block0np1-1)',
'Twall'     :   os.getenv('TWALL_VALUE', '1.0'),
'SuthT'     :   os.getenv('SUTHT_VALUE', '110.4'),
'RefT'      :   os.getenv('REFT_VALUE', '273.15'),
'inv_rfact0_block0'     :   '1.0/Delta0block0',
'inv_rfact1_block0'     :   '1.0/Delta1block0',
'shock_factor'      :   os.getenv('SHOCK_FACTOR_VALUE', '1.0'),
}

# Problem dimension
ndim = 2
# # Constants that are used
constants = ["Re", "Pr", "gama", "Minf", "RefT", "SuthT"]
# # symbol for the coordinate system in the equations
coordinate_symbol = "x"
metriceq = MetricsEquation()
metriceq.generate_transformations(ndim, coordinate_symbol, [(True, True), (True, True)], 2)
#Create an optional substitutions dictionary, this will be used to modify the equations when parsed
optional_subs_dict = metriceq.metric_subs
# symbol for the coordinate system in the equations
conservative = True
NS = NS_Split('KGP', ndim, constants, coordinate_symbol=coordinate_symbol, conservative=conservative, viscosity='dynamic', energy_formulation='enthalpy')
mass, momentum, energy = NS.mass, NS.momentum, NS.energy
# Expand the simulation equations, for this create a simulation equations class
simulation_eq = SimulationEquations()
simulation_eq.add_equations(mass)
simulation_eq.add_equations(momentum)
simulation_eq.add_equations(energy)

# Constituent relations used in the system
if conservative:
    pressure = "Eq(p, (gama-1)*(rhoE - (1/2)*rho*(KD(_i,_j)*u_i*u_j)))"
    velocity = "Eq(u_i, rhou_i/rho)"
    enthalpy = "Eq(H, (rhoE + p) / rho)"
else:
    pressure = "Eq(p, rho*(gama-1)*(Et - (1/2)*(KD(_i,_j)*u_i*u_j)))"
    enthalpy = "Eq(H, Et + p / rho)"

temperature = "Eq(T, p*gama*Minf*Minf/(rho))"
viscosity = "Eq(mu, (T**(1.5)*(1.0+SuthT/RefT)/(T+SuthT/RefT)))"

# Expand the constituent relations and them to the constituent relations class
constituent = ConstituentRelations()  # Instantiate constituent relations object
einstein_eq = EinsteinEquation()
einstein_eq.optional_subs_dict = optional_subs_dict
metric_vel = "Eq(U_i, D_i_j*u_j)"
eqns = einstein_eq.expand(metric_vel, ndim, coordinate_symbol, [], constants)
for eq in eqns:
    einstein_eq.optional_subs_dict[eq.lhs] = eq.rhs

# Expand momentum add the expanded equations to the constituent relations
if conservative:
    eqns = einstein_eq.expand(velocity, ndim, coordinate_symbol, [], constants)
    constituent.add_equations(eqns)
# Expand pressure add the expanded equations to the constituent relations
eqns = einstein_eq.expand(pressure, ndim, coordinate_symbol, [], constants)
constituent.add_equations(eqns)
# Expand enthalpy add the expanded equations to the constituent relations
eqns = einstein_eq.expand(enthalpy, ndim, coordinate_symbol, [], constants)
constituent.add_equations(eqns)
# Expand temperature add the expanded equations to the constituent relations
eqns = einstein_eq.expand(temperature, ndim, coordinate_symbol, [], constants)
constituent.add_equations(eqns)
# Expand viscosity add the expanded equations to the constituent relations
eqns = einstein_eq.expand(viscosity, ndim, coordinate_symbol, [], constants)
constituent.add_equations(eqns)
# Create a simulation block
block = SimulationBlock(ndim, block_number=0, conservative=conservative)
# Transform to curvilinear
simulation_eq.apply_metrics(metriceq)

# Local dictionary for parsing the expressions
local_dict = {"block": block, "GridVariable": GridVariable, "DataObject": DataObject}

# Initial conditions as strings
u0 = "Eq(GridVariable(u0),1.0)"
u1 = "Eq(GridVariable(u1), 0.0,)"
p = "Eq(GridVariable(p), 1/(gama*Minf*Minf))"
r = "Eq(GridVariable(r), gama*Minf*Minf*p)"

if conservative:
    rho = "Eq(DataObject(rho), r)"
    rhou0 = "Eq(DataObject(rhou0), r*u0)"
    rhou1 = "Eq(DataObject(rhou1), r*u1)"
    rhoE = "Eq(DataObject(rhoE), p/(gama-1) + 0.5* r *(u0**2+ u1**2))"
else:
    rho = "Eq(DataObject(rho), r)"
    rhou0 = "Eq(DataObject(u0), u0)"
    rhou1 = "Eq(DataObject(u1), u1)"
    rhoE = "Eq(DataObject(Et), p/(rho*(gama-1)) + 0.5*(u0**2+ u1**2))"
eqns = [u0, u1, p, r, rho, rhou0, rhou1, rhoE]

# parse the initial conditions
initial_equations = [parse_expr(eq, local_dict=local_dict) for eq in eqns]
initial = GridBasedInitialisation()
initial.add_equations(initial_equations)

# Create a schemes dictionary to be used for discretisation
schemes = {}
# Central scheme for spatial discretisation and add to the schemes dictionary
# Low storage optimisation for the central scheme
fns = 'u0 u1 T'
cent = StoreSome(4, fns)
schemes[cent.name] = cent
# RungeKutta scheme for temporal discretisation and add to the schemes dictionary
rk = RungeKuttaLS(3, formulation='SSP', stages=3)
schemes[rk.name] = rk

# Create boundaries, one for each side per dimension
q_vector = flatten(simulation_eq.time_advance_arrays)
boundaries = []
direction = 0
# Apply a periodic boundary over the shared mesh line
boundaries += [PeriodicBC(direction, 0, halos=[-2,2], corners=True)]
boundaries += [PeriodicBC(direction, 1, halos=[-2,2], corners=True)]
# Isothermal wall in x1 direction
gama, Minf, Twall = symbols('gama Minf Twall', **{'cls': ConstantObject})
# Energy on the wall is set
if conservative:
    wall_energy = [Eq(q_vector[-1], Twall*q_vector[0] / (gama * Minf**2.0 * (gama - S.One)))]
else:
    wall_energy = [Eq(q_vector[-1], Twall / (gama * Minf**2.0 * (gama - S.One)))]

direction = 1
lower_wall_eq = wall_energy[:]
boundaries += [IsothermalWallBC(direction, 0, lower_wall_eq)]
# Far field boundary
direction, side = 1,1
boundaries += [DirichletBC(direction, side, initial_equations)]
# set the boundaries for the block
block.set_block_boundaries(boundaries)

# Set the IO class to write out arrays
kwargs = {'iotype': "Write", "write_constants" : True}
h5 = iohdf5(save_every=5000, **kwargs)
h5.add_arrays(simulation_eq.time_advance_arrays)
h5.add_arrays([DataObject('x0'), DataObject('x1')])#, DataObject('kappa')])
kwargs = {'iotype': "Read"}
h5_read = iohdf5(**kwargs)
h5_read.add_arrays([DataObject('x0'), DataObject('x1')])
block.setio([h5, h5_read])

# Set equations 
block.set_equations([constituent, simulation_eq, initial, metriceq])

# Add SFD filtering
if os.getenv("ENABLE_SFD_FILTER", "1") == "1":
    SFD = SFD(block, chifilt=0.1, omegafilt=1.0/0.75)
    block.set_equations(SFD.equation_classes)

j = block.grid_indexes[1]
grid_condition = j >= (int(simulation_parameters['block0np1']) - 4)
if os.getenv("ENABLE_BINOMIAL_FILTER", "1") == "1":
    BF = BinomialFilter(block, order=6, directions=[0,1], grid_condition=grid_condition, sigma=0.2)
    block.set_equations(BF.equation_classes)

if os.getenv("ENABLE_DRP_FILTER", "0") == "1":
    DRP = ExplicitFilter(block, [0,1], width=9, filter_type='DRP', optimized=False, sigma=0.1, multi_block=None)
    block.set_equations(DRP.equation_classes)

# WENO filter for shock-capturing
WF = None
if os.getenv("ENABLE_WENO_FILTER", "0") == "1":
    WF = WENOFilter(block, order=7, metrics=metriceq, flux_type='LLF', airfoil=False)
    block.set_equations(WF.equation_classes)

# set the discretisation schemes
block.set_discretisation_schemes(schemes)

# Monitor residuals within the domain
RM = ResidualMonitor(block, frequency=100)
block.set_equations(RM.equation_classes)

# Discretise the equations on the block
block.discretise()

if WF is not None:
    WF.update_periodic_boundary(block, halos=[-5,5])

# Full 5 swaps for the filter over the interface
# Add some full [-5,5] halo swaps over the periodic directions only when the filter is called
def create_exchange_calls_codes(block, dsets):
    kernels = []
    arrays = [block.location_dataset(a) for a in flatten(dsets)]
    for direction in [0]:
        for side in [0,1]:
            BC = PeriodicBC(direction, side, halos=[-5,5], corners=True)
            kernels += [BC.apply(arrays, block)]
    return kernels

# Make some full swaps for interfaces before filtering
if conservative:
    filter_swaps = create_exchange_calls_codes(block, ['rho', 'rhou0', 'rhou1', 'rhoE'])
else:
    filter_swaps = create_exchange_calls_codes(block, ['rho', 'u0', 'u1', 'Et'])

if WF is not None:
    for no, eq in enumerate(block.list_of_equation_classes):
        if isinstance(eq, UserDefinedEquations):
            if eq.full_swap:
                eq.Kernels += filter_swaps

# Monitor the residuals
# Simulation monitor
# arrays = ['L2_R0', 'L2_R1', 'L2_R2', 'L2_R3', 'u1_B0']
# probe_locations = ['residual', 'residual', 'residual', 'residual', (0, 100)]
# SM = SimulationMonitor(arrays, probe_locations, block, print_frequency=100, output_file='residuals.log')

# Create algorithm
alg = TraditionalAlgorithmRK(block)
# set the simulation data type, for more information on the datatypes see opensbli.core.datatypes
SimulationDataType.set_datatype(Double)
# Write the code for the algorithm
OPSC(alg)
# Add the simulation constants to the OPS C code
substitute_simulation_parameters(simulation_parameters.keys(), simulation_parameters.values())
if os.getenv("ENABLE_NAN_CHECK", "1") == "1":
    print_iteration_ops(NaN_check='rho')
else:
    print_iteration_ops()
