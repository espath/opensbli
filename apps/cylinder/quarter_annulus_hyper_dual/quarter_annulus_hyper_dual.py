#!/usr/bin/env python3
"""Quarter-annulus cylinder case with hyper dual-velocity equations.

Domain: x<=0, y>=0 quarter annulus around a circular wall.

Boundary map for generated quarter mesh (theta in [pi/2, pi]):
- dir0, side0: theta=pi/2 radial line (x=0, y>0) [pressure outlet]
- dir0, side1: theta=pi   radial line (x<0, y=0)  [symmetry]
- dir1, side0: inner arc (cylinder wall)
- dir1, side1: outer arc (farfield)
"""

import os
from sympy import Eq, S
from opensbli import *
from opensbli.physical_models.hyper_dual_velocity_physics import HyperDualVelocitySplit
from opensbli.core.boundary_conditions.inviscid_wall import InviscidWallBC

simulation_parameters = {
    'gama': os.getenv('GAMMA_VALUE', '1.4'),
    'Minf': os.getenv('MINF_VALUE', '2.0'),
    'Pr': os.getenv('PR_VALUE', '0.71'),
    # Lower than airfoil default due much coarser mesh.
    'Re': os.getenv('RE_VALUE', '1.0e4'),
    'Ml': os.getenv('ML_VALUE', '0.0'),
    'dt': os.getenv('DT_VALUE', '1.0e-5'),
    'niter': os.getenv('NITER_VALUE', '50000'),
    'Twall': os.getenv('TWALL_VALUE', '1.0'),
    'SuthT': os.getenv('SUTHT_VALUE', '110.4'),
    'RefT': os.getenv('REFT_VALUE', '273.15'),
    # Mesh constants (quarter mesh generated in apps/cylinder/grid_generation)
    'block0np0': os.getenv('NP0_VALUE', '241'),
    'block0np1': os.getenv('NP1_VALUE', '241'),
    # IMPORTANT: for this mesh layout, dir0 is angular and dir1 is radial.
    'Delta0block0': '((M_PI/2.0)/(block0np0-1.0))',
    'Delta1block0': '((20.0-0.5)/(block0np1-1.0))',
    'inv_rfact0_block0': '1.0/Delta0block0',
    'inv_rfact1_block0': '1.0/Delta1block0',
}

ndim = 2
conservative = True
coordinate_symbol = 'x'
constants = ['Re', 'Pr', 'gama', 'Minf', 'Ml', 'RefT', 'SuthT']

block = SimulationBlock(ndim, block_number=0)
SimulationDataType.set_datatype(Double)

metriceq = MetricsEquation()
metriceq.generate_transformations(ndim, coordinate_symbol, [(True, True), (True, True)], 2)

optional_subs = metriceq.metric_subs
EE = EinsteinEquation()
EE.optional_subs_dict = optional_subs

# Keep lagrangian velocity symbol consistent with standard velocity in this setup.
for eq in EE.expand('Eq(vl_i, u_i)', ndim, coordinate_symbol, [], constants):
    EE.optional_subs_dict[eq.lhs] = eq.rhs

HDV = HyperDualVelocitySplit(
    ndim,
    constants,
    coordinate_symbol=coordinate_symbol,
    conservative=conservative,
    viscosity='dynamic',
    include_skew_work=False,
    use_reynolds=True,
    split_type='KGP',
    energy_formulation='enthalpy',
)

simulation_eq = SimulationEquations()
simulation_eq.add_equations(HDV.mass)
simulation_eq.add_equations(HDV.momentum)
simulation_eq.add_equations(HDV.energy)

constituent = ConstituentRelations()
if conservative:
    for eq in EE.expand('Eq(u_i, rhou_i/rho)', ndim, coordinate_symbol, [], constants):
        constituent.add_equations(eq)
    pressure = 'Eq(p, (gama-1)*(rhoE - (1/2)*rho*(KD(_i,_j)*u_i*u_j)))'
    enthalpy = 'Eq(H, (rhoE + p) / rho)'
else:
    pressure = 'Eq(p, (gama-1)*(Et - (1/2)*(KD(_i,_j)*u_i*u_j)))'
    enthalpy = 'Eq(H, Et + p/rho)'

temperature = 'Eq(T, p*gama*Minf*Minf/(rho))'
viscosity = 'Eq(mu, (T**(1.5)*(1.0+SuthT/RefT)/(T+SuthT/RefT)))'

for expr in [pressure, temperature, viscosity, enthalpy]:
    constituent.add_equations(EE.expand(expr, ndim, coordinate_symbol, [], constants))

simulation_eq.apply_metrics(metriceq)

# Initial freestream state.
q_vector = flatten(simulation_eq.time_advance_arrays)
d, u0, u1, p = symbols('d u0 u1 p', **{'cls': GridVariable})
gama, Minf = symbols('gama Minf', **{'cls': ConstantObject})

init_eq = [
    Eq(d, 1.0),
    Eq(u0, 1.0),
    Eq(u1, 0.0),
    Eq(p, 1.0/(gama*Minf**2.0)),
]
if conservative:
    init_eq += [
        Eq(q_vector[0], d),
        Eq(q_vector[1], d*u0),
        Eq(q_vector[2], d*u1),
        Eq(q_vector[3], p/(gama-1.0) + 0.5*d*(u0**2 + u1**2)),
    ]
else:
    init_eq += [
        Eq(q_vector[0], d),
        Eq(q_vector[1], u0),
        Eq(q_vector[2], u1),
        Eq(q_vector[3], p/(d*(gama-1.0)) + 0.5*(u0**2 + u1**2)),
    ]

initial = GridBasedInitialisation()
initial.add_equations(init_eq)

# Schemes
schemes = {}
rk = RungeKuttaLS(4)
cent = StoreSome(4, 'u0 u1 T', merged=True)
schemes[rk.name] = rk
schemes[cent.name] = cent
block.set_discretisation_schemes(schemes)

# Boundary conditions.
boundaries = []

# dir0/side0 (theta=pi/2): pressure outlet
back_pressure = float(os.getenv('BACK_PRESSURE_VALUE', '1.0')) / (float(os.getenv('GAMMA_VALUE', '1.4')) * float(os.getenv('MINF_VALUE', '2.0'))**2.0)
boundaries += [PressureOutletBC(direction=0, side=0, back_pressure=back_pressure)]
# dir0/side1 (theta=pi): symmetry/slip plane
# Use InviscidWallBC (not SymmetryBC) so the boundary plane itself is set
# to zero normal velocity, not only reflected in halo points.
boundaries += [InviscidWallBC(direction=0, side=1)]

# dir1/side0 inner arc: isothermal no-slip wall
gama, Minf, Twall = symbols('gama Minf Twall', **{'cls': ConstantObject})
if conservative:
    wall_energy = [Eq(q_vector[-1], q_vector[0]*Twall / (gama * Minf**2.0 * (gama - S.One)))]
else:
    wall_energy = [Eq(q_vector[-1], Twall / (gama * Minf**2.0 * (gama - S.One)))]
boundaries += [IsothermalWallBC(direction=1, side=0, equations=wall_energy)]

# dir1/side1 outer arc: farfield
boundaries += [DirichletBC(direction=1, side=1, equations=init_eq)]

block.set_block_boundaries(boundaries)

# Optional shock capturing filters
eq_classes = [constituent, simulation_eq, initial, metriceq]
if os.getenv('ENABLE_WENO', '0') == '1':
    eq_classes += WENOFilter(block, order=5, formulation='Z', flux_type='LLF', airfoil=False, metrics=metriceq).equation_classes
if os.getenv('ENABLE_TVD', '0') == '1':
    eq_classes += TVDFilter(block, airfoil=False).equation_classes

block.set_equations(eq_classes)

# IO: write conservative fields, read grid from data.h5
write_h5 = iohdf5(save_every=int(os.getenv('SAVE_EVERY_VALUE', '1000')), **{'iotype': 'Write'})
write_h5.add_arrays(simulation_eq.time_advance_arrays)
read_h5 = iohdf5(**{'iotype': 'Read'})
x, y = symbols('x0 x1', **{'cls': DataObject})
read_h5.add_arrays([x, y])
block.setio([write_h5, read_h5])

block.discretise()

SM = SimulationMonitor(
    [block.location_dataset('rho'), block.location_dataset('u0'), block.location_dataset('u1')],
    [(10, 10), (120, 120), (230, 230)],
    block,
    print_frequency=100,
    fp_precision=12,
    output_file='quarter_annulus_monitor.log',
)

alg = TraditionalAlgorithmRK(block, simulation_monitor=SM)
OPSC(alg, OPS_diagnostics=1)
if os.getenv('ENABLE_NAN_CHECK', '1') == '1':
    print_iteration_ops(NaN_check='rho', every=100, nblocks=1)
else:
    print_iteration_ops(every=100, nblocks=1)
substitute_simulation_parameters(simulation_parameters.keys(), simulation_parameters.values())
