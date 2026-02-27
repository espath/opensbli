#!/usr/bin/env python3
"""Full-annulus cylinder smoke test using baseline NS_Split equations."""

import os
from sympy import Eq, S
from opensbli import *

simulation_parameters = {
    'gama': os.getenv('GAMMA_VALUE', '1.4'),
    'Minf': os.getenv('MINF_VALUE', '2.0'),
    'Pr': os.getenv('PR_VALUE', '0.71'),
    'Re': os.getenv('RE_VALUE', '1.0e4'),
    'dt': os.getenv('DT_VALUE', '1.0e-5'),
    'niter': os.getenv('NITER_VALUE', '100'),
    'Twall': os.getenv('TWALL_VALUE', '1.0'),
    'SuthT': os.getenv('SUTHT_VALUE', '110.4'),
    'RefT': os.getenv('REFT_VALUE', '273.15'),
    'block0np0': os.getenv('NP0_VALUE', '241'),
    'block0np1': os.getenv('NP1_VALUE', '721'),
    'Delta0block0': '((20.0-0.5)/(block0np0-1.0))',
    'Delta1block0': '((2.0*M_PI)/(block0np1))',
    'inv_rfact0_block0': '1.0/Delta0block0',
    'inv_rfact1_block0': '1.0/Delta1block0',
}

ndim = 2
constants = ['Re', 'Pr', 'gama', 'Minf', 'RefT', 'SuthT']
coordinate_symbol = 'x'
conservative = True

block = SimulationBlock(ndim, block_number=0, conservative=conservative)
SimulationDataType.set_datatype(Double)

metriceq = MetricsEquation()
metriceq.generate_transformations(ndim, coordinate_symbol, [(True, True), (True, True)], 2)

EE = EinsteinEquation()
EE.optional_subs_dict = metriceq.metric_subs
for eq in EE.expand('Eq(U_i, D_i_j*u_j)', ndim, coordinate_symbol, [], constants):
    EE.optional_subs_dict[eq.lhs] = eq.rhs

NS = NS_Split(
    'KGP', ndim, constants,
    coordinate_symbol=coordinate_symbol,
    conservative=conservative,
    viscosity='dynamic',
    energy_formulation='enthalpy',
)

simulation_eq = SimulationEquations()
simulation_eq.add_equations(NS.mass)
simulation_eq.add_equations(NS.momentum)
simulation_eq.add_equations(NS.energy)

constituent = ConstituentRelations()
for eq in EE.expand('Eq(u_i, rhou_i/rho)', ndim, coordinate_symbol, [], constants):
    constituent.add_equations(eq)
for expr in [
    'Eq(p, (gama-1)*(rhoE - (1/2)*rho*(KD(_i,_j)*u_i*u_j)))',
    'Eq(T, p*gama*Minf*Minf/(rho))',
    'Eq(mu, (T**(1.5)*(1.0+SuthT/RefT)/(T+SuthT/RefT)))',
    'Eq(H, (rhoE + p) / rho)',
]:
    constituent.add_equations(EE.expand(expr, ndim, coordinate_symbol, [], constants))

simulation_eq.apply_metrics(metriceq)

q_vector = flatten(simulation_eq.time_advance_arrays)
d, u0, u1, p = symbols('d u0 u1 p', **{'cls': GridVariable})
gama, Minf = symbols('gama Minf', **{'cls': ConstantObject})
init_eq = [
    Eq(d, 1.0),
    Eq(u0, 1.0),
    Eq(u1, 0.0),
    Eq(p, 1.0/(gama*Minf**2.0)),
    Eq(q_vector[0], d),
    Eq(q_vector[1], d*u0),
    Eq(q_vector[2], d*u1),
    Eq(q_vector[3], p/(gama-1.0) + 0.5*d*(u0**2 + u1**2)),
]

initial = GridBasedInitialisation()
initial.add_equations(init_eq)

rk = RungeKuttaLS(4)
cent = StoreSome(4, 'u0 u1 T', merged=True)
block.set_discretisation_schemes({rk.name: rk, cent.name: cent})

boundaries = []
wall_gamma, wall_minf, Twall = symbols('gama Minf Twall', **{'cls': ConstantObject})
wall_energy = [Eq(q_vector[-1], q_vector[0] * Twall / (wall_gamma * wall_minf**2.0 * (wall_gamma - S.One)))]
boundaries += [IsothermalWallBC(direction=0, side=0, equations=wall_energy)]
boundaries += [DirichletBC(direction=0, side=1, equations=init_eq)]
boundaries += [PeriodicBC(direction=1, side=0, halos=[-2, 2], corners=True)]
boundaries += [PeriodicBC(direction=1, side=1, halos=[-2, 2], corners=True)]
block.set_block_boundaries(boundaries)

block.set_equations([constituent, simulation_eq, initial, metriceq])

if os.getenv('ENABLE_WENO', '1') == '1':
    block.set_equations(WENOFilter(block, order=5, formulation='Z', flux_type='LLF', airfoil=False, metrics=metriceq).equation_classes)
if os.getenv('ENABLE_SFD_FILTER', '1') == '1':
    block.set_equations(SFD(block, chifilt=0.1, omegafilt=1.0 / 0.75).equation_classes)
if os.getenv('ENABLE_BINOMIAL_FILTER', '1') == '1':
    j = block.grid_indexes[0]
    grid_condition = j >= (int(simulation_parameters['block0np0']) - 4)
    block.set_equations(BinomialFilter(block, order=6, directions=[0, 1], grid_condition=grid_condition, sigma=0.2).equation_classes)

write_h5 = iohdf5(save_every=int(os.getenv('SAVE_EVERY_VALUE', '100')), **{'iotype': 'Write'})
write_h5.add_arrays(simulation_eq.time_advance_arrays)
read_h5 = iohdf5(**{'iotype': 'Read'})
read_h5.add_arrays([DataObject('x0'), DataObject('x1')])
block.setio([write_h5, read_h5])

block.discretise()
alg = TraditionalAlgorithmRK(block)
OPSC(alg, OPS_diagnostics=1)
if os.getenv('ENABLE_NAN_CHECK', '1') == '1':
    print_iteration_ops(NaN_check='rho', every=100, nblocks=1)
else:
    print_iteration_ops(every=100, nblocks=1)
substitute_simulation_parameters(simulation_parameters.keys(), simulation_parameters.values())
