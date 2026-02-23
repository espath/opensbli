#!/usr/bin/env python3
"""Half-annulus cylinder case with hyper dual-velocity equations.

Domain: y>=0 half annulus around a circular wall.
Generated mesh uses theta in [0, pi].

For OpenSBLI HDF5 layout in this case:
- computational dir0 = radial (size[0] = nr)
- computational dir1 = angular (size[1] = ntheta)

Boundary map:
- dir0, side0: inner arc (cylinder wall)
- dir0, side1: outer arc (farfield)
- dir1, side0: theta=0 radial line (x>0, y=0) [outflow]
- dir1, side1: theta=pi radial line (x<0, y=0) [inflow]
"""

import os
from pathlib import Path
import h5py
from sympy import Eq, S
from opensbli import *
from opensbli.physical_models.hyper_dual_velocity_physics import HyperDualVelocitySplit

grid_path = Path(__file__).with_name('data.h5')
with h5py.File(grid_path, 'r') as _f:
    _dset = _f['opensbliblock00']['x0_B0']
    _size = _dset.attrs.get('size')
if _size is None or len(_size) != 2:
    raise RuntimeError(f'Could not read mesh size metadata from {grid_path}')
mesh_np0_default = str(int(_size[0]))
mesh_np1_default = str(int(_size[1]))

simulation_parameters = {
    'gama': os.getenv('GAMMA_VALUE', '1.4'),
    'Minf': os.getenv('MINF_VALUE', '2.0'),
    'Pr': os.getenv('PR_VALUE', '0.71'),
    'Re': os.getenv('RE_VALUE', '1.0e4'),
    'Ml': os.getenv('ML_VALUE', '0.0'),
    'dt': os.getenv('DT_VALUE', '5.0e-6'),
    'niter': os.getenv('NITER_VALUE', '50000'),
    'Twall': os.getenv('TWALL_VALUE', '1.0'),
    'SuthT': os.getenv('SUTHT_VALUE', '110.4'),
    'RefT': os.getenv('REFT_VALUE', '273.15'),
    # Mesh constants (half mesh generated in apps/cylinder/grid_generation)
    'block0np0': os.getenv('NP0_VALUE', mesh_np0_default),   # radial points
    'block0np1': os.getenv('NP1_VALUE', mesh_np1_default),   # angular points
    # IMPORTANT: computational direction 0 = radial, 1 = angular for this grid file.
    'Delta0block0': '((20.0-0.5)/(block0np0-1.0))',
    'Delta1block0': '(M_PI/(block0np1-1.0))',
    'inv_rfact0_block0': '1.0/Delta0block0',
    'inv_rfact1_block0': '1.0/Delta1block0',
}

ndim = 2
conservative = True
coordinate_symbol = 'x'
constants = ['Re', 'Pr', 'gama', 'Minf', 'Ml', 'RefT', 'SuthT']

block = SimulationBlock(ndim, block_number=0, conservative=conservative)
SimulationDataType.set_datatype(Double)

metriceq = MetricsEquation()
metriceq.generate_transformations(ndim, coordinate_symbol, [(True, True), (True, True)], 2)

EE = EinsteinEquation()
EE.optional_subs_dict = metriceq.metric_subs
for eq in EE.expand('Eq(vl_i, u_i)', ndim, coordinate_symbol, [], constants):
    EE.optional_subs_dict[eq.lhs] = eq.rhs
for eq in EE.expand('Eq(U_i, D_i_j*u_j)', ndim, coordinate_symbol, [], constants):
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
for eq in EE.expand('Eq(u_i, rhou_i/rho)', ndim, coordinate_symbol, [], constants):
    constituent.add_equations(eq)
for expr in [
    'Eq(p, (gama-1)*(rhoE - (1/2)*rho*(KD(_i,_j)*u_i*u_j)))',
    'Eq(T, p*gama*Minf*Minf/(rho))',
    'Eq(mu, (T**(1.5)*(1.0+SuthT/RefT)/(T+SuthT/RefT)))',
    'Eq(H, (rhoE + p)/rho)',
]:
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
    Eq(q_vector[0], d),
    Eq(q_vector[1], d*u0),
    Eq(q_vector[2], d*u1),
    Eq(q_vector[3], p/(gama-1.0) + 0.5*d*(u0**2 + u1**2)),
]

initial = GridBasedInitialisation()
initial.add_equations(init_eq)

schemes = {}
rk = RungeKuttaLS(4)
cent = StoreSome(4, 'u0 u1 T', merged=True)
schemes[rk.name] = rk
schemes[cent.name] = cent
block.set_discretisation_schemes(schemes)

# Boundary conditions.
boundaries = []
# Inner arc: isothermal no-slip wall
gama, Minf, Twall = symbols('gama Minf Twall', **{'cls': ConstantObject})
wall_energy = [Eq(q_vector[-1], q_vector[0]*Twall / (gama * Minf**2.0 * (gama - S.One)))]
boundaries += [IsothermalWallBC(direction=0, side=0, equations=wall_energy)]

# Outer arc: farfield
boundaries += [DirichletBC(direction=0, side=1, equations=init_eq)]

# Radial cut lines (y=0): symmetry boundaries.
boundaries += [SymmetryBC(direction=1, side=0)]
boundaries += [SymmetryBC(direction=1, side=1)]

block.set_block_boundaries(boundaries)

# Optional shock capturing filters
eq_classes = [constituent, simulation_eq, initial, metriceq]
if os.getenv('ENABLE_WENO', '0') == '1':
    eq_classes += WENOFilter(block, order=5, formulation='Z', flux_type='LLF', airfoil=False, metrics=metriceq).equation_classes
if os.getenv('ENABLE_TVD', '0') == '1':
    eq_classes += TVDFilter(block, airfoil=False).equation_classes

block.set_equations(eq_classes)

write_h5 = iohdf5(save_every=int(os.getenv('SAVE_EVERY_VALUE', '1000')), **{'iotype': 'Write'})
write_h5.add_arrays(simulation_eq.time_advance_arrays)
read_h5 = iohdf5(**{'iotype': 'Read'})
x, y = symbols('x0 x1', **{'cls': DataObject})
read_h5.add_arrays([x, y])
block.setio([write_h5, read_h5])

block.discretise()

np0 = int(simulation_parameters['block0np0'])
np1 = int(simulation_parameters['block0np1'])
SM = SimulationMonitor(
    [block.location_dataset('rho'), block.location_dataset('u0'), block.location_dataset('u1')],
    [
        (max(2, int(0.05*np0)), max(2, int(0.05*np1))),
        (max(2, int(0.5*np0)), max(2, int(0.5*np1))),
        (max(2, int(0.95*np0)), max(2, int(0.95*np1))),
    ],
    block,
    print_frequency=100,
    fp_precision=12,
    output_file='half_annulus_monitor.log',
)

alg = TraditionalAlgorithmRK(block, simulation_monitor=SM)
OPSC(alg, OPS_diagnostics=1)
if os.getenv('ENABLE_NAN_CHECK', '1') == '1':
    print_iteration_ops(NaN_check='rho', every=100, nblocks=1)
else:
    print_iteration_ops(every=100, nblocks=1)
substitute_simulation_parameters(simulation_parameters.keys(), simulation_parameters.values())
