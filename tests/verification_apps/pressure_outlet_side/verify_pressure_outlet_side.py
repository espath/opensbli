#!/usr/bin/env python3
"""Tiny side0/side1 regression case for PressureOutletBC.

Goal:
- same tiny domain and equations
- switch outlet side in x-direction via env var
- keep all else identical
"""

from __future__ import annotations

import os

from opensbli import *
from opensbli.utilities.helperfunctions import substitute_simulation_parameters


outlet_side = int(os.getenv("OUTLET_SIDE", "1"))
if outlet_side not in (0, 1):
    raise ValueError("OUTLET_SIDE must be 0 or 1")

# Make flow point from inlet -> outlet depending on chosen side.
u0_inf = -1.0 if outlet_side == 0 else 1.0

simulation_parameters = {
    "gama": os.getenv("GAMMA_VALUE", "1.4"),
    "dt": os.getenv("DT_VALUE", "2.0e-4"),
    "niter": os.getenv("NITER_VALUE", "200"),
    "block0np0": os.getenv("NP0_VALUE", "81"),
    "block0np1": os.getenv("NP1_VALUE", "41"),
    "Delta0block0": "2.0/(block0np0)",
    "Delta1block0": "1.0/(block0np1)",
}

ndim = 2
coordinate_symbol = "x"
constants = ["gama"]
substitutions = []
sc1 = "**{'scheme':'Weno'}"

mass = "Eq(Der(rho,t), - Conservative(rhou_j,x_j,%s))" % sc1
momentum = "Eq(Der(rhou_i,t), -Conservative(rhou_i*u_j + KD(_i,_j)*p,x_j,%s))" % sc1
energy = "Eq(Der(rhoE,t), -Conservative((p+rhoE)*u_j,x_j,%s))" % sc1

velocity = "Eq(u_i, rhou_i/rho)"
pressure = "Eq(p, (gama-1)*(rhoE - rho*(1/2)*(KD(_i,_j)*u_i*u_j)))"
speed_of_sound = "Eq(a, (gama*p/rho)**0.5)"

eq = EinsteinEquation()
simulation_eq = SimulationEquations()
for base in [mass, momentum, energy]:
    simulation_eq.add_equations(eq.expand(base, ndim, coordinate_symbol, substitutions, constants))

constituent = ConstituentRelations()
for cr in [velocity, pressure, speed_of_sound]:
    constituent.add_equations(eq.expand(cr, ndim, coordinate_symbol, substitutions, constants))

block = SimulationBlock(ndim, block_number=0)

schemes = {}
rk = RungeKuttaLS(4)
schemes[rk.name] = rk
cent = StoreSome(4, "u0 u1", merged=True)
schemes[cent.name] = cent
block.set_discretisation_schemes(schemes)

q_vector = flatten(simulation_eq.time_advance_arrays)
d, u0, u1, p = symbols("d u0 u1 p", **{"cls": GridVariable})
gama = symbols("gama", **{"cls": ConstantObject})
init_eq = [
    Eq(d, 1.0),
    Eq(u0, u0_inf),
    Eq(u1, 0.0),
    Eq(p, 1.0),
    Eq(q_vector[0], d),
    Eq(q_vector[1], d * u0),
    Eq(q_vector[2], d * u1),
    Eq(q_vector[3], p / (gama - 1.0) + 0.5 * d * (u0**2 + u1**2)),
]
initial = GridBasedInitialisation()
initial.add_equations(init_eq)

boundaries = [[0, 0] for _ in range(ndim)]
inlet_side = 1 - outlet_side
boundaries[0][outlet_side] = PressureOutletBC(direction=0, side=outlet_side, back_pressure=1.0)
boundaries[0][inlet_side] = DirichletBC(direction=0, side=inlet_side, equations=init_eq)
# Remove corner complications: periodic in y
boundaries[1][0] = PeriodicBC(1, 0)
boundaries[1][1] = PeriodicBC(1, 1)
block.set_block_boundaries(boundaries)

kwargs = {"iotype": "Write"}
h5 = iohdf5(save_every=200, **kwargs)
h5.add_arrays(simulation_eq.time_advance_arrays)
block.setio([h5])

block.set_equations([constituent, simulation_eq, initial])
block.discretise()

monitor = SimulationMonitor(
    [block.location_dataset("rho"), block.location_dataset("u0"), block.location_dataset("u1")],
    [(5, 5), (40, 20), (75, 35)],
    block,
    print_frequency=20,
    fp_precision=12,
    output_file=f"pressure_outlet_side{outlet_side}_monitor.log",
)

alg = TraditionalAlgorithmRK(block, simulation_monitor=monitor)
SimulationDataType.set_datatype(Double)
OPSC(alg)
substitute_simulation_parameters(simulation_parameters.keys(), simulation_parameters.values())
print_iteration_ops(NaN_check="rho_B0", every=20)
