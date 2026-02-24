#!/usr/bin/env python3
import os
from sympy import Eq, S
from opensbli import *
from opensbli.physical_models.split_forms import NS_Split

simulation_parameters = {
    'gama': os.getenv('GAMMA_VALUE', '1.4'),
    'Minf': os.getenv('MINF_VALUE', '2.0'),
    'Pr': os.getenv('PR_VALUE', '0.71'),
    'Re': os.getenv('RE_VALUE', '1.0e4'),
    'dt': os.getenv('DT_VALUE', '1.0e-5'),
    'niter': os.getenv('NITER_VALUE', '50000'),
    'Twall': os.getenv('TWALL_VALUE', '1.0'),
    'SuthT': os.getenv('SUTHT_VALUE', '110.4'),
    'RefT': os.getenv('REFT_VALUE', '273.15'),
    'block0np0': os.getenv('NP0_VALUE', '241'),
    'block0np1': os.getenv('NP1_VALUE', '241'),
    # IMPORTANT: for this mesh layout, dir0 is angular and dir1 is radial.
    'Delta0block0': '((M_PI/2.0)/(block0np0-1.0))',
    'Delta1block0': '((20.0-0.5)/(block0np1-1.0))',
    'inv_rfact0_block0': '1.0/Delta0block0',
    'inv_rfact1_block0': '1.0/Delta1block0',
}

ndim=2
constants=['Re','Pr','gama','Minf','RefT','SuthT']
coordinate_symbol='x'
conservative=True
block=SimulationBlock(ndim, block_number=0)
SimulationDataType.set_datatype(Double)

metriceq=MetricsEquation()
metriceq.generate_transformations(ndim, coordinate_symbol, [(True,True),(True,True)], 2)

NS=NS_Split('KGP', ndim, constants, coordinate_symbol=coordinate_symbol, conservative=conservative, viscosity='dynamic', energy_formulation='enthalpy')
simulation_eq=SimulationEquations(); simulation_eq.add_equations(NS.mass); simulation_eq.add_equations(NS.momentum); simulation_eq.add_equations(NS.energy)

EE=EinsteinEquation(); EE.optional_subs_dict=metriceq.metric_subs
metric_vel='Eq(U_i, D_i_j*u_j)'
for eq in EE.expand(metric_vel, ndim, coordinate_symbol, [], constants):
    EE.optional_subs_dict[eq.lhs]=eq.rhs

constituent=ConstituentRelations()
for eq in EE.expand('Eq(u_i, rhou_i/rho)', ndim, coordinate_symbol, [], constants):
    constituent.add_equations(eq)
for expr in ['Eq(p, (gama-1)*(rhoE - (1/2)*rho*(KD(_i,_j)*u_i*u_j)))','Eq(T, p*gama*Minf*Minf/(rho))','Eq(mu, (T**(1.5)*(1.0+SuthT/RefT)/(T+SuthT/RefT)))','Eq(H, (rhoE + p) / rho)']:
    constituent.add_equations(EE.expand(expr, ndim, coordinate_symbol, [], constants))

simulation_eq.apply_metrics(metriceq)
q_vector=flatten(simulation_eq.time_advance_arrays)
d,u0,u1,p=symbols('d u0 u1 p', **{'cls':GridVariable})
gama,Minf=symbols('gama Minf', **{'cls':ConstantObject})
init_eq=[Eq(d,1.0),Eq(u0,1.0),Eq(u1,0.0),Eq(p,1.0/(gama*Minf**2.0)),Eq(q_vector[0],d),Eq(q_vector[1],d*u0),Eq(q_vector[2],d*u1),Eq(q_vector[3],p/(gama-1.0)+0.5*d*(u0**2+u1**2))]
initial=GridBasedInitialisation(); initial.add_equations(init_eq)

rk=RungeKuttaLS(4); cent=StoreSome(4,'u0 u1 T',merged=True); block.set_discretisation_schemes({rk.name:rk, cent.name:cent})

boundaries=[]
back_pressure = float(os.getenv('BACK_PRESSURE_VALUE', '1.0')) / (float(os.getenv('GAMMA_VALUE', '1.4')) * float(os.getenv('MINF_VALUE', '2.0'))**2.0)
boundaries += [PressureOutletBC(direction=0, side=0, back_pressure=back_pressure)]
boundaries += [SymmetryBC(direction=0, side=1)]
gama,Minf,Twall=symbols('gama Minf Twall', **{'cls':ConstantObject})
wall_energy=[Eq(q_vector[-1], q_vector[0]*Twall/(gama*Minf**2.0*(gama-S.One)))]
boundaries += [IsothermalWallBC(direction=1, side=0, equations=wall_energy)]
boundaries += [DirichletBC(direction=1, side=1, equations=init_eq)]
block.set_block_boundaries(boundaries)

block.set_equations([constituent, simulation_eq, initial, metriceq])
write_h5=iohdf5(save_every=int(os.getenv('SAVE_EVERY_VALUE','1000')), **{'iotype':'Write'}); write_h5.add_arrays(simulation_eq.time_advance_arrays)
read_h5=iohdf5(**{'iotype':'Read'}); x,y=symbols('x0 x1', **{'cls':DataObject}); read_h5.add_arrays([x,y]); block.setio([write_h5, read_h5])
block.discretise()

SM=SimulationMonitor([block.location_dataset('rho')],[(10,10)],block,print_frequency=100,fp_precision=12,output_file='quarter_ns_monitor.log')
alg=TraditionalAlgorithmRK(block, simulation_monitor=SM)
OPSC(alg, OPS_diagnostics=1)
if os.getenv('ENABLE_NAN_CHECK','1')=='1':
    print_iteration_ops(NaN_check='rho', every=100, nblocks=1)
else:
    print_iteration_ops(every=100, nblocks=1)
substitute_simulation_parameters(simulation_parameters.keys(), simulation_parameters.values())
