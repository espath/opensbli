#!/usr/bin/env python
from opensbli import *
import copy
from opensbli.utilities.helperfunctions import substitute_simulation_parameters
from TGV_setup import * # Multi-block configuration of blocks and initial condition
###################################################################################################
simulation_parameters = {
'Re'        :   '1600.0',
'gama'      :   '1.4',
'Minf'      :   '1.25',
'Pr'        :   '0.71',
'dt'        :   '0.0005',
'niter'     :   '40000',
'shock_factor'      :   '1',
'TENO_CT'       :   '1e-6',
'eps'       :   '1.0e-16',
'teno_a1'       :   '10.5',
'teno_a2'       :   '4.5',
}
##################################################################################################
# Number of dimensions of the system to be solved
ndim = 3
nblocks = 2
multi_block = MultiBlock(ndim, nblocks)
# Add block-specific simulation parameters
simulation_parameters = generate_simulation_parameters(simulation_parameters, multi_block)
# # Constants that are used
constants = ["Re", "Pr", "gama", "Minf", "SuthT", "RefT"]
teno = True
einstein_eq = EinsteinEquation()
# Central scheme plus WENO filtering, otherwise pure TENO
if not teno:
    NS = NS_Split('KGP', ndim, constants, coordinate_symbol="x", conservative=True, viscosity='dynamic', energy_formulation='enthalpy', debug=False)
    mass, momentum, energy = NS.mass, NS.momentum, NS.energy
    # Expand the simulation equations, for this create a simulation equations class
    simulation_eq = SimulationEquations()
    simulation_eq.add_equations([mass, momentum, energy])
else:
    sc1 = "**{\'scheme\':\'Teno\'}"
    # Define the compresible Navier-Stokes equations in Einstein notation.
    mass = "Eq(Der(rho,t), - Conservative(rhou_j,x_j,%s))" % sc1
    momentum = "Eq(Der(rhou_i,t) , -Conservative(rhou_i*u_j + KD(_i,_j)*p,x_j , %s) + Der(tau_i_j,x_j) )" % sc1
    energy = "Eq(Der(rhoE,t), - Conservative((p+rhoE)*u_j,x_j, %s) - Der(q_j,x_j) + Der(u_i*tau_i_j ,x_j) )" % sc1
    stress_tensor = "Eq(tau_i_j, (mu/Re)*(Der(u_i,x_j)+ Der(u_j,x_i) - (2/3)* KD(_i,_j)* Der(u_k,x_k)))"
    heat_flux = "Eq(q_j, (-mu/((gama-1)*Minf*Minf*Pr*Re))*Der(T,x_j))"
    # Substitutions
    substitutions = [stress_tensor, heat_flux]
    base_eqns = [mass, momentum, energy]
    # Expand the base equations
    for i, base in enumerate(base_eqns):
        base_eqns[i] = einstein_eq.expand(base, ndim, "x", substitutions, constants)
    simulation_eq = SimulationEquations()
    for eqn in base_eqns:
        simulation_eq.add_equations(eqn)
##################################################################################################
# Constituent relations
pressure = "Eq(p, (gama-1)*(rhoE - (1/2)*rho*(KD(_i,_j)*u_i*u_j)))"
velocity = "Eq(u_i, rhou_i/rho)"
enthalpy = "Eq(H, (rhoE + p) / rho)"
temperature = "Eq(T, p*gama*Minf*Minf/(rho))"
viscosity = "Eq(mu, (T**(1.5)*(1.4042)/(T+0.40417)))" ## Modified sutherland law
speed_of_sound = "Eq(a, (gama*p/rho)**0.5)"
# Expand the constituent relations and them to the constituent relations class
constituent = ConstituentRelations()  # Instantiate constituent relations object
for input_eqn in [velocity, pressure, temperature, enthalpy, viscosity, speed_of_sound]:
    eqns = einstein_eq.expand(input_eqn, ndim, "x", [], constants)
    constituent.add_equations(eqns)
##################################################################################################
# Create a schemes dictionary to be used for discretisation
schemes = {}
# Central scheme for spatial discretisation and add to the schemes dictionary
fns = 'u0 u1 u2'
cent = StoreSome(4, fns)
schemes[cent.name] = cent
# RungeKutta scheme for temporal discretisation and add to the schemes dictionary
rk = RungeKuttaLS(3, formulation='SSP')
schemes[rk.name] = rk
if teno:
    halos = [-4, 4]
    teno_order = 6
    Avg = RoeAverage([0, 1])
    LF = LFTeno(teno_order, averaging=Avg, flux_type='LLF')
    schemes[LF.name] = LF
else:
    halos = [-2, 2]
    exit()
    # WENO filter for shock-capturing
    # WF = WENOFilter(block, order=7, dissipation_sensor='Ducros', flux_type='HLLC', airfoil=False, store_filter=True)
    # block.set_equations(WF.equation_classes)
# Set the discretisation schemes
multi_block.set_discretisation_schemes(schemes)
##################################################################################################
# Set the initial conditions on each of the blocks
mb_initial_conditions = TGV_initial_condition(multi_block)
multi_block.set_initial_conditions(mb_initial_conditions)
##################################################################################################
mb_bcs = TGV_boundaries(multi_block, halos)
multi_block.set_block_boundaries(mb_bcs)
##################################################################################################
# # Add DRP filtering on each block
# filter_list = []
# for no, block in enumerate(multi_block.blocks):
#     filter_list += [ExplicitFilter(block, [0,1,2], width=9, filter_type='DRP', optimized=False, sigma=0.2, multi_block=multi_block).equation_classes]
# multi_block.set_filters(filter_list)
##################################################################################################
## Post-processing for TGV case, kinetic energy and enstrophy reductions
# Velocity in 3D
vel = symbols("u0:%d"%ndim,  **{'cls':DataObject})
# Vorticity-z
wx, wy, wz = symbols("wx wy wz",  **{'cls':GridVariable})
# coordinates
coord = symbols("x0:%d"%ndim,  **{'cls':CoordinateObject})
# Matrix of derivatives
der_matrix = Matrix(ndim,ndim,[CentralDerivative(u,x) for u in vel for x in coord])
post = UserDefinedEquations()
post.kernel_merge = True
post.algorithm_place = InTheSimulation(frequency=100)
post.computation_name = 'Taylor-Green vortex post-processing'
post.order = 10000000 # appear at the end of the kernels at the end of the time-loop
# # X vorticity
vortx = Eq(wx, der_matrix[2,1] - der_matrix[1,2])
# Y vorticity
vorty = Eq(wy, der_matrix[0,2] - der_matrix[2,0])
# Z vorticity
vortz = Eq(wz, der_matrix[1,0] - der_matrix[0,1])
# Dilatation
divV = symbols("divV", **{'cls':DataObject})
dil = Eq(divV, der_matrix[0,0] + der_matrix[1,1] + der_matrix[2,2])
post.add_equations([vortx, vorty, vortz, dil])
# Evaluate reduction quantities required for dissipation measures
rhom, KE, eps_D, eps_S = ReductionSum('rhom'), ReductionSum('KE'), ReductionSum('dilatation_dissipation'), ReductionSum('enstrophy_dissipation')
rho_eqn = OpenSBLIEq(rhom, rhom + DataObject('rho'))
ke_eqn = OpenSBLIEq(KE, KE + 0.5*DataObject('rho')*sum([u**2 for u in vel]))
dilatation_eqn = OpenSBLIEq(eps_D, eps_D + Rational(4,3)*DataObject('mu')*divV**2)
enstrophy_eqn = OpenSBLIEq(eps_S, eps_S + DataObject('mu')*(wx**2 + wy**2 + wz**2))
post.add_equations([rho_eqn, ke_eqn, dilatation_eqn, enstrophy_eqn])
##################################################################################################
multi_block.set_equations([simulation_eq, constituent, post])
##################################################################################################
q_vector = simulation_eq.time_advance_arrays
x,y,z = symbols("x0, x1, x2", **{'cls':DataObject})
kwargs = {'iotype': "Write"}
h5 = iohdf5(save_every=10000, **kwargs)
h5.add_arrays(q_vector + [x, y, z])
multi_block.setio([h5])
##################################################################################################
# Slicing output
for bn in range(nblocks):
    coordinate_slice = iohdf5_slices(blocknumber=bn, **{'iotype': "Init"})
    coordinate_slice.add_slices([([DataObject('x1'), DataObject('x2')], 0, 'block%dnp0/2' % bn)])
    q_slice = iohdf5_slices(save_every=1000, blocknumber=bn, **{'iotype': "Write"})
    q_slice.add_slices([(q_vector, 0, 'block%dnp0/2' % bn)])
    multi_block.setio([coordinate_slice, q_slice])
##################################################################################################
# Discretise the equations on the block
multi_block.discretise()
# Apply a periodic BC for WENO filter
if not teno:
    WF.update_periodic_boundary(block, halos=[-5,5])
##################################################################################################
# Simulation monitor
arrays = [['KE', 'dilatation_dissipation', 'enstrophy_dissipation', 'rhom', 'p'] for _ in range(nblocks)] # monitor the same quantities on every block
probe_locations = [[(None), (None), (None), (None), (3, 4, 5)] for _ in range(nblocks)]
SM = SimulationMonitor(arrays, probe_locations, multi_block, output_file='TGV.log', print_frequency=100)
##################################################################################################
alg = TraditionalAlgorithmRKMB(multi_block, SM)
SimulationDataType.set_datatype(Double)
##################################################################################################
# Write the code for the algorithm
OPSC(alg, OPS_diagnostics=5, OPS_V2=True)
# Add the simulation constants to the OPS C code
substitute_simulation_parameters(simulation_parameters.keys(), simulation_parameters.values())
print_iteration_ops(NaN_check='rho', nblocks=nblocks)
##################################################################################################