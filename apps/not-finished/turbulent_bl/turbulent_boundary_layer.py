from opensbli import *
#from turbulent_inflow_bc import *
import numpy as np
import copy
import time_averaging
#import grid
from opensbli.core.boundary_conditions.turbulent_inflow.core import *
from opensbli.core.boundary_conditions.turbulent_inflow.bc import TurbulentInflowBC
from opensbli.equation_types.opensbliequations import OpenSBLIEquation as Eq
from sympy import pprint

#########################################################################################################################
#
# Inputs, constituent equations and scheme set-up
#
#########################################################################################################################

# Declare inputs
constants = ['gama', 'Minf', 'Pr', 'Re', 'Twall', 'dt', 'niter', 'block0np0', 'block0np1', 'block0np2', 'Lx0',
     'Lx1', 'Lx2', 'Delta0block0', 'Delta1block0', 'Delta2block0', 'SuthT', 'RefT', 'eps', 'TENO_CT', 'by', 'bz', 
     'restart_iteration_no','referenceRe','configuration']
values = ['1.4', '2.0', '0.72', '500.0', '1.676', '0.016', '10000', '800', '300', '90', '640.0', '40.0', '40.0', 
    'Lx0/(block0np0-1)', 'Lx1/(block0np1-1)', 'Lx2/(block0np2)', '110.4', '288.0', '1e-15', '1e-6','3.0', '0.0', 
    '0', '670','bottom']

inputs = get_all_vars(constants, values)

# Build dictionary of the inputs
input_dictionary = dict(zip(constants, values))

# define spacial scheme type, order and formulation (where applicable)
schm = 'teno'
schm_order = 6
formulation = 'Z'
curvilinear = True

# problem dimensions and coordinate symbol
ndim = 3
coordinate_symbol = "x"

# set statistics gathering on/off
stats = True
classes_stats = []
if stats:
    classes_stats = time_averaging.get_stats_classes()

Lx0 = ConstantObject('Lx0')
Lx1 = ConstantObject('Lx1')
Lx2 = ConstantObject('Lx2')
CTD.add_constant([Lx0, Lx1, Lx2])

# Create SimulationEquations and Constituent relations, add the expanded equations
metriceq = MetricsEquation()
metriceq.generate_transformations(ndim, coordinate_symbol, [(False, False), (True, False), (False, False)], 2)

# Define the compresible Navier-Stokes equations in Einstein notation.
sc1 = "**{\'scheme\':\'%s\'}" % schm.capitalize()
# governing equations for NS
if curvilinear:
    a = "Conservative(detJ * rho*U_j,xi_j,%s)" % sc1
    mass = "Eq(Der(rho,t), - %s/detJ)" % (a)
    a = "Conservative(detJ * (rhou_i*U_j + p*D_j_i), xi_j , %s)" % sc1
    momentum = "Eq(Der(rhou_i,t) , -  %s/detJ)" % (a)
    a = "Conservative(detJ * (p+rhoE)*U_j,xi_j, %s)" % sc1
    energy = "Eq(Der(rhoE,t), - %s/detJ)" % (a)
else:
    mass = "Eq(Der(rho,t), - Conservative(rho*u_j,x_j,%s))" % sc1 ## This equation changed from rhou_j to rho*u_j
    momentum = "Eq(Der(rhou_i,t) , - Conservative(rhou_i*u_j + KD(_i,_j)*p,x_j , %s) + Der(tau_i_j,x_j) )" % sc1
    energy = "Eq(Der(rhoE,t), - Conservative((p+rhoE)*u_j,x_j, %s) - Der(q_j,x_j) + Der(u_i*tau_i_j ,x_j) )" % sc1

# auxilliary equations for NS
stress_tensor = "Eq(tau_i_j, (mu/Re)*(Der(u_i,x_j)+ Der(u_j,x_i) - (2/3)* KD(_i,_j)* Der(u_k,x_k)))"
heat_flux = "Eq(q_j, (-mu/((gama-1)*Minf*Minf*Pr*Re))*Der(T,x_j))"

# Substitutions
substitutions = [stress_tensor, heat_flux]
constants_eqns = ["Re", "Pr", "gama", "Minf", "SuthT", "RefT"]
# Formulas for the variables used in the equations
velocity = "Eq(u_i, rhou_i/rho)"
pressure = "Eq(p, (gama-1)*(rhoE - rho*(1/2)*(KD(_i,_j)*u_i*u_j)))"
speed_of_sound = "Eq(a, (gama*p/rho)**0.5)"
temperature = "Eq(T, p*gama*Minf*Minf/(rho))"
viscosity = "Eq(mu, (T**(1.5)*(1.0+SuthT/RefT)/(T+SuthT/RefT)))"

einstein = EinsteinEquation()

# Viscous momentum and energy components
visc_momentum = "Eq(Der(rhou_i, t), Der(tau_i_j,x_j))"
visc_momentum = einstein.expand(visc_momentum, ndim, coordinate_symbol, substitutions, constants)
visc_momentum = [metriceq.apply_transformation(v) for v in visc_momentum]

visc_energy = "Eq(Der(rhoE, t), -Der(q_j,x_j) + Der(u_i*tau_i_j ,x_j))"
visc_energy = einstein.expand(visc_energy, ndim, coordinate_symbol, substitutions, constants)
visc_energy = metriceq.apply_transformation(visc_energy)

# Create an optional substitutions dictionary, this will be used to modify the equations when parsed
optional_subs_dict = metriceq.metric_subs
einstein.optional_subs_dict = optional_subs_dict

metric_vel = "Eq(U_i, D_i_j*u_j)"
simulation_eq = SimulationEquations()
constituent = ConstituentRelations()

# change coordinate symbol to curvilinear
if curvilinear:
    coordinate_symbol = "xi"
base_eqns = [mass, momentum, energy]
constituent_eqns = [velocity, pressure, speed_of_sound, temperature, viscosity, metric_vel]
# Expand the base equations
for i, base in enumerate(base_eqns):
    base_eqns[i] = einstein.expand(base, ndim, coordinate_symbol, substitutions, constants_eqns)
    if base==momentum:
        for no, b in enumerate(base_eqns[i]):
            base_eqns[i][no] = Eq(base_eqns[i][no].lhs, base_eqns[i][no].rhs + visc_momentum[no].rhs)
    else:
        if base==energy:
            base_eqns[i] = Eq(base_eqns[i].lhs, base_eqns[i].rhs + visc_energy.rhs)

# Expand the constituent relations
for i, CR in enumerate(constituent_eqns):
    constituent_eqns[i] = einstein.expand(CR, ndim, coordinate_symbol, substitutions, constants_eqns)

for eqn in base_eqns:
    simulation_eq.add_equations(eqn)

for eqn in constituent_eqns:
    constituent.add_equations(eqn)

# reset coordinate symbol
coordinate_symbol = "x"
 
# Initialise simulation block
block = SimulationBlock(ndim, block_number=0)

if schm == 'weno':
    # Set scheme parameters
    weno_order = schm_order
    Avg = RoeAverage([0, 1])
    LF = LFWeno(weno_order, formulation=formulation, averaging=Avg)
    schemes = {}
    schemes[LF.name] = LF
    cent = Central(4)
    schemes[cent.name] = cent

elif schm == 'teno':
    teno_order = schm_order
    Avg = RoeAverage([0, 1])
    LF = LFTeno(teno_order, averaging=Avg)
    schemes = {}
    schemes[LF.name] = LF
    cent = Central(4)
    schemes[cent.name] = cent

rk = RungeKuttaLS(3)
schemes[rk.name] = rk
block.set_discretisation_schemes(schemes)

#########################################################################################################################
#
# Mean profile and fluctuation amplitudes
#
#########################################################################################################################

Nt = int(input_dictionary['niter'])
Nx = int(input_dictionary['block0np0'])
Ny = int(input_dictionary['block0np1'])
Nz = int(input_dictionary['block0np2'])

# Get mean profiles, RMS scaling and turbulent length scales from turbulent_inflow
# wall_normal_v option turns on/off the generation of a wall-normal velocity profile (supersonic cases)
# fixed_seeds option allows for replicable random number generation
Umean_2D, Vmean_2D, Wmean_2D, Tmean_2D, Dmean_2D, A, seeds = inflow_init(constants, values, wall_normal_v=False, fixed_seed=False)

# Arrays for inflow boundary condition
halos = [-5, 5]
shape = tuple([Nz+2*halos[1], Ny+2*halos[1], 1])
Umean_in, Vmean_in, Wmean_in, Tmean_in, Dmean_in = [np.zeros(shape) for _ in range(5)]
A11_in, A21_in, A22_in, A31_in, A32_in, A33_in = [np.zeros(shape) for _ in range(6)]
seed1_in, seed2_in, seed3_in, seed4_in, seed5_in, seed6_in = [np.zeros(shape) for _ in range(6)]

# Convert from 1D array to 3D data sets for the inflow BC and the initial condition
for k in range(Nz+2*halos[1]):
    for j in range(Ny+2*halos[1]):
        for i in range(1):
            if j<halos[1]:
                j1 = 0
            elif j>=Ny + 5:
                j1 = Ny-1
            else:
                j1 = j-5
            if k<halos[1]:
                k1 = 0
            elif k>=Nz + 5:
                k1 = Nz-1
            else:
                k1 = k-5
            Umean_in[k,j,i] = Umean_2D[j1,k1]
            Vmean_in[k,j,i] = Vmean_2D[j1,k1]
            Wmean_in[k,j,i] = Wmean_2D[j1,k1]
            Tmean_in[k,j,i] = Tmean_2D[j1,k1]
            Dmean_in[k,j,i] = Dmean_2D[j1,k1]
            A11_in[k,j,i], A21_in[k,j,i], A22_in[k,j,i], A31_in[k,j,i], A32_in[k,j,i], A33_in[k,j,i] = [A[j1,k1,no] for no in range(6)]
            seed1_in[k,j,i], seed2_in[k,j,i], seed3_in[k,j,i], seed4_in[k,j,i], seed5_in[k,j,i], seed6_in[k,j,i] = [seeds[j1,k1,no] for no in range(6)]

# Declare constant values
M, gama, dt, restart_iteration_no = symbols("Minf gama dt restart_iteration_no", **{'cls':ConstantObject})
restart_iteration_no.datatype = Int()
CTD.add_constant([M, gama, dt, restart_iteration_no])

Mach = float(input_dictionary['Minf'])
g = float(input_dictionary['gama'])
total_x_points = Nx+2*halos[1]
division = 1.0/((g - 1.0)*g*Mach**2)

## data arrays for initial data
rho = np.repeat(Dmean_in, total_x_points, axis=2)
rhou = np.repeat(Umean_in, total_x_points, axis=2) * rho
rhov = np.repeat(Vmean_in, total_x_points, axis=2) * rho
rhow = np.repeat(Wmean_in, total_x_points, axis=2) * rho
rhoE = np.repeat(Tmean_in, total_x_points, axis=2) * rho * division + 0.5*(rhou**2 + rhov**2 + rhow**2)/rho

# Write data to HDF5
npoints = [1, Ny, Nz]
nhalos = [[0,0],[-5,5],[-5,5]]
to_write = [Umean_in, Vmean_in, Wmean_in, Dmean_in, A11_in, A21_in, A22_in, A31_in, A32_in, A33_in, seed1_in, seed2_in, seed3_in, seed4_in, seed5_in, seed6_in]
names = ['Umean', 'Vmean', 'Wmean', "Dmean", "A11", "A21", "A22", "A31", "A32", "A33", 'seed1', 'seed2', 'seed3', 'seed4', 'seed5', 'seed6',]
output_hdf5(to_write, names, nhalos, npoints, block, **{'filename':'correlationdata.h5'})

## write to initial.h5
npoints = [Nx, Ny, Nz]
nhalos = [[-5,5],[-5,5],[-5,5]]
to_write = [rho, rhou, rhov, rhow, rhoE]
names = ['rho', 'rhou0', 'rhou1', "rhou2", "rhoE"]
output_hdf5(to_write, names, nhalos, npoints, block, **{'filename':'initial.h5'})

# write zeros to tbc.h5
npoints = [Nx, Ny, Nz]
nhalos = [[-5,5],[-5,5],[-5,5]]
shape = tuple([Nz+2*halos[1], Ny+2*halos[1], Nx+2*halos[1]])
ustar, vstar, wstar =  [np.zeros(shape) for _ in range(3)]
phi_u_old, phi_v_old, phi_w_old = [np.zeros(shape) for _ in range(3)]
psi_u_old, psi_v_old, psi_w_old = [np.zeros(shape) for _ in range(3)]
to_write = [ustar, vstar, wstar, phi_u_old, phi_v_old, phi_w_old, psi_u_old, psi_v_old, psi_w_old]
names = ["ustar", "vstar", "wstar", "phi_u_old", "phi_v_old", "phi_w_old", "psi_u_old", "psi_v_old", "psi_w_old"]
output_hdf5(to_write, names, nhalos, npoints, block, **{'filename':'tbc.h5'})

# declare data objects for data read from hdf5
seed1, seed2, seed3, seed4, seed5, seed6 = DataObject('seed1'), DataObject('seed2'), DataObject('seed3'), DataObject('seed4'), DataObject('seed5'), DataObject('seed6')
A11, A21, A22, A31, A32, A33 = DataObject('A11'), DataObject('A21'), DataObject('A22'), DataObject('A31'), DataObject('A32'), DataObject('A33')
Umean, Vmean, Wmean, Dmean, Pmean = DataObject('Umean'), DataObject('Vmean'), DataObject('Wmean'), DataObject('Dmean'), DataObject('Pmean')
Uinit, Vinit, Winit, Tinit, Dinit, Pinit = DataObject('Umean'), DataObject('Vinit'), DataObject('Winit'), DataObject('Tinit'), DataObject('Dinit'), DataObject('Pinit')

#########################################################################################################################
#
# Grid
#
#########################################################################################################################

# variable dictionaries
local_dict = {"block": block, "GridVariable": GridVariable, "DataObject": DataObject}
local_dict['by'], local_dict['bz'] = ConstantObject('by'), ConstantObject('bz')
local_dict['block0np0'], local_dict['block0np1'], local_dict['block0np2'] = ConstantObject('block0np0'), ConstantObject('block0np1'), ConstantObject('block0np2')
# grid equations
x0 = parse_expr("Eq(DataObject(x0), block.deltas[0]*block.grid_indexes[0])", local_dict=local_dict)
x1 = parse_expr("Eq(DataObject(x1), Lx1*sinh(by*block.deltas[1]*block.grid_indexes[1]/Lx1)/sinh(by))", local_dict=local_dict)
x2 = parse_expr("Eq(DataObject(x2), block.deltas[2]*block.grid_indexes[2])", local_dict=local_dict)
# set grid initialisation
initial = GridBasedInitialisation()
initial_eqns = [x0, x1, x2]
initial.add_equations(initial_eqns)


#########################################################################################################################
#
# Boundary conditions
#
#########################################################################################################################

to_read = [Umean, Vmean, Wmean, Dmean, A11, A21, A22, A31, A32, A33, seed1, seed2, seed3, seed4, seed5, seed6]
kwargs = {'iotype': "Read", 'filename':'correlationdata.h5'}
h5_read1 = iohdf5(**kwargs)
h5_read1.add_arrays(to_read)
block.add_io([h5_read1])

# # Reading correlation data
to_read = [DataObject('ustar'), DataObject('vstar'), DataObject('wstar'), 
    DataObject('phi_u_old'), DataObject('phi_v_old'), DataObject('phi_w_old'), 
    DataObject('psi_u_old'), DataObject('psi_v_old'), DataObject('psi_w_old')]
kwargs = {'iotype': "Read", 'filename':'tbc.h5'}
h5_read2 = iohdf5(**kwargs)
h5_read2.add_arrays(to_read)
block.add_io([h5_read2])
boundaries = [[0, 0] for t in range(ndim)]

## turbulent inflow BC at inlet plane
direction = 0
side = 0
bc = TurbulentInflowBC(direction, side)

## apply integral length scale data
lx = [10.0, 4.0, 4.0]
ly = [1.5, 1.75, 1.0]
lz = [1.5, 1.0, 1.75]
lscale_y = Matrix([lx, ly, lz]) # bottom wall
lscale_z = Matrix([lx, [lz[0], lz[2], lz[1]], [ly[0], ly[2], ly[1]]]) # side wall
bc.lengthscales = lscale_y, lscale_z
boundaries[direction][side] = bc

# Right extrapolation at outlet
direction, side = 0, 1
boundaries[direction][side] = ExtrapolationBC(direction, side, order=0)

# Isothermal condition on bottom wall
direction = 1
side = 0
wall_const = ["Minf", "Twall"]
for con in wall_const:
    local_dict[con] = ConstantObject(con)

rhoE_wall = parse_expr("Eq(DataObject(rhoE), DataObject(rho)*Twall/(gama*(gama-1.0)*Minf**2.0))", local_dict=local_dict)
wall_eqns = [rhoE_wall]
boundaries[direction][side] = IsothermalWallBC(direction, side, wall_eqns)

# Extrapolation at top boundary
direction = 1
side = 1
boundaries[direction][side] = ExtrapolationBC(direction, side, order=0)

# Periodic condition on left and right walls
direction = 2
boundaries[direction][0] = PeriodicBC(direction, 0)
boundaries[direction][1] = PeriodicBC(direction, 1)

# set boundary conditions
block.set_block_boundaries(boundaries)

#########################################################################################################################
#
# setting equations, data output and generating ops code
#
#########################################################################################################################

## Interval for writing outputs
save_every = 10000

# Read restart file
kwargs = {'iotype': "Read", 'filename':'initial.h5'}
h5_read = iohdf5(**kwargs)
h5_read.add_arrays(simulation_eq.time_advance_arrays)
block.add_io([h5_read])

# Output q variables
kwargs = {'iotype': "Write"}
h5 = iohdf5(save_every=save_every, **kwargs)
h5.add_arrays(simulation_eq.time_advance_arrays)
block.setio(copy.deepcopy(h5))

# Output grid coordinates
kwargs = {'iotype': "Write", "name":"grid.h5"}
to_write = [DataObject('x0'), DataObject('x1'), DataObject('x2')]
h5 = iohdf5(**kwargs)
h5.add_arrays(to_write)
block.setio(copy.deepcopy(h5))

# Output correlation data
kwargs = {'iotype': "Write", "name":"tbc_output"}
to_write = [DataObject('ustar'), DataObject('vstar'), DataObject('wstar'), 
    DataObject('phi_u_old'), DataObject('phi_v_old'), DataObject('phi_w_old'), 
    DataObject('psi_u_old'), DataObject('psi_v_old'), DataObject('psi_w_old')]
h5 = iohdf5(save_every=save_every, **kwargs)
h5.add_arrays(to_write)
block.setio(copy.deepcopy(h5))

# Write statistics
if stats:
    stats_arrays = time_averaging.get_arrays()
    kwargs = {'iotype': "Write", "name":"statistics_output.h5"}
    h5 = iohdf5(**kwargs)
    h5.add_arrays(stats_arrays)
    block.setio(copy.deepcopy(h5))

# Set equations on the block and discretise
sim_eq = copy.copy(simulation_eq)
CR = copy.copy(constituent)
block.set_equations([CR, sim_eq, initial, metriceq] + classes_stats)
block.discretise()

# Force shock sensor to 1 in the halos, so that WENO is used at the boundary
kernel1 = Kernel(block, computation_name="sensor1")
kernel1.set_grid_range(block)
schemes = block.discretisation_schemes
for d in range(block.ndim):
    for sc in schemes:
        if schemes[sc].schemetype == "Spatial":
            kernel1.set_halo_range(d, 0, schemes[sc].halotype)
            kernel1.set_halo_range(d, 1, schemes[sc].halotype)
eq = block.dataobjects_to_datasets_on_block([OpenSBLIEq(DataObject('sensor'), 1.0)])
kernel1.add_equation(eq)
kernel1.update_block_datasets(block)
initial.Kernels += [kernel1]

# generate OPS code and substitute params
alg = TraditionalAlgorithmRK(block)
SimulationDataType.set_datatype(Double)
OPSC(alg, OPS_diagnostics=5)
substitute_simulation_parameters(constants, values)
print_iteration_ops(every=1, NaN_check='rho_B0')
