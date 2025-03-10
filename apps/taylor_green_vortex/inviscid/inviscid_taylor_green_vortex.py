#!/usr/bin/env python
# Import all the functions from opensbli
from opensbli import *
import copy
from opensbli.utilities.helperfunctions import substitute_simulation_parameters
import os, shutil, glob

simulation_parameters = {
'Re'        :   '1600.0',
'gama'      :   '1.4',
'Minf'      :   '0.4',
'Pr'        :   '0.71',
'dt'        :   '0.004',
'niter'     :   '20000',
'block0np0'     :   '64',
'block0np1'     :   '64',
'block0np2'     :   '64',
'Delta0block0'      :   '2*M_PI/block0np0',
'Delta1block0'      :   '2*M_PI/block0np1',
'Delta2block0'      :   '2*M_PI/block0np2',
}

split_schemes = ['Divergence', 'Blaisdell', 'Jameson', 'Kok', 'Feiereisen', 'KGP', 'KEEP']

for split_scheme in split_schemes:
    try:
        os.mkdir(split_scheme)
    except:
        pass
    # Number of dimensions of the system to be solved
    ndim = 3
    # # Constants that are used
    constants = ["Re", "Pr", "gama", "Minf", "mu"]
    # # symbol for the coordinate system in the equations
    coordinate_symbol = "x"
    # symbol for the coordinate system in the equations
    conservative = True
    NS = NS_Split(split_scheme, ndim, constants, coordinate_symbol=coordinate_symbol, conservative=conservative, viscosity='inviscid', energy_formulation='enthalpy')

    mass, momentum, energy = NS.mass, NS.momentum, NS.energy
    # Expand the simulation equations, for this create a simulation equations class
    simulation_eq = SimulationEquations()
    simulation_eq.add_equations(mass)
    simulation_eq.add_equations(momentum)
    simulation_eq.add_equations(energy)

    # Constituent relations used in the system
    velocity = "Eq(u_i, rhou_i/rho)"
    if conservative:
        pressure = "Eq(p, (gama-1)*(rhoE - (1/2)*rho*(KD(_i,_j)*u_i*u_j)))"
        velocity = "Eq(u_i, rhou_i/rho)"
    else:
        pressure = "Eq(p, rho*(gama-1)*(Et - (1/2)*(KD(_i,_j)*u_i*u_j)))"

    internal_energy = "Eq(e, p / (rho*(gama-1)))"
    enthalpy = "Eq(H, (rhoE + p) / rho)"
    temperature = "Eq(T, p*gama*Minf*Minf/(rho))"
    # divV = "Eq(divV, Der(u_j, x_j))"

    # Expand the constituent relations and them to the constituent relations class
    constituent = ConstituentRelations()  # Instantiate constituent relations object
    einstein_eq = EinsteinEquation()

    # Expand momentum add the expanded equations to the constituent relations
    if conservative:
        eqns = einstein_eq.expand(velocity, ndim, coordinate_symbol, [], constants)
        constituent.add_equations(eqns)
    # Expand pressure add the expanded equations to the constituent relations
    eqns = einstein_eq.expand(pressure, ndim, coordinate_symbol, [], constants)
    constituent.add_equations(eqns)
    # Expand temperature add the expanded equations to the constituent relations
    eqns = einstein_eq.expand(temperature, ndim, coordinate_symbol, [], constants)
    constituent.add_equations(eqns)
    # Expand temperature add the expanded equations to the constituent relations
    eqns = einstein_eq.expand(internal_energy, ndim, coordinate_symbol, [], constants)
    constituent.add_equations(eqns)
    # # Expand temperature add the expanded equations to the constituent relations
    # eqns = einstein_eq.expand(divV, ndim, coordinate_symbol, [], constants)
    # constituent.add_equations(eqns)
    # Expand temperature add the expanded equations to the constituent relations
    eqns = einstein_eq.expand(enthalpy, ndim, coordinate_symbol, [], constants)
    constituent.add_equations(eqns)


    # Write the expanded equations to a Latex file with a given name and titile
    latex = LatexWriter()
    latex.open('equations.tex', "Einstein Expansion of the simulation equations")
    latex.write_string("Simulation equations\n")
    for index, eq in enumerate(flatten(simulation_eq.equations)):
        latex.write_expression(eq)

    latex.write_string("Constituent relations\n")
    for index, eq in enumerate(flatten(constituent.equations)):
        latex.write_expression(eq)

    latex.close()

    # Create a simulation block
    block = SimulationBlock(ndim, block_number=0, conservative=conservative)

    # Local dictionary for parsing the expressions
    local_dict = {"block": block, "GridVariable": GridVariable, "DataObject": DataObject}

    # Initial conditions as strings
    x0 = "Eq(GridVariable(x0), block.deltas[0]*block.grid_indexes[0])"
    x1 = "Eq(GridVariable(x1), block.deltas[1]*block.grid_indexes[1])"
    x2 = "Eq(GridVariable(x2), block.deltas[2]*block.grid_indexes[2])"

    u0 = "Eq(GridVariable(u0),sin(x0)*cos(x1)*cos(x2))"
    u1 = "Eq(GridVariable(u1),-cos(x0)*sin(x1)*cos(x2))"
    u2 = "Eq(GridVariable(u2), 0.0)"
    p = "Eq(GridVariable(p), 1.0/(gama*Minf*Minf)+ (1.0/16.0) * (cos(2.0*x0)+cos(2.0*x1))*(2.0 + cos(2.0*x2)))"
    r = "Eq(GridVariable(r), gama*Minf*Minf*p)"


    if conservative:
        rho = "Eq(DataObject(rho), r)"
        rhou0 = "Eq(DataObject(rhou0), r*u0)"
        rhou1 = "Eq(DataObject(rhou1), r*u1)"
        rhou2 = "Eq(DataObject(rhou2), r*u2)"
        rhoE = "Eq(DataObject(rhoE), p/(gama-1) + 0.5* r *(u0**2+ u1**2 + u2**2))"
    else:
        rho = "Eq(DataObject(rho), r)"
        rhou0 = "Eq(DataObject(u0), u0)"
        rhou1 = "Eq(DataObject(u1), u1)"
        rhou2 = "Eq(DataObject(u2), u2)"
        rhoE = "Eq(DataObject(Et), p/(r*(gama-1)) + 0.5*(u0**2+ u1**2 + u2**2))"    

    eqns = [x0, x1, x2, u0, u1, u2, p, r, rho, rhou0, rhou1, rhou2, rhoE]

    # parse the initial conditions
    initial_equations = [parse_expr(eq, local_dict=local_dict) for eq in eqns]
    initial = GridBasedInitialisation()
    initial.add_equations(initial_equations)

    # Create a schemes dictionary to be used for discretisation
    schemes = {}
    # Central scheme for spatial discretisation and add to the schemes dictionary
    fns = 'u0 u1 u2'
    cent = StoreSome(4, fns)
    schemes[cent.name] = cent
    # RungeKutta scheme for temporal discretisation and add to the schemes dictionary
    rk = RungeKuttaLS(4)
    schemes[rk.name] = rk

    boundaries = []
    # Create boundaries, one for each side per dimension, so in total 6 BC's for 3D'
    for direction in range(ndim):
        boundaries += [PeriodicBC(direction, 0, full_depth=True)]
        boundaries += [PeriodicBC(direction, 1, full_depth=True)]

    # set the boundaries for the block
    block.set_block_boundaries(boundaries)
    # set the IO class to write out arrays
    kwargs = {'iotype': "Write"}
    h5 = iohdf5(save_every=100000000, **kwargs)
    h5.add_arrays(simulation_eq.time_advance_arrays)
    block.setio(copy.deepcopy(h5))
    # set the equations to be solved on the block

    # # Dispersion relation preserving filters
    # DRP = ExplicitFilter(block, [0,1,2], width=11, filter_type='DRP', optimized=True, sigma=0.2, multi_block=None)
    # block.set_equations(DRP.equation_classes)

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


    block.set_equations([copy.deepcopy(constituent), copy.deepcopy(simulation_eq), initial, post])
    # set the discretisation schemes
    block.set_discretisation_schemes(schemes)

    # Discretise the equations on the block
    block.discretise()

    # Simulation monitor
    arrays = ['KE', 'dilatation_dissipation', 'enstrophy_dissipation', 'rhom']
    probe_locations = [(None), (None), (None), (None)]
    # probe_locations = [None, None, None, None]
    SM = SimulationMonitor(arrays, probe_locations, block, output_file='TGV.log', print_frequency=100)

    # set the simulation data type, for more information on the datatypes see opensbli.core.datatypes
    SimulationDataType.set_datatype(Double)
    # Add the simulation monitor to the algorithm
    alg = TraditionalAlgorithmRK(block, simulation_monitor=SM)

    # Write the code for the algorithm
    OPSC(alg, OPS_diagnostics=2, OPS_V2=True)
    # Add the simulation constants to the OPS C code
    substitute_simulation_parameters(simulation_parameters.keys(), simulation_parameters.values())
    print_iteration_ops(NaN_check='rho')
    print("Finished generating code with: {}".format(split_scheme))
    # Move to the correct directory
    for file in glob.glob(r'./*.cpp'):
        shutil.copy(file, split_scheme)
    for file in glob.glob(r'./*.h'):
        shutil.copy(file, split_scheme)
    shutil.copy("Makefile", split_scheme)
