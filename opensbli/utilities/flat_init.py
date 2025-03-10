from opensbli.utilities.numerical_functions import spline, splint
from sympy import Piecewise
from scipy.integrate import solve_ivp
import numpy as np
import numpy.polynomial.polynomial as poly
import matplotlib.pyplot as plt
from opensbli.initialisation import GridBasedInitialisation
from opensbli.core.opensbliobjects import DataObject, ConstantObject
from opensbli.core.grid import GridVariable
from opensbli.core.kernel import Kernel
import warnings
from opensbli.equation_types.opensbliequations import OpenSBLIEq

plt.style.use('classic')

class Newton_Main():

    # init method or constructor
    def __init__(self, req_tol, tol, delta):
        self.req_tol = req_tol
        self.tol = tol
        self.delta = delta

    def constants(self, eta_max, eta_points, sigma, gw, fwall):
        self.eta_max = eta_max
        self.eta_points = eta_points
        self.time_array = np.linspace(0, eta_max, eta_points)

        self.sigma = sigma
        self.gw = gw
        self.fwall = fwall

    def define_inputs(self, ode, ini_con, guess, pos_ini, boundary, pos_bd):
        self.ode = ode
        self.ini_con = ini_con
        self.guess = guess
        self.pos_ini = pos_ini
        self.boundary = boundary
        self.pos_bd = pos_bd

    # solve with newton method
    def solve_newton(self):
        # print('Solving for initial solutions')
        tolerance_stop = 1

        while tolerance_stop > self.req_tol: 
            # create a p matrix
            p_mat_size = (len(self.pos_bd),len(self.pos_bd))
            p_matrix = np.zeros(p_mat_size)

            boundary_array = []
            r_array = np.zeros(len(self.pos_bd))

            # assign variables
            for i in range(0, len(self.pos_ini)):
                self.ini_con[self.pos_ini[i]] = self.guess[i]

            # create delta matrix
            del_mat_size = (len(self.pos_ini)+1,len(self.ini_con))
            delta_matrix = np.zeros(del_mat_size)

            for g in range(0,len(self.pos_ini)):
                delta_matrix[g+1, self.pos_ini[g]] = self.delta
            
            y_initial = self.ini_con + delta_matrix

            for yi in range(0, len(y_initial)):
                # print(y_initial[yi, :])
                extract_cond = y_initial[yi, :]
                ini_conditions = extract_cond.tolist()

                ode_solution = solve_ivp(self.ode, [self.time_array[0], self.time_array[-1]], ini_conditions, t_eval = self.time_array, dense_output = 'true', rtol = 1e-10, atol = 1e-10, max_step = 0.01)

                boundary_array.append(ode_solution.y[:, -1])

            for p in range(0, len(self.pos_bd)):
                for  o in range(0, len(self.pos_bd)):
                    p_matrix[p, o] = (boundary_array[o+1][self.pos_bd[p]] - boundary_array[0][self.pos_bd[p]])/self.delta
                
                r_array[p] = self.boundary[p] - boundary_array[0][self.pos_bd[p]]

            resolve = np.linalg.solve(p_matrix, r_array.T)
            self.guess = self.guess + resolve.T
            
            # monitor tolerance
            tolerance_stop = abs(r_array[0])

    def new_ini(self):
        # print('------------------------------------ Returing calculated initial conditions in their respective positions')
        return self.ini_con 

class Boundary_layer_profile(object):
    """ Performs a similarity solution (Viscous fluid flow, F.White 1974),
    to obtain u and T profiles for a laminar compressible boundary-layer.

    :arg float xmach: Free-stream Mach numnber.
    :arg float Pr: Prandtl number.
    :arg float gama: Ratio of specific heats.
    :arg float Tw: Wall temperature, use -1 for adibatic wall conditions.
    :arg float Re: Free-stream Reynolds number.
    :arg float Tinf: Dimensional free-stream temperature for Sutherland's law."""

    def __init__(self, xmach, Pr, gama, Tw, Re, Tinf):
        self.y, self.u, self.T, self.scale = self.generate_boundary_layer_profile(xmach, Pr, gama, Tw, Tinf, Re)
        self.Re = Re
        self.n = np.size(self.y)
        return

    def compbl(self, t, y):
        
        suth = self.suth
        c = np.sqrt(y[3])*(1.0+suth)/(y[3]+suth)
        dcdg = 1.0/(2.0*np.sqrt(y[3])) - np.sqrt(y[3])/(y[3]+suth)
        dcdg *= (1.0+suth) / (y[3]+suth)
        cp = dcdg*y[4]

        return [y[1], y[2], -y[2]*(cp+y[0])/c, y[4], -y[4]*(cp+self.pr*y[0])/c - self.pr*(self.gama-1)*self.xmach**2 * y[2]**2]


    def generate_boundary_layer_profile(self, xmach, pr, gama, Tw, Tinf, Re):
        """ Generates a boundary layer initial profile. Solves the mean flow
        in a compressible boundary layer. (Equations 7.32) in White (1974).

        :arg float xmach: Mach number.
        :arg float pr: Prandtl number.
        :arg float gama: Ratio of specific heats.
        :arg float Tw: Wall temperature Tw/Tinf (< 0 for adiabatic)
        :arg float Tinf: Freestream reference temperature (Kelvin)
        :arg float Re: Reynolds number."""

        # self.nvisc = 1  # 1 Sutherlands law, 2 Power law, 3 Chapman-Rubesin approximation.
        self.pr, self.gama, self.xmach, self.Re, Tw = pr, gama, xmach, Re, Tw
        
        if Tw == -1:
            self.adiabatic = True
            self.Tw = 1.0 
            print('Setting conditions for an adiabatic wall')

            etamax, jmax = 10.0, 1041
            nstep = jmax - 1
            self.suth = 110.4/Tinf # sutherland's constant
            self.time_array = np.linspace(0, etamax, jmax)

        else:
            self.adiabatic = False
            self.Tw = Tw
            print('Setting conditions for an isothermal wall with wall temperature, Tw =',Tw)

            etamax, jmax = 20.0, 1041
            nstep = jmax - 1
            self.suth = 110.4/Tinf # sutherland's constant
            self.time_array = np.linspace(0, etamax, jmax)

        sigma, gw, fwall = 0.72, self.Tw, 0.0

        # newton method to get initial conditions
        fn_newton = Newton_Main(1e-9, 1e-10, 1e-3)
        fn_newton.constants(etamax, jmax, sigma, gw, fwall)

        if self.adiabatic == True:
            fn_newton.define_inputs(self.compbl, [0, 0, 0, gw, 0], [0.2, 10.0], [2, 3], [1, 1], [1, 3])
        else:
            fn_newton.define_inputs(self.compbl, [0, 0, 0, gw, 0], [0.2, 10.0], [2, 4], [1, 1], [1, 3])

        fn_newton.solve_newton()
        initial_conditions_newton = fn_newton.new_ini()
        
        etamax, jmax = 10.0, 1001
        nstep = jmax - 1
        self.suth = 110.4/Tinf # sutherland's constant
        self.time_array = np.linspace(0, etamax, jmax)

        ode_mtd1 = solve_ivp(self.compbl, [self.time_array[0], self.time_array[-1]], initial_conditions_newton, t_eval = self.time_array, dense_output = 'true', rtol = 1e-8, atol = 1e-8, max_step = 0.02)

        self.eta = ode_mtd1.t
        self.f, self.df, self.ddf = ode_mtd1.y[0], ode_mtd1.y[1], ode_mtd1.y[2]
        self.g, self.dg = ode_mtd1.y[3], ode_mtd1.y[4]

        # Save the wall temperature value
        self.Twall = self.g[0]
        print("The wall temperature is :", self.Twall)
        y, u, T, scale_factor = self.integrate_boundary_layer(nstep)
        self.scale_factor = scale_factor
        print("The scale factor is :", self.scale_factor)
        print("Wall normal derivative of velocity at the wall is :", self.dudy)
        return y, u, T, scale_factor

    def integrate_boundary_layer(self, n):
        """ Integrates the boundary-layer and calculates the scale factor from displacement thickness.

        :arg int n: Iteration number from the iterative solver.
        :returns: ndarray: y: Wall normal coordinates.
        :returns: ndarray: u: Streamwise velocity component profile.
        :returns: ndarray: T: Temperature profile.
        :returns: float: scale: Scale factor of the boundary-layer."""
        sumd, record_z = 0, 0
        z = np.zeros(n+1)
        d_eta = self.eta[1]*0.5
        # self.soln[1,:] is the u velocity, should be 1 in free stream
        for i in range(1, n+1):
            z[i] = z[i-1] + d_eta*(self.g[i] + self.g[i-1])
            dm = self.g[i-1]- self.df[i-1]
            dd = self.g[i] - self.df[i]
            sumd += d_eta*(dd+dm)
            if(self.df[i] > 0.999 and record_z < 1.0):
                # print "recording at iteration: ", i
                # dlta = z[i]
                record_z = 2.0
            scale = sumd
        # print("delta is :", dlta)
        # print("conversion factor is: ", scale)
        # print("scaled delta is: ", dlta/scale)
        # Rescale with displacement thickness and convert to FLOWER variable normalisation

        y, u, T = z/scale, self.df[:], self.g[:]
        # Calculate du/dy at the wall
        dy = y[1]
        # self.dudy = (-3*u[0]+4*u[1]-u[2])/(2.0*dy)
        self.dudy = (-1.83333333333334*u[0]+3.00000000000002*u[1]-1.50000000000003*u[2]+0.333333333333356*u[3]-8.34657956545823e-15*u[4]+1.06910315192207e-15*u[5])/dy
        self.dTdy = (-1.83333333333334*T[0]+3.00000000000002*T[1]-1.50000000000003*T[2]+0.333333333333356*T[3]-8.34657956545823e-15*T[4]+1.06910315192207e-15*T[5])/dy
        return y, u, T, scale

class Initialise_Flatplate(GridBasedInitialisation):
    """ Generates the initialiastion equations for the boundary-layer profile.

    :arg list npoints: Numerical values of the number of points in each direction.
    :arg list lengths: Numerical values of the problem dimensions.
    :arg list directions: Integer values of the problem directions.
    :arg list betas: Stretching factors for stretched grids.
    :arg int n_coeffs: Desired number of coefficients for the polynomial fit.
    :arg float Re: Reynolds number.
    :arg float xMach: Free-stream Mach number"""
    def __new__(cls, bl_directions, n_coeffs, Re, xMach, Tinf, Twall=True, coordinate_evaluations=None):
        ret = super(Initialise_Flatplate, cls).__new__(cls)

        if ret.find_constant_values([Twall])[0] == True:
            ret.find_constant_values([Twall])[0] == 'adiabtic'
            twallprint = 'adiabtic'
            print("Polynomial boundary-layer initialiastion called with Re = %f, Mach = %f, T_inf = %f. T_wall = %s. " % (Re, xMach, Tinf, twallprint))
        else:
            print("Polynomial boundary-layer initialiastion called with Re = %f, Mach = %f, T_inf = %f. T_wall = %f. " % (Re, xMach, Tinf, Twall))
        
        ret.coordinates = [x[1] for x in bl_directions]
        ret.bl_directions = bl_directions
        ret.n_coeffs = n_coeffs
        ret.coordinate_evaluations = coordinate_evaluations
        ret.Re = ret.find_constant_values([Re])[0]
        ret.Tinf = ret.find_constant_values([Tinf])[0]
        ret.equations = []
        ret.xMach = ret.find_constant_values([xMach])[0]
        ret.Twall = ret.find_constant_values([Twall])[0]
        ret.Adiabatic_condition = ret.find_constant_values([Twall])[0]
        return ret

    def find_constant_values(self, input):
        outlist = []
        for l in input:
            if isinstance(l, ConstantObject):
                if isinstance(l.value, str):
                    raise ValueError("")
                else:
                    outlist += [l.value]
            else:
                outlist += [l]
        return outlist

    def check_inputs(self, block):
        bl_directions = [x[0] for x in self.bl_directions]
        if sum(bl_directions) < 1:
            raise ValueError("Provide the directions to apply a boundary layer profile in.")
        if len(bl_directions) != block.ndim:
            raise ValueError("The list of polynomial directions must match the dimensions of the problem.")
        if self.n_coeffs < 10:
            raise ValueError("Higher number of polynomial coefficients are required for a good polynomial fit.")
        return

    def spatial_discretisation(self, block):
        self.equations = []
        self.block = block
        self.idxs = block.grid_indexes
        self.check_inputs(block)

        # Check if user has passed equations to evaluate coordinates, and add them to the kernel
        if self.coordinate_evaluations:
            self.equations += self.coordinate_evaluations

        # addition for initialisation around bumps
        x1b0 = [OpenSBLIEq(GridVariable('x1b0'), self.coordinate_evaluations[1].rhs.xreplace({self.block.grid_indexes[1] : 0.0}))]
        self.equations += x1b0

        self.initial = self.generate_initial_condition()
        # Add polynomial equations to initialise the solution
        self.equations += self.eqns

        # Convert inputs to DataSets
        self.equations = block.dataobjects_to_datasets_on_block(self.equations)
        # Ensure coordinate arrays are also restarted when required
        self.check_coordinate_evaluation(block)
        # Create the Katzer kernel
        from opensbli.core.kernel import Kernel
        flat_kernel = Kernel(block, computation_name="Similiarity solution laminar boundary-layer initialisation%d" % self.order)
        flat_kernel.set_grid_range(block)
        # Set halo range
        from opensbli.schemes.spatial.scheme import CentralHalos_defdec
        for d in range(block.ndim):
            # Initialize all five halos
            flat_kernel.set_halo_range(d, 0, CentralHalos_defdec())
            flat_kernel.set_halo_range(d, 1, CentralHalos_defdec())
        flat_kernel.add_equation(self.equations)
        flat_kernel.update_block_datasets(block)
        self.Kernels = [flat_kernel]
        return

    def generate_initial_condition(self):
        n_coeffs = self.n_coeffs
        # Load from similarity solution class
        y, u, T, rho, n = self.load_similarity()
        # Tolerance for finding the edge of the boundary layer.
        tolerance = 1e-10
        bl_directions = [x[0] for x in self.bl_directions]
        # print(bl_directions)

        if sum(bl_directions) == 1:  # 2D Flatplate and 3D spanwise periodic Flatplate, boundary layer in one direction
            # Create a large array of coordinates for this direction to interpolate the profile onto
            dire = [i for i, x in enumerate(self.bl_directions) if x[0]][0]
            poly_coordinates = self.uniform_1d_coordinate()
            # Interpolate u, T, rho onto the grid
            u_new = self.interpolate_onto_grid(y, poly_coordinates, u, self.dudy, 0)
            T_new = self.interpolate_onto_grid(y, poly_coordinates, T, 0, 0)
            rho_new = 1.0/T_new
            rhou_new = rho_new*u_new
            # Solve continuity equation to obtain rhov
            rhov_new = self.solve_continuity(poly_coordinates, u_new, rho_new)
            edge = self.find_edge_of_bl(u_new, tolerance)

            # Obtain polynomial fit coefficients
            rhou_coeffs = self.fit_polynomial(poly_coordinates, rhou_new, edge, n_coeffs)
            rhov_coeffs = self.fit_polynomial(poly_coordinates, rhov_new, edge, n_coeffs)
            T_coeffs = self.fit_polynomial(poly_coordinates, T_new, edge, n_coeffs)
            self.generate_one_wall_equations([rhou_new, rhov_new, T_new], [rhou_coeffs, rhov_coeffs, T_coeffs], dire, edge, poly_coordinates)

        elif sum(bl_directions) == 2:  # 3D with one side wall in x2 # WARNING: This currently does the same thing twice,
            edges, coeffs, profiles, normal_coeffs, normal_profiles = [], [], [], [], []
            directions = [i for i, x in enumerate(self.bl_directions) if x[0]]
            for dire in directions:  # For different bl thicknesses on different walls this still needs to be a loop
                poly_coordinates = self.uniform_1d_coordinate()
                # Interpolate u, and T onto the grid
                u_new = self.interpolate_onto_grid(y, poly_coordinates, u, self.dudy, 0)
                T_new = self.interpolate_onto_grid(y, poly_coordinates, T, 0, 0)
                # Temperature scaling function in region [0, 1]
                Tw, Tinf = self.Twall, 1.0  # Free-stream normalised temperature of 1.0, Wall temp taken from similarity solution
                g = self.temperature_scaling(T_new, Tw, Tinf)
                rho_new = 1.0/T_new
                rhou_new = rho_new*u_new
                # Solve continuity equation to obtain rho*wall_normal_velocity
                rho_vel_normal = self.solve_continuity(poly_coordinates, u_new, rho_new)
                profiles.append([rhou_new, T_new])
                normal_profiles.append(rho_vel_normal)
                edge = self.find_edge_of_bl(u_new, tolerance)
                edges.append(edge)
                # Obtain polynomial fit coefficients
                rhou_coeffs = self.fit_polynomial(poly_coordinates, rhou_new, edge, n_coeffs)
                g_coeffs = self.fit_polynomial(poly_coordinates, g, edge, n_coeffs)
                rho_vel_normal_coeffs = self.fit_polynomial(poly_coordinates, rho_vel_normal, edge, n_coeffs)
                coeffs.append([rhou_coeffs, g_coeffs])
                normal_coeffs.append(rho_vel_normal_coeffs)
            self.generate_two_wall_equations(profiles, coeffs, directions, edges, normal_profiles, normal_coeffs, poly_coordinates)
        else:
            raise NotImplementedError("Boundary layer initialisation is not implemented for walls in 3 dimensions.")
        return

    def uniform_1d_coordinate(self):
        n_elem = 10000
        return np.linspace(0, 20.0, n_elem)

    def temperature_scaling(self, temp_profile, Tw, Tinf):
        """ Computes the temperature profile between [0,1]

        :arg ndarray temp_profile: Temperature profile values between wall temperature and freestream.
        :arg float Tw: Wall temperature.
        :arg float Tinf: Free-stream temperature.
        :returns: ndarray: g: Temperature profile ranging from 0 to 1. """
        g = (temp_profile - Tw)/(Tinf - Tw)
        return g

    def form_equation(self, variable, name, coefficients, direction, edge, poly_coordinates):
        """ Creates the piecewise equations for the cases of 2D and 3D span-periodic boundary-layer profiles.

        :arg ndarray variable: Array of values for a given flow variable, used to obtain the free-stream value.
        :arg string name: Name of the variable.
        :arg ndarray coefficients: Coefficients for the polynomial fit.
        :arg int direction: Spatial direction to apply the equation to.
        :arg int edge: Grid index for the edge of the boundary-layer.
        returns: Eq: eqn: OpenSBLI equation to add to the initialisation kernel."""
        bl_edge_coordinate = poly_coordinates[edge]
        
        powers = [i for i in range(np.size(coefficients))][::-1]
        eqn = sum([coeff*(self.coordinates[direction] - GridVariable('x1b0'))**power for (coeff, power) in zip(coefficients, powers)])  # TODO set to exactl 1.0 if required

        # -----------------------------------------------------------------------------------------
        # here: potentiall addition of y+y0
        eqn = OpenSBLIEq(GridVariable('%s' % name), Piecewise((eqn, self.coordinates[direction] - GridVariable('x1b0') < bl_edge_coordinate), (variable[edge], True)))
        return eqn

    def form_mixed_equation(self, profiles, names, coefficients, directions, edges, normal_profiles, normal_coeffs, poly_coordinates):
        """ Generates the equations for the 3D SBLI sidewall case.

        :arg list profiles: Arrays of values for the rhou and [0,1] temperature profiles.
        :arg list names: Variable names as strings.
        :arg list coefficients: Coefficients for the polynomial fit for rhou and temperature profiles.
        :arg list directions: Directions that contain a wall.
        :arg list edges: Indices for the boundary layer edges in each direction.
        :arg list normal_profiles: Arrays of values for the rhov and rhow profiles.
        :arg list normal_coeffs: Coefficients for the polynomial fit for rhov and rhow.
        :returns: list: piecewise_eqns: Piecewise initialisation equations to be added to the initialisation class."""

        # Assuming we have the same number of poly coefficients in each direction, change later if required
        powers = [i for i in range(np.size(coefficients[0][0]))][::-1]
        # Loop over rhou, and T profiles
        piecewise_eqns = []
        direction1 = directions[0]
        direction2 = directions[1]
        # x variables and coordinates at the edge of the boundary layer in each
        coord1, coord2 = self.coordinates[direction1], self.coordinates[direction2]
        bl_coord1, bl_coord2 = poly_coordinates[edges[0]], poly_coordinates[edges[1]]
        # Create rhou profiles
        eqn1 = sum([coeff*coord1**power for (coeff, power) in zip(coefficients[0][0], powers)])
        eqn2 = sum([coeff*self.local_coordinate**power for (coeff, power) in zip(coefficients[1][0], powers)])  # profiles[0][max(edges)]
        u_var1, u_var2 = GridVariable('%s_1' % names[0]), GridVariable('%s_2' % names[0])
        # freestream_value = np.max([profiles[0][0][edges[0]], profiles[1][0][edges[1]]])
        freestream_value = 1.0
        piecewise_eqns.append(OpenSBLIEq(u_var1, Piecewise((eqn1, coord1 < bl_coord1), (freestream_value, True))))
        piecewise_eqns.append(OpenSBLIEq(u_var2, Piecewise((eqn2, self.local_coordinate < bl_coord2), (freestream_value, True))))
        piecewise_eqns.append(OpenSBLIEq(GridVariable('%s' % names[0]), u_var1*u_var2))
        # Create T profiles, g = (T-Tw)/(Tinf - Tw) ---> T = g*(Tinf-Tw) + Tw
        Tw = ConstantObject('Twall')
        eqn1 = sum([coeff*coord1**power for (coeff, power) in zip(coefficients[0][1], powers)])
        eqn2 = sum([coeff*self.local_coordinate**power for (coeff, power) in zip(coefficients[1][1], powers)])  # profiles[0][max(edges)]
        T_var1, T_var2 = GridVariable('%s_1' % names[1]), GridVariable('%s_2' % names[1])
        freestream_value = np.max([profiles[0][1][edges[0]], profiles[1][1][edges[1]]])
        freestream_value = 1.0  # CHECK THIS VALUE
        piecewise_eqns.append(OpenSBLIEq(T_var1, Piecewise((eqn1, coord1 < bl_coord1), (freestream_value, True))))
        piecewise_eqns.append(OpenSBLIEq(T_var2, Piecewise((eqn2, self.local_coordinate < bl_coord2), (freestream_value, True))))
        piecewise_eqns.append(OpenSBLIEq(GridVariable('%s' % names[1]), T_var1*T_var2*(freestream_value - Tw) + Tw))
        # Create normal velocity component profiles, rhov:
        eqn1 = sum([coeff*coord1**power for (coeff, power) in zip(normal_coeffs[0], powers)])
        temp1 = GridVariable('%s' % names[2])
        rhov_inf = normal_profiles[0][edges[0]]*u_var2  # Multiplying by rhou from the other direction
        piecewise_eqns.append(OpenSBLIEq(temp1, Piecewise((eqn1*u_var2, coord1 < bl_coord1), (rhov_inf, True))))
        # rhow:
        eqn2 = sum([coeff*self.local_coordinate**power for (coeff, power) in zip(normal_coeffs[1], powers)])
        temp2 = GridVariable('%s' % names[3])
        # rhow should reduce to zero at the symmetry plane
        rhow_inf = normal_profiles[1][edges[1]]*(1 - (self.local_coordinate-bl_coord2)/(0.5*self.domain_length-bl_coord2))*u_var1
        piecewise_eqns.append(OpenSBLIEq(temp2, self.side_fact*Piecewise((eqn2*u_var1, self.local_coordinate < bl_coord2), (rhow_inf, True))))
        return piecewise_eqns

    def generate_one_wall_equations(self, data, coeffs, direction, edge, poly_coordinates):
        """ Generates the equations for 2D SBLI and 3D span-periodic cases.

        :arg list data: Profile arrays for rhou0, rhou1 and temperature.
        :arg list coeffs: Coefficients for the polynomial fits.
        :arg list direction: Direction normal to the wall.
        :arg int edge: Grid index for the edge of the boundary-layer."""
        self.eqns = []
        rhou0_eqn = self.form_equation(data[0], 'rhou0', coeffs[0], direction, edge, poly_coordinates)
        rhou1_eqn = self.form_equation(data[1], 'rhou1', coeffs[1], direction, edge, poly_coordinates)
        T_eqn = self.form_equation(data[2], 'T', coeffs[2], direction, edge, poly_coordinates)

        # Set conservative values
        rho, rhou0, rhou1, T = GridVariable('rho'), GridVariable('rhou0'), GridVariable('rhou1'), GridVariable('T')
        rho_eqn = OpenSBLIEq(rho, 1.0/T)
        rho_store = OpenSBLIEq(DataObject('rho'), rho)
        rhou0_store = OpenSBLIEq(DataObject('rhou0'), rhou0)
        rhou1_store = OpenSBLIEq(DataObject('rhou1'), rhou1)
        gama, Minf = ConstantObject("gama"), ConstantObject("Minf")
        rhoE_store = OpenSBLIEq(DataObject('rhoE'), rho*T/(gama*(gama-1)*Minf**2) + 0.5*(rhou0**2 + rhou1**2)/rho)
        self.eqns += [rhou0_eqn, rhou1_eqn, T_eqn, rho_eqn, rho_store, rhou0_store, rhou1_store, rhoE_store]
        if self.block.ndim == 3:  # Periodic case, rhow = 0
            self.eqns += [OpenSBLIEq(DataObject('rhou2'), 0.0)]
        return

    def generate_two_wall_equations(self, profiles, coeffs, directions, edges, normal_profile, normal_coeffs, poly_coordinates):
        """ Generates the equations for 3D case with sidewalls.

        :arg list profiles: Profile arrays for rhou0 and temperature.
        :arg list coeffs: Coefficients for the polynomial fits.
        :arg list directions: Directions normal to the wall.
        :arg list edges: Grid indexes for the edge of the boundary-layer.
        :arg list normal_profile: Profile for the wall normal velocity components.
        :arg list normal_coeffs: Coefficients for the wall normal polynomial fit."""
        self.eqns = []
        names = ['rhou0', 'T', 'rhou1', 'rhou2']
        # Hard code to spanwise direction for 2 walls
        direction = 2
        array_coord = self.coordinates[direction]
        self.domain_length = self.block.deltas[direction]*(self.block.shape[direction]-1)
        self.side_fact = GridVariable('side_fact')
        side_fact_eqn = OpenSBLIEq(self.side_fact, Piecewise((1, array_coord <= self.domain_length/2.0), (-1, True)))

        self.local_coordinate = GridVariable('local_coord')
        local_coord_eqn = OpenSBLIEq(self.local_coordinate, Piecewise((array_coord, array_coord <= self.domain_length/2.0), (self.domain_length-array_coord, True)))
        self.eqns += [local_coord_eqn]
        self.eqns += [side_fact_eqn]
        # Create the piecewise equations formed from the boundary layers
        bl_equations = self.form_mixed_equation(profiles, names, coeffs, directions, edges, normal_profile, normal_coeffs, poly_coordinates)
        # Set conservative values
        rho, rhou0, rhou1, rhou2, T = GridVariable('rho'), GridVariable('rhou0'), GridVariable('rhou1'), GridVariable('rhou2'), GridVariable('T')
        rho_eqn = OpenSBLIEq(rho, 1.0/T)
        rho_store = OpenSBLIEq(DataObject('rho'), rho)
        rhou0_store = OpenSBLIEq(DataObject('rhou0'), rhou0)
        rhou1_store = OpenSBLIEq(DataObject('rhou1'), rhou1)
        rhou2_store = OpenSBLIEq(DataObject('rhou2'), rhou2)
        gama, Minf = ConstantObject("gama"), ConstantObject("Minf")
        rhoE_store = OpenSBLIEq(DataObject('rhoE'), rho*T/(gama*(gama-1)*Minf**2) + 0.5*(rhou0**2 + rhou1**2 + rhou2**2)/rho)
        self.eqns += bl_equations + [rho_eqn, rho_store, rhou0_store, rhou1_store, rhou2_store, rhoE_store]
        return

    def fit_polynomial(self, coords, variable, bl_edge, n_coeffs):
        """ Fits a polynomial to the input data, coefficients are returned.

        :arg ndarray coords: Independent variable of the input data.
        :arg ndarray variable: Dependent variable of the input data.
        :arg int bl_edge: Array index at the edge of the boundary-layer.
        :arg int n_coeffs: Desired number of coefficients for the polynomial.
        :returns: ndarray: coeffs: Coefficients of the polynomial fit."""
        coords = coords[0:bl_edge]
        variable = variable[0:bl_edge]
        # with warnings.catch_warnings():
        #     warnings.filterwarnings('error')
        #     try:
        warnings.filterwarnings("ignore")
        coeffs = poly.polyfit(coords, variable, n_coeffs)
        # except np.RankWarning:
        # print "Poorly conditioned fit"
        # ffit = poly.polyval(coords, coeffs)
        # plt.plot(coords, ffit, label='fit')
        # plt.plot(coords, variable, label='original_data')
        # plt.legend(loc="best")
        # plt.show()
        # Reverse coefficients so they are in descending order
        return coeffs[::-1]

    def solve_continuity(self, y, u, rho):
        """ Solves the continuity equation to obtain the wall normal velocity profile.

        :arg ndarray y: Dependent coordinate values.
        :arg ndarray u: Streamwise velocity component values.
        :arg ndarray rho: Density values.
        :returns: ndarray: rhov: Array of values for the wall normal velocity components. """
        # Grid offset delta to form derivative approximation
        n = np.size(y)
        ya2 = y[:]
        delta, scale, re = 0.001, self.scale_factor, self.Re
        rex0 = 0.5*(re/scale)**2
        x0 = 0.5*re/scale**2

        print('starting location of simulation x', x0)
        # print "Domain inlet is when the boundary-layer has developed for a length of x0 = %.10f" % x0
        drudx, rhov = np.zeros_like(y), np.zeros_like(y)
        # Local Reynolds number scaling to obtain a v profile
        sqrex = np.sqrt(rex0)
        delsx = np.sqrt(2.0)*scale*(x0)/sqrex
        ya2 = delsx*ya2
        d2y_u = spline(ya2, u, n, self.dudy, 0)
        d2y_rho = spline(ya2, rho, n, 0, 0)
        for j in range(0, n):
            ya2[j] = delsx*ya2[j]
            dstarp = delsx*np.sqrt((x0+delta)/(x0))
            dstarm = delsx*np.sqrt((x0-delta)/(x0))
            yp, ym = ya2[j]/dstarp, ya2[j]/dstarm
            uxp = splint(ya2, u, d2y_u, n, yp)
            uxm = splint(ya2, u, d2y_u, n, ym)
            rhoxp = splint(ya2, rho, d2y_rho, n, yp)
            rhoxm = splint(ya2, rho, d2y_rho, n, ym)
            drudx[j] = (rhoxp*uxp-rhoxm*uxm)/(2.0*delta)
            rhov[j] = rhov[j-1]-0.5*(ya2[j]-ya2[j-1])*(drudx[j]+drudx[j-1])
        return rhov

    def load_similarity(self):
        """ Solves the compressible boundary-layer equations via similarity solution."""
        Re, xMach, Tinf, Twall = self.Re, self.xMach, self.Tinf, self.Twall
        Adiabatic_condition = self.Adiabatic_condition
        Pr, gama = 0.72, 1.4  # Prandtl number, ratio of specific heats
        # bl = Boundary_layer_profile(xMach, Pr, gama, -1, Re, Tinf)  # -1 for Tw sets an adiabatic wall
        # bl = Boundary_layer_profile(xMach, Pr, gama, Twall, Re, Tinf)  # -1 for Tw sets an adiabatic wall

        if Adiabatic_condition == True:
            bl = Boundary_layer_profile(xMach, Pr, gama, -1, Re, Tinf)  # -1 for Tw sets an adiabatic wall
        else:
            bl = Boundary_layer_profile(xMach, Pr, gama, Twall, Re, Tinf)  # -1 for Tw sets an adiabatic wall

        # bl = Boundary_layer_profile(xMach, Pr, gama, 2.0, Re, Tinf)
        y, u, T, rho, n = bl.y, bl.u, bl.T, 1.0/bl.T, np.size(bl.y)
        self.Twall, self.scale_factor = bl.Twall, bl.scale_factor  # Wall temperature and scale factor from the similarity solution
        self.dudy = bl.dudy  # du/dy at the wall
        return y, u, T, rho, n

    def interpolate_onto_grid(self, y_in, y_out, var_in, y0, yn):
        # Create interpolating second derivative spline
        # n = size of the original data
        n = np.size(y_in)
        n_out = np.size(y_out)
        d2y = spline(y_in, var_in, n, y0, yn)
        # Array for variable interpolated onto the grid
        var_out = np.zeros_like(y_out)
        for i in range(n_out):
            var_out[i] = splint(y_in, var_in, d2y, n, y_out[i])
        return var_out

    def find_edge_of_bl(self, variable, tolerance):
        """ Finds the edge of the boundary layer and returns the index of that grid point.

        :arg ndarray variable: Array of values for a given flow variable.
        :arg float tolerance: Stopping tolerance for the difference between two successive grid points.
        :returns: int: index: Index of the boundary-layer edge."""
        index = 1
        while np.abs(variable[index]-variable[index-1]) > tolerance:
            index += 1
        return index

    def freestream_value(self, variable, index):
        """ Returns the value of a flow variable for a given grid index.

        :arg ndarray variable: Array of values for a given flow variable.
        :arg index Array index to use.
        :returns: float: variable[index]: Value of the flow variable at that index."""
        return variable[index]
