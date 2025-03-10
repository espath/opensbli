""" Max Walker & David J Lusher (09/2023): TVD non-linear filter for shock-capturing."""
from opensbli import *
from opensbli.equation_types.opensbliequations import NonSimulationEquations
from opensbli.postprocess.post_process_eq import *
from opensbli.schemes.spatial.TVD import *
from opensbli.filters.WENO_filter import NonLinearFilterBase

class TVDFilter(NonSimulationEquations, NonLinearFilterBase):
    """ Class to apply a TVD-based non-linear filter after a full time-step of a non-dissipative high order base scheme. The governing
    equations in the user script should be central derivatives in a skew-symmetric formulation to improve numerical stability."""
    def __init__(self, block, metrics=None, airfoil=False, species=None, optimize=False):
        print("Using non-linear TVD filtering on block {:}.".format(block.blocknumber))
        self.species = species
        # Get the shared functionality between TVD/WENO non-linear filters
        NonLinearFilterBase.__init__(self, airfoil, block, metrics, optimize=optimize)
        # Main class to generate the filter
        self.optimize = optimize
        self.main(block)
        return

    def main(self, block):
        """ Main calling function to generate the kernels for the TVD filter."""
        # Counter to order the kernels. Put the TVD filtering kernels at the very end of the time loop
        self.component_counter = 1000 + block.blocknumber*1000
        # Create the equations for TVD
        if self.species == 'passive_scalar':
            eqn = self.Euler_equations_passive_scalar(block, 'TVD')
        elif self.species == 'N N2':
            eqn = self.Euler_equations_N_N2(block, 'TVD')
        else:
            eqn = self.Euler_equations(block, 'TVD')
        # Convert the equations to datasets on this block
        self.equations = self.convert_to_datasets(block, eqn)
        # Create a TVD scheme
        self.SF = TVDFlux(averaging=SimpleAverage([0, 1]), shock_filter=True, conservative=block.conservative, species=self.species)
        self.halo_type = set()
        self.halo_type.add(self.SF.halotype)
        # Start the discretisation and create residual arrays for the equations
        self.Kernels = []
        self.create_residual_arrays(block)
        CR, solution_vector, reductions = self.SF.discretise(self, block)
        # Q vector
        self.solution_vector = flatten(self.time_advance_arrays)
        # Swap over the TVD stencil if periodic boundaries
        # bc_kernels = self.update_periodic_boundary(block, [-1, 1])
        # Shock sensor evaluation to find which points to evaluate the WENO scheme on
        if self.optimize:
            self.kappa = self.evaluate_shock_sensor(block)
        else:
            self.kappa = 1
        # Constituent relations evaluations on the Q vector at the end of the full RK time-step
        if self.species == 'N N2':
            self.constituent_relations_N_N2(block)
        else:
            self.constituent_relations(block)

        # Zero the work arrays
        self.zero_work_arrays(block)
        # Create the TVD reconstruction kernels
        reconstruction_kernels = []
        for direction, ker in enumerate(self.reconstruction_kernels):
            # Hybrid mode
            if self.optimize:
                ker = self.hybrid_condition(ker, block, direction)
            halo_ranges = ker.halo_ranges
            reconstruction_kernels.append(self.create_kernel('TVD reconstruction direction %d' % direction, ker.equations, halo_ranges, block))
            self.component_counter += 1

        self.add_kernel(reconstruction_kernels)
        # Check if there any wall boundary conditions or interfaces defined on the block.
        self.detect_wall_boundaries()
        self.detect_interface_boundaries()
        # # Create the residual kernel
        self.TVD_filter_application(block)
        return

    def evaluate_shock_sensor(self, block):
        # Add a shock sensor for the TVD filter
        SS = ShockSensor()
        if block.ndim > 1: # no shock sensor defined for ndim=1 currently
            # Ducros dilatation part
            sensor_evaluations, kappa = SS.ducros_equations(block, "x", metrics=self.metric_class, name='kappa')
        else:
            raise ValueError("No hybrid method for ndim=1.")
        # Update the CRs needed, and set the shock sensor to 1 on the outer boundaries
        self.constituent_relations(block, kappa)
        # Halo points for the sensor kernel
        sensor_halos = []
        for _ in range(self.ndim):
            sensor_halos.append([self.halo_type, self.halo_type])
        sensor_kernel = self.create_kernel('Shock sensor', flatten(sensor_evaluations), sensor_halos, block)
        self.add_kernel(sensor_kernel)
        self.component_counter += 1
        return kappa

    def TVD_filter_application(self, block):
        """ Applies the non-linear TVD filter by subtracting from the q vector after a full RK time-step."""
        resid_kernel = self.residual_kernels[0]
        filter_equations = []

        nvars = len(self.solution_vector)
        # Turn off the sensor at the walls
        wall_detection, wall_equations = self.wall_control(depth=5)
        filter_equations += wall_equations
        kappa_fact = self.kappa
        if isinstance(kappa_fact, DataSet):
            check = self.kappa
            for direction in range(self.ndim):
                for loc in [-1, 0, 1]:
                    check = Max(check, increment_dataset(self.kappa, direction, loc))        # for direction in range(self.ndim):
            kappa_fact = block.location_dataset('TVD_filter')
            # Selection of which kappa points to apply the filter to
            DS = ConstantObject('Ducros_select')
            DS.value = 0.05
            CTD.add_constant(DS)
            check = check >= DS
            cond1 = ExprCondPair(1, check)
            cond2 = ExprCondPair(0.0, True)
            filter_equations += [OpenSBLIEq(kappa_fact, Piecewise(*[cond1, cond2]))]
        else:
            kappa_fact = 1
        # detJ if needed, need to improve these scaling
        if self.curvilinear and self.airfoil:
            if self.ndim == 3:
                filter_equations += [OpenSBLIEq(gv('inv_detJ'), 1 / (Abs(block.location_dataset('detJ')) /  self.block.deltas[2])) ] ## Assumes span-periodic for now, for scaling
            else:
                filter_equations += [OpenSBLIEq(gv('inv_detJ'), 1 / Abs(block.location_dataset('detJ')))]
            detJ_term = gv('inv_detJ')
        else:
            detJ_term = 1

        dt = ConstantObject('dt')
        for i, eqn in enumerate(resid_kernel.equations):
            print(eqn)
            tvd_eqn = eqn.rhs.xreplace({ConstantObject('Delta%dblock%d' % (1, block.blocknumber)) : ConstantObject('Delta%dblock%d' % (1, block.blocknumber))*1/wall_detection})
            rhs = kappa_fact*ConstantObject('dt')*tvd_eqn*detJ_term
            filter_equations.append(OpenSBLIEq(self.solution_vector[i], self.solution_vector[i] + rhs))

   
        # Finish creating the kernel
        resid_kernel.equations = filter_equations
        residual_kernel = self.create_kernel('Non-linear TVD Filter application', filter_equations, resid_kernel.halo_ranges, block)
        self.component_counter += 1
        self.add_kernel(residual_kernel)
        return

    def hybrid_condition(self, kernel, block, direction):
        """ Checks the Ducros sensor, if it is a shock we perform the TVD reconstruction, else do nothing."""
        from sympy import And, Or
        input_equations = flatten(kernel.equations)
        kernel.equations = []

        if self.optimize:
            """ Only evaluate the TVD kernels at certain points, based on the shock sensor result. Improves performance."""
            DC = ConstantObject('Ducros_check')
            DC.value = 0.05
            CTD.add_constant(DC)
            locations = [-3,-2,-1,1,2]
            term = self.kappa
            # for dire in range(block.ndim):
            dire = direction # Only check 1D kappa
            for loc in locations:
                term = Max(term, increment_dataset(self.kappa, dire, loc))
            check = term > DC
            cond1 = ExprCondPair(input_equations, check)
            zeroed_equations = flatten([OpenSBLIEq(dset, 0.0) for dset in self.SF.temp_wk_arrays[direction]])
            cond2 = ExprCondPair(zeroed_equations, True)
            kernel.add_equation([GroupedPiecewise(cond1, cond2)])
        else:
            kernel.add_equation(input_equations)
        return kernel