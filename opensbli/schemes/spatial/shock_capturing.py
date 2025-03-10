"""@brief
   @author Satya Pramod Jammy, David J Lusher
   @contributors Max Walker
   @details
"""

from sympy import Symbol, Rational, zeros, Abs, Matrix, flatten, Max, diag, Function, count_ops, simplify, factor, sign, Min, symbols, sqrt, And, Piecewise, pi, sin, Equality
from sympy.core.numbers import Zero
from opensbli.core.opensbliobjects import EinsteinTerm, DataSetBase, ConstantObject, DataSet, DataObject, ReductionMax
from opensbli.equation_types.opensbliequations import OpenSBLIEq
from opensbli.core.kernel import Kernel, ConstantsToDeclare
from opensbli.core.grid import GridVariable as gv
from opensbli.utilities.helperfunctions import increment_dataset
from opensbli.physical_models.euler_eigensystem import EulerEquations
from sympy import factor, pprint
from opensbli.schemes.spatial.averaging import SimpleAverage, RoeAverage
from opensbli.core.grid import GridVariable


class ShockCapturing(object):
    """ Class that contains functions common to all shock capturing schemes."""

    def evaluate_residuals(self, block, eqns, local_ders):
        residual_eqn = []
        for eq in eqns:
            substitutions = {}
            for d in eq.rhs.atoms(Function):
                if d in local_ders:
                    substitutions[d] = d._discretise_derivative(block)
                else:
                    substitutions[d] = 0
            residual_eqn += [OpenSBLIEq(eq.residual, eq.rhs.subs(substitutions))]
        residual_kernel = Kernel(block, computation_name="%s Residual" % self.__class__.__name__)
        residual_kernel.set_grid_range(block)
        residual_kernel.add_equation(residual_eqn)
        return residual_kernel

    def generate_left_reconstruction_variables(self, expression_matrix, derivatives):
        if isinstance(expression_matrix, Matrix):
            stencil = expression_matrix.stencil_points
            for i in range(expression_matrix.shape[0]):
                settings = derivatives[i].settings
                if "single_reconstruction_variable" in settings and settings["single_reconstruction_variable"]:
                    name = "Recon_%d" % i
                else:
                    name = 'L_X%d_%d' % (self.direction, i)
                lv = type(self.reconstruction_classes[1])(name)
                lv.settings = settings
                for p in sorted(set(self.reconstruction_classes[1].func_points)):
                    lv.function_stencil_dictionary[p] = expression_matrix[i, stencil.index(p)]
                derivatives[i].add_reconstruction_classes([lv])
        else:
            raise TypeError("Input should be a matrix.")
        return

    def generate_right_reconstruction_variables(self, expression_matrix, derivatives):
        if isinstance(expression_matrix, Matrix):
            stencil = expression_matrix.stencil_points
            for i in range(expression_matrix.shape[0]):
                settings = derivatives[i].settings
                if "single_reconstruction_variable" in settings and settings["single_reconstruction_variable"]:
                    name = "Recon_%d" % i
                else:
                    name = 'R_X%d_%d' % (self.direction, i)
                rv = type(self.reconstruction_classes[0])(name)
                for p in sorted(set(self.reconstruction_classes[0].func_points)):
                    rv.function_stencil_dictionary[p] = expression_matrix[i, stencil.index(p)]
                rv.settings = settings
                derivatives[i].add_reconstruction_classes([rv])
        else:
            raise TypeError("Input should be a matrix.")
        return

    def interpolate_reconstruction_variables(self, derivatives, block, MP_limiter=False, positivity_preservation=False):
        """ Perform the WENO/TENO interpolation on the reconstruction variables.

        :arg list derivatives: A list of the TENO derivatives to be computed.
        :arg object kernel: The current computational kernel."""
        output_eqns = []
        # Extract info from the RVs
        # self.process_rv(derivatives)
        for no, d in enumerate(derivatives):
            for rv in d.reconstructions:
                if isinstance(rv, type(self.reconstruction_classes[1])):
                    original_rv = self.reconstruction_classes[1]
                elif isinstance(rv, type(self.reconstruction_classes[0])):
                    original_rv = self.reconstruction_classes[0]
                else:
                    raise ValueError("Reconstruction must be left or right")
                rv.update_quantities(original_rv)
                rv.evaluate_quantities()
                # Apply a sensor to each characteristic wave if using filtering methods
                if self.sensor_evaluation is not None:
                    output_eqns += [rv.final_equations[0:-1]]
                    if isinstance(rv, type(self.reconstruction_classes[0])):
                        output_eqns += [OpenSBLIEq(gv('rj%d' % no), self.sensor_evaluation[0].rhs)]
                    elif isinstance(rv, type(self.reconstruction_classes[1])):
                        output_eqns += [OpenSBLIEq(gv('rj%d' % no), Max(gv('rj%d' % no),self.sensor_evaluation[1].rhs))]
                        # output_eqns += [OpenSBLIEq(gv('rj%d' % no), self.sensor_evaluation[1].rhs)]
                    output_eqns += [rv.final_equations[-1]]
                else:
                    output_eqns += [rv.final_equations]
                # if MP_limiter:
                #     output_eqns += self.monotonicity_limiter(rv)
        return output_eqns


    def DMM(self, x, y):
        """ MinMod limiter function for 2 arguments."""
        return 0.5*(sign(x, evaluate=False) + sign(y, evaluate=False))*Min(Abs(x), Abs(y))

    def DM4(self, w, x, y, z):
        """ MinMod limiter function for 4 arguments."""
        return 0.125*(sgn(w, evaluate=False) + sgn(x, evaluate=False))*Abs((sgn(w, evaluate=False) + sgn(y, evaluate=False))*(sgn(w, evaluate=False) + sgn(z, evaluate=False)))*Min(Abs(w), Min(Abs(x), Min(Abs(y),Abs(z))))

    def median(self, x, y, z):
        return x + self.DMM(y-x, z-x)

    def process_rv(self, derivatives):
        # Collects the flux terms used for the reconstructions, per derivative
        for no, d in enumerate(derivatives):
            fn = {}
            for rv in d.reconstructions:
                for key, value in rv.function_stencil_dictionary.items():
                    for term in value.atoms(gv):
                        if 'CF' in str(term):
                            fn[key] = term
            # Add to both RV
            for rv in d.reconstructions:
                rv.fluxes = fn
        return


    def monotonicity_limiter(self, rv):
        from opensbli.schemes.spatial.weno import RightWenoReconstructionVariable, LeftWenoReconstructionVariable
        from opensbli.schemes.spatial.teno import RightTenoReconstructionVariable, LeftTenoReconstructionVariable
        output_eqns = []
        # Original reconstructed value
        original = rv.reconstructed_symbol
        # Variable to control CFL number that can be used, default alpha=2
        alpha = ConstantObject('alpha_MP')
        alpha.value = 2.0
        ConstantsToDeclare.add_constant(alpha)
        # Function values required to take the limits
        fn = rv.fluxes
        print(fn)
        # dj = u_(j+1) - 2*u_(j) + u_(j-1) curvature measure
        if isinstance(rv, RightWenoReconstructionVariable) or isinstance(rv, RightTenoReconstructionVariable): # e.g. [-2,2] stencil, 5th order WENO
            # # f(i), f(i+1), f(i+1/2)^MD
            # djm1 = fn[0] - 2*fn[-1] + fn[-2]
            # dj = fn[1]-2*fn[0]+fn[-1]
            # djp1 = fn[2]-2*fn[1]+fn[0]
            # MM1, MM2 = gv('MM1'), gv('MM2')
            # output_eqns += [OpenSBLIEq(MM1, self.DMM(dj, djp1))]
            # output_eqns += [OpenSBLIEq(MM2, self.DMM(djm1, dj))]
            # # f_j, f_(j+1), 0.5*(u_j + u_(j+1)) - 0.5*MM(d_j, d_(j+1))
            # c_1 = [fn[0], fn[1], 0.5*(fn[0] + fn[1]) - 0.5*MM1]
            # # f_j, f_j + alpha*(f_j - f_(j+1)), f_j + 0.5*(f_j - f_(j+1)) + 4/3 * MM(dj, d_(j-1))
            # c_2 = [fn[0], fn[0] + alpha*(fn[0] - fn[1]), fn[0] + 0.5*(fn[0] - fn[-1]) + Rational(4,3)*MM2]
            # # The two lists of arguments
            # u_min = Max(Min(c_1[0], Min(c_1[1], c_1[2])), Min(c_2[0], Min(c_2[1], c_2[2])))
            # u_max = Min(Max(c_1[0], Max(c_1[1], c_1[2])), Max(c_2[0], Max(c_2[1], c_2[2])))
            # # rv.limiter = gv('limiter')
            # output_eqns += [OpenSBLIEq(original, self.median(original, u_min, u_max))]
            pass
        elif isinstance(rv, LeftWenoReconstructionVariable) or isinstance(rv, LeftTenoReconstructionVariable):
            # f(i), f(i+1), f(i+1/2)^MD
            djm1 = fn[1] - 2*fn[2] + fn[3]
            dj = fn[0]-2*fn[1]+fn[2]
            djp1 = fn[-1]-2*fn[0]+fn[1]
            MM1, MM2 = gv('MM1'), gv('MM2')
            output_eqns += [OpenSBLIEq(MM1, self.DMM(dj, djp1))]
            output_eqns += [OpenSBLIEq(MM2, self.DMM(djm1, dj))]
            # f_j, f_(j+1), 0.5*(u_j + u_(j+1)) - 0.5*MM(d_j, d_(j+1))
            c_1 = [fn[1], fn[0], 0.5*(fn[1] + fn[0]) - 0.5*MM1]
            # f_j, f_j + alpha*(f_j - f_(j+1)), f_j + 0.5*(f_j - f_(j+1)) + 4/3 * MM(dj, d_(j-1))
            c_2 = [fn[1], fn[1] + alpha*(fn[1] - fn[0]), fn[1] + 0.5*(fn[1] - fn[2]) + Rational(4,3)*MM2]
            # The two lists of arguments
            u_min = Max(Min(c_1[0], Min(c_1[1], c_1[2])), Min(c_2[0], Min(c_2[1], c_2[2])))
            u_max = Min(Max(c_1[0], Max(c_1[1], c_1[2])), Max(c_2[0], Max(c_2[1], c_2[2])))
            # rv.limiter = gv('limiter')
            output_eqns += [OpenSBLIEq(original, self.median(original, u_min, u_max))]
            pass
        else:
            raise ValueError("Input to MP should be left or right biased reconstruction variable.")
        # Final corrected flux
        # for eqn in output_eqns:
        #     pprint(eqn)
        return output_eqns

    def update_constituent_relation_symbols(self, sym, direction):
        """ Function to take the set of required quantities from the constituent relations in symbolic form
        and update the directions in which they are used.

        :arg set sym: Set of required symbols.
        :arg int direction: The axis on which WENO is being applied to (x0, x1 ..)."""
        if isinstance(sym, Symbol):
            sym = [sym]
        elif isinstance(sym, list) or isinstance(sym, set):
            pass
        else:
            raise ValueError("The symbol provided should be either a list or symbols")
        for s in sym:
            if isinstance(s, DataSetBase):
                s = s.noblockname
            if s in self.required_constituent_relations_symbols.keys():
                self.required_constituent_relations_symbols[s] += [direction]
            else:
                self.required_constituent_relations_symbols[s] = [direction]
        return

    def generate_constituent_relations_kernels(self, block):
        """ Generates constituent relation kernels on the block.

        :arg object block: The current block."""
        crs = {}
        for key in self.required_constituent_relations_symbols:
            kernel = Kernel(block, computation_name="CR%s" % key)
            kernel.set_grid_range(block)
            for direction in self.required_constituent_relations_symbols[key]:
                kernel.set_halo_range(direction, 0, self.halotype)
                kernel.set_halo_range(direction, 1, self.halotype)
            crs[block.location_dataset(key)] = kernel
        return crs

    def check_constituent_relations(self, block, list_of_eq, current_constituents):
        """ Checks all the datasets in equations provided are evaluated in constituent relations."""
        arrays = []
        for eq in flatten(list_of_eq):
            arrays += list(eq.atoms(DataSet))
        arrays = set(arrays)
        undefined = arrays.difference(current_constituents.keys())

        for dset in undefined:
            current_constituents[dset] = Kernel(block, computation_name="CR%s" % dset)
            current_constituents[dset].set_grid_range(block)
        return current_constituents

    def create_reconstruction_kernel(self, direction, halos, block):
        """ Creates a kernel for the WENO/TENO reconstruction, currently done per direction."""
        kernel = Kernel(block, computation_name="%s_reconstruction_%d_direction" % (self.__class__.__name__, direction))
        kernel.set_grid_range(block)
        # WENO reconstruction should be evaluated for extra point on each side
        kernel.set_halo_range(direction, 0, halos)
        kernel.set_halo_range(direction, 1, halos)
        block.set_block_boundary_halos(direction, 0, self.halotype)
        block.set_block_boundary_halos(direction, 1, self.halotype)
        return kernel


class EigenSystem(object):
    """ Class to hold the routines required by the characteristic decomposition of the Euler equations. The input
    is the eigensystems used to diagonalise the Jacobian in each direction.

    :arg object physics: Physics object, defaults to NSPhysics."""

    def __init__(self, physics):
        self.physics = physics
        return

    def instantiate_eigensystem(self, block, species):
        if self.physics is None:
            Euler_eq = EulerEquations(block.ndim, species)
            Euler_eq.generate_eig_system(block)
        else:
            self.physics.generate_eig_system(block)
        self.euler = Euler_eq
        self.eigen_value = {}
        self.left_eigen_vector = {}
        self.right_eigen_vector = {}
        return

    def get_symbols_in_ev(self, direction):
        """ Retrieve the unique symbols present in the eigenvalue matrix for a given direction."""
        return self.eigen_value[direction].atoms(EinsteinTerm).difference(self.eigen_value[direction].atoms(ConstantObject))

    def get_symbols_in_LEV(self, direction):
        """ Retrieve the unique symbols present in the left eigenvector matrix for a given direction."""
        return self.left_eigen_vector[direction].atoms(EinsteinTerm).difference(self.left_eigen_vector[direction].atoms(ConstantObject))

    def get_symbols_in_REV(self, direction):
        """ Retrieve the unique symbols present in the right rigenvector matrix for a given direction."""
        return self.right_eigen_vector[direction].atoms(EinsteinTerm).difference(self.right_eigen_vector[direction].atoms(ConstantObject))

    def generate_grid_variable_ev(self, direction, name):
        """ Create a matrix of eigenvalue GridVariable elements. """
        name = '%s_lambda' % (name)
        return self.symbol_matrix(self.eigen_value[direction], name)

    def generate_grid_variable_REV(self, direction, name):
        """ Create a matrix of right eigenvector GridVariable elements. """
        name = '%s_%d_REV' % (name, direction)
        return self.symbol_matrix(self.right_eigen_vector[direction], name)

    def generate_grid_variable_LEV(self, direction, name):
        """ Create a matrix of left eigenvector GridVariable elements. """
        name = '%s_%d_LEV' % (name, direction)
        return self.symbol_matrix(self.left_eigen_vector[direction], name)

    def symbol_matrix(self, mat, name):
        """ Generic function to populate a matrix of GridVariables for a given name and shape."""
        shape = mat.shape
        symbolic_matrix = zeros(*shape)
        for i in range(shape[0]):
            for j in range(shape[1]):
                if mat[i, j]:
                    symbolic_matrix[i, j] = gv('%s_%d%d' % (name, i, j))
        return symbolic_matrix

    def convert_matrix_to_grid_variable(self, mat, name):
        """ Converts the given symbolic matrix to grid variable equivalent.

        :arg Matrix mat: Symbolic matrix to convert to GridVariable elements.
        :arg str name: Base name to use for the GridVariables
        :returns: mat: The matrix updated to contain GridVariables."""
        syms = list(mat.atoms(EinsteinTerm).difference(mat.atoms(ConstantObject)))
        new_syms = [gv('%s_%s' % (name, str(sym))) for sym in syms]
        substitutions = dict(zip(syms, new_syms))
        # Find Metric terms which are datasets
        dsets = list(mat.atoms(DataSet).difference(mat.atoms(ConstantObject)))
        new_dsets = [gv('%s_%s' % (name, d.base.simplelabel())) for d in dsets]
        substitutions.update(dict(zip(dsets, new_dsets)))
        mat = mat.subs(substitutions)
        return mat

    def generate_equations_from_matrices(self, lhs_matrix, rhs_matrix):
        """ Forms a matrix containing Sympy equations at each element, given two input matrices.

        :arg Matrix lhs_matrix: The elements of lhs_matrix form the LHS of the output equations.
        :arg Matrix rhs_matrix: The elements alpha rhs_matrix form the RHS of the output equations.
        returns: Matrix: equations: A matrix containing an Eq(lhs, rhs) pair in each element."""
        if lhs_matrix.shape != rhs_matrix.shape:
            raise ValueError("Matrices should have the same dimension.")
        equations = zeros(*lhs_matrix.shape)
        for no, v in enumerate(lhs_matrix):
            if rhs_matrix[no] != 0:
                equations[no] = OpenSBLIEq(v, factor(rhs_matrix[no]))
        return equations


class Characteristic(EigenSystem):
    """ Class containing the routines required to perform the characteristic decomposition.

        :arg object physics: Physics object, defaults to NSPhysics."""

    def __init__(self, physics, flux_split, shock_filter=None, TVD=False):
        self.flux_split = flux_split
        self.shock_filter = shock_filter
        self.TVD = TVD
        self.temp_wk_arrays = {}
        EigenSystem.__init__(self, physics)
        return

    def get_characteristic_equations(self, direction, derivatives, solution_vector, block):
        """ Performs the three stages required for a characteristic based reconstruction."""
        print("Flux split: {}".format(self.flux_split))
        print("Shock filter: {}".format(self.shock_filter))
        print("TVD: {}".format(self.TVD))
        if 'HLLC' in self.flux_type:
            recon_setting = False
        else:
            recon_setting = True
        settings = {"single_reconstruction_variable": recon_setting, "shock_filter": self.shock_filter}
        for i in range(len(derivatives)):
            derivatives[i].update_settings(**settings)
        pre_process_eqns, reduction_eqns = self.pre_process(direction, derivatives, solution_vector, block)

        if not self.TVD:
            interpolated_eqns = self.interpolate_reconstruction_variables(derivatives, block)
        else:
            interpolated_eqns = []

        if self.flux_split:
            post_process_eqns = self.post_process(direction, derivatives, block)
        else:
            post_process_eqns = self.post_process_Q(direction, derivatives, block)
        return [pre_process_eqns, reduction_eqns, interpolated_eqns, post_process_eqns]

    def remove_zero_equations(self, equations):
        for eqn in equations[:]:
            if isinstance(eqn, Zero) is True:
                equations.remove(eqn)
        return equations

    def direction_flux_vector(self, derivatives, direction):
        flux_vector = []
        for d in derivatives:
            if d.get_direction[0] != direction:
                raise ValueError("Derivatives provided for flux vector are not homogeneous in direction.")
            flux_vector += [d.args[0]]
        return flux_vector

    def convert_symbolic_to_dataset(self, symbolics, location, direction, block):
        """ Converts symbolic terms to DataSets.

        :arg object symbolics: Expression containing Symbols to be converted into DataSets.
        :arg int location: The integer location to apply to the DataSet.
        :arg int direction: The integer direction (axis) of the DataSet to apply the location to.
        returns: object: symbolics: The original expression updated to be in terms of DataSets rather than Symbols."""
        dsets = symbolics.atoms(DataSet).difference(symbolics.atoms(ConstantObject))
        symbols = symbolics.atoms(EinsteinTerm).difference(symbolics.atoms(ConstantObject))
        substitutions = {}
        # Increment flow variables
        for s in symbols:
            dest = block.create_datasetbase(s)
            loc = dest.location
            loc[direction] = loc[direction] + location
            substitutions[s] = dest[loc]
        # Increment metric datasets
        for d in dsets:
            substitutions[d] = increment_dataset(d, direction, location)
        return symbolics.subs(substitutions)

    def set_global_eigenvalues(self, block):
        """ Reduction variables for the global flux-splitting."""
        self.global_eigenvalue_reductions, self.global_eigenvalues = {}, {}
        block_no = block.blocknumber
        # Get the eigenvalues for all directions from the characteristic system
        for dire in range(block.ndim):
            ev_dict, LEV_dict, REV_dict, required_metrics, inv_metric = self.euler.apply_direction(dire)
            self.global_eigenvalues[dire] = ev_dict[dire]
            # name = str(ev_dict[dire][0,0])#
            name = str('u%d' % dire) # hard-coded eigenvalue names for now
            if block.ndim == 1:
                reduction_names = [name+'_minus'+'_max'+'_B%d' % block_no] + [name+'_max'+'_B%d' % block_no] + [name+'_plus'+'_max'+'_B%d' % block_no] 
            else:
                reduction_names = [name+'_max'+'_B%d' % block_no for _ in range(block.ndim)] + [name+'_plus'+'_max'+'_B%d' % block_no] + [name+'_minus'+'_max'+'_B%d' % block_no]
            reduction_vars = [ReductionMax(x) for x in reduction_names]

            symbolic_matrix = zeros(*(block.ndim+2, block.ndim+2))
            for i in range(block.ndim+2):
                for j in range(block.ndim+2):
                    if i == j:
                        symbolic_matrix[i, j] = reduction_vars[i]
            self.global_eigenvalue_reductions[dire] = symbolic_matrix
        return

    def solution_vector_to_characteristic(self, solution_vector, direction, name):
        stencil_points = sorted(list(set(self.reconstruction_classes[0].func_points + self.reconstruction_classes[1].func_points)))
        solution_vector_stencil = zeros(len(solution_vector), len(stencil_points))
        CS_matrix = zeros(len(solution_vector), len(stencil_points))
        for j, val in enumerate(stencil_points):  # j in fv stencil matrix
            for i, flux in enumerate(solution_vector):
                solution_vector_stencil[i, j] = increment_dataset(flux, direction, val)
                CS_matrix[i, j] = gv('CS_%d%d' % (i, j))
        grid_LEV = self.generate_grid_variable_LEV(direction, name)
        characteristic_solution_stencil = grid_LEV*solution_vector_stencil
        characteristic_solution_stencil.stencil_points = stencil_points
        CS_matrix.stencil_points = stencil_points
        return characteristic_solution_stencil, CS_matrix

    def characteristic_setup(self, direction, name, derivatives, block, index=0):
        """ Perform the initial characteristic steps used in both LLF and RF."""
        ev_dict, LEV_dict, REV_dict, required_metrics, inv_metric = self.euler.apply_direction(direction)
        self.eigen_value.update(ev_dict)
        self.left_eigen_vector.update(LEV_dict)
        self.right_eigen_vector.update(REV_dict)
        averaged_suffix_name = 'AVG_%d' % direction
        self.averaged_suffix_name = averaged_suffix_name
        # Finding flow variables to average
        required_symbols = self.get_symbols_in_ev(direction).union(self.get_symbols_in_LEV(direction)).union(self.get_symbols_in_REV(direction))
        required_terms = required_symbols.union(required_metrics)
        # Get the averaged eigenvalue variables
        locations = [index, index+1]
        averaged_equations = self.average(required_terms, direction, averaged_suffix_name, block, locations)
        # Add symbols from the derivatives: e.g. pressure is required
        for d in derivatives:
            required_symbols = required_symbols.union(d.atoms(DataSetBase))
        # Dictionary for reduction operatons
        self.set_global_eigenvalues(block)
        return inv_metric, averaged_equations, required_symbols

    def create_LEV_inverses(self, direction, avg_LEV_values):
        """ Optimizations to avoid repeated divides."""
        inverses = [gv('inv_AVG_a'), gv('inv_AVG_rho')]
        to_be_replaced = [1/gv('AVG_%d_a' % direction), 1/gv('AVG_%d_rho' % direction)]
        if len(self.inv_metric.atoms(gv)) > 0:  # Don't substitute if there are no metrics
            inverses += [gv('inv_AVG_met_fact')]
            to_be_replaced += [1/self.inv_metric]
        inverse_evals = [OpenSBLIEq(a, b) for (a, b) in zip(inverses, to_be_replaced)]
        for i, term in enumerate(avg_LEV_values):
            for old, new in zip(to_be_replaced, inverses):
                term = term.subs({old: new})
            avg_LEV_values[i] = term
        return inverse_evals, avg_LEV_values

    def flux_vector_to_characteristic(self, derivatives, direction, name):
        fv = self.direction_flux_vector(derivatives, direction)
        stencil_points = sorted(list(set(self.reconstruction_classes[0].func_points + self.reconstruction_classes[1].func_points)))
        flux_stencil = zeros(len(fv), len(stencil_points))
        CF_matrix = zeros(len(fv), len(stencil_points))
        for j, val in enumerate(stencil_points):  # j in fv stencil matrix
            for i, flux in enumerate(fv):
                flux_stencil[i, j] = increment_dataset(flux, direction, val)
                CF_matrix[i, j] = gv('CF_%d%d' % (i, j))
        grid_LEV = self.generate_grid_variable_LEV(direction, name)
        characteristic_flux_stencil = grid_LEV*flux_stencil
        characteristic_flux_stencil.stencil_points = stencil_points
        CF_matrix.stencil_points = stencil_points
        return characteristic_flux_stencil, CF_matrix

    def replace_gamma_factor(self, eqns):
        ConstantsToDeclare.add_constant(ConstantObject('gama'))
        gamma_minus_one = ConstantObject('gamma_m1')
        gamma_minus_one.value = (ConstantObject('gama') - 1)
        ConstantsToDeclare.add_constant(gamma_minus_one)
        for i, eqn in enumerate(eqns):
            eqns[i] = eqn.subs({gamma_minus_one.value: gamma_minus_one})
        return eqns

    def create_REV_inverses(self, direction):
        averaged_suffix_name = self.averaged_suffix_name
        avg_REV_values = self.convert_matrix_to_grid_variable(self.right_eigen_vector[self.direction], averaged_suffix_name)
        # Manually remove the divides
        inverses = [gv('inv_AVG_a'), gv('inv_AVG_rho')]
        to_be_replaced = [1/gv('AVG_%d_a' % direction), 1/gv('AVG_%d_rho' % direction)]
        if len(self.inv_metric.atoms(gv)) > 0:  # Don't substitute if there are no metrics
            inverses += [gv('inv_AVG_met_fact')]
            to_be_replaced += [1/self.inv_metric]

        for i, term in enumerate(avg_REV_values):
            for old, new in zip(to_be_replaced, inverses):
                term = term.replace(old, new)
            avg_REV_values[i] = term
        return avg_REV_values


    def create_output_wk_arrays(self, direction, derivatives, block):
        # Re use arrays to reduce memory usage if WENO is being applied as a shock filter step
        if block.shock_filter:
            if direction == 0:
                reconstructed_work = [d.reconstruction_work for d in derivatives]
            elif direction == 1:
                reconstructed_work = [block.location_dataset('Residual%d' % i) for i in range(len(derivatives))]
            else:
                reconstructed_work = [block.location_dataset('%s_RKold' % fn.base.label) for fn in self.input_solution_vector]
            # Update the work arrays
            for i, d in enumerate(derivatives):
                d.reconstruction_work = reconstructed_work[i]

                if direction in self.temp_wk_arrays.keys():
                    self.temp_wk_arrays[direction].append(reconstructed_work[i])
                else:
                    self.temp_wk_arrays[direction] = [reconstructed_work[i]]
        else:
            reconstructed_work = [d.reconstruction_work for d in derivatives]
        return reconstructed_work

    def central_diff_formula(self, component, reconstruction_variable, flux_type, derivative=None):
        """ Central difference formula based on the f_i = 0.5*(f_(i+1/2) - f_(i-1/2)) half-node locations. Applied in characteristic space
        using the CF local variables."""
        if (self.order+1) == 4:
            weights = [Rational(-1,12), Rational(7,12), Rational(7,12), Rational(-1,12)]
            locations = [-1, 0, 1, 2]
        elif (self.order+1) == 6:
            weights = [Rational(1,60), Rational(-8,60), Rational(37,60), Rational(37,60), Rational(-8,60), Rational(1,60)]
            locations = [-2, -1, 0, 1, 2, 3]
        elif (self.order+1) == 8:
            weights = [Rational(-1,280), Rational(29, 840), Rational(-139,840), Rational(533, 840), Rational(533,840), Rational(-139,840), Rational(29,840), Rational(-1, 280)]
            locations = [-3, -2, -1, 0, 1, 2, 3, 4]
        else:
            raise NotImplementedError("Only 4th, 6th, and 8th order Central are implemented for the WENO-Filter.")
        # Take a central difference of the characteristic fluxes
        terms = []
        if flux_type == 'LLF' or flux_type == 'GLF':
            # Derivatives of the transformed fluxes
            terms += [gv('CF_%d%d' % (component, j)) for j in range(len(weights))]
        else:
            # HLLC solver, build from the flux vector directly
            for location in locations:
                dsets = [x for x in [derivative.args[0]]]
                combined = sum([increment_dataset(dsets[i], derivative.args[1].direction, location) for i in range(len(dsets))])
                terms.append(combined)
        formula = factor(sum([x*y for (x,y) in zip(weights, terms)]))
        output_equation = [OpenSBLIEq(reconstruction_variable, gv('rj%d' % component)*(reconstruction_variable - formula))]
        return output_equation

    def define_substitution_dictionaries(self, direction, derivatives, block, CS, velocities, pressure):
        """ Dictionaries to link the datasets present in the flux vector and the left/right evaluated states after the WENO reconstruction."""
        # Conservative variables
        q_dsets = self.input_solution_vector
        left, right = [], []
        left += [(q_dsets[0], CS[0][0])]
        left += [(q_dsets[i+1], CS[0][i+1]) for i in range(block.ndim)]
        left += [(q_dsets[-1], CS[0][-1])]
        right += [(q_dsets[0], CS[1][0])]
        right += [(q_dsets[i+1], CS[1][i+1]) for i in range(block.ndim)]
        right += [(q_dsets[-1], CS[1][-1])]
        # Pressure
        left += [(block.location_dataset('p'), pressure[0])]
        right += [(block.location_dataset('p'), pressure[1])]
        # Velocities
        left += [(block.location_dataset('u%d' % i), velocities[0][i]) for i in range(block.ndim)]
        right += [(block.location_dataset('u%d' % i), velocities[1][i]) for i in range(block.ndim)]
        SD_L = dict(left)
        SD_R = dict(right)
        return SD_L, SD_R

class TVDCharacteristic(Characteristic):
    """ Implements the TVD equations

    :arg object physics: Physics object, defaults to NSPhysics.
    :arg object averaging: The averaging procedure to be applied for characteristics, defaults to Simple averaging."""

    def __init__(self, physics, averaging=None, flux_split=False):
        Characteristic.__init__(self, physics, flux_split)
        if averaging is None:
            self.average = SimpleAverage([0, 1]).average
        else:
            self.average = averaging.average
        self.flux_split = flux_split
        return

    def pre_process(self, direction, derivatives, solution_vector, block):
        """ TVD implementation over j-1/2, j+1/2, j+3/2, Max Walker & David J. Lusher (09/2023)."""
        self.direction = direction
        self.input_solution_vector = solution_vector
        # Number of variables to solve
        nvars = len(solution_vector)
        # Output equations
        pre_process_equations, reduction_equations = [], []
        # Update the ev, LEV and REV dicts and perform averaging
        avg_name = 'AVG_%d' % direction
        # Loop over 3 indices to get the eigensystem at j-1/2, j+1/2, j+3/2
        offsets = [-1, 0, 1]
        # Create the alpha variables
        alphas = Matrix([[GridVariable('alpha_%d%d' % (i, j)) for i in range(nvars)] for j in range(3)])
        for recon_index, jbase in enumerate(offsets):
            # Step 1: obtain the averaged eigenstate
            inv_metric, averaged_equations, required_symbols = self.characteristic_setup(direction, avg_name, derivatives, block, index=jbase)
            pre_process_equations += averaged_equations
            # Update required CRs
            if recon_index == 0:
                self.update_constituent_relation_symbols(required_symbols, direction)
            # Step 2: Store the R^-1 matrix components (LEV) based on this averaged state
            # # Inverse metric term: WARNING: Metrics need to be looked at for curvilinear cases
            self.inv_metric = self.convert_matrix_to_grid_variable(inv_metric, avg_name)
            # # Eigensystem based on averaged quantities
            avg_LEV_values = self.convert_matrix_to_grid_variable(self.left_eigen_vector[direction], avg_name)
            # # Manually replace re-used divides by inverses
            inverse_evals, avg_LEV_values = self.create_LEV_inverses(direction, avg_LEV_values)
            pre_process_equations += inverse_evals
            # # Grid variables to store averaged eigensystems
            grid_LEV_matrix = self.generate_grid_variable_LEV(direction, avg_name)
            pre_process_equations += flatten(self.generate_equations_from_matrices(grid_LEV_matrix, avg_LEV_values))
            # Step 3: Transform the du terms, alpha = R^-1 *(u_{i+1} - u_i)
            du = Matrix([increment_dataset(v, direction, jbase+1) - increment_dataset(v, direction, jbase) for v in solution_vector])
            alpha_vector = grid_LEV_matrix*du
            # Store the alphas for later use
            pre_process_equations += [OpenSBLIEq(alphas[recon_index, a], alpha_vector[a]) for a in range(len(alpha_vector))]
            # Step 5: If j+1/2 location, store the wavespeeds ws_0, ws_1, ws_2, and REV for later use
            if jbase == 0:
                # Step 5a: Wavespeed calculation
                ev = self.eigen_value[direction]
                out = zeros(*ev.shape)
                # Replace the eigenvalues with their averaged state instead
                subs_dict = {}
                for eqn in averaged_equations:
                    dset_base = list(eqn.rhs.atoms(DataSetBase))[0]
                    subs_dict[EinsteinTerm(str(dset_base.noblockname))] = eqn.lhs
                for i in range(ev.shape[0]):
                    for j in range(ev.shape[1]):
                        ev[i,j] = ev[i,j].subs(subs_dict)
                wavespeeds = self.generate_grid_variable_ev(direction, 'ws')
                ev_equations = self.generate_equations_from_matrices(wavespeeds, ev)
                ev_equations = [x for x in ev_equations if x != 0]
                pre_process_equations += ev_equations
                # Step 5b: REV calculation
                avg_REV_values = self.create_REV_inverses(direction)
                self.grid_REV = self.generate_grid_variable_REV(direction, 'AVG')
                REV_equations = [OpenSBLIEq(x, y, evaluate=False) for (x,y) in zip(self.grid_REV, avg_REV_values)]
                # Remove any 0=0 evaluations
                REV_equations = [x for x in REV_equations if x.lhs != 0 and x.rhs != 0]
                pre_process_equations += REV_equations

        # Step 6: Calculate sigma values from Yee et al (1999)
        sigmas_0 = [GridVariable('sigma_0_%d' % i) for i in range(nvars)]
        sigmas_1 = [GridVariable('sigma_1_%d' % i) for i in range(nvars)]
        delta = ConstantObject('delta_TVD')
        delta.value = 0.5
        ConstantsToDeclare.add_constant(delta)
        for i, sigma in enumerate(sigmas_0):
            if_expr, else_expr = Abs(wavespeeds[i,i]), (wavespeeds[i,i]**2 + delta**2) / (2*delta)
            pre_process_equations += [OpenSBLIEq(sigmas_0[i], Piecewise((if_expr, Abs(wavespeeds[i,i]) >= delta), (else_expr, True)))]
        # Step 7: Calculate the 'g' functions, which are the upwind limiter functions
        S = GridVariable('S')
        g_equations = []
        # g terms at j and j+1
        g_terms = [[GridVariable('g_%d%d' % (i,j)) for i in range(nvars)] for j in range(2)]
        # Temporary terms
        t1, t2 = GridVariable('t1'), GridVariable('t2')
        for j in [1, 2]:
            for i in range(nvars):
                g_equations += [OpenSBLIEq(S, sign(alphas[j,i]))]
                a1, a2 = alphas[j,i], alphas[j-1,i]
                g_equations += [OpenSBLIEq(t1, Max(0.0, Min(2*Abs(a1), S*a2))), OpenSBLIEq(t2, Min(Abs(a1), 2*S*a2))]
                g_eval = S*Max(t1, t2)
                g_equations += [OpenSBLIEq(g_terms[j-1][i], g_eval)]
        pre_process_equations += g_equations
        # Step 8: Calculate the phi_{j+1/2} terms
        gamma_terms = [GridVariable('gamma_%d' % i) for i in range(nvars)]
        eps = ConstantObject('eps_TVD')
        eps.value = 0.00000001
        ConstantsToDeclare.add_constant(eps)
        for i, gamma in enumerate(gamma_terms):
            if_expr, else_expr = 0, sigmas_0[i]*alphas[1,i]*(g_terms[1][i] - g_terms[0][i]) / (alphas[1,i]**2 + eps)
            pre_process_equations += [OpenSBLIEq(gamma_terms[i], Piecewise((if_expr, Equality(alphas[1,i], 0)), (else_expr, True)))]
        # new step: sigma of gamma terms
        for i, sigma in enumerate(sigmas_1):
            if_expr, else_expr = Abs(gamma_terms[i]), (gamma_terms[i]**2 + delta**2) / (2*delta)
            pre_process_equations += [OpenSBLIEq(sigmas_1[i], Piecewise((if_expr, Abs(gamma_terms[i]) >= delta), (else_expr, True)))]
        # Step 9: Calculate the phi terms
        phi_terms = [GridVariable('phi_%d' % i) for i in range(nvars)]
        # Lambda needs to be one per direction
        TVD_lambda = ConstantObject('lambda%d_TVD' % direction)
        for i, phi in enumerate(phi_terms):
            rhs = 0.5*sigmas_0[i]*(g_terms[1][i] + g_terms[0][i]) - Abs(sigmas_1[i]+sigmas_0[i])*alphas[1,i]
            pre_process_equations += [OpenSBLIEq(phi, rhs)]
        # Step 10: Calculate switch to control the dissipation of the limiter (Harten switch, (eqns 2.22 and 2.23 in Yee et al. 1999))
        theta_hat_terms = [[GridVariable('theta_hat_%d%d' % (i,j)) for i in range(nvars)] for j in range(2)]
        for j in [1, 2]:
            for i in range(nvars):
                a1, a2 = alphas[j,i], alphas[j-1,i]
                pre_process_equations += [OpenSBLIEq(t1, Abs(a1) - Abs(a2))]
                pre_process_equations += [OpenSBLIEq(t2, Abs(a1) + Abs(a2) + eps)]
                rhs = Abs(t1/t2)**2
                pre_process_equations += [OpenSBLIEq(theta_hat_terms[j-1][i], rhs)]
        # Step 11: Calculate theta terms at j+1/2 interface
        theta_terms = [GridVariable('theta_%d' % i) for i in range(nvars)]
        pre_process_equations += [OpenSBLIEq(theta_terms[i], Max(theta_hat_terms[0][i], theta_hat_terms[1][i])) for i in range(nvars)]
        # Step 12: Calculate phi star terms
        self.phi_star_terms = [GridVariable('phi_star_%d' % i) for i in range(nvars)]
        kappa = ConstantObject('kappa_TVD')
        ConstantsToDeclare.add_constant(kappa)
        pre_process_equations += [OpenSBLIEq(self.phi_star_terms[i], kappa*theta_terms[i]*phi_terms[i]) for i in range(nvars)]
        # NOTE: Setting the work arrays is performed in the post_process function to be consistent with WENO
        # Remove '0' entries and gamma - 1 factors from pre_process_equations
        pre_process_equations = self.remove_zero_equations(pre_process_equations)
        pre_process_equations = self.replace_gamma_factor(pre_process_equations)
        # for eqn in pre_process_equations:
        #     pprint(eqn)
        # # exit()
        return pre_process_equations, reduction_equations

    def post_process_Q(self, dire, derivatives, block):
        """ Transforms the characteristic WENO interpolated fluxes back into real space by multiplying by the right
        eigenvector matrix.

        :arg list derivatives: The derivatives to perform the characteristic decomposition and WENO on.
        :arg object kernel: The current computational kernel."""
        pp_equations = []
        ndim = block.ndim
        # Create output arrays to store the final flux reconstructions
        reconstructed_work = self.create_output_wk_arrays(dire, derivatives, block)
        # Transformation matrix back to physical (non-characteristic) space
        # Step 13: Create the numerical filter flux: (eq. 2.2 in Yee et al. 1999)
        phi_star_terms = Matrix(self.phi_star_terms)
        reconstructions = factor(0.5*self.grid_REV*phi_star_terms)
        # Assign to work arrays
        pp_equations += [OpenSBLIEq(wk, reconstructions[i]) for i, wk in enumerate(reconstructed_work)]
        return pp_equations

class LFCharacteristic(Characteristic):
    """ This class contains the base Local Lax-Fedrich scheme performed in characteristic space.

    :arg object physics: Physics object, defaults to NSPhysics.
    :arg object averaging: The averaging procedure to be applied for characteristics, defaults to Simple averaging."""

    def __init__(self, physics, flux_split, flux_type=None, averaging=None, shock_filter=None):
        Characteristic.__init__(self, physics, flux_split, shock_filter)
        if averaging is None:
            self.average = SimpleAverage([0, 1]).average
        else:
            self.average = averaging.average
        if flux_type is None:
            flux_type = 'LLF' # default to local Lax Friedrichs
        self.flux_type = flux_type
        self.flux_split = flux_split
        return

    def pre_process(self, direction, derivatives, solution_vector, block):
        """ Performs the transformation of the derivatives into characteristic space using the eigensystems provided to Characteristic. Flux splitting is then applied
        to the characteristic variables in preparation for the WENO interpolation. Required quantities are added to pre_process_equations.

        :arg int direction: Integer direction to apply the characteristic decomposition and WENO (x0, x1, ...).
        :arg list derivatives: The derivatives to perform the characteristic decomposition and WENO on.
        :arg list solution_vector: Solution vector from the Euler equations (rho, rhou0, rhou1, rhou2, rhoE) in vector form."""
        self.direction = direction
        self.input_solution_vector = solution_vector
        # Output equations
        pre_process_equations = []
        # Update the ev, LEV and REV dicts and perform averaging
        avg_name = 'AVG_%d' % direction
        inv_metric, averaged_equations, required_symbols = self.characteristic_setup(direction, avg_name, derivatives, block)
        pre_process_equations += averaged_equations
        # Update required CRs
        self.update_constituent_relation_symbols(required_symbols, direction)

        # Inverse metric term
        self.inv_metric = self.convert_matrix_to_grid_variable(inv_metric, avg_name)
        # Eigensystem based on averaged quantities
        avg_LEV_values = self.convert_matrix_to_grid_variable(self.left_eigen_vector[direction], avg_name)
        # Manually replace re-used divides by inverses
        inverse_evals, avg_LEV_values = self.create_LEV_inverses(direction, avg_LEV_values)
        pre_process_equations += inverse_evals

        # Grid variables to store averaged eigensystems
        grid_LEV = self.generate_grid_variable_LEV(direction, avg_name)
        pre_process_equations += flatten(self.generate_equations_from_matrices(grid_LEV, avg_LEV_values))

        # Create characteristic matrices and add their evaluations to pre_process
        # Check whether the time advanced quantities are in conservative form or not
        if not self.conservative:
            rho = solution_vector[0]
            assert str(rho.base.simplelabel()) == 'rho'
            solution_vector = [rho] + [rho*x for x in solution_vector[1:]]
        evaluations, CS_matrix, CF_matrix = self.create_characteristic_matrices(direction, derivatives, solution_vector, avg_name)
        pre_process_equations += evaluations
        # Get max wavespeeds and their evaluations, eigenvalues evaluated either local or globally
        if self.flux_type == 'LLF':
            self.grid_EV, pre_process_equations = self.create_max_characteristic_wave_speed(pre_process_equations, direction, block)
            reduction_equations = []
        else:
            self.grid_EV, reduction_equations, pre_process_equations = self.calculate_eigenvalue_reductions(pre_process_equations, direction, block)
        # Transform both the flux vector and the solution vector to characteristic space, flux splitting directly into the WENO/TENO procedure
        self.characteristic_flux_splitting(self.grid_EV, CS_matrix, CF_matrix, derivatives)

        # Remove '0' entries and gamma - 1 factors from pre_process_equations
        pre_process_equations = self.remove_zero_equations(pre_process_equations)
        pre_process_equations = self.replace_gamma_factor(pre_process_equations)
        return pre_process_equations, reduction_equations

    def calculate_eigenvalue_reductions(self, pre_process_equations, direction, block):
        """ Performs a reduction of the eigenvalues over the entire domain."""
        reductions = []
        ndim = block.ndim
        # Create a reduction kernel to compute the eigenvalues globally over the domain for all directions
        if direction == 0:
            for dire in range(block.ndim):
                global_EV_reductions = self.global_eigenvalue_reductions[dire]
                global_EVs = self.global_eigenvalues[dire]
                # u, u+a, u-a eigenvalue ordering assumed
                u = self.convert_symbolic_to_dataset(global_EVs[0,0], 0, 0, block)
                upa = self.convert_symbolic_to_dataset(global_EVs[ndim,ndim], 0, 0, block)
                uma = self.convert_symbolic_to_dataset(global_EVs[ndim+1,ndim+1], 0, 0, block)
                reductions += [OpenSBLIEq(global_EV_reductions[0,0], Abs(u))]
                reductions += [OpenSBLIEq(global_EV_reductions[ndim,ndim], Abs(upa))]
                reductions += [OpenSBLIEq(global_EV_reductions[ndim+1,ndim+1], Abs(uma))]
        # Assign the max wave speed to the correct reduced variables for this direction
        grid_vars, reduction_vars = self.generate_grid_variable_ev(direction, 'max'), self.global_eigenvalue_reductions[direction]
        pre_process_equations += [x for x in self.generate_equations_from_matrices(grid_vars, reduction_vars) if x != 0]
        return grid_vars, reductions, pre_process_equations

    def post_process(self, dire, derivatives, block):
        """ Transforms the characteristic WENO interpolated fluxes back into real space by multiplying by the right
        eigenvector matrix.

        :arg list derivatives: The derivatives to perform the characteristic decomposition and WENO on.
        :arg object kernel: The current computational kernel."""
        pp_equations = []
        ndim = block.ndim
        # Create output arrays to store the final flux reconstructions
        reconstructed_work = self.create_output_wk_arrays(dire, derivatives, block)
        # Transformation matrix back to physical (non-characteristic) space
        avg_REV_values = self.create_REV_inverses(dire)
        reconstructed_characteristics = Matrix([d.evaluate_reconstruction for d in derivatives])
        # Single reconstruction variable?
        if derivatives[0].settings["single_reconstruction_variable"]:
            # Apply a shock sensor if the WENO is being applied as a filter step
            if block.shock_filter:
                for i, recon in enumerate(reconstructed_characteristics):
                    pp_equations += flatten([self.central_diff_formula(i, recon, self.flux_type)])
            reconstructed_flux = avg_REV_values*reconstructed_characteristics
            pp_equations += [OpenSBLIEq(x, y) for x, y in zip(reconstructed_work, reconstructed_flux)]
        else:
            print("Using multiple reconstruction variables in WENO/TENO post-process.")
            left_F, right_F = symbols("fL:%d" % (ndim+2), **{'cls':GridVariable}), symbols("fR:%d" % (ndim+2), **{'cls':GridVariable})
            wL, wR = Matrix(reconstructed_characteristics[:,1]), Matrix(reconstructed_characteristics[:,0])
            pp_equations += [OpenSBLIEq(x, y) for x, y in zip(right_F, avg_REV_values*wR)]
            pp_equations += [OpenSBLIEq(x, y) for x, y in zip(left_F, avg_REV_values*wL)]

            for i, component in enumerate(reconstructed_work):
                pp_equations += [OpenSBLIEq(component, left_F[i] + right_F[i])]

        pp_equations = self.replace_gamma_factor(pp_equations)
        return pp_equations

    def post_process_Q(self, dire, derivatives, block):
        """ Transforms the characteristic WENO interpolated fluxes back into real space by multiplying by the right
        eigenvector matrix.

        :arg list derivatives: The derivatives to perform the characteristic decomposition and WENO on.
        :arg object kernel: The current computational kernel."""
        pp_equations = []
        ndim = block.ndim
        # Create output arrays to store the final flux reconstructions
        reconstructed_work = self.create_output_wk_arrays(dire, derivatives, block)

        reconstructed_characteristics = Matrix([d.evaluate_reconstruction for d in derivatives])
        # Transformation matrix back to physical (non-characteristic) space
        avg_REV_values = self.create_REV_inverses(dire)
        # LLF: Reconstruct the left and right states separately and then transform each back to physical space
        left_q, right_q = symbols("qL:%d" % (ndim+2), **{'cls':GridVariable}), symbols("qR:%d" % (ndim+2), **{'cls':GridVariable})
        wL, wR = Matrix(reconstructed_characteristics[:,1]), Matrix(reconstructed_characteristics[:,0])
        pp_equations += [OpenSBLIEq(x, y) for x, y in zip(right_q, avg_REV_values*wR)]
        pp_equations += [OpenSBLIEq(x, y) for x, y in zip(left_q, avg_REV_values*wL)]
        # LLF: Calculate the speed of sound and pressure at each interface using the reconstructed values
        pL, pR = symbols("pL", **{'cls':GridVariable}), symbols("pR", **{'cls':GridVariable})
        rhoL, rhoR = left_q[0], right_q[0]
        rhoEL, rhoER = left_q[-1], right_q[-1]
        vel_L, vel_R = symbols("uL:%d" % (ndim), **{'cls':GridVariable}), symbols("uR:%d" % (ndim), **{'cls':GridVariable})
        # Density inversion factor
        inv_rhoL, inv_rhoR = symbols("inv_rhoL", **{'cls':GridVariable}), symbols("inv_rhoR", **{'cls':GridVariable})
        pp_equations += [OpenSBLIEq(inv_rhoL, 1.0/rhoL), OpenSBLIEq(inv_rhoR, 1.0/rhoR)]
        # Velocity component in the dire of reconstruction
        uL, uR = vel_L[dire], vel_R[dire]
        for i, u in enumerate(vel_L):
            pp_equations += [OpenSBLIEq(vel_L[i], left_q[i+1]*inv_rhoL)]
            pp_equations += [OpenSBLIEq(vel_R[i], right_q[i+1]*inv_rhoR)]
        # WARNING: ideal gas law assumed
        gama = ConstantObject('gama')
        pp_equations += [OpenSBLIEq(pL, (gama- 1)*(left_q[-1] - 0.5*rhoL*sum([x**2 for x in vel_L])))]
        pp_equations += [OpenSBLIEq(pR, (gama- 1)*(right_q[-1] - 0.5*rhoR*sum([x**2 for x in vel_R])))]

        # Build the system of flux components
        # Create a substitution dictionary from the flux components
        SD_L, SD_R = self.define_substitution_dictionaries(dire, derivatives, block, [left_q, right_q], [vel_L, vel_R], [pL, pR])
        # Vectors for the right/left states
        U_L, U_R = Matrix([0 for _ in range(ndim+2)]), Matrix([0 for _ in range(ndim+2)])
        F_L, F_R = Matrix([0 for _ in range(ndim+2)]), Matrix([0 for _ in range(ndim+2)])
        # Substitute left/right states into the input Q and flux vector
        for i, d in enumerate(derivatives):
            input_args = d.args[0]
            U_L[i] = self.input_solution_vector[i].subs(SD_L)
            U_R[i] = self.input_solution_vector[i].subs(SD_R)
            F_L[i] = input_args.subs(SD_L)
            F_R[i] = input_args.subs(SD_R)
        # Flux split with wave-speed selection
        flux_split = factor(Rational(1,2)*(F_L+F_R - self.grid_EV*(U_R - U_L)))

        # Assign the fluxes to the storage arrays
        for i, component in enumerate(reconstructed_work):
            pp_equations += [OpenSBLIEq(component, flux_split[i])]
        # Apply a shock sensor if the WENO is being applied as a filter step
        if block.shock_filter:
            for i in range(len(reconstructed_work)):
                pp_equations += flatten([self.central_diff_formula(i, reconstructed_work[i], self.flux_type, derivatives[i])])
        # Replace gamma factors if required
        pp_equations = self.replace_gamma_factor(pp_equations)
        return pp_equations

    def create_characteristic_matrices(self, direction, derivatives, solution_vector, name):
        # Reconstruct the combined flux or only the Q vector
        if self.flux_split:
            characteristic_flux_vector, CF_matrix = self.flux_vector_to_characteristic(derivatives, direction, name)
            CF_evaluations = flatten(self.generate_equations_from_matrices(CF_matrix, characteristic_flux_vector))
        else:
            CF_matrix = None

        characteristic_solution_vector, CS_matrix = self.solution_vector_to_characteristic(solution_vector, direction, name)
        CS_evaluations = flatten(self.generate_equations_from_matrices(CS_matrix, characteristic_solution_vector))
        # Optimise by grouping evaluations by stencil location
        n_rows, n_cols = CS_matrix.shape[0], CS_matrix.shape[1]
        reordered = []
        for j in range(n_cols):
            reordered_CS, reordered_CF = [], []
            if self.flux_split:
                for i in range(n_rows):
                    reordered.append(CF_evaluations[j+i*n_cols])
            for i in range(n_rows):
                reordered.append(CS_evaluations[j+i*n_cols])
        return reordered, CS_matrix, CF_matrix

    def characteristic_flux_splitting(self, ev_matrix, CS_matrix, CF_matrix, derivatives):
        if self.flux_split:
            positive = factor(Rational(1, 2)*(CF_matrix + ev_matrix*CS_matrix))
            negative = factor(Rational(1, 2)*(CF_matrix - ev_matrix*CS_matrix))
        else:
            positive = CS_matrix
            negative = CS_matrix

        positive_flux = zeros(*positive.shape)
        negative_flux = zeros(*negative.shape)
        # Assign the values
        for i in range(positive_flux.shape[0]):
            for j in range(positive_flux.shape[1]):
                positive_flux[i, j] = factor(positive[i, j])
                negative_flux[i, j] = factor(negative[i, j])
        positive_flux.stencil_points = CS_matrix.stencil_points
        negative_flux.stencil_points = CS_matrix.stencil_points
        self.generate_right_reconstruction_variables(positive_flux, derivatives)
        self.generate_left_reconstruction_variables(negative_flux, derivatives)
        return

    def create_max_characteristic_wave_speed(self, pre_process_equations, direction, block):
        """ Creates the equations for local Lax-Friedrich wave speeds, maximum eigenvalues over the local
        WENO/TENO stencils are found."""
        stencil_points = sorted(list(set(self.reconstruction_classes[0].func_points + self.reconstruction_classes[1].func_points)))
        ev = self.eigen_value[direction]
        out = zeros(*ev.shape)
        stencil_points = [0,1]
        for p in stencil_points:
            location_ev = self.convert_symbolic_to_dataset(ev, p, direction, block)
            for no, val in enumerate(location_ev):
                out[no] = Max(out[no], Abs(val))
        max_wave_speed = self.generate_grid_variable_ev(direction, 'max')
        ev_equations = self.generate_equations_from_matrices(max_wave_speed, out)
        # pprint(ev_equations)
        FC = ConstantObject('shock_filter_control')
        FC.value = 1.0 # Default condition has no scaling
        ConstantsToDeclare.add_constant(FC)
        ev_equations = [x for x in ev_equations if x != 0]
        ev_lhs = [x.lhs for x in ev_equations]
        ev_rhs = [FC*x.rhs for x in ev_equations]
        # If there are repeated eigenvalues we don't compute them multiple times
        for no, eqn in enumerate(ev_rhs):
            if no > 0:
                if eqn == ev_rhs[0]:
                    ev_rhs[no] = ev_lhs[0]  # u, u, u, u+a, u-a ordering
        ev_equations = [OpenSBLIEq(left, right) for (left, right) in zip(ev_lhs, ev_rhs)]
        pre_process_equations += ev_equations
        max_wave_speed = diag(*([max_wave_speed[i, i] for i in range(ev.shape[0])]))
        return max_wave_speed, pre_process_equations


class HLLCCharacteristic(Characteristic):
    """ This class contains the base Local Lax-Fedrich scheme performed in characteristic space.

    :arg object physics: Physics object, defaults to NSPhysics.
    :arg object averaging: The averaging procedure to be applied for characteristics, defaults to Simple averaging."""

    def __init__(self, physics, flux_type='HLLC-LM', shock_filter=False, averaging=None):
        flux_split = False # HLLC has flux splitting disabled
        Characteristic.__init__(self, physics, flux_split, shock_filter)
        if averaging is None:
            self.average = RoeAverage([0, 1]).average
        else:
            self.average = averaging.average
        # if-else version of the HLLC flux construction
        self.conditional = True
        # self.flux_split = True
        self.flux_type = flux_type
        return

    def pre_process(self, direction, derivatives, solution_vector, block):
        """ Performs the transformation of the derivatives into characteristic space using the eigensystems provided to Characteristic. Flux splitting is then applied
        to the characteristic variables in preparation for the WENO interpolation. Required quantities are added to pre_process_equations.

        :arg int direction: Integer direction to apply the characteristic decomposition and WENO (x0, x1, ...).
        :arg list derivatives: The derivatives to perform the characteristic decomposition and WENO on.
        :arg list solution_vector: Solution vector from the Euler equations (rho, rhou0, rhou1, rhou2, rhoE) in vector form."""
        self.direction = direction
        self.input_solution_vector = solution_vector
        pre_process_equations = []
        # Update the ev, LEV and REV dicts and perform averaging
        avg_name = 'AVG_%d' % direction
        inv_metric, averaged_equations, required_symbols = self.characteristic_setup(direction, avg_name, derivatives, block)
        pre_process_equations += averaged_equations
        # Update required CRs
        self.update_constituent_relation_symbols(required_symbols, direction)

        # Inverse metric term
        self.inv_metric = self.convert_matrix_to_grid_variable(inv_metric, avg_name)
        # Eigensystem based on averaged quantities
        avg_LEV_values = self.convert_matrix_to_grid_variable(self.left_eigen_vector[direction], avg_name)
        # Manually replace re-used divides by inverses
        inverse_evals, avg_LEV_values = self.create_LEV_inverses(direction, avg_LEV_values)
        pre_process_equations += inverse_evals

        # Grid variables to store averaged eigensystems
        grid_LEV = self.generate_grid_variable_LEV(direction, avg_name)
        pre_process_equations += flatten(self.generate_equations_from_matrices(grid_LEV, avg_LEV_values))

        # Create characteristic matrices and add their evaluations to pre_process
        # Check whether the time advanced quantities are in conservative form or not
        if not self.conservative:
            rho = solution_vector[0]
            assert str(rho.base.simplelabel()) == 'rho'
            solution_vector = [rho] + [rho*x for x in solution_vector[1:]]
        evaluations, CS_matrix = self.create_characteristic_matrices(direction, derivatives, solution_vector, avg_name)

        pre_process_equations += evaluations
        self.generate_reconstruction_variables(CS_matrix, derivatives)
        # Remove '0' entries and gamma - 1 factors from pre_process_equations
        pre_process_equations = self.remove_zero_equations(pre_process_equations)
        pre_process_equations = self.replace_gamma_factor(pre_process_equations)
        reduction_equations = []
        return pre_process_equations, reduction_equations


    def create_characteristic_matrices(self, direction, derivatives, solution_vector, name):
        characteristic_solution_vector, CS_matrix = self.solution_vector_to_characteristic(solution_vector, direction, name)
        CS_evaluations = flatten(self.generate_equations_from_matrices(CS_matrix, characteristic_solution_vector))
        # Optimise by grouping evaluations by stencil location
        n_rows, n_cols = CS_matrix.shape[0], CS_matrix.shape[1]
        reordered = []
        for j in range(n_cols):
            reordered_CS = []
            for i in range(n_rows):
                reordered.append(CS_evaluations[j+i*n_cols])
        return reordered, CS_matrix

    def generate_reconstruction_variables(self, CS_matrix, derivatives):
        positive = CS_matrix
        positive_flux = zeros(*positive.shape)
        negative = CS_matrix
        negative_flux = zeros(*negative.shape)
        for i in range(positive_flux.shape[0]):
            for j in range(positive_flux.shape[1]):
                positive_flux[i, j] = factor(positive[i, j])
                negative_flux[i, j] = factor(negative[i, j])
        positive_flux.stencil_points = CS_matrix.stencil_points
        negative_flux.stencil_points = CS_matrix.stencil_points
        self.generate_right_reconstruction_variables(positive_flux, derivatives)
        self.generate_left_reconstruction_variables(negative_flux, derivatives)
        return

    def post_process_Q(self, dire, derivatives, block):
        """ Transforms the characteristic WENO interpolated fluxes back into real space by multiplying by the right
        eigenvector matrix.

        :arg list derivatives: The derivatives to perform the characteristic decomposition and WENO on.
        :arg object kernel: The current computational kernel."""
        pp_equations = []
        ndim = block.ndim

        reconstructed_characteristics = Matrix([d.evaluate_reconstruction for d in derivatives])
        # Transformation matrix back to physical (non-characteristic) space
        avg_REV_values = self.create_REV_inverses(dire)
        # HLLC: Reconstruct the left and right states separately and then transform each back to physical space
        left_q, right_q = symbols("qL:%d" % (ndim+2), **{'cls':GridVariable}), symbols("qR:%d" % (ndim+2), **{'cls':GridVariable})
        wL, wR = Matrix(reconstructed_characteristics[:,1]), Matrix(reconstructed_characteristics[:,0])
        pp_equations += [OpenSBLIEq(x, y) for x, y in zip(right_q, avg_REV_values*wR)]
        pp_equations += [OpenSBLIEq(x, y) for x, y in zip(left_q, avg_REV_values*wL)]
        # HLLC: Calculate the speed of sound and pressure at each interface using the reconstructed values
        aL, aR = symbols("aL", **{'cls':GridVariable}), symbols("aR", **{'cls':GridVariable})
        pL, pR = symbols("pL", **{'cls':GridVariable}), symbols("pR", **{'cls':GridVariable})
        rhoL, rhoR = left_q[0], right_q[0]
        rhoEL, rhoER = left_q[-1], right_q[-1]
        vel_L, vel_R = symbols("uL:%d" % (ndim), **{'cls':GridVariable}), symbols("uR:%d" % (ndim), **{'cls':GridVariable})
        # Density inversion factor
        inv_rhoL, inv_rhoR = symbols("inv_rhoL", **{'cls':GridVariable}), symbols("inv_rhoR", **{'cls':GridVariable})
        pp_equations += [OpenSBLIEq(inv_rhoL, 1.0/rhoL), OpenSBLIEq(inv_rhoR, 1.0/rhoR)]
        # Velocity component in the dire of reconstruction
        uL, uR = vel_L[dire], vel_R[dire]
        for i, u in enumerate(vel_L):
            pp_equations += [OpenSBLIEq(vel_L[i], left_q[i+1]*inv_rhoL)]
            pp_equations += [OpenSBLIEq(vel_R[i], right_q[i+1]*inv_rhoR)]
        # WARNING: ideal gas law assumed
        gama = ConstantObject('gama')
        pp_equations += [OpenSBLIEq(pL, (gama- 1)*(left_q[-1] - 0.5*rhoL*sum([x**2 for x in vel_L])))]
        pp_equations += [OpenSBLIEq(pR, (gama- 1)*(right_q[-1] - 0.5*rhoR*sum([x**2 for x in vel_R])))]
        pp_equations += [OpenSBLIEq(aL, sqrt(gama*pL*inv_rhoL))]
        pp_equations += [OpenSBLIEq(aR, sqrt(gama*pR*inv_rhoR))]
        # Compute wave speeds
        sL, sR = symbols("sL", **{'cls':GridVariable}), symbols("sR", **{'cls':GridVariable})
        smin, smax = symbols("smin", **{'cls':GridVariable}), symbols("smax", **{'cls':GridVariable})
        pp_equations += [OpenSBLIEq(sL, Min(GridVariable('AVG_%d_u%d' % (dire,dire)) -  GridVariable('AVG_%d_a' % (dire)), uL - aL))]
        pp_equations += [OpenSBLIEq(sR, Max(GridVariable('AVG_%d_u%d' % (dire,dire)) +  GridVariable('AVG_%d_a' % (dire)), uR + aR))]
        # Intermediate star speed, need to modify these for metrics?
        s_star = symbols('s_star', **{'cls':GridVariable})
        pp_equations += [OpenSBLIEq(s_star, (pR - pL + rhoL*uL*(sL - uL) - rhoR*uR*(sR - uR)) / (rhoL*(sL - uL) - rhoR*(sR - uR)))]
        # Build the system of flux components
        # Create a substitution dictionary from the flux components
        SD_L, SD_R = self.define_substitution_dictionaries(dire, derivatives, block, [left_q, right_q], [vel_L, vel_R], [pL, pR])
        # Create output arrays to store the final flux reconstructions
        reconstructed_work = self.create_output_wk_arrays(dire, derivatives, block)
        # Vectors for the right/left states
        U_L, U_R = Matrix([0 for _ in range(ndim+2)]), Matrix([0 for _ in range(ndim+2)])
        F_L, F_R = Matrix([0 for _ in range(ndim+2)]), Matrix([0 for _ in range(ndim+2)])
        # Substitute left/right states into the input Q and flux vector
        for i, d in enumerate(derivatives):
            input_args = d.args[0]
            U_L[i] = self.input_solution_vector[i].subs(SD_L)
            U_R[i] = self.input_solution_vector[i].subs(SD_R)
            F_L[i] = input_args.subs(SD_L)
            F_R[i] = input_args.subs(SD_R)


        # Build the star states, left and right
        if ndim == 1:
            USTAR_L = Matrix([1, s_star, rhoEL*inv_rhoL + (s_star - vel_L[dire])*(s_star + pL/(rhoL*(sL-vel_L[dire])))])
            USTAR_R = Matrix([1, s_star, rhoER*inv_rhoR + (s_star - vel_R[dire])*(s_star + pR/(rhoR*(sL-vel_R[dire])))])
        elif ndim == 2:
            if dire == 0:
                USTAR_L = Matrix([1, s_star, vel_L[1], rhoEL*inv_rhoL + (s_star - vel_L[dire])*(s_star + pL/(rhoL*(sL-vel_L[dire])))])
                USTAR_R = Matrix([1, s_star, vel_R[1], rhoER*inv_rhoR + (s_star - vel_R[dire])*(s_star + pR/(rhoR*(sL-vel_R[dire])))])
            elif dire == 1:
                USTAR_L = Matrix([1, vel_L[0], s_star, rhoEL*inv_rhoL + (s_star - vel_L[dire])*(s_star + pL/(rhoL*(sL-vel_L[dire])))])
                USTAR_R = Matrix([1, vel_R[0], s_star, rhoER*inv_rhoR + (s_star - vel_R[dire])*(s_star + pR/(rhoR*(sL-vel_R[dire])))])
        elif ndim == 3:
            if dire == 0:
                USTAR_L = Matrix([1, s_star, vel_L[1], vel_L[2], rhoEL*inv_rhoL + (s_star - vel_L[dire])*(s_star + pL/(rhoL*(sL-vel_L[dire])))])
                USTAR_R = Matrix([1, s_star, vel_R[1], vel_R[2], rhoER*inv_rhoR + (s_star - vel_R[dire])*(s_star + pR/(rhoR*(sL-vel_R[dire])))])
            elif dire == 1:
                USTAR_L = Matrix([1, vel_L[0], s_star, vel_L[2], rhoEL*inv_rhoL + (s_star - vel_L[dire])*(s_star + pL/(rhoL*(sL-vel_L[dire])))])
                USTAR_R = Matrix([1, vel_R[0], s_star, vel_R[2], rhoER*inv_rhoR + (s_star - vel_R[dire])*(s_star + pR/(rhoR*(sL-vel_R[dire])))])
            elif dire == 2:
                USTAR_L = Matrix([1, vel_L[0], vel_L[1], s_star, rhoEL*inv_rhoL + (s_star - vel_L[dire])*(s_star + pL/(rhoL*(sL-vel_L[dire])))])
                USTAR_R = Matrix([1, vel_R[0], vel_R[1], s_star, rhoER*inv_rhoR + (s_star - vel_R[dire])*(s_star + pR/(rhoR*(sL-vel_R[dire])))])

        # Add extra dissipation via scaling term if required on the signal speeds
        FC = ConstantObject('shock_filter_control')
        FC.value = 1.0 # Default condition has no scaling
        ConstantsToDeclare.add_constant(FC)

        # Outside density and wave-speed factor
        USTAR_L *= rhoL*((sL - vel_L[dire])/(sL - s_star))
        USTAR_R *= rhoR*((sR - vel_R[dire])/(sR - s_star))

        # Flux variables
        flux_vars = [gv('F%d' % i) for i in range(ndim+2)]
        if self.flux_type == 'HLLC-LM': # low-Mach correction formulation of HLLC solver
            # Mach number reduction of nonlinear signal speeds (HLLC-LM variant)
            phi, M_local = GridVariable('phi'), GridVariable('M_local')
            pp_equations += [OpenSBLIEq(M_local, Max(Abs(uL/aL), Abs(uR/aR)))]
            Mach_limit = 0.1 # Ensure the correction only applies when uL is less than 10% of the local sound speed
            pp_equations += [OpenSBLIEq(phi, sin(Min(1.0, M_local/Mach_limit)*Rational(1,2)*pi))]
            # Modify the wave-speeds
            pp_equations += [OpenSBLIEq(sL, sL*phi), OpenSBLIEq(sR, sR*phi)]
            # Build the conditional states
            condition1 = (F_L, sL >= 0)
            condition2 = (F_R, sR <= 0)
            F_STAR = Rational(1,2)*(F_L + F_R) + Rational(1,2)*(sL*(USTAR_L - U_L) + Abs(s_star)*(USTAR_L - USTAR_R) + sR*(USTAR_R - U_R))
            condition3 = (F_STAR, True)
            # Assign the fluxes to the storage arrays
            for i, component in enumerate(flux_vars):
                pp_equations += [OpenSBLIEq(flux_vars[i], Piecewise(*[(condition1[0][i], condition1[1]), (condition2[0][i], condition2[1]), (condition3[0][i], condition3[1])]))]

        else: # standard HLLC formulation
            # Build the conditional states
            condition1 = (F_L, sL >= 0)
            condition2 = (F_L + sL*(USTAR_L - U_L), And(sL <= 0, 0 <= s_star))
            condition3 = (F_R + sR*(USTAR_R - U_R), And(s_star <= 0, 0 <= sR))
            condition4 = (F_R, sR <= 0)
            condition5 = ([0 for _ in range(block.ndim+2)], True)
            # Assign the fluxes to the storage arrays
            for i, component in enumerate(flux_vars):
                pp_equations += [OpenSBLIEq(flux_vars[i], Piecewise(*[(condition1[0][i], condition1[1]), (condition2[0][i], condition2[1]), (condition3[0][i], condition3[1]), (condition4[0][i], condition4[1]), (condition5[0][i], condition5[1])]))]

        # Check for positive density/pressures
        # pp_equations += self.positivity_limiter(dire, block, derivatives, pL, pR, rhoL, rhoR, flux_vars)
        # Assign to global storage work arrays
        pp_equations += [OpenSBLIEq(reconstructed_work[i], flux_vars[i]) for i in range(ndim+2)]
        # Apply a shock sensor if the WENO is being applied as a filter step
        if block.shock_filter:
            for i in range(len(reconstructed_work)):
                pp_equations += flatten([self.central_diff_formula(i, reconstructed_work[i], self.flux_type, derivatives[i])])
        # Replace gamma factors if required
        pp_equations = self.replace_gamma_factor(pp_equations)
        return pp_equations


    def positivity_limiter(self, dire, block, derivatives, pL, pR, rhoL, rhoR, flux_vars):
        from sympy.functions.elementary.piecewise import ExprCondPair
        from opensbli.core.opensbliobjects import GroupedPiecewise
        from opensbli.schemes.spatial.weno import RightWenoReconstructionVariable, LeftWenoReconstructionVariable
        from opensbli.schemes.spatial.teno import RightTenoReconstructionVariable, LeftTenoReconstructionVariable
        ndim = block.ndim
        # Initialise variables
        thm, thp, theta_rho, theta_pressure = gv('theta_m'), gv('theta_p'), gv('theta_rho'), gv('theta_pressure')
        output_eqns = []
        # Minimum density and pressure
        # Variable to control CFL number that can be used, default alpha=2
        eps_rho, eps_p = ConstantObject('eps_rho'), ConstantObject('eps_p')
        eps_rho.value, eps_p.value = 1.0e-13, 1.0e-13
        ConstantsToDeclare.add_constant(eps_rho)
        ConstantsToDeclare.add_constant(eps_p)

        # Check for low density and pressures
        rho_i, rho_i1 = increment_dataset(block.location_dataset('rho'), dire, 0), increment_dataset(block.location_dataset('rho'), dire, 1)
        # Right (positive) state
        check = rhoR < eps_rho
        solve_theta = OpenSBLIEq(thp, (eps_rho - rho_i)/(rhoR - rho_i))
        cond1 = ExprCondPair(solve_theta, check)
        cond2 = ExprCondPair(OpenSBLIEq(thp,1), True)
        output_eqns += [GroupedPiecewise(cond1, cond2)]
        # Left (negative) state
        check = rhoL < eps_rho
        solve_theta = OpenSBLIEq(thm, (eps_rho - rho_i1)/(rhoL - rho_i1))
        cond1 = ExprCondPair(solve_theta, check)
        cond2 = ExprCondPair(OpenSBLIEq(thm,1), True)
        output_eqns += [GroupedPiecewise(cond1, cond2)]
        # Find minimum
        output_eqns += [OpenSBLIEq(theta_rho, Min(thm, thp))]
        # Check for pressures
        p_i, p_i1 = increment_dataset(block.location_dataset('p'), dire, 0), increment_dataset(block.location_dataset('p'), dire, 1)
        # Right (positive) state
        check = pR < eps_p
        solve_theta = OpenSBLIEq(thp, (eps_p - p_i)/(pR - p_i))
        cond1 = ExprCondPair(solve_theta, check)
        cond2 = ExprCondPair(OpenSBLIEq(thp,1), True)
        output_eqns += [GroupedPiecewise(cond1, cond2)]
        # Left (negative) state
        check = pR < eps_p
        solve_theta = OpenSBLIEq(thm, (eps_p - p_i1)/(pL - p_i1))
        cond1 = ExprCondPair(solve_theta, check)
        cond2 = ExprCondPair(OpenSBLIEq(thm,1), True)
        output_eqns += [GroupedPiecewise(cond1, cond2)]
        # Find minimum
        output_eqns += [OpenSBLIEq(theta_pressure, Min(thm, thp))]
        # Correct the fluxes as a convex combination of the first order LF flux
        # F_(i+1/2) = 0.5*(F_i + F_(i+1) - (u+a)*(U_i - U_(i+1))
        F_fixed = Matrix([0 for _ in range(ndim+2)])
        # if ndim == 1:
        #     alpha = self.grid_EV[-1,-1]
        # else:
        #     alpha = self.grid_EV[-2,-2]
        # Max wave-speed
        alpha = Abs(block.location_dataset('u%d' % dire)) + block.location_dataset('a')
        for i in range(ndim+2):
            q, q1 = self.input_solution_vector[i], increment_dataset(self.input_solution_vector[i], dire, 1)
            F, F1 = derivatives[i].args[0], increment_dataset(derivatives[i].args[0], dire, 1)
            LF = Rational(1,2)*(F+F1 + alpha*(q - q1))
            F_fixed[i] = (1 - theta_rho*theta_pressure)*LF + theta_rho*theta_pressure*flux_vars[i]
            # pprint(F_fixed)

        # Set the corrected flux
        output_eqns += [OpenSBLIEq(flux_vars[i], F_fixed[i]) for i in range(ndim+2)]
        return output_eqns
