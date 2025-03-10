"""@brief
   @authors David J Lusher
   @contributors
   @details
"""

from sympy import flatten, Idx, sqrt, Rational, pprint
from opensbli.core.opensbliobjects import ConstantObject, ConstantIndexed, Globalvariable
from opensbli.core.grid import GridVariable
from opensbli.equation_types.opensbliequations import OpenSBLIEq
from opensbli.core.kernel import Kernel
from opensbli.core.datatypes import Int
from opensbli.schemes.spatial.scheme import Scheme, TemporalSolution
from opensbli.core.kernel import ConstantsToDeclare as CTD


class RungeKuttaLS(Scheme):
    """ Applies a low storage (2 Register, Williamson form) Runge-Kutta scheme from: Fourth-Order Kutta Schemes 2N-Storage (Carpenter, 1994).
        Optimal SSP coefficients for the 3rd order scheme are taken from "Gottlieb, Shu, Tadmor (2001): Strong Stability-Preserving High-Order Time
        Discretization Methods, SIAM Review Vol. 43, No.1, pp 89-112.

        For m stages and current time 'n' with u^0 = u^n, du^0 = 0:
        do i = 1..m:
            du^i = A[i]*du^(i-1) + dt*Residual
            u^i = u^(i-1) + B[i]*du^i
        u^n+1 = u^m."

        :arg int order: The order of accuracy of the scheme."""

    def __init__(cls, order, stages=None, formulation=None):
        Scheme.__init__(cls, "RungeKutta", order)
        cls.solution = {}
        cls.schemetype = "Temporal"
        cls.formulation = formulation
        cls.stages, cls.order = stages, order
        # Create constants
        cls.create_constants(order, stages)
        cls.add_constants()
        # Update the RK coefficients
        cls.get_coefficients
        if formulation == 'SSP':
            if order != 3:
                raise ValueError("RK-SSP is only available for the 3rd order scheme.")
            else:
                print("An SSP Runge-Kutta scheme of order %d is being used for time-stepping." % order)
        else:
            print("A Runge-Kutta scheme of order %d is being used for time-stepping." % order)
        return

    def set_stages(cls, order, stages):
        if order == 3:  # 3rd order schemes are 3-stage
            if cls.formulation == 'SSP':
                if stages == None:
                    n_stages = 3 # default is (3,3) for RK-SSP
                else:
                    n_stages = stages # User input stages
            else:
                n_stages = order # Regular (3,3) RK scheme, non-SSP
        elif order == 4:  # 4th order scheme is 5-stage
            n_stages = order + 1
        return n_stages

    def create_constants(cls, order, stages):
        # Set the number of stages if provided
        cls.n_stages = cls.set_stages(order, stages)
        cls.stage = Idx('stage', cls.n_stages)
        cls.stage.restart = None
        cls.solution_coeffs = ConstantIndexed('rkB', cls.stage)
        cls.stage_coeffs = ConstantIndexed('rkA', cls.stage)
        cls.niter_symbol = ConstantObject('niter', integer=True)
        cls.niter_symbol.datatype = Int()
        cls.iteration_number = Globalvariable("iter", integer=True)
        cls.iteration_number._value = None
        cls.iteration_number.datatype = Int()
        # As iteration number is used in a for loop we dont add them to constants to declare
        cls.temporal_iteration = Idx(cls.iteration_number, cls.niter_symbol)
        cls.constant_time_step = True
        cls.time_step = ConstantObject("dt")
        # Variables to hold the simulation time and starting iteration
        cls.start_time = ConstantObject('simulation_time', restart=True)
        cls.start_time.value = 0.0
        cls.start_iter = ConstantObject('start_iter', restart=True, integer=True)
        cls.start_iter.value = 0
        cls.start_iter.datatype = Int()
        cls.temporal_iteration.restart = cls.start_iter
        return

    def add_constants(cls):
        CTD.add_constant(cls.niter_symbol)
        CTD.add_constant(cls.solution_coeffs)
        CTD.add_constant(cls.stage_coeffs)
        CTD.add_constant(cls.time_step)
        CTD.add_constant(cls.start_time)
        CTD.add_constant(cls.start_iter)
        return

    @property
    def get_coefficients(cls):
        """ Create A (intermediate update) and B (solution advance) coefficients for the RK scheme.
        From: S. Ruuth. GLOBAL OPTIMIZATION OF EXPLICIT STRONG-STABILITY-PRESERVING RUNGE-KUTTA METHODS (2005)."""
        if cls.order == 3:
            if cls.formulation == 'SSP':
                if cls.n_stages == 3: # CFL coefficient 0.322349301195940
                    B1, B2, B3 = 0.924574112262461, 0.287712943868770, 0.626538293270800
                    A1, A2, A3 = 0.0, -2.915493957701923, 0.0
                    cls.solution_coeffs.value = [B1, B2, B3]
                    cls.stage_coeffs.value = [A1, A2, A3]
                elif cls.n_stages == 4: # CFL coefficient 0.634274456962008
                    B1, B2, B3, B4 = 1.086620745813428, 0.854115548251602, -1.576604558206099, -0.278475500113052
                    A1, A2, A3, A4 = 0.0, -0.449336503268844, 0.0, -4.661555711601366
                    cls.solution_coeffs.value = [B1, B2, B3, B4]
                    cls.stage_coeffs.value = [A1, A2, A3, A4]
                elif cls.n_stages == 5: # CFL coefficient 1.40154693827206
                    B1, B2, B3, B4, B5 = 0.713497331193829, 0.133505249805329, 0.713497331193829, 0.149579395628565, 0.384471116121269
                    A1, A2, A3, A4, A5 = 0.0, -4.344339134485095, 0.0, -3.770024161386381, -0.046347284573284
                    cls.solution_coeffs.value = [B1, B2, B3, B4, B5]
                    cls.stage_coeffs.value = [A1, A2, A3, A4, A5]               
                else:
                    raise ValueError("The 3rd order RK-SSP is defined for 3, 4, or 5 stages.")
            else:
                A1, A2, A3 = 0, Rational(-5, 9), Rational(-153, 128)
                B1, B2, B3 = Rational(1, 3), Rational(15, 16), Rational(8, 15)
                cls.solution_coeffs.value = [B1, B2, B3]
                cls.stage_coeffs.value = [A1, A2, A3]
        elif cls.order == 4:
            A1, A2, A3, A4, A5 = 0, -0.4178904745, -1.192151694643, -1.697784692471, -1.514183444257
            B1, B2, B3, B4, B5 = 0.1496590219993, 0.3792103129999, 0.8229550293869, 0.6994504559488, 0.1530572479681
            cls.solution_coeffs.value = [B1, B2, B3, B4, B5]
            cls.stage_coeffs.value = [A1, A2, A3, A4, A5]
        else:
            raise NotImplementedError("Only 3rd and 4th order RK schemes are currently implemented.")
        return

    def __str__(cls):
        return "%s" % (cls.__class__.__name__)

    def get_local_function(cls, list_of_components):
        """ Finds the time derivatives to be advanced."""

        from opensbli.core.opensblifunctions import TemporalDerivative
        CD_fns = []
        for c in flatten(list_of_components):
            CD_fns += list(c.atoms(TemporalDerivative))
        return CD_fns

    def discretise(cls, type_of_eq, block):
        """ Main discretise function for the temporal advancement."""
        # We need only the equations as they contain residual residual_arrays
        if type_of_eq in cls.solution.keys():
            pass
        else:
            cls.solution[type_of_eq] = TemporalSolution()
        td_fns = cls.get_local_function(type_of_eq.equations)
        if td_fns:
            # Create a Kernel for the update ()
            temp_data_sets = cls.create_temp_data_sets(td_fns, block)
            new_data_sets = [eq.time_advance_array for eq in td_fns]
            cls.var_solved = new_data_sets
            cls.temp_RK_arrays = temp_data_sets
            # Create the stage and solution updates
            residuals = [eq.residual for eq in flatten(type_of_eq.equations)]
            zipped = zip(temp_data_sets, new_data_sets, residuals)
            kernels = cls.create_discretisation_kernel(zipped, block)
            cls.solution[type_of_eq].kernels += kernels
            # New testing
            type_of_eq.temporalsolution = TemporalSolution()
            type_of_eq.temporalsolution.kernels += kernels
            type_of_eq.temporalsolution.start_kernels += cls.solution[type_of_eq].start_kernels
        # Re apply the constants for multi-block, deepcopy was clearing them
        cls.create_constants(cls.order, cls.stages)
        cls.add_constants()
        return

    def convert_to_conservative(cls, equations, block):
        # Convert between conservative/primitive form before and after the time update
        rho = block.location_dataset('rho')
        rho_inv = GridVariable('rho_inv')
        output = []
        primitive_to_conservative = [OpenSBLIEq(var, rho*var) for var in cls.var_solved[1:]]
        conservative_to_primitive = [OpenSBLIEq(var, rho_inv*var) for var in cls.var_solved[1:]]
        inv = [OpenSBLIEq(rho_inv, 1./rho)]
        output += primitive_to_conservative + equations + inv + conservative_to_primitive
        return output

    def create_discretisation_kernel(cls, zipped, block):
        """ Creates the kernels for the intermediate step and time update.

        :arg list zipped: List of tuples containing the intermediate, solution and residual arrays for each equation being solved.
        :arg object block: OpenSBLI SimulationBlock.
        :returns: list: List of the two discretised Kernels required for the RK scheme."""
        solution_update_kernel = Kernel(block, computation_name="Temporal solution advancement")
        # Update the range of evaluation
        solution_update_kernel.set_grid_range(block)
        # Update the solution and stages
        if cls.constant_time_step:
            solution_update = cls.constant_time_step_solution(zipped)
        if not block.conservative:
            solution_update = cls.convert_to_conservative(solution_update, block)

        solution_update_kernel.add_equation(solution_update)
        solution_update_kernel.update_block_datasets(block)
        return [solution_update_kernel]

    def constant_time_step_solution(cls, zipped):
        """ Creates the equations for the intermediate step and solution update kernels.

        :arg list zipped: List of tuples containing the intermediate, solution and residual arrays for each equation being solved.
        :returns: list solution_update: Equations for the time advancement of the solution.
        :returns: list intermediate_update: Equations for the intermediate update step."""
        dt = cls.time_step
        solution_update = []
        for z in zipped:
            solution_update += [OpenSBLIEq(z[0], cls.stage_coeffs*z[0] + dt*z[2], evaluate=False)]
            solution_update += [OpenSBLIEq(z[1], z[1] + cls.solution_coeffs*z[0], evaluate=False)]
        return solution_update

    def create_temp_data_sets(cls, equations, block):
        """Creates the arrays to store the intermediate update.

        :arg list equations: Equations to be advanced in time.
        :arg object block: OpenSBLI SimulationBlock.
        :returns: list temp_datasets: List of the temporary storage DataSets."""
        temp_data_sets = []
        for no, eq in enumerate(flatten(equations)):
            fn = eq.time_advance_array
            temp_data_sets += [block.work_array('%s_RKold' % fn.base.label)]
        return temp_data_sets

    def generate_inner_loop(cls, kernels):
        from opensbli.code_generation.algorithm.algorithm import DoLoop
        rkloop = DoLoop(cls.stage)
        rkloop.add_components(kernels)
        return rkloop
