from opensbli.equation_types.opensbliequations import NonSimulationEquations, Discretisation, Solution
from opensbli.schemes.spatial.scheme import CentralDerivative
from opensbli.core.kernel import Kernel
from sympy import flatten, pprint, Equality
from opensbli.core.opensbliobjects import GroupedPiecewise

class UserDefinedEquations(NonSimulationEquations, Discretisation, Solution):
    """User defined equations. No checking is performed.
    Just forms a kernel on the range and places the kernel in the algorithm place passed as an
    input to the class. """

    def __new__(cls, **kwargs):
        ret = super(UserDefinedEquations, cls).__new__(cls)
        ret.order = 0
        ret.equations = []
        ret.kwargs = kwargs
        ret._place = []
        ret.computation_name = None
        # Optional halo type
        ret.halos = None
        ret.kernel_merge = False # by default don't merge the kernels
        cls._full_swap = False
        ret.custom_grid_range = None
        return ret

    @property
    def evaluated_datasets(cls):
        evaluated = set()
        for eq in flatten(cls.equations):
            if isinstance(eq, Equality):
                evaluated = evaluated.union(eq.lhs_datasetbases)
            elif isinstance(eq, GroupedPiecewise):
                evaluated = evaluated.union(eq.lhs_datasetbases)
        return evaluated

    @property
    def algorithm_place(cls):
        return cls._place

    @algorithm_place.setter
    def algorithm_place(cls, place):
        cls._place += [place]
        return

    @property
    def full_swap(cls):
        return cls._full_swap

    @full_swap.setter
    def full_swap(cls, full_swap):
        cls._full_swap = full_swap
        return

    def merge_kernels(cls, kernel_list, no_derivatives, block):
        """ Merges the multiple kernels provided by discretize into a single evaluation to evaluate the entire UDF in a single kernel.
        Currently assumes the grid range on all computations is the same. May need to change later."""
        UDF_equations = cls.equations
        derivative_evaluations = []
        for ker in kernel_list:
            for eqn in ker.equations:
                derivative_evaluations.append(eqn)
        
        merged_kernel = Kernel(block, computation_name=cls.computation_name)
        # Derivative evaluations
        merged_kernel.add_equation(derivative_evaluations)
        # Equations that depend on these derivative evaluations
        merged_kernel.add_equation(no_derivatives)
        merged_kernel.ranges = block.ranges[:]
        return [merged_kernel]


    def spatial_discretisation(cls, block):
        """ Applies the spatial discretisation of the equations by calling the discretisation of each spatial scheme provided on the block

        :param SimulationBlock block: the block on which the equations are solved
        :return: None """

        # Instantiate the solution class
        # cls.solution = Solution()

        # Discretize any derivatives in the equations
        spatialschemes = []
        # Get the schemes on the block
        schemes = block.discretisation_schemes
        for sc in schemes:
            if schemes[sc].schemetype == "Spatial":
                spatialschemes += [sc]
        # Perform spatial Discretisation if any in constituent relations evaluation
        # Input equations are saved here, before discretisation
        equations = cls.equations

        UDF_derivative_kernels = []
        evaluations = [] # these evaluations are never used here
        no_derivatives = []

        for eq in flatten(equations):
            if isinstance(eq, GroupedPiecewise):
                pass
            elif len(eq.rhs.atoms(CentralDerivative)) == 0: # checking for equations requiring derivative evaluation
                no_derivatives += [eq]
            else: # Need to compute the derivative
                cls.equations = [eq]
            # for sc in spatialschemes:
            #     # Constituent relations are returned
            evaluations.append(schemes[str(CentralDerivative)].discretise(cls, block)) # only Central should be used
            UDF_derivative_kernels.append(cls.Kernels[:])
            # Reset discretized Kernels and equations for this equation
            cls.Kernels = []
            cls.equations = []
        # Original input equations are restored here
        cls.equations = equations

        # No discretisation required
        if not flatten(UDF_derivative_kernels):
            user_defined_kernel = Kernel(block)
            user_defined_kernel.set_computation_name("User kernel: %s" % str(cls.computation_name))
            user_defined_kernel.set_grid_range(block)
            # Evaluate the kernel into the halos as required
            if cls.halos:
                for direction in range(block.ndim):
                    user_defined_kernel.set_halo_range(direction, 0, cls.halos[direction][0])
                    user_defined_kernel.set_halo_range(direction, 1, cls.halos[direction][1])
            user_defined_kernel.add_equation(cls.equations)
            # user_defined_kernel.update_block_datasets(block)
            cls.Kernels += [user_defined_kernel]
        else:
            cls.Kernels += flatten(UDF_derivative_kernels)

        if cls.kernel_merge: # merge kernels to evaluate all within one kernel
            cls.Kernels = cls.merge_kernels(cls.Kernels, no_derivatives, block)

        # Process the kernels to update parameters on the block
        cls.process_kernels(block)
        # Apply a custom grid range if defined
        if cls.custom_grid_range is not None:
            for ker in cls.Kernels:
                ker.ranges = cls.custom_grid_range
        return

    def process_kernels(cls, block):
        """A function to update some dependant parameters of each kernel

        :param SimulationBlock block: the block on which the equations are solved
        :return: None """
        for kernel in cls.Kernels:
            kernel.update_block_datasets(block)
        return

    def apply_boundary_conditions(cls, block):
        return

# MBCHANGE
    def apply_interface_bc(cls, block, multiblock_descriptor):
        return
# MBCHANGE
