"""@brief
   @author Satya Pramod Jammy
   @contributors David Lusher
   @details
"""

from opensbli.equation_types.opensbliequations import NonSimulationEquations
from opensbli.core.kernel import Kernel
from opensbli.core.kernel import ConstantsToDeclare as CTD
from opensbli.core.opensbliobjects import GroupedPiecewise, ConstantObject, DataObject, DataSet
from sympy import Equality, flatten, pprint
from opensbli.code_generation.algorithm.common import BeforeSimulationStarts
from opensbli.schemes.spatial.scheme import CentralHalos_defdec


class GridBasedInitialisation(NonSimulationEquations):
    def __new__(cls, order=None, **kwargs):
        ret = super(GridBasedInitialisation, cls).__new__(cls)
        if order:  # Local order if multiple instances of the class are declared on the block more for future work
            ret.order = order
        else:
            ret.order = 0
        ret.equations = []
        ret.kwargs = kwargs
        # Variable to control restarting
        cls.restart = ConstantObject('restart', integer=True)
        cls.restart.value = 0
        # cls.restart.datatype = Int()
        CTD.add_constant(cls.restart)
        condition = Equality(cls.restart, 0)
        ret.algorithm_place = [BeforeSimulationStarts(start_condition=condition)]
        return ret

    def __hash__(self):
        h = hash(self._hashable_content())
        self._mhash = h
        return h

    def _hashable_content(self):
        return "GridBasedInitialisation"

    @property
    def evaluated_datasets(cls):
        evaluated = set()
        for eq in flatten(cls.equations):
            if isinstance(eq, Equality):
                evaluated = evaluated.union(eq.lhs_datasetbases)
            elif isinstance(eq, GroupedPiecewise):
                evaluated = evaluated.union(eq.lhs_datasetbases)
        return evaluated

    def spatial_discretisation(cls, block):
        """ Compute the equations for the initial condition, evaluated only once at the start of a simulation if restart == 0."""
        # Check for coordinate dependencies
        cls.check_coordinate_evaluation(block)
        kernel = Kernel(block, computation_name="Grid_based_initialisation%d" % cls.order)
        kernel.set_grid_range(block)
        schemes = block.discretisation_schemes
        for d in range(block.ndim):
            # Initialize all five halos
            kernel.set_halo_range(d, 0, CentralHalos_defdec())
            kernel.set_halo_range(d, 1, CentralHalos_defdec())
        kernel.add_equation(cls.equations)
        kernel.update_block_datasets(block)
        cls.Kernels = [kernel]
        return

    def check_coordinate_evaluation(cls, block):
        """ If the grid was created inside the initialisation kernel (not read in from a separate grid file), 
        the array declarations must be changed to also restart the coordinate arrays."""
        coordinate_symbol = "x"
        coordinate_names = [coordinate_symbol + '%d' % i + '_B%d' % block.blocknumber for i in range(block.ndim)]
        LHS = [x.lhs for x in flatten(cls.equations) if isinstance(x.lhs, DataObject)]
        LHS += [x.lhs for x in flatten(cls.equations) if isinstance(x.lhs, DataSet)]
        block.coordinate_arrays_to_restart = [x for x in LHS if str(x) in coordinate_names]
        return

    def apply_boundary_conditions(cls, block):
        """No boundary conditions in the Initialisation currently, any logic needed should be implemented here
        """
        return

    def apply_interface_bc(cls, block, multiblock_descriptor):
        """ No BC for initialisation."""
        return

