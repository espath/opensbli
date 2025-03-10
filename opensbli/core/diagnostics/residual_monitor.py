""" Calculates the maximum of the residuals over the domain, to check convergence to a steady state."""

from sympy import symbols, exp, pprint
from opensbli.core.opensbliobjects import DataObject, ConstantObject, ReductionSum
from opensbli.equation_types.opensbliequations import OpenSBLIEquation, SimulationEquations
from opensbli.postprocess.post_process_eq import *
from opensbli.core.kernel import ConstantsToDeclare as CTD
from opensbli.code_generation.algorithm.common import *
from opensbli.utilities.user_defined_kernels import UserDefinedEquations
from opensbli.equation_types.opensbliequations import OpenSBLIEq


class ResidualMonitor(object):
    """ Calculates the L2 norm of each of the residuals (dQ). The scaling is done before the SimulationMonitoring printing, separately."""
    def __init__(self, block, frequency=100):
        self.block = block
        self.ndim = block.ndim
        self.frequency = frequency
        self.equation_classes = []
        # Arrays for the filtered solution
        self.create_equations(block)
        return

    def create_equations(self, block):
        # Create a kernel at the end of the time loop, every iteration (no frequency)
        filter_class = UserDefinedEquations()
        filter_class.algorithm_place = InTheSimulation(frequency=self.frequency)
        filter_class.computation_name = 'Reduction of residuals'
        # Get the residuals defined on the block
        residuals = []
        for eqn_class in block.list_of_equation_classes:
            if isinstance(eqn_class, SimulationEquations):
                no_residuals = len(eqn_class.equations) + 1
        reduction_dsets = [block.location_dataset('Residual%d' % (i)) for i in range(no_residuals)]
        reduction_names = ['L2_R%d' % i for i in range(no_residuals)]
        reduction_vars = [ReductionSum(x) for x in reduction_names]
        # Create reduction equations
        output_equations = [OpenSBLIEq(x, reduction_dsets[i]**2) for i, x in enumerate(reduction_vars)]
        filter_class.add_equations(output_equations)
        self.equation_classes.append(filter_class)
        return
