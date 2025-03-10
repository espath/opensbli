"""@brief
   @authors Max Walker, David J. Lusher (09/2023)
   @contributors
   @details
"""

from sympy import IndexedBase, Symbol, Rational, solve, interpolating_poly, integrate, Abs, Float, flatten, S
from opensbli.core.opensblifunctions import TVDDerivative
from opensbli.core.opensbliobjects import ConstantObject
from opensbli.core.kernel import Kernel
from opensbli.equation_types.opensbliequations import SimulationEquations, OpenSBLIEq, NonSimulationEquations
from opensbli.core.grid import GridVariable as gv
from .scheme import Scheme
from sympy import horner, pprint
from opensbli.schemes.spatial.shock_capturing import ShockCapturing, TVDCharacteristic

class TVDHalos(object):
    """ Object for TVD halos.

    :arg int order: Order of the TVD scheme.
    :arg bool reconstruction: True if halos for a reconstruction. """

    def __init__(self, order, reconstruction=None):
        k = int(0.5*(order+1))
        if not reconstruction:
            self.halos = [-k, k+1]
        else:
            self.halos = [-1, 1]
        return

    def get_halos(self, side):
        return self.halos[side]

class TVD(Scheme, ShockCapturing):
    """ Base class for TVD implementation."""

    def __init__(self):
        Scheme.__init__(self, "TVDDerivative", order=1)
        self.schemetype = "Spatial"
        self.halotype = TVDHalos(2)
        self.required_constituent_relations_symbols = {}
        return []

    def reconstruction_halotype(self, order, reconstruction=True):
        return TVDHalos(order, reconstruction)

    def group_by_direction(self, eqs):
        """ Groups the input equations by the direction (x0, x1, ...) they depend upon.

        :arg list eqs: List of equations to group by direction.
        :returns: dict: grouped: Dictionary of {direction: equations} key, value pairs for equations grouped by direction."""
        all_WDS = []
        for eq in eqs:
            all_WDS += list(eq.atoms(TVDDerivative))
        grouped = {}
        for cd in all_WDS:
            direction = cd.get_direction[0]
            if direction in grouped.keys():
                grouped[direction] += [cd]
            else:
                grouped[direction] = [cd]
        return grouped

class TVDFlux(TVDCharacteristic, TVD):
    """ TVDFlux class, called within the TVDFilter file as a non-linear filter step.""

    :arg object physics: Physics object, defaults to NSPhysics.
    :arg object averaging: The averaging procedure to be applied for characteristics. """

    def __init__(self, physics=None, averaging=None, shock_filter=None, conservative=True, flux_split=False, species=None):
        self.flux_split = flux_split
        self.species = species
        self.flux_type = 'TVD'
        self.temp_wk_arrays = []
        TVDCharacteristic.__init__(self, physics, averaging, flux_split)
        self.conservative = conservative
        if shock_filter is not None:
            self.shock_filter = shock_filter
            self.sensor_evaluation = TVD.__init__(self)
        return

    def discretise(self, type_of_eq, block):
        """ Discretise function for calling the TVD routines."""
        if isinstance(type_of_eq, SimulationEquations):
            pass # Currently only added TVD as a filter step option
        # Apply WENO as a non-linear filter step instead
        elif isinstance(type_of_eq, NonSimulationEquations):
            eqs = flatten(type_of_eq.equations)
            grouped = self.group_by_direction(eqs)
            all_derivatives_evaluated_locally = []
            reconstruction_halos = self.reconstruction_halotype(self.order, reconstruction=True)
            solution_vector = flatten(type_of_eq.time_advance_arrays)

            # Instantiate eigensystems with block, but don't add metrics yet
            self.instantiate_eigensystem(block, self.species)
            for direction, derivatives in sorted(grouped.items()):
                all_derivatives_evaluated_locally += derivatives
                for no, deriv in enumerate(derivatives):
                    deriv.create_reconstruction_work_array(block)
                # Kernel for the reconstruction in this direction
                kernel = self.create_reconstruction_kernel(direction, reconstruction_halos, block)
                # Get the pre, interpolations and post equations for characteristic reconstruction
                pre_process, reductions, interpolated, post_process = self.get_characteristic_equations(direction, derivatives, solution_vector, block)
                if direction == 0:
                    reduction_output = reductions
                # Add the equations to the kernel and add the kernel to SimulationEquations
                kernel.add_equation(pre_process + interpolated + post_process)
                type_of_eq.reconstruction_kernels += [kernel]
            # Generate kernels for the constituent relations
            if grouped:
                constituent_relations = None
                type_of_eq.residual_kernels += [self.evaluate_residuals(block, eqs, all_derivatives_evaluated_locally)]
                # constituent_relations = self.check_constituent_relations(block, eqs, constituent_relations)
            return constituent_relations, solution_vector, reduction_output
