from sympy import Rational, Min, Abs, sqrt, tanh, pprint, Max, exp, Rational
from opensbli.core.opensblifunctions import CentralDerivative as CD
from opensbli.core.parsing import EinsteinEquation as EE
from opensbli.core.opensbliobjects import ConstantObject, CoordinateObject, DataObject
from opensbli.equation_types.opensbliequations import OpenSBLIEq
from opensbli.core.kernel import ConstantsToDeclare as CTD
from opensbli.utilities.helperfunctions import increment_dataset
from opensbli.core.grid import GridVariable as gv


class ShockSensor(object):
    def __init__(self):
        self.epsilon = 1.0e-40
        return

    def ducros_equations(self, block, coordinate_symbol, metrics=None, name='kappa', Mach=None):
        """ Create the non-discretized equations for the modified Ducros shock sensor and applies a metric transformation if required.
        :arg object block: OpenSBLI simulation block.
        :arg string coordinate_symbol: Coordinate symbol to perform the derivatives with.
        :arg object metrics: OpenSBLI metric class to apply curvilinear coordinates to the sensor if required."""

        if Mach == None:
            Mach = 1

        ndim = block.ndim
        substitutions, constants, output_eqns = [], [], []
        cart = CoordinateObject('x_i')
        cartesian_coordinates = [cart.apply_index(cart.indices[0], dim) for dim in range(ndim)]
        sensor_array = block.location_dataset('%s' % name)

        # Calculate vorticity
        if block.ndim == 2:
            dx, dy = block.deltas
            x0, x1 = cartesian_coordinates[0], cartesian_coordinates[1]
            u0, u1 = DataObject('u0'), DataObject('u1')
            vorticity_sq = (CD(u1, x0) - CD(u0, x1))**2
        elif block.ndim == 3:
            dx, dy, dz = block.deltas
            x0, x1, x2 = cartesian_coordinates[0], cartesian_coordinates[1], cartesian_coordinates[2]
            u0, u1, u2 = DataObject('u0'), DataObject('u1'), DataObject('u2')
            vorticity_sq = (CD(u2, x1) - CD(u1, x2))**2 + (CD(u0, x2) - CD(u2, x0))**2 + (CD(u1, x0) - CD(u0, x1))**2
        else:
            raise ValueError("vorticity is not defined for one-dimensional cases.")

        divergence = "Eq(divergence, (KD(_i,_j)*Der(u_i, x_j)))"
        divergence = EE().expand(divergence, ndim, coordinate_symbol, substitutions, constants)

        # Apply metric transformation if required
        if metrics:
            vorticity_sq = metrics.apply_transformation(vorticity_sq)
            divergence = metrics.apply_transformation(divergence)

        # Tanh function doesn't vectorise, use exponentials instead
        a = 200 / sqrt(sum([block.deltas[i]**2 for i in range(ndim)]))
        tanh_filter = 0.5*(1 - tanh(2.5*(1+a*divergence.rhs)))
        # tanh_inner = 2.5*(1 + c*divergence.rhs)
        # tanh_filter = Rational(1, 2)*(1 - (exp(2*tanh_inner, evaluate=False) - 1) / (exp(2*tanh_inner, evaluate=False) + 1))
        output_eqns += [OpenSBLIEq(sensor_array, tanh_filter*divergence.rhs**2 / (divergence.rhs**2 + vorticity_sq + self.epsilon))]
        return output_eqns, sensor_array

    def Ren_sensor(self, block, name='kappa'):
        eps = 0.001
        # r_j
        base_loc = 0
        pm, p, pp = increment_dataset(block.location_dataset('p'), 0, base_loc -1), increment_dataset(block.location_dataset('p'), 0, base_loc), increment_dataset(block.location_dataset('p'), 0, base_loc + 1)
        ph, mh = pp - p, p - pm
        rj1 = (Abs(2*ph*mh) + eps) / (ph**2 + mh**2 + eps)
        # r_(j+1)
        base_loc = 1
        pm, p, pp = increment_dataset(block.location_dataset('p'), 0, base_loc -1), increment_dataset(block.location_dataset('p'), 0, base_loc), increment_dataset(block.location_dataset('p'), 0, base_loc + 1)
        ph, mh = pp - p, p - pm
        rj2 = (Abs(2*ph*mh) + eps) / (ph**2 + mh**2 + eps)
        output = 1 - Min(rj1, rj2)
        for dire in range(1, block.ndim):
            base_loc = 0
            pm, p, pp = increment_dataset(block.location_dataset('p'), dire, base_loc -1), increment_dataset(block.location_dataset('p'), dire, base_loc), increment_dataset(block.location_dataset('p'), dire, base_loc + 1)
            ph, mh = pp - p, p - pm
            rj1 = (Abs(2*ph*mh) + eps) / (ph**2 + mh**2 + eps)
            # r_(j+1)
            base_loc = 1
            pm, p, pp = increment_dataset(block.location_dataset('p'), dire, base_loc -1), increment_dataset(block.location_dataset('p'), dire, base_loc), increment_dataset(block.location_dataset('p'), dire, base_loc + 1)
            ph, mh = pp - p, p - pm
            rj2 = (Abs(2*ph*mh) + eps) / (ph**2 + mh**2 + eps)
            output = Max(output, 1 - Min(rj1, rj2))

        sensor_array = block.location_dataset('%s' % name)
        output = [OpenSBLIEq(sensor_array, output)]
        return output, sensor_array

    def WENO_1D_sensor(self, block, name='kappa'):
        # Evaluate based on pressure
        pm, p, pp = increment_dataset(block.location_dataset('p'), 0, -1), block.location_dataset('p'), increment_dataset(block.location_dataset('p'), 0, 1)
        # Weighting coefficients
        a, b = Rational(1,4), Rational(13,12)
        output = (a*(pp - pm)**2 + b*(pp - 2*p + pm)**2)**2
        sensor_array = block.location_dataset('%s' % name)
        output = [OpenSBLIEq(sensor_array, output)]
        return output, sensor_array


    def Jameson_sensor(self, block, name='kappa'):
        pm, p, pp = increment_dataset(block.location_dataset('p'), 0, -1), block.location_dataset('p'), increment_dataset(block.location_dataset('p'), 0, 1)
        output = Abs((pp - 2*p + pm) / (pp + 2*p + pm))
        for dire in range(1, block.ndim):
            pm, p, pp = increment_dataset(block.location_dataset('p'), dire, -1), block.location_dataset('p'), increment_dataset(block.location_dataset('p'), dire, 1)
            output = Max(output, Abs((pp - 2*p + pm) / (pp + 2*p + pm)))
        output = [OpenSBLIEq(sensor_array, output)]
        return output, sensor_array


