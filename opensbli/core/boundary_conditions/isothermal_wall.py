from opensbli.core.boundary_conditions.bc_core import BoundaryConditionBase, ModifyCentralDerivative, WallBC
from opensbli.core.boundary_conditions.Carpenter_scheme import Carpenter
from opensbli.utilities.helperfunctions import increment_dataset
from opensbli.equation_types.opensbliequations import OpenSBLIEq
from opensbli.physical_models.ns_physics import NSphysics
from opensbli.core.opensbliobjects import ConstantObject
from sympy import Rational, S, Float, pprint
from opensbli.core.grid import GridVariable
from opensbli.schemes.spatial.weno import ShockCapturing


class IsothermalWall_ZeroPressureGradBC(ModifyCentralDerivative, BoundaryConditionBase, WallBC):
    """ Navier-Stokes specific boundary condition. Applies a no-slip viscous wall condition,
    velocity components are zero on the wall. Temperature is fixed with a prescribed wall temperature,
    given as the rhoE equation passed to this BC.

    :arg int direction: Spatial direction to apply boundary condition to.
    :arg int side: Side 0 or 1 to apply the boundary condition for a given direction.
    :arg Eq equations: Equation for conservative variable rhoE to set on the wall with constant temperature.
    :arg object scheme: Boundary scheme if required, defaults to Carpenter boundary treatment.
    :arg bool plane: True/False: Apply boundary condition to full range/split range only."""

    def __init__(self, direction, side, equations, scheme=None, corners=True, plane=True, multi_block=False):
        BoundaryConditionBase.__init__(self, direction, side, plane)
        self.bc_name = 'IsothermalWall_ZeroPressureGrad'
        self.corners = corners
        self.equations = equations
        self.multi_block = multi_block
        if not scheme:
            self.modification_scheme = Carpenter()
        else:
            self.modification_scheme = scheme
        return

    def apply(self, arrays, block):
        halos, kernel = self.generate_boundary_kernel(block, self.bc_name, corners=self.corners)
        n_halos = abs(halos[self.direction][self.side])
        # Using Navier Stokes physics object, create conservative variables
        NS = NSphysics(block)
        from_side_factor, to_side_factor = self.set_side_factor()
        wall_eqns = []
        # Update halos if a shock capturing scheme is being used.
        # if any(isinstance(sc, ShockCapturing) for sc in block.discretisation_schemes.values()) and self.shock:
        pw, gamma, Minf, Twall = GridVariable('P_above'), NS.specific_heat_ratio(), NS.mach_number(), ConstantObject('Twall')
        from opensbli.core.kernel import ConstantsToDeclare as CTD
        CTD.add_constant(Twall)
        # Wall pressure
        wall_eqns += [OpenSBLIEq(pw, increment_dataset(NS.pressure(relation=True, conservative=True), self.direction, to_side_factor*1))]
        # Set wall conditions, momentum zero
        wall_eqns += [OpenSBLIEq(x, Float(S.Zero)) for x in NS.momentum()]
        # Set density and energy on the wall
        wall_eqns += [OpenSBLIEq(NS.density(), pw*Twall*gamma*Minf*Minf)] # p_inf = 1/(gamma*Minf^2)
        wall_eqns += [OpenSBLIEq(NS.total_energy(), pw / (gamma - 1))]
        kernel.add_equation(wall_eqns)
        # Direct copy indices over the wall
        kernel.update_block_datasets(block)
        print("Using the zero pressure gradient version of the isothermal wall. Please set the residual arrays to zero on the wall plane for this BC.")
        return kernel

class IsothermalWallBC(ModifyCentralDerivative, BoundaryConditionBase, WallBC):
    """ Navier-Stokes specific boundary condition. Applies a no-slip viscous wall condition,
    velocity components are zero on the wall. Temperature is fixed with a prescribed wall temperature,
    given as the rhoE equation passed to this BC.

    :arg int direction: Spatial direction to apply boundary condition to.
    :arg int side: Side 0 or 1 to apply the boundary condition for a given direction.
    :arg Eq equations: Equation for conservative variable rhoE to set on the wall with constant temperature.
    :arg object scheme: Boundary scheme if required, defaults to Carpenter boundary treatment.
    :arg bool plane: True/False: Apply boundary condition to full range/split range only."""

    def __init__(self, direction, side, equations, scheme=None, corners=True, plane=True, multi_block=False):
        BoundaryConditionBase.__init__(self, direction, side, plane)
        self.bc_name = 'IsothermalWall'
        self.corners = corners
        self.equations = equations
        self.multi_block = multi_block
        if not scheme:
            self.modification_scheme = Carpenter()
        else:
            self.modification_scheme = scheme
        return

    def apply(self, arrays, block):
        halos, kernel = self.generate_boundary_kernel(block, self.bc_name, corners=self.corners)
        n_halos = abs(halos[self.direction][self.side])
        # Using Navier Stokes physics object, create conservative variables
        NS = NSphysics(block)
        from_side_factor, to_side_factor = self.set_side_factor()
        wall_eqns = []
        # Set wall conditions, momentum zero
        if block.conservative:
            wall_eqns += [OpenSBLIEq(x, Float(S.Zero)) for x in NS.momentum()]
        else:
            wall_eqns += [OpenSBLIEq(x, Float(S.Zero)) for x in NS.velocity()]
        # kernel.add_equation(OpenSBLIEq(GridVariable('x0'), block.location_dataset('x0')))  # Used for sidewall case to set the ramp
        # Set wall energy, density is left to be evaluated
        wall_eqns += self.equations[:]
        kernel.add_equation(wall_eqns)
        # Update halos if a shock capturing scheme is being used.
        # if any(isinstance(sc, ShockCapturing) for sc in block.discretisation_schemes.values()) and self.shock:
        # if not self.multi_block: # Don't evaluate halos if using multi-block. Causing issues at the corners
        if block.conservative:
            pw, gamma, Minf, Twall = GridVariable('Pwall'), NS.specific_heat_ratio(), NS.mach_number(), ConstantObject('Twall')
            # Wall pressure
            wall_eqns = [OpenSBLIEq(pw, increment_dataset(NS.pressure(relation=True, conservative=True), self.direction, 0))]
            kernel.add_equation(wall_eqns)
            # Direct copy indices over the wall
            transfer_indices = [tuple([from_side_factor*t, to_side_factor*t]) for t in range(1, n_halos + 1)]
            # Evaluate velocities in the halos
            halo_velocities = []
            for index, momentum_comp in enumerate(NS.momentum()):
                for pair_index, index_pair in enumerate(transfer_indices):
                    momentum, density = increment_dataset(momentum_comp, self.direction, index_pair[1]), increment_dataset(NS.density(), self.direction, index_pair[1])
                    halo_velocities += [OpenSBLIEq(GridVariable('u%d%d' % (index, pair_index+1)), momentum/density)]
            kernel.add_equation(halo_velocities)
            # Evaluate the temperature one point above the wall, use it to interpolate temperature in the halos
            T_above = GridVariable('T_above')
            kernel.add_equation(OpenSBLIEq(GridVariable('T_above'), increment_dataset(NS.temperature(relation=True, conservative=True), self.direction, to_side_factor)))
            halo_densities = []
            halo_momentum = []
            halo_energy = []
            # Evaluate halo density using wall pressure and interpolated temperatures

            for i in range(1, n_halos+1):
                # Interpolate temperature into the halos and calculate halo density
                rho_halo = GridVariable('rho_halo_%d' % i)
                T_halo = GridVariable('T%d' % i)
                kernel.add_equation([OpenSBLIEq(T_halo, (i+1)*Twall - i*T_above)])
                halo_densities += [OpenSBLIEq(rho_halo, pw*gamma*Minf*Minf/T_halo)]
                # Set momentum in the halos, reverse direction
                velocities = []
                for index, momentum_comp in enumerate(NS.momentum()):
                    velocity = GridVariable('u%d%d' % (index, i))
                    halo_momentum += [OpenSBLIEq(increment_dataset(momentum_comp, self.direction, from_side_factor*i), -rho_halo*velocity)]
                    velocities += [velocity**2]
                # Set energy in halos based on wall pressure, interpolated temperature and the calculated density
                energy_rhs = pw/(gamma-1) + Rational(1, 2)*rho_halo*(sum(velocities))
                halo_energy += [OpenSBLIEq(increment_dataset(NS.total_energy(), self.direction, from_side_factor*i), energy_rhs)]
                halo_densities += [OpenSBLIEq(increment_dataset(NS.density(), self.direction, from_side_factor*i), rho_halo)]
                kernel.add_equation(halo_densities)
                kernel.add_equation(halo_momentum)
                kernel.add_equation(halo_energy)
        else:
            pw, gamma, Minf, Twall = GridVariable('Pwall'), NS.specific_heat_ratio(), NS.mach_number(), ConstantObject('Twall')
            rho, Et = block.location_dataset('rho'), block.location_dataset('Et')
            # Wall pressure
            # rho*(gama-1)*(Et - (1/2)*(KD(_i,_j)*u_i*u_j))
            vel_sq = sum([x**2 for x in NS.velocity()])
            wall_eqns = OpenSBLIEq(pw, rho*(gamma - 1)*(Et - 0.5*vel_sq))
            kernel.add_equation(wall_eqns)
            # Direct copy indices over the wall
            transfer_indices = [tuple([from_side_factor*t, to_side_factor*t]) for t in range(1, n_halos + 1)]
            # Evaluate the temperature one point above the wall, use it to interpolate temperature in the halos
            T_above = GridVariable('T_above')
            # Et = p/(rho*(gamma-1)) + 0.5*U**2
            # p = rho*R*T, T=p*rho*gamma*Minf**2
            # Et = T/(Minf**2 * (gamma - 1)) + 0.5*U**2
            T = (Et - 0.5*vel_sq)*Minf**2 * (gamma - 1)
            kernel.add_equation(OpenSBLIEq(GridVariable('T_above'), increment_dataset(T, self.direction, to_side_factor)))
            # Evaluate halo density using wall pressure and interpolated temperatures
            for i in range(1, n_halos+1):
                halo_densities, halo_velocities, halo_energy = [], [], []
                # Interpolate temperature into the halos and calculate halo density
                rho_halo = GridVariable('rho_halo_%d' % i)
                T_halo = GridVariable('T%d' % i)
                kernel.add_equation([OpenSBLIEq(T_halo, (i+1)*Twall - i*T_above)])
                halo_densities += [OpenSBLIEq(rho_halo, pw*gamma*Minf*Minf/T_halo)]
                # Set velocity in the halos, reverse direction
                velocities = []
                for index, velocity_comp in enumerate(NS.velocity()):
                    halo_velocities += [OpenSBLIEq(increment_dataset(velocity_comp, self.direction, from_side_factor*i), -increment_dataset(velocity_comp, self.direction, to_side_factor*i))]
                # Set energy in halos based on wall pressure, interpolated temperature and the calculated density
                energy_rhs = pw/(rho_halo*(gamma-1)) + Rational(1, 2)*increment_dataset(vel_sq, self.direction, from_side_factor*i)
                halo_energy += [OpenSBLIEq(increment_dataset(Et, self.direction, from_side_factor*i), energy_rhs)]
                halo_densities += [OpenSBLIEq(increment_dataset(rho, self.direction, from_side_factor*i), rho_halo)]
                kernel.add_equation(halo_densities)
                kernel.add_equation(halo_velocities)
                kernel.add_equation(halo_energy)
        kernel.update_block_datasets(block)
        return kernel
