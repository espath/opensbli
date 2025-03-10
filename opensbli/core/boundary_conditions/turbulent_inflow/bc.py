""" Boundary Condition for turbulent inflow"""

from sympy import symbols, exp
from opensbli.equation_types.opensbliequations import OpenSBLIEquation as Eq
from opensbli.core.kernel import ConstantsToDeclare as CTD
from opensbli.code_generation.algorithm.common import *
from opensbli.core.boundary_conditions.Carpenter_scheme import Carpenter
from opensbli.core.io_hdf5 import iohdf5
from sympy import Function
from opensbli.core.grid import GridVariable
from opensbli.core.opensbliobjects import DataObject, ConstantObject
from sympy.functions.elementary.piecewise import ExprCondPair, Piecewise
from opensbli.utilities.helperfunctions import increment_dataset
from opensbli.core.boundary_conditions.bc_core import BoundaryConditionBase, ModifyCentralDerivative, WallBC
from opensbli.schemes.spatial.shock_capturing import ShockCapturing
from opensbli.core.boundary_conditions.periodic import PeriodicBC

class UserFunction(Function):
    def __new__(cls, *args):
        ret = super(UserFunction, cls).__new__(cls, *args, evaluate=False)
        return ret
    
    def name(self):
        return self.args[0]


class TurbulentInflowBC(ModifyCentralDerivative, BoundaryConditionBase):
    """ Applies a synthetic turbulence signal at inflow.

    :arg int direction: Spatial direction to apply boundary condition to.
    :arg int side: Side 0 or 1 to apply the boundary condition for a given direction."""
    def __init__(self, boundary_direction, side, scheme=None, plane=True):
        BoundaryConditionBase.__init__(self, boundary_direction, side, plane)
        self.bc_name = 'turbulent inflow'
        if not scheme:
            self.modification_scheme = Carpenter()
        else:
            self.modification_scheme = scheme
        return
    
    @property
    def lengthscales(self):
        return self._lscale_y, self._lscale_z
    
    @lengthscales.setter
    def lengthscales(self, values):
        self._lscale_y = values[0]
        self._lscale_z = values[1]

    def apply(self, arrays, block):
        # Integral length scales and filter variables
        lscale_y, lscale_z = self.lengthscales
        # Lscale = Matrix(3, 3, lambda i, j: GridVariable('Lscale%d%d' % (i, j)))
        Lscale = lscale_y
        
        # Fluctuating components
        Upp, Vpp, Wpp = GridVariable('Upp'), GridVariable('Vpp'), GridVariable('Wpp')

        # Random numbers
        seed1, seed2, seed3, seed4, seed5, seed6 = DataObject('seed1'), DataObject('seed2'), DataObject('seed3'), DataObject('seed4'), DataObject('seed5'), DataObject('seed6')
        rand_u, rand_v, rand_w = GridVariable('rand_u'), GridVariable('rand_v'), GridVariable('rand_w')
        # rand_u, rand_v, rand_w = DataObject('rand_u'), DataObject('rand_v'), DataObject('rand_w')

        # Total velocity why data objects these can be grid variables SPJ
        Utotal, Vtotal, Wtotal, Dtotal, Ptotal = symbols("Utotal, Vtotal, Wtotal, Dtotal, Ptotal", **{'cls':GridVariable})
        rho, rhou0, rhou1, rhou2, rhoE = DataObject('rho'), DataObject('rhou0'), DataObject('rhou1'), DataObject('rhou2'), DataObject('rhoE')
        
        # Phi, psi: intermediate variables for correlation mehod
        phi_u, phi_v, phi_w = DataObject('phi_u'), DataObject('phi_v'), DataObject('phi_w')
        # phi_u_ij, phi_v_ij, phi_w_ij = block.location_dataset(phi_u), block.location_dataset(phi_v), block.location_dataset(phi_w)
        # phi_u_next, phi_v_next, phi_w_next = increment_dataset(phi_u_ij, 1, 1), increment_dataset(phi_v_ij, 1, 1), increment_dataset(phi_w_ij, 1, 1)

        psi_u, psi_v, psi_w = DataObject('psi_u'), DataObject('psi_v'), DataObject('psi_w')
        # psi_u_ij, psi_v_ij, psi_w_ij = block.location_dataset(psi_u), block.location_dataset(psi_v), block.location_dataset(psi_w)
        # psi_u_next, psi_v_next, psi_w_next = increment_dataset(psi_u_ij, 2, 1), increment_dataset(psi_v_ij, 2, 1), increment_dataset(psi_w_ij, 2, 1)

        phi_u_old, phi_v_old, phi_w_old = DataObject('phi_u_old'), DataObject('phi_v_old'), DataObject('phi_w_old')
        psi_u_old, psi_v_old, psi_w_old = DataObject('psi_u_old'), DataObject('psi_v_old'), DataObject('psi_w_old')

        phi_u_ij, phi_v_ij, phi_w_ij = block.location_dataset(phi_u_old), block.location_dataset(phi_v_old), block.location_dataset(phi_w_old)
        phi_u_prev, phi_v_prev, phi_w_prev = increment_dataset(phi_u_ij, 1, -1), increment_dataset(phi_v_ij, 1, -1), increment_dataset(phi_w_ij, 1, -1)
        psi_u_ij, psi_v_ij, psi_w_ij = block.location_dataset(psi_u_old), block.location_dataset(psi_v_old), block.location_dataset(psi_w_old)
        psi_u_prev, psi_v_prev, psi_w_prev = increment_dataset(psi_u_ij, 2, -1), increment_dataset(psi_v_ij, 2, -1), increment_dataset(psi_w_ij, 2, -1)

        # RMS amplitudes & mean velocity profiles 
        # These are to be read in 
        ustar, vstar, wstar = DataObject('ustar'), DataObject('vstar'), DataObject('wstar')
        A11, A21, A22, A31, A32, A33 = DataObject('A11'), DataObject('A21'), DataObject('A22'), DataObject('A31'), DataObject('A32'), DataObject('A33')
        Umean, Vmean, Wmean, Dmean = DataObject('Umean'), DataObject('Vmean'), DataObject('Wmean'), DataObject('Dmean')
        
        dummy0, dummy1, dummy2, dummy3, dummy4 = GridVariable('dummy0'), GridVariable('dummy1'), GridVariable('dummy2'), GridVariable('dummy3'), GridVariable('dummy4')
        '''Set values to be read from hdf5. undeclared values will default to zeros'''
        
        self._read_from_hdf5_arrays = [Umean, Vmean, Wmean, Dmean, A11, A21, A22, A31, A32, A33, seed1, seed2, seed3, seed4, seed5, seed6]
        # Unscaled fluctations # These will be updated every time step

        '''
        Equations
        '''
        eqns = []
        eqns_phi = []
        eqns_psi = []
        eqns_ustar = []

        # Inlet conditions
        M, gama, dt, restart_iteration_no, block0np1 = symbols("Minf gama dt restart_iteration_no block0np1", **{'cls':ConstantObject})
        # Add the constants to CTD
        CTD.add_constant([M, gama, dt, restart_iteration_no, block0np1])
        dx, dz = block.deltas[0], block.deltas[2]
        nx, ny, nz = block.shape

        y = DataObject('x1')
        y_current = block.location_dataset(y)
        y_next = increment_dataset(y_current, 1, 1)
        dy = GridVariable('dy')
        eqns_phi += [Eq(dy, y_next - y_current)]

        # turn on/off fluctuation scaling and time/space correlation
        scaling = True
        spatial_correlation = True
        temporal_correlation = True

        # Equations for defining velocity, density, pressure
        eqns += [Eq(Utotal, Umean + ustar*A11)]
        eqns += [Eq(Vtotal, Vmean + (ustar*A21 + vstar*A22))]
        eqns += [Eq(Wtotal, Wmean + (ustar*A31 + vstar*A32 +wstar*A33))]
        eqns += [Eq(Dtotal, Dmean)]
        eqns += [Eq(Ptotal, 1.0 /(gama*M**2))]

        # Normal distibution random numbers
        i, j, k = block.grid_indexes
        current_iter = block.get_temporal_schemes[0].iteration_number + restart_iteration_no
        eqns_phi += [Eq(rand_u, UserFunction('random_generator', current_iter, seed1, seed2))]
        eqns_phi += [Eq(rand_v, UserFunction('random_generator', current_iter, seed3, seed4))]
        eqns_phi += [Eq(rand_w, UserFunction('random_generator', current_iter, seed5, seed6))]
        
        # Correlation function equations - u component
        eqns_phi += [Eq(phi_u, phi_u_prev*exp(-dy/Lscale[1, 0]) + rand_u*(1 - exp(-2*dy/Lscale[1, 0]))**0.5)]
        eqns_psi += [Eq(psi_u, psi_u_prev*exp(-dz/Lscale[2, 0]) + phi_u*(1 - exp(-2*dz/Lscale[2, 0]))**0.5)]
        eqns_ustar += [Eq(ustar, ustar*exp(-Umean*dt/Lscale[0, 0]) + psi_u*(1 - exp(-2*Umean*dt/Lscale[0, 0]))**0.5)]

        # Correlation method - v component
        eqns_phi += [Eq(phi_v, phi_v_prev*exp(-dy/Lscale[1, 1]) + rand_v*(1 - exp(-2*dy/Lscale[1, 1]))**0.5)]
        eqns_psi += [Eq(psi_v, psi_v_prev*exp(-dz/Lscale[2, 1]) + phi_v*(1 - exp(-2*dz/Lscale[2, 1]))**0.5)]
        eqns_ustar += [Eq(vstar, vstar*exp(-Umean*dt/Lscale[0, 1]) + psi_v*(1 - exp(-2*Umean*dt/Lscale[0, 1]))**0.5)]

        # Correlation method - w component
        eqns_phi += [Eq(phi_w, phi_w_prev*exp(-dy/Lscale[1, 2]) + rand_w*(1 - exp(-2*dy/Lscale[1, 2]))**0.5)]
        eqns_psi += [Eq(psi_w, psi_w_prev*exp(-dz/Lscale[2, 2]) + phi_w*(1 - exp(-2*dz/Lscale[2, 2]))**0.5)]
        eqns_ustar += [Eq(wstar, wstar*exp(-Umean*dt/Lscale[0, 2]) + psi_w*(1 - exp(-2*Umean*dt/Lscale[0, 2]))**0.5)]

        # set limits for very large fluctuations
        eqns_ustar += [Eq(ustar, Piecewise((5.0, ustar>5.0), (-5.0, ustar<-5.0), (ustar, True)))]
        eqns_ustar += [Eq(vstar, Piecewise((5.0, vstar>5.0), (-5.0, vstar<-5.0), (vstar, True)))]
        eqns_ustar += [Eq(wstar, Piecewise((5.0, wstar>5.0), (-5.0, wstar<-5.0), (wstar, True)))]

        # update old values
        eqns_ustar += [Eq(phi_u_old, phi_u)]
        eqns_ustar += [Eq(psi_u_old, psi_u)]
        eqns_ustar += [Eq(phi_v_old, phi_v)]
        eqns_ustar += [Eq(psi_v_old, psi_v)]
        eqns_ustar += [Eq(phi_w_old, phi_w)]
        eqns_ustar += [Eq(psi_w_old, psi_w)]

        # compute state variables
        eqns += [Eq(rho, Dtotal)]
        eqns += [Eq(rhou0, Dtotal*Utotal)]
        eqns += [Eq(rhou1, Dtotal*Vtotal)]
        eqns += [Eq(rhou2, Dtotal*Wtotal)]
        eqns += [Eq(rhoE, Ptotal/(gama-1) + 0.5* Dtotal *(Utotal**2 + Vtotal**2 + Wtotal**2))]

        # Dummy variables used as temporary fix - ops R/W issue with weno
        eqns += [Eq(dummy0, rho)]
        eqns += [Eq(dummy1, rhou0)]
        eqns += [Eq(dummy2, rhou1)]
        eqns += [Eq(dummy3, rhou2)]
        eqns += [Eq(dummy4, rhoE)]

        # set kernels
        halos, kernel_phi = self.generate_boundary_kernel(block, self.bc_name)
        halos, kernel_psi = self.generate_boundary_kernel(block, self.bc_name)
        halos, kernel_ustar = self.generate_boundary_kernel(block, self.bc_name)
        halos, kernel_apply = self.generate_boundary_kernel(block, self.bc_name)

        # Apply equations to kernels      
        for eq in eqns_phi:
            eq1 = eq.convert_to_datasets(block)
            kernel_phi.add_equation(eq1)

        for eq in eqns_psi:
            eq1 = eq.convert_to_datasets(block)
            kernel_psi.add_equation(eq1)

        for eq in eqns_ustar:
            eq1 = eq.convert_to_datasets(block)
            kernel_ustar.add_equation(eq1)

        for eq in eqns:
            eq1 = eq.convert_to_datasets(block)
            kernel_apply.add_equation(eq1)

        # Apply dirichlet bc if weno/teno
        if any(isinstance(sc, ShockCapturing) for sc in block.discretisation_schemes.values()):
            # TODO optimise
            n_halos = abs(halos[self.direction][self.side])
            cons_vars = [block.location_dataset(d) for d in [rho, rhou0, rhou1, rhou2, rhoE]]
            base_loc = list(cons_vars[0].indices)
            rhs_value = [Dtotal, Dtotal*Utotal, Dtotal*Vtotal, Dtotal*Wtotal, Ptotal/(gama-1) + 0.5* Dtotal *(Utotal**2 + Vtotal**2 + Wtotal**2)]
            rhs_values = []
            from_side_factor, to_side_factor = self.set_side_factor()
            for i in range(1, n_halos+1):
                loc_lhs= base_loc[:]
                loc_lhs[self.direction] += from_side_factor*i
                for left, right in zip(cons_vars, rhs_value):
                    left = self.convert_dataset_base_expr_to_datasets(left, loc_lhs)
                    kernel_apply.add_equation(Eq(left, right))
        
        for kernel in [kernel_phi, kernel_psi, kernel_ustar, kernel_apply]:
            kernel.update_block_datasets(block)

        self.set_read_from_hdf5(block)
        # Periodic exchange for psi values
        # for span periodic cases - ensures z-wise correlation between boundaries
        all_kernels = [kernel_phi, kernel_psi]

        # do halo exchange in y direction
        exchange_array = [block.location_dataset(e) for e in [phi_u, phi_v, phi_w]]
        bc1 = PeriodicBC(1,0)
        bc2 = PeriodicBC(1,1)
        k1 = bc1.apply(exchange_array, block)
        k2 = bc2.apply(exchange_array, block)
        all_kernels += [k1, k2]

        # do halo exchange in z direction
        exchange_array = [block.location_dataset(e) for e in [psi_u, psi_v, psi_w]]
        bc3 = PeriodicBC(2,0)
        bc4 = PeriodicBC(2,1)
        k3 = bc3.apply(exchange_array, block)
        k4 = bc4.apply(exchange_array, block)
        all_kernels += [k3, k4]

        all_kernels += [kernel_ustar]
        self.iteration_kernels = all_kernels

        return kernel_apply
    
    def set_read_from_hdf5(self, block):
        kwargs = {'iotype': "Read", 'filename':'correlationdata.h5'}
        h5_read = iohdf5(**kwargs)
        h5_read.add_arrays(self._read_from_hdf5_arrays)
        block.add_io([h5_read])
        return