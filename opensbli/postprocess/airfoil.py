from sympy import sin, exp, pi, tan, cos, symbols, Eq, pprint, Abs, flatten
from sympy.functions.elementary.piecewise import Piecewise, ExprCondPair
from opensbli.core.grid import GridVariable
from opensbli.core.opensbliobjects import ConstantObject, DataObject
from opensbli.core.kernel import Kernel
from opensbli.core.boundary_conditions import DirichletBC
from opensbli.core.boundary_conditions.multi_block import InterfaceBC
import copy

def generate_outlet_sponge(q_vector, block, Lx, npoints):
    """ Applies a sponge boundary on the outflow of the domain to damp oscillations."""
    block_number = block.blocknumber
    sponge_length, lc, sigma = symbols("spongel, lc, sigma", **{'cls':GridVariable})
    gama, Minf = symbols("gama Minf", **{'cls':ConstantObject})
    residual = symbols("Residual0:%d" % (block.ndim+2), **{'cls':DataObject})
    x0 = symbols("x0", **{'cls':DataObject})
    equations = []
    equations += [Eq(sponge_length, 0.85)]
    # Characteristic length
    equations += [Eq(lc, x0 - (Lx - sponge_length)), Eq(sigma,  0.5*(1.0 + cos(pi* lc/sponge_length)))]
    # Freestream values to enforce
    if block.ndim == 2:
        rho_inf, u_inf, v_inf, p_inf = 1.0, 1.0, 0.0, 1.0/(gama*Minf*Minf)
        rhoE_inf = p_inf / (gama - 1.0) + 0.5*rho_inf*(u_inf**2 + v_inf**2)
        boundary_values = [rho_inf, rho_inf*u_inf, rho_inf*v_inf, rhoE_inf]
    else:
        rho_inf, u_inf, v_inf, w_inf, p_inf = 1.0, 1.0, 0.0, 0.0, 1.0/(gama*Minf*Minf)
        rhoE_inf = p_inf / (gama - 1.0) + 0.5*rho_inf*(u_inf**2 + v_inf**2 + w_inf**2)
        boundary_values = [rho_inf, rho_inf*u_inf, rho_inf*v_inf, rho_inf*w_inf, rhoE_inf]

    # Application
    for residual_array, qvec, freestream in zip(residual, q_vector, boundary_values):
        equations += [Eq(residual_array, residual_array - sigma * (qvec - freestream))]
    eqns = block.dataobjects_to_datasets_on_block(equations)
    ker = Kernel(block, computation_name="Outlet sponge zone block%d" % block_number)
    ker.kernelname = "outlet_sponge_block%d" % block_number
    ker.add_equation(eqns)
    ranges = copy.deepcopy(block.ranges)
    ker.ranges = ranges
    # Reduce the evaluation range
    ker.ranges[0][0] =  ranges[0][1] - npoints
    ker.update_block_datasets(block)
    return ker

def generate_farfield_sponge(q_vector, block, Ly, npoints):
    """ Applies a sponge boundary on the farfield of the domain to damp oscillations."""
    # WARNING: Need to update these for swept cases
    block_number = block.blocknumber
    sponge_length, lc, sigma = symbols("spongel, lc, sigma", **{'cls':GridVariable})
    gama, Minf = symbols("gama Minf", **{'cls':ConstantObject})
    residual = symbols("Residual0:%d" % (block.ndim+2), **{'cls':DataObject})
    equations = []
    equations += [Eq(sponge_length, 0.85)]
    # Characteristic length
    if block_number == 0 or block_number == 2:
        x1 = symbols("x1", **{'cls':DataObject})
        equations += [Eq(lc, Abs(x1) - (Ly - sponge_length)), Eq(sigma,  0.5*(1.0 + cos(pi* lc/sponge_length)))]
    elif block_number == 1:
        x0 = symbols("x0", **{'cls':DataObject})
        x1 = symbols("x1", **{'cls':DataObject})
        equations += [Eq(lc, Abs(x0**2 + x1**2) - (Ly - sponge_length)), Eq(sigma,  0.5*(1.0 + cos(pi* lc/sponge_length)))]
    else:
        raise ValueError("Sponge zones are configured for three blocks.")
    # Freestream values to enforce
    if block.ndim == 2:
        rho_inf, u_inf, v_inf, p_inf = 1.0, 1.0, 0.0, 1.0/(gama*Minf*Minf)
        rhoE_inf = p_inf / (gama - 1.0) + 0.5*rho_inf*(u_inf**2 + v_inf**2)
        boundary_values = [rho_inf, rho_inf*u_inf, rho_inf*v_inf, rhoE_inf]
    else:
        rho_inf, u_inf, v_inf, w_inf, p_inf = 1.0, 1.0, 0.0, 0.0, 1.0/(gama*Minf*Minf)
        rhoE_inf = p_inf / (gama - 1.0) + 0.5*rho_inf*(u_inf**2 + v_inf**2 + w_inf**2)
        boundary_values = [rho_inf, rho_inf*u_inf, rho_inf*v_inf, rho_inf*w_inf, rhoE_inf]
    # Application
    for residual_array, qvec, freestream in zip(residual, q_vector, boundary_values):
        equations += [Eq(residual_array, residual_array - sigma * (qvec - freestream))]
    eqns = block.dataobjects_to_datasets_on_block(equations)
    ker = Kernel(block, computation_name="Farfield sponge zone block%d" % block.blocknumber)
    ker.kernelname = "farfield_sponge_block%d" % block.blocknumber
    ker.add_equation(eqns)
    ranges = copy.deepcopy(block.ranges)
    for eqn in ker.equations:
        pprint(eqn)
    ker.ranges = ranges
    # Reduce the evaluation range
    ker.ranges[1][0] =  ranges[1][1] - npoints
    ker.update_block_datasets(block)
    return ker


def generate_wake_kernel(q_vector, multi_block, wall_energy):
    """ Function to generate a kernel that averages over the wake line between block 0 and block 2."""
    # Apply to block 0 (lower wake block)
    block = multi_block.get_block(0)
    # Temporary work arrays
    nvars = len(q_vector)
    wk = symbols("wk0:%d" % nvars, **{'cls':DataObject})
    # x index
    idx = block.grid_indexes[0]
    # Set the density
    equations = [Eq(q_vector[0], 0.5 * (q_vector[0] + wk[0]))]
    # for rhou,v,w
    for b0, b1 in zip(q_vector[1:-1], wk[1:-1]):
        pairs = [ExprCondPair(0.0, idx == 1000000000), ExprCondPair(0.5* (b0 + b1), True)]
        equations += [Eq(b0, Piecewise(*pairs, evaluate=False))]
    pairs = [ExprCondPair(wall_energy.rhs, idx == 1000000000), ExprCondPair(0.5* (q_vector[-1] + wk[-1]), True)]
    equations += [Eq(q_vector[-1], Piecewise(*pairs, evaluate=False))]
    
    equations = block.dataobjects_to_datasets_on_block(equations)
    direction = 1
    side = 0
    bc = DirichletBC(direction, side, equations)
    application_kernel = bc.apply([], block)
    application_kernel.kernelname = "wake_treatment"
    application_kernel.computation_name = "Wake treatment"
    application_kernel.halo_ranges[1][0] = set()
    # Wake exchanges from block2 wakeline (q_vector) to blokck0 work_arrays
    block2 = multi_block.get_block(2)
    bc = InterfaceBC(direction, side,  halos=[-4, 4], match=(0, 1, 0, False))
    arrays = [block2.work_array(str(a)) for a in flatten(q_vector)]
    other_arrays = [block.work_array(str(a)) for a in flatten(wk)]
    wake_transfer1 = bc.apply_interface(arrays, block2, multi_block, other_arrays=other_arrays)
    wake_transfer1.transfer_size[1] = 1
    wake_transfer1.transfer_from[1] = 0
    wake_transfer1.transfer_to[1] = 0
    wake_transfer1.computation_name = "wake_block0_to_block2_"
    
    bc = InterfaceBC(direction, side,  halos=[-4, 4], match=(2, 1, 0, False))
    arrays = [block.work_array(str(a)) for a in flatten(q_vector)]
    wake_transfer2 = bc.apply_interface(arrays, block, multi_block)
    wake_transfer2.transfer_size[1] = 1
    wake_transfer2.transfer_from[1] = 0
    wake_transfer2.transfer_to[1] = 0
    wake_transfer2.computation_name = "wake_block2_to_block0_"
    return [wake_transfer1, application_kernel, wake_transfer2]

