#!/usr/bin/env python
# Import all the functions from opensbli
from opensbli import *
from sympy import pi

def generate_simulation_parameters(sim_params, multi_block):
    nblocks = multi_block.nblocks
    ndim = multi_block.blocks[0].ndim
    # Number of points per direction total
    TNx, TNy, TNz = 64, 64, 64
    for block_no, block in enumerate(multi_block.blocks):
        if nblocks == 2:
            nbx, nby, nbz = 2, 1, 1
            Lx, Ly, Lz = 'M_PI', '2*M_PI', '2*M_PI'
            Nx, Ny, Nz = int(TNx/nbx), int(TNy/nby), int(TNz/nbz)
        elif nblocks == 4:
            nbx, nby, nbz = 2, 2, 1
            Lx, Ly, Lz = 'M_PI', 'M_PI', '2*M_PI'
            Nx, Ny, Nz = int(TNx/nbx), int(TNy/nby), int(TNz/nbz)
        elif nblocks == 8:
            nbx, nby, nbz = 2, 2, 2
            Lx, Ly, Lz = 'M_PI', 'M_PI', 'M_PI'
            Nx, Ny, Nz = int(TNx/nbx), int(TNy/nby), int(TNz/nbz)
        else:
            raise ValueError("Only implemented simulation parameters for 2 < nblocks < 8.")
        block_params = {
        # Block parameters
        'block%dnp0' % block_no    :   '%d' % Nx,
        'block%dnp1' % block_no    :   '%d' % Ny,
        'block%dnp2' % block_no    :   '%d' % Nz,
        'Delta0block%d'  % block_no    :   '%s/block%dnp0' % (Lx, block_no),
        'Delta1block%d'  % block_no    :   '%s/block%dnp1' % (Ly, block_no),
        'Delta2block%d'  % block_no    :   '%s/block%dnp2' % (Lz, block_no),
        }
        sim_params = {**sim_params, **block_params}
    return sim_params

# Initial condition on each block
def TGV_initial_condition(multi_block):
    """ Create the TGV initial condition based on a multi-block decomposition of equally sized blocks."""
    nblocks = multi_block.nblocks
    ndim = multi_block.blocks[0].ndim
    mb_initial_conditions = [None for _ in range(nblocks)]
    mb_initial_conditions = dict(list(enumerate(mb_initial_conditions)))
    if nblocks < 2:
        raise ValueError("This problem is set up for nblocks 2 <= 8.")
    if nblocks % 2 != 0:
        raise ValueError("Number of blocks should be a multiple of 2, for symmetry.")

    # Create block decomposition depending on number of blocks (2, 4, 8)
    if nblocks == 2:
        nbx, nby, nbz = 2, 1, 1
    elif nblocks == 4:
        nbx, nby, nbz = 2, 1, 1
    elif nblocks == 8:
        nbx, nby, nbz = 2, 2, 2
    else:
        raise ValueError("This problem is set up for nblocks 2 <= 8.")

    x, y, z = symbols("x0:%d" % ndim,  **{'cls':DataObject})
    # Loop over all the blocks to set the coordinates
    for block_no, block in enumerate(multi_block.blocks):
        # Create an initial condition kernel
        init_class = GridBasedInitialisation()
        # Local dictionary for parsing the expressions
        local_dict = {"block": block, "GridVariable": GridVariable, "DataObject": DataObject}
        # Grid spacing and indexing on this block
        dx, dy, dz = block.deltas
        i, j, k = block.grid_indexes
        if nblocks == 2:
            if block_no == 0:
                x0 = i*dx
            elif block_no == 1:
                x0 = pi + i*dx
            x1, x2 = j*dy, k*dz
        elif nblocks == 4:
            if block_no == 0:
                x0, x1 = i*dx, j*dy
            elif block_no == 1:
                x0, x1 = pi + i*dx, j*dy
            elif block_no == 2:
                x0, x1 = i*dx, pi + j*dy
            elif block_no == 3:
                x0, x1 = pi + i*dx, pi + j*dy
            x2 = k*dz
        elif nblocks == 8:
            if block_no == 0:
                x0, x1, x2 = i*dx, j*dy, k*dz
            elif block_no == 1:
                x0, x1, x2 = pi + i*dx, j*dy, k*dz
            elif block_no == 2:
                x0, x1, x2 = i*dx, pi + j*dy, k*dz
            elif block_no == 3:
                x0, x1, x2 = pi + i*dx, pi + j*dy, k*dz
            elif block_no == 4:
                x0, x1, x2 = i*dx, j*dy, pi + k*dz
            elif block_no == 5:
                x0, x1, x2 = pi + i*dx, j*dy, pi + k*dz
            elif block_no == 6:
                x0, x1, x2 = i*dx, pi + j*dy, pi + k*dz
            elif block_no == 7:
                x0, x1, x2 = pi + i*dx, pi + j*dy, pi + k*dz
        else:
            raise ValueError("This problem is set up for nblocks 2 <= 8.")

        coordinates = [OpenSBLIEq(x, x0), OpenSBLIEq(y, x1), OpenSBLIEq(z, x2)]
        # Initial conditions as strings
        u0 = "Eq(GridVariable(u0),sin(DataObject(x0))*cos(DataObject(x1))*cos(DataObject(x2)))"
        u1 = "Eq(GridVariable(u1),-cos(DataObject(x0))*sin(DataObject(x1))*cos(DataObject(x2)))"
        u2 = "Eq(GridVariable(u2), 0.0)"
        p = "Eq(GridVariable(p), 1.0/(gama*Minf*Minf)+ (1.0/16.0) * (cos(2.0*DataObject(x0))+cos(2.0*DataObject(x1)))*(2.0 + cos(2.0*DataObject(x2))))"
        r = "Eq(GridVariable(r), gama*Minf*Minf*p)"
        # Set the conservative arrays
        rho = "Eq(DataObject(rho), r)"
        rhou0 = "Eq(DataObject(rhou0), r*u0)"
        rhou1 = "Eq(DataObject(rhou1), r*u1)"
        rhou2 = "Eq(DataObject(rhou2), r*u2)"
        rhoE = "Eq(DataObject(rhoE), p/(gama-1) + 0.5* r *(u0**2+ u1**2 + u2**2))"
        # Parse the initial conditions
        vortex_condition = [parse_expr(eq, local_dict=local_dict) for eq in [u0, u1, u2, p, r, rho, rhou0, rhou1, rhou2, rhoE]]
        init_class.add_equations(coordinates + vortex_condition)
        # Store the condition
        mb_initial_conditions[block_no] = init_class
    return mb_initial_conditions

def block_matches(multi_block):
    # Create boundaries, one for each side per dimension, so in total 6 BC's for 3D'
    nblocks = multi_block.nblocks
    ndim = multi_block.blocks[0].ndim
    # Matching conditions
    match_conditions = [None for _ in range(nblocks)]
    match_conditions = dict(list(enumerate(match_conditions)))
    if nblocks == 2:
        # block number.        # direction       # side
        match_conditions[0] = [[(1, 0, 1, False), (1, 0, 0, False)]]
        match_conditions[1] = [[(0, 0, 1, False), (0, 0, 0, False)]]
    elif nblocks == 4:
        match_conditions[0] = [[(1, 0, 1, False), (1, 0, 0, False)], [(2, 1, 1, False), (2, 1, 0, False)]]
        match_conditions[1] = [[(0, 0, 1, False), (0, 0, 0, False)], [(3, 1, 1, False), (3, 1, 0, False)]]
        match_conditions[2] = [[(3, 0, 1, False), (3, 0, 0, False)], [(0, 1, 1, False), (0, 1, 0, False)]]
        match_conditions[3] = [[(2, 0, 1, False), (2, 0, 0, False)], [(1, 1, 1, False), (1, 1, 0, False)]]
    return match_conditions

def TGV_boundaries(multi_block, halos):
    """ Create the TGV boundary conditions for multi-block decomposition of equally sized blocks."""
    nblocks = multi_block.nblocks
    ndim = multi_block.blocks[0].ndim
    mb_boundaries = [None for _ in range(nblocks)]
    mb_boundaries = dict(list(enumerate(mb_boundaries)))
    # Get the matching configuration depending on the number of blocks
    match_conditions = block_matches(multi_block)
    for block_no, block in enumerate(multi_block.blocks):
        if nblocks == 2:
            xm = SharedInterfaceBC(direction=0, side=0,  halos=halos, match=match_conditions[block_no][0][0])
            xp = SharedInterfaceBC(direction=0, side=1,  halos=halos, match=match_conditions[block_no][0][1])
            ym = PeriodicBC(direction=1, side=0,  halos=halos)
            yp = PeriodicBC(direction=1, side=1,  halos=halos)
            zm = PeriodicBC(direction=2, side=0,  halos=halos)
            zp = PeriodicBC(direction=2, side=1,  halos=halos)
        elif nblocks == 4:
            xm = SharedInterfaceBC(direction=0, side=0,  halos=halos, match=match_conditions[block_no][0][0])
            xp = SharedInterfaceBC(direction=0, side=1,  halos=halos, match=match_conditions[block_no][0][1])
            ym = SharedInterfaceBC(direction=1, side=0,  halos=halos, match=match_conditions[block_no][1][0])
            yp = SharedInterfaceBC(direction=1, side=1,  halos=halos, match=match_conditions[block_no][1][1])
            zm = PeriodicBC(direction=2, side=0,  halos=halos)
            zp = PeriodicBC(direction=2, side=1,  halos=halos)
        elif nblocks == 8:
            pass
        mb_boundaries[block_no] = [xm, xp, ym, yp, zm, zp]
    return mb_boundaries
