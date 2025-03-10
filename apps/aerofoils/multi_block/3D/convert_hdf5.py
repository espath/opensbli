#] Converts the SBLI 2d block grid into HDF5 format. For reading into OpenSBLI
# Structured mesh, sharp trailing edge version.

import numpy as np
import h5py
from opensbli import *
import matplotlib.pylab as plt

def apply_group_attributes(group, block):
    group.attrs.create("dims", [block.ndim], dtype="int32")
    group.attrs.create("ops_type", u"ops_block",dtype="S8")
    group.attrs.create("index", [block.blocknumber], dtype="int32")
    return

def set_hdf5_metadata(dset, halos, npoints, block):
    """ Function to set hdf5 metadata required by OPS to a dataset. """
    d_m = halos[0]
    d_p = halos[1]

    dset.attrs.create("d_p", d_p, dtype="int32")
    dset.attrs.create("d_m", d_m, dtype="int32")
    dset.attrs.create("dim", [1], dtype="int32")
    dset.attrs.create("ops_type", u"ops_dat",dtype="S10")
    dset.attrs.create("block_index", [block.blocknumber], dtype="int32")
    dset.attrs.create("base", [0 for i in range(block.ndim)], dtype="int32")
    dset.attrs.create("type", u"double",dtype="S15")
    dset.attrs.create("block", u"%s" % block.blockname,dtype="S25")
    dset.attrs.create("size", npoints, dtype="int32")
    return

def output_hdf5(array, array_name, halos, npoints, block):
    """ Creates an HDF5 file for reading in data to a simulation, 
    sets the metadata required by the OPS library. """
    if not isinstance(array, list):
        array = [array]
    if not isinstance(array_name, list):
        array_name = [array_name]
    assert len(array) == len(array_name)
    with h5py.File('data.h5', 'w') as hf:
        # Create a group
        if (isinstance(block, MultiBlock)):
            all_blocks = block.blocks
        else:
            all_blocks = [block]
        for b in all_blocks:
            g1 = hf.create_group(b.blockname)
            # Loop over all the dataset inputs and write to the hdf5 file
            for ar, name in zip(array, array_name):
                g1.attrs.create("dims", [b.ndim], dtype="int32")
                g1.attrs.create("ops_type", u"ops_block",dtype="S9")
                g1.attrs.create("index", [b.blocknumber], dtype="int32")
                block_dset_name = b.location_dataset(name).base
                dset = g1.create_dataset('%s' % (block_dset_name), data=ar)
                set_hdf5_metadata(dset, halos, npoints, b)
    return

def create_z_coordinates(block_data, full_z):
    inner_shape = block_data[block_number]['x'].shape[::-1]
    full_shape = full_z.shape[::-1]
    # Uniform spacing in the span
    dz = Lz / (float(nz))
    print("Grid spacing in z is: {:.5f}".format(dz))
    z_coordinates = [k*dz for k in range(full_z.shape[0] - 2*nhalo)]
    zm = [z_coordinates[0] - k*dz for k in reversed(range(1, nhalo+1))]
    zp = [z_coordinates[-1] + k*dz for k in range(1, nhalo+1)]
    z_coordinates = np.around(np.array(zm + z_coordinates + zp), decimals=10)
    # print("Z coordinates including halo points are:", z_coordinates)
    for k in range(full_z.shape[0]):
        z = np.full((full_shape[0], full_shape[1]), z_coordinates[k])
        full_z[k, :, :] = np.transpose(z)
    return full_z


def fill_halo_coordinates(block_data, block_number):
    x, y = block_data[block_number]['x'], block_data[block_number]['y']
    # Create an array with zeros padded around the data
    shape = [nz] + list(x.shape) 
    new_shape = tuple([shape[i]+ 2*nhalo for i in range(ndim)])
    print("Reversed shape for C-style indexing", new_shape)
    # Full arrays with halo points on the outside
    full_x = np.zeros(new_shape)
    full_y = np.zeros(new_shape)
    full_z = np.zeros(new_shape) 

    # Fill out the interior data
    # Create slices of the interior points to reuse (Nz, Ny, Nz) (they have been transposed into C style)
    x_slice = np.s_[nhalo:new_shape[2] -nhalo]
    y_slice = np.s_[nhalo:new_shape[1] -nhalo]
    # z_slice = np.s_[nhalo:new_shape[0] -nhalo]

    # Filling out the full data
    for k in range(new_shape[0]):        
        full_x[k, y_slice, x_slice] = x
        full_y[k, y_slice, x_slice] = y

    full_z = create_z_coordinates(block_data, full_z)

    for k in range(new_shape[0]):
        # Bottom right wake block
        if block_number == 0:
            # Negative x halos in block 0 are the first coordinate values in block 1 (Nz, Ny, Nx)
            ## Starting from index 0 because we have already taken off one line of coordinates from block 1
            full_x[k, y_slice, 4] = block_data[1]['x'][:,0]
            full_x[k, y_slice, 3] = block_data[1]['x'][:,1]
            full_x[k, y_slice, 2] = block_data[1]['x'][:,2]
            full_x[k, y_slice, 1] = block_data[1]['x'][:,3]
            full_x[k, y_slice, 0] = block_data[1]['x'][:,4]

            full_y[k, y_slice, 4] = block_data[1]['y'][:,0]
            full_y[k, y_slice, 3] = block_data[1]['y'][:,1]
            full_y[k, y_slice, 2] = block_data[1]['y'][:,2]
            full_y[k, y_slice, 1] = block_data[1]['y'][:,3]
            full_y[k, y_slice, 0] = block_data[1]['y'][:,4]

            # Positive x (farfield) add with constant spacing
            dx = np.abs((full_x[k, y_slice, -7] - full_x[k, y_slice, -6]))
            full_x[k, y_slice, -5] = full_x[k, y_slice, -6] + dx
            full_x[k, y_slice, -4] = full_x[k, y_slice, -5] + dx
            full_x[k, y_slice, -3] = full_x[k, y_slice, -4] + dx
            full_x[k, y_slice, -2] = full_x[k, y_slice, -3] + dx
            full_x[k, y_slice, -1] = full_x[k, y_slice, -2] + dx
            # Copy the y coordinates into the halos
            full_y[k, y_slice, -5] = full_y[k, y_slice, -6]
            full_y[k, y_slice, -4] = full_y[k, y_slice, -6]
            full_y[k, y_slice, -3] = full_y[k, y_slice, -6]
            full_y[k, y_slice, -2] = full_y[k, y_slice, -6]
            full_y[k, y_slice, -1] = full_y[k, y_slice, -6]
            # Add the extended positive y (farfield)
            full_x[k, -5, x_slice] = full_x[k, -6, x_slice]
            full_x[k, -4, x_slice] = full_x[k, -6, x_slice]
            full_x[k, -3, x_slice] = full_x[k, -6, x_slice]
            full_x[k, -2, x_slice] = full_x[k, -6, x_slice]
            full_x[k, -1, x_slice] = full_x[k, -6, x_slice]

            dy = np.abs((full_y[k, -7, x_slice] - full_y[k, -6, x_slice]))
            full_y[k, -5, x_slice] = full_y[k, -6, x_slice] + dy
            full_y[k, -4, x_slice] = full_y[k, -5, x_slice] + dy
            full_y[k, -3, x_slice] = full_y[k, -4, x_slice] + dy
            full_y[k, -2, x_slice] = full_y[k, -3, x_slice] + dy
            full_y[k, -1, x_slice] = full_y[k, -2, x_slice] + dy

            # Negative y halos in block 0 are the first coordinate values in block 2
            full_x[k, 4, x_slice] = block_data[2]['x'][1,:]
            full_x[k, 3, x_slice] = block_data[2]['x'][2,:]
            full_x[k, 2, x_slice] = block_data[2]['x'][3,:]
            full_x[k, 1, x_slice] = block_data[2]['x'][4,:]
            full_x[k, 0, x_slice] = block_data[2]['x'][5,:]

            full_y[k, 4, x_slice] = block_data[2]['y'][1,:]
            full_y[k, 3, x_slice] = block_data[2]['y'][2,:]
            full_y[k, 2, x_slice] = block_data[2]['y'][3,:]
            full_y[k, 1, x_slice] = block_data[2]['y'][4,:]
            full_y[k, 0, x_slice] = block_data[2]['y'][5,:]

        # Aerofoil block
        elif block_number == 1:
            # Negative x halos in block 1 are the first coordinate values in block 1
            full_x[k, y_slice, 4] = block_data[0]['x'][:,0]
            full_x[k, y_slice, 3] = block_data[0]['x'][:,1]
            full_x[k, y_slice, 2] = block_data[0]['x'][:,2]
            full_x[k, y_slice, 1] = block_data[0]['x'][:,3]
            full_x[k, y_slice, 0] = block_data[0]['x'][:,4]

            full_y[k, y_slice, 4] = block_data[0]['y'][:,0]
            full_y[k, y_slice, 3] = block_data[0]['y'][:,1]
            full_y[k, y_slice, 2] = block_data[0]['y'][:,2]
            full_y[k, y_slice, 1] = block_data[0]['y'][:,3]
            full_y[k, y_slice, 0] = block_data[0]['y'][:,4]

        #   # Positive x halos in block 1 are the first coordinate values in block 2
            full_x[k, y_slice, -5] = block_data[2]['x'][:,0]
            full_x[k, y_slice, -4] = block_data[2]['x'][:,1]
            full_x[k, y_slice, -3] = block_data[2]['x'][:,2]
            full_x[k, y_slice, -2] = block_data[2]['x'][:,3]
            full_x[k, y_slice, -1] = block_data[2]['x'][:,4]

            full_y[k, y_slice, -5] = block_data[2]['y'][:,0]
            full_y[k, y_slice, -4] = block_data[2]['y'][:,1]
            full_y[k, y_slice, -3] = block_data[2]['y'][:,2]
            full_y[k, y_slice, -2] = block_data[2]['y'][:,3]
            full_y[k, y_slice, -1] = block_data[2]['y'][:,4]

            # Below the wall coordinates ## CHECK again later
            full_x[k, 4, x_slice] = full_x[k, 5, x_slice] 
            full_x[k, 3, x_slice] = full_x[k, 5, x_slice] 
            full_x[k, 2, x_slice] = full_x[k, 5, x_slice] 
            full_x[k, 1, x_slice] = full_x[k, 5, x_slice] 
            full_x[k, 0, x_slice] = full_x[k, 5, x_slice] 
            dy = np.abs((full_y[k, 6, x_slice] - full_y[k, 5, x_slice])) 
            full_y[k, 4, x_slice] = full_y[k, 1, x_slice] + full_y[k, 6, x_slice]*dy
            full_y[k, 3, x_slice] = full_y[k, 2, x_slice] + full_y[k, 6, x_slice]*dy
            full_y[k, 2, x_slice] = full_y[k, 3, x_slice] + full_y[k, 6, x_slice]*dy
            full_y[k, 1, x_slice] = full_y[k, 4, x_slice] + full_y[k, 6, x_slice]*dy
            full_y[k, 0, x_slice] = full_y[k, 5, x_slice] + full_y[k, 6, x_slice]*dy

            # Farfield, positive y
            dx = np.abs((full_x[k, y_slice, -7] - full_x[k, y_slice, -6]))
            full_x[k, -5, x_slice] = full_x[k, -6, x_slice] 
            full_x[k, -4, x_slice] = full_x[k, -6, x_slice] 
            full_x[k, -3, x_slice] = full_x[k, -6, x_slice] 
            full_x[k, -2, x_slice] = full_x[k, -6, x_slice] 
            full_x[k, -1, x_slice] = full_x[k, -6, x_slice] 
            dy = np.abs((full_y[k, -7, x_slice] - full_y[k, -6, x_slice])) 
            full_y[k, -5, x_slice] = full_y[k, -6, x_slice] + full_y[k, 6, x_slice]*dy
            full_y[k, -4, x_slice] = full_y[k, -5, x_slice] + full_y[k, 6, x_slice]*dy
            full_y[k, -3, x_slice] = full_y[k, -4, x_slice] + full_y[k, 6, x_slice]*dy
            full_y[k, -2, x_slice] = full_y[k, -3, x_slice] + full_y[k, 6, x_slice]*dy
            full_y[k, -1, x_slice] = full_y[k, -2, x_slice] + full_y[k, 6, x_slice]*dy

        elif block_number == 2:
            # Negative x halos in block 2 are the last coordinate values in block 1
            full_x[k, y_slice, 4] = block_data[1]['x'][:,-1]
            full_x[k, y_slice, 3] = block_data[1]['x'][:,-2]
            full_x[k, y_slice, 2] = block_data[1]['x'][:,-3]
            full_x[k, y_slice, 1] = block_data[1]['x'][:,-4]
            full_x[k, y_slice, 0] = block_data[1]['x'][:,-5]

            full_y[k, y_slice, 4] = block_data[1]['y'][:,-5]
            full_y[k, y_slice, 3] = block_data[1]['y'][:,-4]
            full_y[k, y_slice, 2] = block_data[1]['y'][:,-3]
            full_y[k, y_slice, 1] = block_data[1]['y'][:,-2]
            full_y[k, y_slice, 0] = block_data[1]['y'][:,-1]

            # Add x farfield outlet halos
            dx = np.abs((full_x[k, y_slice, -7] - full_x[k, y_slice, -6]))
            full_x[k, y_slice, -5] = full_x[k, y_slice, -6] + dx
            full_x[k, y_slice, -4] = full_x[k, y_slice, -5] + dx
            full_x[k, y_slice, -3] = full_x[k, y_slice, -4] + dx
            full_x[k, y_slice, -2] = full_x[k, y_slice, -3] + dx
            full_x[k, y_slice, -1] = full_x[k, y_slice, -2] + dx
            # Copy the y coordinates into the halos
            full_y[k, y_slice, -5] = full_y[k, y_slice, -6]
            full_y[k, y_slice, -4] = full_y[k, y_slice, -6]
            full_y[k, y_slice, -3] = full_y[k, y_slice, -6]
            full_y[k, y_slice, -2] = full_y[k, y_slice, -6]
            full_y[k, y_slice, -1] = full_y[k, y_slice, -6]

            # Negative y halos in block 2 are the first coordinate values in block 0
            full_x[k, 4, x_slice] = block_data[0]['x'][1,:]
            full_x[k, 3, x_slice] = block_data[0]['x'][2,:]
            full_x[k, 2, x_slice] = block_data[0]['x'][3,:]
            full_x[k, 1, x_slice] = block_data[0]['x'][4,:]
            full_x[k, 0, x_slice] = block_data[0]['x'][5,:]

            full_y[k, 4, x_slice] = block_data[0]['y'][1,:]
            full_y[k, 3, x_slice] = block_data[0]['y'][2,:]
            full_y[k, 2, x_slice] = block_data[0]['y'][3,:]
            full_y[k, 1, x_slice] = block_data[0]['y'][4,:]
            full_y[k, 0, x_slice] = block_data[0]['y'][5,:]

            # Add y farfield outlet halos
            # Add the extended positive y (farfield)
            full_x[k, -5, x_slice] = full_x[k, -6, x_slice]
            full_x[k, -4, x_slice] = full_x[k, -6, x_slice]
            full_x[k, -3, x_slice] = full_x[k, -6, x_slice]
            full_x[k, -2, x_slice] = full_x[k, -6, x_slice]
            full_x[k, -1, x_slice] = full_x[k, -6, x_slice]

            dy = np.abs((full_y[k, -7, x_slice] - full_y[k, -6, x_slice]))
            full_y[k, -5, x_slice] = full_y[k, -6, x_slice] + dy
            full_y[k, -4, x_slice] = full_y[k, -5, x_slice] + dy
            full_y[k, -3, x_slice] = full_y[k, -4, x_slice] + dy
            full_y[k, -2, x_slice] = full_y[k, -3, x_slice] + dy
            full_y[k, -1, x_slice] = full_y[k, -2, x_slice] + dy

    return full_x,  full_y, full_z

def plot_mesh(x_reduce, y_reduce, block_number, x, y):

    colors = ['g', 'b', 'r']
    nlines = 100
    total_grid_points = 0
    grid_line_thickness = 0.2
    edge_thickness = 0.6

    if block_number == 1:
        # Plot the airfoil
        plt.plot(x[:,0], y[:,0], color='k', lw=edge_thickness)
        # Plot the grid lines
        plt.plot(x_reduce[:,1:nlines], y_reduce[:,1:nlines], color=colors[block_number], lw=grid_line_thickness)
        plt.plot(x_reduce[:,1:nlines].T, y_reduce[:,1:nlines].T, color=colors[block_number], lw=grid_line_thickness)
        # Plot the edges on this block
        plt.plot(x_reduce[0,1:nlines].T, y_reduce[0,1:nlines].T, color='k', lw=edge_thickness)
        plt.plot(x_reduce[-1,1:nlines].T, y_reduce[-1,1:nlines].T, color='k', lw=edge_thickness)
        # Add trip location
        xPS, xSS = x[0:550,0], x[550:,0]
        yPS, ySS = y[0:550,0], y[550:,0]
        ps = np.argmin(np.abs(xPS - 0.1))
        ss = np.argmin(np.abs(xSS - 0.1))
        print("Indices: {}, {}".format(ps, ss))
        print("Trip location at: {}".format(x[ps,0]))
        print("Trip location at: {}".format(x[ss,0]))

        plt.scatter(xPS[ps+1], yPS[ps+1]-0.0075, s=15, marker=(4, 0, 45-10), color='magenta', zorder=10000)
        plt.scatter(xSS[ss], ySS[ss]+0.0075, s=15, marker=(4, 0, 45+5), color='magenta', zorder=10000, label='Trip Location')
    else:
        # Plot the grid lines
        plt.plot(x_reduce[0:,1:nlines], y_reduce[0:,1:nlines], color=colors[block_number], lw=grid_line_thickness)
        plt.plot(x_reduce[0:,1:nlines].T, y_reduce[0:,1:nlines].T, color=colors[block_number], lw=grid_line_thickness)
        # Plot the block edges
        if block_number == 2:
            # plt.plot(x_reduce[0,1:nlines], y_reduce[0,1:nlines], color='k', lw=edge_thickness)
            plt.plot(x_reduce[:,0], y_reduce[:,0], color='k', lw=edge_thickness)
            

    # Annotate block number
    ax = plt.gca()
    if block_number == 0:
        labelx, labely = 0.825, 0.03
        text_to_write = 'Wake \nBlock: %d' % block_number
    elif block_number == 1:
        labelx, labely = 0.05, 0.86
        text_to_write = 'Aerofoil \nBlock: %d' % block_number
    else:
        labelx, labely = 0.825, 0.86
        text_to_write = 'Wake \nBlock: %d' % block_number

    t = ax.text(labelx, labely, text_to_write,
    verticalalignment='bottom', horizontalalignment='left',
    transform=ax.transAxes,
    color='black', fontsize=fontsize-2)
    t.set_bbox(dict(facecolor='white', alpha=0.85, edgecolor='black', linewidth=0.75))
    return


def read_blocks(block_data):
    # Plotting mesh
    fig, ax = plt.subplots()
    # Loop over the blocks
    total_grid_points = 0
    for block_number, block in enumerate(input_files):
        print("\n\nReading from %s." % block)
        f = open(block)
        nx,ny = map(int, f.readlines()[0].split())
        print("Nx, Ny from the input file for block %d" % block_number, nx, ny)
        f.close()
        # Read the data
        x,z,y = np.loadtxt(block, skiprows =1, unpack=True)
        # For grid use
        x = x.reshape(nx, ny)
        y = y.reshape(nx, ny)
        # Reduced grid lines for plotting
        nplotx = 8 # plot only every 10th grid line
        nploty = 10 # plot only every 10th grid line
        nplot = 7
        x_reduce = x[::nplot, ::nplot]
        y_reduce = y[::nplot, ::nplot]
        
        # Plot the aerofoil
        plot_mesh(x_reduce, y_reduce, block_number, x, y)
        # exit()
        print("Min/max x = ({}, {})".format(np.min(np.abs(x)), np.max(np.abs(x))))
        print("Min/max y = ({}, {})".format(np.min(np.abs(y)), np.max(np.abs(y))))

        if block_number == 2:
            ax.axis('equal')
            ax.set_xlim([-0.125, 1.6])
            ax.set_ylim([-0.05, 0.05])
            leg = ax.legend(bbox_to_anchor=(0.425, 0.14), prop={'size': fontsize-2}, facecolor='white', framealpha=0.85, edgecolor='black')
            leg.get_frame().set_linewidth(0.75)
            plt.xlabel(r"$x'$")
            plt.ylabel(r"$y'$")
            plt.savefig('CRM_Mesh.pdf', bbox_inches='tight')
    # Sharp trailing edge -> take away one point at the start and end of the grid in x, from block 2
        if block_number == 1:
            if sharp_TE:
                print("Taking off 2 columns in x direction for the sharp trailing edge.")
                x = x[1:-1,:]
                y = y[1:-1,:]
        shape = list(x.shape)
        total = shape[0]*shape[1]
        print("Block %d has %e grid points." % (block_number, int(total)))
        total_grid_points +=  total
        print("Original 2D shape: %s" % shape)
        # Transpose into OPS format
        x, y = np.transpose(x), np.transpose(y)
        # Save this block data to the dictionary
        block_data[block_number] = {'x': x, 'y': y}
    # Finish
    print("Total grid points: %g" % total_grid_points)
    return


block_data = {}

fontsize = 16
plt.rcParams.update({'font.size': fontsize})
plt.rcParams.update({"text.usetex": True,"font.family": "serif","font.serif": ["Palatino"]})

    # Specify the input grid files
input_path = 'generated-grid/'
input_files = [input_path + "crm.c-grid.block-0.dat", input_path + "crm.c-grid.block-1.dat", input_path + "crm.c-grid.block-2.dat"]
# Number of halo points to add on each side of each direction (default 5)
nhalo = 5
ndim = 3
nblocks = len(input_files)
# Output grid file name
fname = "data.h5"
h5f = h5py.File(fname, 'w')

# Number of points in the periodic span.
nz = 500
# Grid spacing
# Span width
Lz = 0.5

sharp_TE = True

total_grid_points = 0

# Loop over all of the grid points
for block_number, block in enumerate(input_files):
    # Read all of the blocks at once
    if block_number == 0:
        print("Reading all of the x,y data for the %d blocks" % nblocks)
        read_blocks(block_data)

    x, y = block_data[block_number]['x'], block_data[block_number]['y']
    # Extend the coordinates into the halo points to avoid issues with the metric calculations
    full_x, full_y, full_z = fill_halo_coordinates(block_data, block_number)

    # Shape without halo points (Nx, Ny, Nz) required for OPS attribute, without halos
    OPS_shape = list(x.shape)[::-1] + [nz]
    # Make an OpenSBLI block
    b = SimulationBlock(3, block_number=block_number)
    g1 = h5f.create_group(b.blockname)
    halo = [[-nhalo, -nhalo, -nhalo], [nhalo, nhalo, nhalo]]  
    apply_group_attributes(g1, b)
    block_dset_name = b.location_dataset("x0").base
    print("OpenSBLI block shape without halo points: %s" % OPS_shape)
    # Create x coordinates
    dset = g1.create_dataset('%s' % (block_dset_name), data=full_x)
    set_hdf5_metadata(dset, halos=halo, npoints=[OPS_shape[0], OPS_shape[1], OPS_shape[2]], block=b)
    # Create y coordinates
    block_dset_name = b.location_dataset("x1").base
    dset = g1.create_dataset('%s' % (block_dset_name), data=full_y)
    set_hdf5_metadata(dset, halos=halo, npoints=[OPS_shape[0], OPS_shape[1], OPS_shape[2]], block=b)
    # Create z coordinates
    block_dset_name = b.location_dataset("x2").base
    dset = g1.create_dataset('%s' % (block_dset_name), data=full_z)
    set_hdf5_metadata(dset, halos=halo, npoints=[OPS_shape[0], OPS_shape[1], OPS_shape[2]], block=b)
    
h5f.close()
