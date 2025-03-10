# Converts the SBLI 2d block grid into HDF5 format. For reading into OpenSBLI
# Structured mesh, sharp trailing edge version.

import numpy as np
import h5py
from opensbli import *

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


# Specify the input grid files
input_files = ["../Bl1.dat", "../Bl2.dat","../Bl3.dat"]
nblocks = len(input_files)
# Number of halo points to add on each side of each direction (default 5)
nhalo = [5, 5] # was 5, 5, 2
# Output grid file name
fname = "data.h5"
h5f = h5py.File(fname, 'w')

ndim = 2
sharp_TE = True

def read_blocks(block_data):
    total_grid_points = 0
    for block_number, block in enumerate(input_files):
        print("\n\nReading from %s." % block)
        f = open(block)
        nx,ny = map(int, f.readlines()[0].split())
        print("Nx, Ny from the input file for block %d" % block_number, nx, ny)
        f.close()
        # Read the data
        x,y,z = np.loadtxt(block, skiprows =1, unpack=True)
        x = x.reshape(nx, ny)
        y = y.reshape(nx, ny)
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

def fill_halo_coordinates(block_data, block_number):
    x, y = block_data[block_number]['x'], block_data[block_number]['y']
    # Create an array with zeros padded around the data
    shape = list(x.shape)
    new_shape = tuple([shape[i]+ 2*nhalo[i] for i in range(ndim)])
    print("Reversed shape for C-style indexing", new_shape)
    # Full arrays with halo points on the outside
    full_x = np.zeros(new_shape)
    full_y = np.zeros(new_shape)    

    # Fill out the interior data
    # Create slices of the interior points to reuse (Nz, Ny, Nz) (they have been transposed into C style)
    x_slice = np.s_[nhalo[1]:new_shape[1] -nhalo[1]]
    y_slice = np.s_[nhalo[0]:new_shape[0] -nhalo[0]]
    # Filling out the full data
    full_x[y_slice, x_slice] = x
    full_y[y_slice, x_slice] = y
    print("Block 0")
    print(full_x.shape)
    print("Block 1")
    print(block_data[1]['x'].shape)
    # Bottom right wake block
    if block_number == 0:
        # Negative x halos in block 0 are the first coordinate values in block 1 (Nz, Ny, Nx)
        ## Starting from index 0 because we have already taken off one line of coordinates from block 1
        full_x[y_slice, 4] = block_data[1]['x'][:,0]
        full_x[y_slice, 3] = block_data[1]['x'][:,1]
        full_x[y_slice, 2] = block_data[1]['x'][:,2]
        full_x[y_slice, 1] = block_data[1]['x'][:,3]
        full_x[y_slice, 0] = block_data[1]['x'][:,4]

        full_y[y_slice, 4] = block_data[1]['y'][:,0]
        full_y[y_slice, 3] = block_data[1]['y'][:,1]
        full_y[y_slice, 2] = block_data[1]['y'][:,2]
        full_y[y_slice, 1] = block_data[1]['y'][:,3]
        full_y[y_slice, 0] = block_data[1]['y'][:,4]

        # Positive x (farfield) add with constant spacing
        dx = np.abs((full_x[y_slice, -7] - full_x[y_slice, -6]))
        full_x[y_slice, -5] = full_x[y_slice, -6] + dx
        full_x[y_slice, -4] = full_x[y_slice, -5] + dx
        full_x[y_slice, -3] = full_x[y_slice, -4] + dx
        full_x[y_slice, -2] = full_x[y_slice, -3] + dx
        full_x[y_slice, -1] = full_x[y_slice, -2] + dx
        # Copy the y coordinates into the halos
        full_y[y_slice, -5] = full_y[y_slice, -6]
        full_y[y_slice, -4] = full_y[y_slice, -6]
        full_y[y_slice, -3] = full_y[y_slice, -6]
        full_y[y_slice, -2] = full_y[y_slice, -6]
        full_y[y_slice, -1] = full_y[y_slice, -6]
        # Add the extended positive y (farfield)
        full_x[-5, x_slice] = full_x[-6, x_slice]
        full_x[-4, x_slice] = full_x[-6, x_slice]
        full_x[-3, x_slice] = full_x[-6, x_slice]
        full_x[-2, x_slice] = full_x[-6, x_slice]
        full_x[-1, x_slice] = full_x[-6, x_slice]

        dy = np.abs((full_y[-7, x_slice] - full_y[-6, x_slice]))
        full_y[-5, x_slice] = full_y[-6, x_slice] + dy
        full_y[-4, x_slice] = full_y[-5, x_slice] + dy
        full_y[-3, x_slice] = full_y[-4, x_slice] + dy
        full_y[-2, x_slice] = full_y[-3, x_slice] + dy
        full_y[-1, x_slice] = full_y[-2, x_slice] + dy

        # Negative y halos in block 0 are the first coordinate values in block 2
        full_x[4, x_slice] = block_data[2]['x'][1,:]
        full_x[3, x_slice] = block_data[2]['x'][2,:]
        full_x[2, x_slice] = block_data[2]['x'][3,:]
        full_x[1, x_slice] = block_data[2]['x'][4,:]
        full_x[0, x_slice] = block_data[2]['x'][5,:]

        full_y[4, x_slice] = block_data[2]['y'][1,:]
        full_y[3, x_slice] = block_data[2]['y'][2,:]
        full_y[2, x_slice] = block_data[2]['y'][3,:]
        full_y[1, x_slice] = block_data[2]['y'][4,:]
        full_y[0, x_slice] = block_data[2]['y'][5,:]

    # Aerofoil block
    elif block_number == 1:
        # Negative x halos in block 1 are the first coordinate values in block 1
        full_x[y_slice, 4] = block_data[0]['x'][:,0]
        full_x[y_slice, 3] = block_data[0]['x'][:,1]
        full_x[y_slice, 2] = block_data[0]['x'][:,2]
        full_x[y_slice, 1] = block_data[0]['x'][:,3]
        full_x[y_slice, 0] = block_data[0]['x'][:,4]

        full_y[y_slice, 4] = block_data[0]['y'][:,0]
        full_y[y_slice, 3] = block_data[0]['y'][:,1]
        full_y[y_slice, 2] = block_data[0]['y'][:,2]
        full_y[y_slice, 1] = block_data[0]['y'][:,3]
        full_y[y_slice, 0] = block_data[0]['y'][:,4]

    #   # Positive x halos in block 1 are the first coordinate values in block 2
        full_x[y_slice, -5] = block_data[2]['x'][:,0]
        full_x[y_slice, -4] = block_data[2]['x'][:,1]
        full_x[y_slice, -3] = block_data[2]['x'][:,2]
        full_x[y_slice, -2] = block_data[2]['x'][:,3]
        full_x[y_slice, -1] = block_data[2]['x'][:,4]

        full_y[y_slice, -5] = block_data[2]['y'][:,0]
        full_y[y_slice, -4] = block_data[2]['y'][:,1]
        full_y[y_slice, -3] = block_data[2]['y'][:,2]
        full_y[y_slice, -2] = block_data[2]['y'][:,3]
        full_y[y_slice, -1] = block_data[2]['y'][:,4]

        # Below the wall coordinates ## CHECK again later
        full_x[4, x_slice] = full_x[5, x_slice] 
        full_x[3, x_slice] = full_x[5, x_slice] 
        full_x[2, x_slice] = full_x[5, x_slice] 
        full_x[1, x_slice] = full_x[5, x_slice] 
        full_x[0, x_slice] = full_x[5, x_slice] 
        dy = np.abs((full_y[6, x_slice] - full_y[5, x_slice])) 
        full_y[4, x_slice] = full_y[1, x_slice] + full_y[6, x_slice]*dy
        full_y[3, x_slice] = full_y[2, x_slice] + full_y[6, x_slice]*dy
        full_y[2, x_slice] = full_y[3, x_slice] + full_y[6, x_slice]*dy
        full_y[1, x_slice] = full_y[4, x_slice] + full_y[6, x_slice]*dy
        full_y[0, x_slice] = full_y[5, x_slice] + full_y[6, x_slice]*dy

        # Farfield, positive y
        dx = np.abs((full_x[y_slice, -7] - full_x[y_slice, -6]))
        full_x[-5, x_slice] = full_x[-6, x_slice] 
        full_x[-4, x_slice] = full_x[-6, x_slice] 
        full_x[-3, x_slice] = full_x[-6, x_slice] 
        full_x[-2, x_slice] = full_x[-6, x_slice] 
        full_x[-1, x_slice] = full_x[-6, x_slice] 
        dy = np.abs((full_y[-7, x_slice] - full_y[-6, x_slice])) 
        full_y[-5, x_slice] = full_y[-6, x_slice] + full_y[6, x_slice]*dy
        full_y[-4, x_slice] = full_y[-5, x_slice] + full_y[6, x_slice]*dy
        full_y[-3, x_slice] = full_y[-4, x_slice] + full_y[6, x_slice]*dy
        full_y[-2, x_slice] = full_y[-3, x_slice] + full_y[6, x_slice]*dy
        full_y[-1, x_slice] = full_y[-2, x_slice] + full_y[6, x_slice]*dy

    elif block_number == 2:
        # Negative x halos in block 2 are the last coordinate values in block 1
        full_x[y_slice, 4] = block_data[1]['x'][:,-1]
        full_x[y_slice, 3] = block_data[1]['x'][:,-2]
        full_x[y_slice, 2] = block_data[1]['x'][:,-3]
        full_x[y_slice, 1] = block_data[1]['x'][:,-4]
        full_x[y_slice, 0] = block_data[1]['x'][:,-5]

        full_y[y_slice, 4] = block_data[1]['y'][:,-5]
        full_y[y_slice, 3] = block_data[1]['y'][:,-4]
        full_y[y_slice, 2] = block_data[1]['y'][:,-3]
        full_y[y_slice, 1] = block_data[1]['y'][:,-2]
        full_y[y_slice, 0] = block_data[1]['y'][:,-1]

        # Add x farfield outlet halos
        dx = np.abs((full_x[y_slice, -7] - full_x[y_slice, -6]))
        full_x[y_slice, -5] = full_x[y_slice, -6] + dx
        full_x[y_slice, -4] = full_x[y_slice, -5] + dx
        full_x[y_slice, -3] = full_x[y_slice, -4] + dx
        full_x[y_slice, -2] = full_x[y_slice, -3] + dx
        full_x[y_slice, -1] = full_x[y_slice, -2] + dx
        # Copy the y coordinates into the halos
        full_y[y_slice, -5] = full_y[y_slice, -6]
        full_y[y_slice, -4] = full_y[y_slice, -6]
        full_y[y_slice, -3] = full_y[y_slice, -6]
        full_y[y_slice, -2] = full_y[y_slice, -6]
        full_y[y_slice, -1] = full_y[y_slice, -6]

        # Negative y halos in block 2 are the first coordinate values in block 0
        full_x[4, x_slice] = block_data[0]['x'][1,:]
        full_x[3, x_slice] = block_data[0]['x'][2,:]
        full_x[2, x_slice] = block_data[0]['x'][3,:]
        full_x[1, x_slice] = block_data[0]['x'][4,:]
        full_x[0, x_slice] = block_data[0]['x'][5,:]

        full_y[4, x_slice] = block_data[0]['y'][1,:]
        full_y[3, x_slice] = block_data[0]['y'][2,:]
        full_y[2, x_slice] = block_data[0]['y'][3,:]
        full_y[1, x_slice] = block_data[0]['y'][4,:]
        full_y[0, x_slice] = block_data[0]['y'][5,:]

        # Add y farfield outlet halos
        # Add the extended positive y (farfield)
        full_x[-5, x_slice] = full_x[-6, x_slice]
        full_x[-4, x_slice] = full_x[-6, x_slice]
        full_x[-3, x_slice] = full_x[-6, x_slice]
        full_x[-2, x_slice] = full_x[-6, x_slice]
        full_x[-1, x_slice] = full_x[-6, x_slice]

        dy = np.abs((full_y[-7, x_slice] - full_y[-6, x_slice]))
        full_y[-5, x_slice] = full_y[-6, x_slice] + dy
        full_y[-4, x_slice] = full_y[-5, x_slice] + dy
        full_y[-3, x_slice] = full_y[-4, x_slice] + dy
        full_y[-2, x_slice] = full_y[-3, x_slice] + dy
        full_y[-1, x_slice] = full_y[-2, x_slice] + dy

    return full_x, full_y


block_data = {}


# Loop over all of the grid points
for block_number, block in enumerate(input_files):
    # Read all of the blocks at once
    if block_number == 0:
        print("Reading all of the x,y data for the %d blocks" % nblocks)
        read_blocks(block_data)

    x, y = block_data[block_number]['x'], block_data[block_number]['y']
    # Shape without halo points (Nx, Ny, Nz) required for OPS attribute, without halos
    OPS_shape = list(x.shape)[::-1]

    # Extend the coordinates into the halo points to avoid issues with the metric calculations
    full_x, full_y = fill_halo_coordinates(block_data, block_number)

    # Make an OpenSBLI block
    b = SimulationBlock(3, block_number=block_number)
    g1 = h5f.create_group(b.blockname)
    halo = [[-i for i in nhalo], nhalo]
    apply_group_attributes(g1, b)
    block_dset_name = b.location_dataset("x0").base
    print("OpenSBLI block shape without halo points: %s" % OPS_shape)
    # Create x coordinates
    dset = g1.create_dataset('%s' % (block_dset_name), data=full_x)
    set_hdf5_metadata(dset, halos=halo, npoints=[OPS_shape[0], OPS_shape[1]], block=b)
    # Create y coordinates
    block_dset_name = b.location_dataset("x1").base
    dset = g1.create_dataset('%s' % (block_dset_name), data=full_y)
    set_hdf5_metadata(dset, halos=halo, npoints=[OPS_shape[0], OPS_shape[1]], block=b)
    
    print("Length in x for block %d:" % block_number, abs(np.amin(x) - np.amax(x)))
    print("Length in y for block %d:" % block_number, abs(np.amin(y) - np.amax(y)))

h5f.close()
