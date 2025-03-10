from opensbli.postprocess.plot_functions import *
from opensbli.core.block import SimulationBlock
import numpy as np
import h5py
import os
""" Extracts a 2D slice from a 3D datafile and writes it back out as an HDF5 file for Tecplot. DJL: 06/22."""

# Step 1 read a file
inPath = './'
outPath = './file_output/'

try:
	os.mkdir(outPath)
except FileExistsError:
	pass
remove = False

if remove:
	# Remove old files
	for f in glob.glob(outPath + '/opensbli_output_*.h5'):
		os.remove(f)

already_processed = sorted(glob.glob(outPath + '/opensbli_output_*.h5'))
already_processed = [re.findall("\d+", s)[0].lstrip('0') for s in already_processed]
# Extract data
PP = OpenSBLIPreProcess()
fnames, iters = PP.find_files(inPath)

blocks = 3
block_numbers = [i for i in range(blocks)]

# Slice location
#check_Mach = True
check_Mach = False

for i, f in enumerate(fnames):
	if check_Mach:
		f = fnames[-1]
	if iters[i] not in already_processed:
		# Read the last block last, to make sure the datafile has all of the required blocks already written to it
		for block_number in block_numbers[::-1]:
			data_2D = {}
			grid_file, block_name, dsets, shape = PP.read_block('./data.h5', block_number)
			# grid_file.close()
			halos = np.abs(grid_file[block_name][dsets[0]].attrs['d_p'])
			nhalo = halos[0]
			zloc = int(shape[0]/2.0)
			#zloc = 45
			#zloc = 0
			surface = np.s_[zloc+nhalo,:,:]
			add_grid = True
			if add_grid:
				data_2D['x0_B%d' % block_number] = PP.read_full_dset(grid_file, block_name, 'x0_B%d' % block_number, partial_slice=surface)
				data_2D['x1_B%d' % block_number] = PP.read_full_dset(grid_file, block_name, 'x1_B%d' % block_number, partial_slice=surface)
			# Read the flow data
			data_file, block_name, dsets, shape = PP.read_block(f, block_number)
			if check_Mach:
				gamma = 1.4
				rho = PP.read_full_dset(data_file, block_name, 'rho_B0')
				u = PP.read_full_dset(data_file, block_name, 'u0_B0')
				v = PP.read_full_dset(data_file, block_name, 'u1_B0')
				w = PP.read_full_dset(data_file, block_name, 'u2_B0')
				Et = PP.read_full_dset(data_file, block_name, 'Et_B0')
				p = (gamma-1)*(Et - 0.5*(u**2 + v**2 + w**2))*rho
				c = np.sqrt(gamma*p/rho)
				Mach = np.sqrt(u**2 + v**2 + w**2) / c
				full_x = PP.read_full_dset(grid_file, block_name, 'x0_B0')
				full_y = PP.read_full_dset(grid_file, block_name, 'x1_B0')
				full_z = PP.read_full_dset(grid_file, block_name, 'x2_B0')

				argmax = np.argmax(Mach)
				Mach = np.max(Mach)
				print("Maximum Mach number is: ", Mach)
				print("Maximum is located at (x,y,z):", np.ravel(full_x)[argmax], np.ravel(full_y)[argmax], np.ravel(full_z)[argmax])
				print(full_z[zloc,10,10])
				exit()

			# Extract all the variables at this slice
			for dset_name in dsets:
				# print(dset_name)
				data_2D[dset_name] = PP.read_full_dset(data_file, block_name, dset_name, partial_slice=surface)

			# Write to a new 2D HDF5 file
			if block_number == block_numbers[-1]:
				h5f = h5py.File(outPath + 'opensbli_output_%s.h5' % iters[i], 'w')
			b = SimulationBlock(2, block_number=block_number)
			g1 = h5f.create_group(b.blockname)
			nhalos = [5, 5]
			halo = [[-i for i in nhalos], nhalos]
			# apply_group_attributes(g1, b)
			for dset_name in data_2D.keys():
				dset = g1.create_dataset('%s' % (dset_name), data=data_2D[dset_name])
				# set_hdf5_metadata(dset, halos=halo, npoints=[OPS_shape[0], OPS_shape[1]], block=b)
		h5f.close()
