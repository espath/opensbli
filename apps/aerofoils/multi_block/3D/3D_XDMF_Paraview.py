import numpy as np
import h5py
import argparse
import os
from opensbli.postprocess.plot_functions import *

inPath = './'
outPath = './paraview_files/'

try:
    os.mkdir(outPath)
except FileExistsError:
    pass



remove = False
if remove:
    # Remove old files
    for f in glob.glob(outPath + '/opensbli_output_*.h5'):
        os.remove(f)
    try:
        os.remove(outPath + 'read.xdmf')
    except:
        pass


class XDMF_Paraview(object):
    def __init__(self, input_files, float_precision):
        self.nblocks = 3
        self.selected_dsets = None
        self.selected_dsets = ['rhou0', 'rhou2']
        self.Mach = True
        self.input_files = input_files
        self.remove_halos = True
        self.tstep_counter = 0
        self.reduce_2D = False
        self.velocities = True
        self.side_panel = False
        self.surface = True
        if float_precision == 'single':
            self.precision = 4
            self.dtype = 'f'
        else:
            self.precision = 8
            self.dtype = 'double'
        # Reduce size of data (for 3 blocks)
        self.reduce_data = True
        self.main(input_files)
        return

    def extract_surface(self, dataset):
        zloc = 1
        surface = np.s_[zloc+5,:,:]
        return dataset[surface]

    def read_dset(self, dataset, halos, reduced_shape):
        size = dataset.shape
        read_start = [abs(d) for d in halos]
        read_end = [s-abs(d) for d, s in zip(halos, (reduced_shape[0]+10, reduced_shape[1]+10, reduced_shape[2]+10))]
        if self.remove_halos:
            if len(read_end) == 1:
                read_data = dataset[read_start[0]:read_end[0]]
            elif len(read_end) == 2:
                self.ndim = 2
                read_data = dataset[read_start[0]:read_end[0], read_start[1]:read_end[1]]
            else:
                self.ndim = 3
                read_data = dataset[read_start[0]:read_end[0], read_start[1]:read_end[1], read_start[2]:read_end[2]]
        else:
            read_data = dataset
            self.ndim = len(read_end)
        return read_data

    def add_Mach_number(self, f, nhalos, blocknumber, block_name, reduced_shape):
        # Q vector
        rho = self.read_dset(f[block_name]['rho_B%d' % blocknumber], nhalos, reduced_shape)
        rhou0 = self.read_dset(f[block_name]['rhou0_B%d' % blocknumber], nhalos, reduced_shape)
        rhou1 = self.read_dset(f[block_name]['rhou1_B%d' % blocknumber], nhalos, reduced_shape)
        rhou2 = self.read_dset(f[block_name]['rhou2_B%d' % blocknumber], nhalos, reduced_shape)
        rhoE = self.read_dset(f[block_name]['rhoE_B%d' % blocknumber], nhalos, reduced_shape)

        u, v, w = rhou0/rho, rhou1/rho, rhou2/rho
        ## Get the constants here later
        gamma = 1.4
        p = (gamma - 1.)*(rhoE - 0.5*rho*(u**2 + v**2 + w**2))
        a = np.sqrt(gamma*p/rho)
        Mach = np.sqrt(u**2 + v**2 + w**2) / a
        return Mach

    def reduce_data_size(self, input_shape, blocknumber):
        # Clip the data to remove outer parts and reduce memory usage
        Nz, Ny, Nx = input_shape
        print(Nz, Ny, Nx)
        if self.side_panel:
            factors = [[1, 1, 1], [1, 1, 1], [1, 1, 1]]
            factor = factors[blocknumber]
            Nz_new, Ny_new, Nx_new = 1, int(factor[1]*Ny), int(factor[2]*Nx)
        elif self.surface:
            factors = [[1, 1, 0.01], [1, 1, 1], [1, 1, 0.01]]
            factor = factors[blocknumber]
            Nz_new, Ny_new, Nx_new = int(factor[0]*Nz), 1, int(factor[2]*Nx)
        else:
            factors = [[1, 0.05, 0.1], [1, 0.05, 1], [1, 0.05, 0.1]]
            factor = factors[blocknumber]
            Nz_new, Ny_new, Nx_new = int(factor[0]*Nz), int(factor[1]*Ny), int(factor[2]*Nx)
        print("Reducing data on Block number: {} from {} to {}.".format(blocknumber, (Nz, Ny, Nx), (Nz_new, Ny_new, Nx_new)))
        print(Nz_new, Ny_new, Nx_new)
        return (Nz_new, Ny_new, Nx_new)


    def main(self, input_files):
        # Strip the halos from the input dataset
        # Loop over all times
        output_names = []
        simulation_times = []
        for i, fname in enumerate(input_files):
            base_name = fname.split('.h5')[0]
            if self.side_panel:
                HDF5_without_halos = base_name + '_side_pp.h5'
            elif self.surface:
                HDF5_without_halos = base_name + '_surface_pp.h5'
            else:
                HDF5_without_halos = base_name + '_pp.h5'
            HDF5_without_halos = HDF5_without_halos.split('./')[-1]
            # Names of the files processed, without halos now
            output_names.append(HDF5_without_halos)
            output_file = h5py.File(outPath + HDF5_without_halos, 'w')
            f = h5py.File(fname, 'r')
            simulation_times.append(f['simulation_time'][0])
            # Loop over all the blocks in the file
            for blocknumber in range(self.nblocks):
                block_name = 'opensbliblock0%d' % blocknumber
                dsets = [x for x in f[block_name].keys()]
                print("Found %d datasets: %s, with dimensions: %s" % (len(dsets), dsets, f[block_name][dsets[0]].shape[::-1]))
                ndim = len(f[block_name][dsets[0]].shape)
                try:
                    nhalos = np.abs(f[block_name][dsets[0]].attrs['d_m'])
                except:
                    nhalos = [5 for _ in range(ndim)]
                input_shape = f[block_name][dsets[0]].shape
                input_shape = (input_shape[0] - 2*nhalos[0], input_shape[1] - 2*nhalos[1], input_shape[2] - 2*nhalos[2])
                reduced_shape = self.reduce_data_size(input_shape, blocknumber)
                if self.velocities:
                    rho = self.read_dset(f[block_name]['rho_B%d' % blocknumber], nhalos, reduced_shape)
                if self.selected_dsets is not None:
                    for dset in self.selected_dsets:
                        print("Processing dataset: {}".format(dset))
                        dset = dset + '_B%d' % blocknumber
                        data = self.read_dset(f[block_name][dset], nhalos, reduced_shape)
                        if self.velocities:
                            data = data/rho
                        output_file.create_dataset("%s" % (dset), data=data, dtype=self.dtype)
                    if self.Mach:
                        data = self.add_Mach_number(f, nhalos, blocknumber, block_name, reduced_shape)
                        output_file.create_dataset("Mach_B%d" % blocknumber, dtype=self.dtype, data=data)

                    if i == 0: # add grid coordinates
                        d = h5py.File('data.h5', 'r')
                        for dim in range(self.ndim):
                            dset = 'x%d_B%d' % (dim, blocknumber)
                            data = self.read_dset(d[block_name][dset], nhalos, reduced_shape)
                            output_file.create_dataset("%s" % (dset), data=data, dtype=self.dtype)
                else:
                    for dset in dsets:
                        data = self.read_dset(f[block_name][dset], nhalos, reduced_shape)
                        output_file.create_dataset("%s" % (dset), data=data, dtype=self.dtype)
        # Write to file
        if self.side_panel:
            self.write_xdmf(output_names, simulation_times, 'read_side.xdmf')
        elif self.surface:
            self.write_xdmf(output_names, simulation_times, 'read_surface.xdmf')
        else:
            self.write_xdmf(output_names, simulation_times, 'read.xdmf')
        output_file.close()
        return


    def mesh_info(self, ndim, array_size):
        self.coord_names = ["X", "Y", "Z"]
        self.node_size_str = " ".join([str(s) for s in array_size])
        return

    def write_xdmf(self, output_files, simulation_times, output_name):
        print("Writing to XDMF file: {}".format(output_name))
        # Header
        code = """<?xml version=\"1.0\" ?>\n<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n<Xdmf Version=\"2.0\">\n<Domain>\n\n"""
        # Loop over all input files
        for iteration, fname in enumerate(self.input_files):

            # First grid entry for multi-block collection
            code += self.block_container(simulation_times[iteration])
            code += """<Grid Name=\"block%d\" GridType=\"Uniform\">""" % 0
            for blocknumber in range(self.nblocks):
                block_name = 'opensbliblock0%d' % blocknumber
                opensbli_file = h5py.File(outPath + output_files[iteration], 'r')
                size =  opensbli_file['rhou0_B%s' % blocknumber].shape
                self.mesh_info(self.ndim, size)
                # Write code
                if blocknumber > 0:
                    code += """<Grid Name=\"block%d\" GridType=\"Uniform\">""" % blocknumber
                # Mesh topology
                code += self.topology()
                # Coordinate domain info
                code += self.coordinate_read(blocknumber, output_files[0])
                # Datasets in the file
                code += ''.join([self.attribute_node(k, output_files[iteration], blocknumber) for k in opensbli_file.keys()])+"""</Grid>\n\n\n"""
                opensbli_file.close()

            # End multiblock
            code += "</Grid>\n"
            # End time
            code += "</Grid>\n"
        # Finished all time-steps
        # End domain and Xdmf
        code += """</Domain>\n</Xdmf>\n"""
        with open(outPath + output_name, 'w') as fout:
            fout.write(code)
        return

    def block_container(self, time):
        code = "<Grid Name='t%d' GridType='Collection' CollectionType='Temporal'>\n" % self.tstep_counter 
        code += "<Grid Name='tstep%d' GridType='Collection' CollectionType='Spatial'>\n" % self.tstep_counter
        code += "<Time Type=\"Single\" Value=\"%s\" />\n\n\n" % time
        self.tstep_counter += 1
        return code

    def topology(self):
        return "<Topology TopologyType=\"%dDSMesh\" NumberOfElements=\"%s\"/>\n\n" % (self.ndim, self.node_size_str)

    @property
    def dataitem_node(self):
        return """<DataItem Dimensions=\"%s\" NumberType=\"Float\" Precision=\"%s\" Format=\"HDF\">\n %s:/%s\n</DataItem>\n"""

    def attribute_node(self, attribute, data_file, blocknumber):
        if 'B%d' % blocknumber in attribute:
            name = attribute.split('_B%d' % blocknumber)[0]
            attr = """<Attribute Name=\"%s\" AttributeType=\"Scalar\" Center=\"Node\">\n""" % (name)
            attr += self.dataitem_node%(self.node_size_str, self.precision, data_file, attribute) + "</Attribute>\n"
        else:
            attr = ''
        return attr

    def coordinate_read(self, blocknumber, data_file):
        coordinate = ["x%d_B%d" % (d, blocknumber) for d in range(self.ndim)]
        geom = "<Geometry GeometryType=\"%s\">\n\n"%("_".join(self.coord_names[0:self.ndim]))
        reading = geom
        for d in range(self.ndim):
            reading += self.dataitem_node%(self.node_size_str, self.precision, data_file, coordinate[d])
        reading += "</Geometry>\n\n\n"
        return reading



# Extract data

def atoi(text):
        return int(text) if text.isdigit() else text

def natural_keys(text):
    return [ atoi(c) for c in re.split(r'(\d+)', text) ]


# Extract data
PP = OpenSBLIPreProcess()
fnames, iters = PP.find_files(inPath)

fnames.sort(key=natural_keys)

iter_no = '0'

fname = [f for f in fnames if iter_no in f]
print(fname)

only_last = False

if only_last:
    fnames = [fnames[-1]]
else:
    fnames = fnames

XD = XDMF_Paraview(fname, 'double')
