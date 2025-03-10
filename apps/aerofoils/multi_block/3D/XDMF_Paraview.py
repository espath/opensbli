import numpy as np
import h5py
import argparse
import os
from opensbli.postprocess.plot_functions import *
import gc
import time
# import numexpr as ne

inPath = './'
outPath = './'

def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    return [ atoi(c) for c in re.split(r'(\d+)', text) ]


class XDMF_Paraview(object):
    def __init__(self, input_files):
        # self.nblocks = 3
        self.input_files = input_files
        self.remove_halos = True
        self.tstep_counter = 0
        self.surface_view = True
        #self.surface_view = False
        self.ndim = 2
        self.reduce_2D = True
        self.main(input_files)
        return

    def calculate_derived_quantities(self, f, snapshot, block_no, block_name):
        # Q vector
        start = time.time()
        rho = f[block_name][snapshot]['rho_B%d' % block_no][()]
        rhou0 = f[block_name][snapshot]['rhou0_B%d' % block_no][()]
        rhou1 = f[block_name][snapshot]['rhou1_B%d' % block_no][()]
        rhou2 = f[block_name][snapshot]['rhou2_B%d' % block_no][()]
        rhoE = f[block_name][snapshot]['rhoE_B%d' % block_no][()]
        end = time.time()
        print("Time taken for reading Q vector: {:.2f}s".format(end - start))

        start = time.time()
        u, v, w = rhou0/rho, rhou1/rho, rhou2/rho
        ## Get the constants here later
        gamma = 1.4
        p = (gamma - 1.)*(rhoE - 0.5*rho*(u**2 + v**2 + w**2))
        a = np.sqrt(gamma*p/rho)
        Mach = np.sqrt(u**2 + v**2 + w**2) / a
        end = time.time()
        #print("Time taken for calculating Mach number: {:.2f}s".format(end - start))
        return Mach, p

    def main(self, input_files):
        # Process all the files
        output_names = []
        simulation_times = []
        # Calculate any new quantities using the Q vector and write it into the same file
        self.blocks_to_process = {}
        for i, fname in enumerate(input_files):
            print("Processing file: %s" % fname)
            base_name = fname.split('.h5')[0]
            # Read the slice file
            f = h5py.File(fname, 'r+')
            self.blocks_to_process[base_name] = [int(str(n)[-1]) for n in f.keys()] # assumes less than 10 blocks
            # Loop over all the blocks in the file
            for blocknumber in self.blocks_to_process[base_name]:
                block_name = 'opensbliblock0%d' % blocknumber
                time_snapshots = [x for x in f[block_name].keys() if x != "0"]
                time_snapshots.sort(key=natural_keys)
                if only_last:
                    time_snapshots = [time_snapshots[-1]]
                print("Found %d time snapshots." % (len(time_snapshots)))
                # Loop over the time instances
                for time_index, time_group in enumerate(time_snapshots):
                    existing = list(f[block_name][time_group].keys())
                    # Compute derived quantities: Mach, pressure, ...
                    if 'Mach_B%d' % blocknumber not in existing and 'p_B%d' not in existing:
                        print("Processing: %s." % time_group)
                        # Calculate derived and write to disk, unless they have already been processed
                        Mach, p = self.calculate_derived_quantities(f, time_group, blocknumber, block_name)
                        f.create_dataset("/%s/%s/Mach_B%d" % (block_name, time_group, blocknumber), data=Mach)
                        f.create_dataset("/%s/%s/p_B%d" % (block_name, time_group, blocknumber), data=p)
                        del Mach, p
                        gc.collect()
                    else:
                        print("Skipping: %s." % time_group)

            print(self.blocks_to_process)
            # Write to file
            start = time.time()
            self.write_xdmf(base_name)
            end = time.time()
            print("Time taken for writing output file: {:.2f}s".format(end - start))
            f.close()
        return


    def mesh_info(self, array_size):
        if not self.reduce_2D:
            self.coord_names = ["X", "Y", "Z"]
        else:
            self.coord_names = ["X", "Y"] # Paraview only works with X-Y, not X-Z for surface
        self.node_size_str = " ".join([str(s) for s in array_size])
        return

    def write_xdmf(self, base_name):
        # Loop over all input files
        output_name = 'read_%s.xdmf' % base_name
        fname = base_name + '.h5'
        print("Writing to XDMF file: {}".format(output_name))
        # Header
        code = """<?xml version=\"1.0\" ?>\n<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n<Xdmf Version=\"2.0\">\n<Domain>\n\n"""
        # Load the file and find the datasets
        print(fname)
        f = h5py.File(fname, 'r')
        # Find all the times (assume all blocks have the same time snapshots)
        time_snapshots = [x for x in f['opensbliblock0%d' % self.blocks_to_process[base_name][0]].keys() if x != "0"]
        time_snapshots.sort(key=natural_keys)

        if only_last:
            time_snapshots = [time_snapshots[-1]]

        # temporary, update with actual simulation time later from constants
        dt = 5e-5
        simulation_times = [dt*float(x) for x in time_snapshots]

        # Loop over all times
        for time_index, time in enumerate(simulation_times):
            print("Processing: %s" % time_snapshots[time_index])
            # First grid entry for multi-block collection
            code += self.block_container(time)
            first_block = int(min(self.blocks_to_process[base_name]))
            code += """<Grid Name=\"block%d\" GridType=\"Uniform\">""" % first_block
            for blocknumber in self.blocks_to_process[base_name]:
                block_name = 'opensbliblock0%d' % blocknumber
                dsets = list(f[block_name][time_snapshots[time_index]].keys())
                size =  f[block_name][time_snapshots[time_index]][dsets[0]].shape
                self.mesh_info(size)
                # Write code
                if blocknumber > first_block:
                    code += """<Grid Name=\"block%d\" GridType=\"Uniform\">""" % blocknumber
                # Mesh topology
                code += self.topology()
                # Coordinate domain info
                code += self.coordinate_read(block_name, base_name + ".h5", list(f[block_name]["0"].keys()))
                # Datasets in the file
                dset_path = base_name + ".h5" + ":/" + block_name + "/" + time_snapshots[time_index]
                code += ''.join([self.attribute_node(k,  blocknumber, dset_path) for k in dsets])+"""</Grid>\n\n\n"""
            # End multiblock
            code += "</Grid>\n"
            # End time
            code += "</Grid>\n"
        f.close()
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
        return """<DataItem Dimensions=\"%s\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n %s/%s\n</DataItem>\n"""

    def attribute_node(self, dset_name, blocknumber, dset_path):
        if 'B%d' % blocknumber in dset_name:
            name = dset_name.split('_B%d' % blocknumber)[0]
            attr = """<Attribute Name=\"%s\" AttributeType=\"Scalar\" Center=\"Node\">\n""" % (name)
            attr += self.dataitem_node%(self.node_size_str, dset_path, dset_name) + "</Attribute>\n"
        else:
            attr = ''
        return attr

    def coordinate_read(self, block_name, data_file, coordinates):
        geom = "<Geometry GeometryType=\"%s\">\n\n"%("_".join(self.coord_names))
        reading = geom
        # Grid is always written to "0" timesnap shot
        path = data_file + ":/" + block_name + "/" + '0' #+ "/"
        for d in range(self.ndim):
            reading += self.dataitem_node%(self.node_size_str, path, coordinates[d])
        reading += "</Geometry>\n\n\n"
        return reading

# Extract data
PP = OpenSBLIPreProcess()

#fname = ['J20.h5']
fname = ['K0.h5']

only_last = False

XD = XDMF_Paraview(fname)
