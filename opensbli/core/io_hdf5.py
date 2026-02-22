"""@brief
   @authors Satya Pramod Jammy, David J Lusher
   @contributors 
   @details
"""
from opensbli.code_generation.algorithm.common import InTheSimulation, AfterSimulationEnds, BeforeSimulationStarts
from opensbli.core.opensbliobjects import Globalvariable, ConstantObject
from opensbli.core.datatypes import Int
from opensbli.core.kernel import ConstantsToDeclare as CTD
from sympy import flatten, pprint
from collections import OrderedDict


class opensbliIO(object):
    group_number = 0

    @staticmethod
    def increase_io_group_number():
        opensbliIO.group_number += 1
        return

class iohdf5(opensbliIO):
    def __new__(cls, arrays=None, save_every=None, **kwargs):
        ret = super(iohdf5, cls).__new__(cls)
        ret.order = 0
        ret.blocknumber = 0
        ret.group_number = cls.group_number
        # Function name to call from the main program file
        ret.func_name = 'HDF5_IO_Write_%d' % ret.group_number
        cls.increase_io_group_number()
        if kwargs:
            ret.kwargs = {}
            for key in kwargs:
                if isinstance(key, str):
                    if isinstance(kwargs[key], str):
                        ret.kwargs[key.lower()] = kwargs[key].lower()
        else:
            # Default IO type is write to hdf5
            ret.kwargs = {'iotype': "write"}
            # Default write placement is the end of the simulation
        # Position of write calls in the output
        if 'position' not in ret.kwargs:
            ret.kwargs['position'] = 'end'
        ret.algorithm_place = []
        # Check if constants should be written to the HDF5 file
        if 'write_constants' in kwargs:
            cls.write_constants = kwargs['write_constants']
        else:
            cls.write_constants = True # by default always write the constants to HDF5 now
        # Constant for file write frequency
        if save_every:
            cls.save_every = ConstantObject('write_output_file', integer=True)
            cls.save_every._value = save_every
            cls.save_every.datatype = Int()
            CTD.add_constant(cls.save_every)
        else:
            cls.save_every = None
        # HDF5 timing switch
        HDF5_timer = ConstantObject('HDF5_timing', integer=True)
        HDF5_timer.value = 0 # no timing of I/O by default
        HDF5_timer.datatype = Int()
        CTD.add_constant(HDF5_timer)

        ret.get_algorithm_location()
        ret.arrays = []
        if arrays:
            ret.add_arrays(arrays)
        return ret

    def get_algorithm_location(cls):
        if cls.kwargs['iotype'] == "write":
            if cls.kwargs['position'] == "init":
                cls.algorithm_place = [BeforeSimulationStarts()]
            else:
                cls.algorithm_place += [AfterSimulationEnds()]
            if cls.save_every:
                cls.algorithm_place += [InTheSimulation(cls.save_every, initial_condition=True)]
        elif cls.kwargs['iotype'] == "read":
            cls.algorithm_place = [BeforeSimulationStarts()]
        else:
            raise ValueError("HDF5 class is missing an iotype.")
        return

    def add_arrays(cls, arrays):
        cls.arrays += flatten(arrays)
        return

    def check_datasets(cls, block):
        """ Checks if the user has added any datasets to the IO class that are not defined within the simulation."""
        simulation_dsets = [str(ar) for ar in block.block_datasets.keys()]
        io_dsets = [str(ar) for ar in cls.arrays]
        missing_dsets = [x for x in io_dsets if x not in simulation_dsets]
        if len(missing_dsets) > 0:
            raise ValueError("The dataset(s): '%s' added to the HDF5 class are not defined in the simulation code. Please check the HDF5 add_arrays input in the problem script." % str(', '.join(missing_dsets)))
        return

    def set_read_from_hdf5_arrays(cls, block):
        """ Tracks which arrays are being read in from HDF5 at the start of the simulation (External grid files, restart files)."""
        if cls.kwargs['iotype'] == "read":
            if 'filename' in cls.kwargs:
                fname = cls.kwargs['filename']
            else:
                fname = 'data.h5'
            for ar in cls.arrays:
                if str(ar) in block.block_datasets.keys():
                    dset = block.block_datasets[str(ar)]
                    dset.read_from_hdf5 = True
                    dset.input_file_name = fname
                    block.block_datasets[str(ar)] = dset
                else:
                    block.block_datasets[str(ar)] = ar
                    block.block_datasets[str(ar)].read_from_hdf5 = True
                    dset.input_file_name = fname
            # Check if restarting from a previous solution or not
            if 'restart_simulation' in cls.kwargs.keys():
                if cls.kwargs['restart_simulation']:
                    for const in CTD.constants:
                        if str(const) == 'restart':
                            const._value = 1
        return

    def set_write_to_hdf5_arrays(cls, block):
        """ Adds an attribute to DataSets which the user has chosen to write out to disk."""
        if cls.kwargs['iotype'] == "write":
           for ar in cls.arrays:
                if str(ar) in block.block_datasets.keys():
                    dset = block.block_datasets[str(ar)]
                    dset.write_to_hdf5 = True
                    block.block_datasets[str(ar)] = dset
                else:
                    block.block_datasets[str(ar)] = ar
                    block.block_datasets[str(ar)].write_to_hdf5 = True
        return

    def write_latex(cls, latex):
        string = ["HDF5 IO type %s on arrays" % (cls.kwargs['iotype'])]
        string += ["%s" % (d) for d in cls.arrays]
        latex.write_string(' '.join(string))
        return

    @property
    def opsc_code(cls):
        code = []
        # Add the constant writing function first
        if cls.kwargs['iotype'] == "write":
            # Write the HDF5 calls to a separate header file
            with open('io.h', 'r') as f:
                lines = f.readlines()
            f.close()
            # Write out a function to add constants to HDF5 files
            if len(lines) == 0:
                cls.constant_writing_opsc_code()
            # Write the body of the HDF5 output calls
            code += cls.hdf5write_opsc_code()
        elif cls.kwargs['iotype'] == "read":
            code += cls.hdf5read_opsc_code()
        else:
            raise ValueError("Cannot classify HDF5 IO object")
        return code

    def set_output_constants(cls, constants):
        cls.constants_to_write += flatten([constants])
        return

    def constant_writing_opsc_code(cls):
        """ Writes the constants defined in the simulation to the HDF5 output files."""
        code = []
        user_constants = sorted([x for x in CTD.constants if isinstance(x, ConstantObject)], key=lambda x: str(x))
        user_constants = [x for x in user_constants if not x.rational]
        # Write a separate function for constant writing
        code += ['void write_constants(const char* %s){' % "filename"]
        for c in user_constants: # Write only user input constants, not rational factors and inverses
            code += ['ops_write_const_hdf5(\"%s\", 1, \"%s\", (char*)&%s, %s);' % (c.name, c.datatype.opsc(), c.name, "filename")]
        # Constants to always write to HDF5
        code += ['ops_write_const_hdf5(\"iter\", 1, \"int\", (char*)&iter, %s);' % ("filename")]
        code += ['}\n\n']
        const_file = open('io.h', 'w')
        const_file.write('\n'.join(flatten(code)))
        const_file.close()
        return

    def hdf5write_opsc_code(cls):
        """ Generates the OPS C code to write data to disk as HDF5 files."""
        code = []
        var_name = 'name%s' % cls.blocknumber
        # Add HDF5 timers to check I/O cost
        code += ['double cpu_start0, elapsed_start0;']
        code += ['if (HDF5_timing == 1){']
        code += ['ops_timers(&cpu_start0, &elapsed_start0);\n}']
        code += ['// Writing OPS datasets']
        if "name" in cls.kwargs:
            if '.h5' in cls.kwargs["name"]:
                name = cls.kwargs["name"]
            elif '.' in cls.kwargs["name"]:
                raise ValueError("")
            elif cls.dynamic_fname:
                name = cls.kwargs["name"]
            else:
                name = cls.kwargs["name"] + '.h5'
            filename = "\"%s\"" % name
        else:
            name = "opensbli_output"
            code += ['char %s[80];' % var_name]
            if cls.dynamic_fname:
                code += ['snprintf(%s, 80, \"%s_%%06d.h5\", %s);' % (var_name, name, cls.control_parameter)]
            else:
                code += ['snprintf(%s, 80, \"%s.h5\");' % (var_name, name)]
            filename = var_name
        dataset_write = []
        if len(cls.arrays) == 0:
            raise ValueError("The IO class: {} does not have any arrays assigned to it for reading/writing from disk.".format(cls.kwargs["name"]))
        else:
            for ar in cls.arrays:
                block_name = ar.base.blockname
                dataset_write += ['ops_fetch_dat_hdf5_file(%s, %s);' % (ar, filename)]
        # Generate the block name
        code += ['ops_fetch_block_hdf5_file(%s, %s);' % (block_name, filename)] + dataset_write
        # Write constants to the HDF5 output file, once per file (not per block)
        if cls.write_constants and cls.blocknumber == 0:
            code += ['// Writing simulation constants']
            code += ['write_constants(%s);' % filename]

        code += ['if (HDF5_timing == 1){']
        code += ['double cpu_end0, elapsed_end0;']
        code += ['ops_timers(&cpu_end0, &elapsed_end0);']
        code += ['ops_printf("-----------------------------------------\\n");']
        code += ['ops_printf("Time to write HDF5 file: %s: %lf\\n", {}, elapsed_end0-elapsed_start0);'.format(filename)]
        code += ['ops_printf("-----------------------------------------\\n");']
        # code += ['fflush(stdout);\n}']
        code += ['\n}']

        # Create a function template
        fname = cls.func_name + '_' + block_name
        # Array list
        ar_list = ', '.join([str(x) for x in cls.arrays])
        header_list = ', '.join(['ops_dat& ' + str(x) for x in cls.arrays])

        if cls.dynamic_fname:
            header = ['void %s_dynamic(ops_block& %s, int iter, %s, int HDF5_timing){' % (fname, block_name, header_list)]
            function_call = ['%s_dynamic(%s, iter, %s, HDF5_timing);' % (fname, block_name, ar_list)]
        else:
            header = ['void %s(ops_block& %s, %s, int HDF5_timing){' % (fname, block_name, header_list)]
            function_call = ['%s(%s, %s, HDF5_timing);' % (fname, block_name, ar_list)]
        # Populate the function
        code = header + code 
        code += ['}\n\n']
        io_file = open('io.h', 'a')
        io_file.write('\n'.join(flatten(code)))
        io_file.close()
        return function_call

    def hdf5read_opsc_code(cls):
        """To keep the abstraction going return nothing for HDF5 OPSC code"""
        return []

    @property
    def evaluated_datasets(cls):
        evaluated = set()
        if cls.kwargs['iotype'] == "read":
            evaluated = evaluated.union(set([a.base for a in cls.arrays]))
        return evaluated

class HDF5_slice(object):
    """ Object to hold information required to generate the HDF5 reduced dimension data output option in OPS."""
    def __init__(cls, arrays, direction, index):
        cls.arrays = arrays
        cls.direction = direction
        cls.index = index
        return

class iohdf5_slices(opensbliIO):
    def __new__(cls, slices=None, save_every=None, blocknumber=None, **kwargs):
        ret = super(iohdf5_slices, cls).__new__(cls)
        ret.order = 0
        if blocknumber == None:
            ret.blocknumber = 0
        else:    
            ret.blocknumber = blocknumber
        ret.group_number = cls.group_number
        cls.increase_io_group_number()
        if kwargs:
            ret.kwargs = {}
            for key in kwargs:
                if isinstance(key, str):
                    if isinstance(kwargs[key], str):
                        ret.kwargs[key.lower()] = kwargs[key].lower()
        else:
            # Default IO type is write to hdf5
            ret.kwargs = {'iotype': "write"}
            # Default write placement is the end of the simulation
        # Position of write calls in the output
        if 'position' not in ret.kwargs:
            ret.kwargs['position'] = 'end'
        ret.algorithm_place = []
        # Check if constants should be written to the HDF5 file
        if 'write_constants' in kwargs:
            cls.write_constants = kwargs['write_constants']
        else:
            cls.write_constants = False # by default always write the constants to HDF5 now
        # Constant for file write frequency
        if save_every:
            cls.save_every = ConstantObject('write_slices', integer=True)
            cls.save_every._value = save_every
            cls.save_every.datatype = Int()
            CTD.add_constant(cls.save_every)
        else:
            cls.save_every = None
        ret.get_algorithm_location()
        ret.arrays = []
        ret.slices = []
        if slices:
            ret.add_slices(slices)
        return ret

    def get_algorithm_location(cls):
        if cls.kwargs['iotype'] == "write": # slices only written in time
            cls.algorithm_place = [InTheSimulation(cls.save_every)]        
        elif cls.kwargs['iotype'] == "init":
            cls.algorithm_place = [BeforeSimulationStarts()]
        else:
            raise ValueError("Please specify the type of IO class (position in the algorithm) for the HDF5 slice writer.")
        return

    def add_slices(cls, slices):
        # add some checks on the input here
        for slc in slices:
            ars, direction, index = flatten(slc[0]), slc[1], slc[2]
            cls.slices += [HDF5_slice(ars, direction, index)]
            cls.arrays += ars
        cls.arrays = list(OrderedDict.fromkeys(cls.arrays)) # Remove duplicates
        return

    def check_datasets(cls, block):
        """ Checks if the user has added any datasets to the IO class that are not defined within the simulation."""
        simulation_dsets = [str(ar) for ar in block.block_datasets.keys()]
        io_dsets = [str(ar) for ar in cls.arrays]
        missing_dsets = [x for x in io_dsets if x not in simulation_dsets]
        if len(missing_dsets) > 0:
            raise ValueError("The dataset(s): '%s' added to the HDF5 class are not defined in the simulation code. Please check the HDF5 add_arrays input in the problem script." % str(', '.join(missing_dsets)))
        return

    def write_latex(cls, latex):
        string = ["HDF5 slice writer IO type %s on arrays" % (cls.kwargs['iotype'])]
        string += ["%s" % (d) for d in cls.arrays]
        latex.write_string(' '.join(string))
        return

    @property
    def opsc_code(cls):
        code = []
        if cls.kwargs['iotype'] == "write":
            code += cls.hdf5write_opsc_code()
        elif cls.kwargs['iotype'] == 'init':
            code += cls.hdf5write_opsc_code(init=True)
        else:
            raise ValueError("Cant classify HDF5io")
        return code

    def set_read_from_hdf5_arrays(cls, block):
        return # hdf5 slicing has no read functionality, output option only

    def set_write_to_hdf5_arrays(cls, block):
        return # for consistency with regular hdf5 io class

    def hdf5write_opsc_code(cls, init=False):
        var_name = 'slice_name%s' % cls.blocknumber
        code = []
        if "name" in cls.kwargs:
            if '.h5' in cls.kwargs["name"]:
                name = cls.kwargs["name"]
            elif '.' in cls.kwargs["name"]:
                raise ValueError("")
            else:
                name = cls.kwargs["name"] + '.h5'
            if cls.dynamic_fname:
                raise ValueError("dynamic fname not allowed ")
            filename = "\"%s\"" % name
        else:
            name = '0'
            code += ['char %s[80];' % var_name]
            if cls.dynamic_fname:
                code += ['snprintf(%s, 80, \"%%d\", %s);' % (var_name, cls.control_parameter)]
            else:
                code += ['snprintf(%s, 80, \"%s\");' % (var_name, "0")]
            filename = var_name
        

        dataset_write = []
        # syntax: ops_write_plane_group_hdf5({{2, 50}}, str, {{rho_B0, rhou0_B0, rhou1_B0, rhou2_B0, rhoE_B0}});
        for slc in cls.slices:
            dataset_write += ['ops_write_plane_group_hdf5({{%s, %s}}, %s, {{%s}});' % (str(slc.direction), str(slc.index), var_name, str(', '.join([str(x)+'_B%d' % cls.blocknumber for x in slc.arrays])))]
        code += dataset_write
        
        # Write constants to the HDF5 output file, once per file (not per block)
        if cls.write_constants and cls.blocknumber == 0:
            # Generate the OPS API calls
            user_constants = [x for x in CTD.constants if isinstance(x, ConstantObject)]
            user_constants = [x for x in user_constants if not x.rational]
            for c in user_constants: # Write only user input constants, not rational factors and inverses
                code += ['ops_write_const_hdf5(\"%s\", 1, \"%s\", (char*)&%s, %s);' % (c.name, c.datatype.opsc(), c.name, filename)]
            # Constants to always write to HDF5
            code += ['ops_write_const_hdf5(\"iter\", 1, \"int\", (char*)&iter, %s);' % (filename)]
        return code

    def hdf5read_opsc_code(cls):
        """To keep the abstraction going return nothing for HDF5 OPSC code"""
        return []

    @property
    def evaluated_datasets(cls):
        evaluated = set() # hdf5 slice output has no read functionality planned
        return evaluated


class IoGroup():
    pass
