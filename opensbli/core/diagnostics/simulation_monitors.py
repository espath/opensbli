"""@brief Simulation monitors to report data values during the simulation.
   @author David J. Lusher
   @details Can be used to generate a large number of flow samples to construct a time signal.
"""
from opensbli.core.datatypes import SimulationDataType
from opensbli.core.block import SimulationBlock
from opensbli.core.opensbliobjects import DataSet, DataObject
from opensbli.multiblock.blockcollection import MultiBlock

class Monitor(object):
    def __init__(self, block, flow_var, probe_loc, numbering):
        if isinstance(flow_var, DataSet):
            self.flow_var = flow_var
        elif isinstance(flow_var, DataObject):
            self.flow_var = str(DataObject) + '_B%d' % block.blocknumber
        elif isinstance(flow_var, str):
            if '_B%d' % block.blocknumber in flow_var:
                self.flow_var = flow_var
            else:
                self.flow_var = flow_var + '_B%d' % block.blocknumber
        else:
            raise ValueError("Unknown simulation monitor input: {}".format(flow_var))
        self.probe_loc = probe_loc
        self.block = block
        if len(probe_loc) != block.ndim:
            raise ValueError("The number of probe indices (i,j,k) must match the dimensions of the problem (block.ndim).")
        self.probe_no = numbering
        return

class ScalarMonitor(object):
    def __init__(self, block, scalar, output=None):
        self.scalar = scalar
        self.block = block
        if output == None:
            self.output = scalar + '_B%d' % block.blocknumber + '_out'
        else:
            self.output = output
        return

class SimulationMonitor(object):
    def __init__(self, arrays, probe_locations, blocks, print_frequency=100, OPS_V2=True, fp_precision=12, NaN_check=None, output_file='output.log'):
        """ Class to enable access of dataset values during the simulation.
        :arg list arrays: A list of DataSets to monitor during the simulation.
        :arg list probe_locations: A list of tuples giving the (i,j,k) grid index location of the probe.
        :arg object block: An OpenSBLI simulation block.
        :arg int print_frequency: The iteration frequency at which to print the output.
        :arg int fp_precision: The number of decimal places to format the output.
        :arg bool NaNcheck: Adds an optional call to ops_NaNcheck().
        :arg str output_file: The option to write the output directly to a log file, otherwise defaults to stdout."""

        # Check number of probes equals the number of input arrays
        if len(arrays) != len(probe_locations):
            raise ValueError("The number of arrays must equal the number of probe locations.")
        # OPS formatting
        if OPS_V2:
            self.ops_headers = {'input': "const ACC<%s> &%s", 'output': 'ACC<%s> &%s', 'inout': 'ACC<%s> &%s'}
        else:
            self.ops_headers = {'input': "const %s *%s", 'output': '%s *%s', 'inout': '%s *%s'}
        # Create the monitor objects
        self.initial_setup(blocks, arrays, probe_locations, output_file)
        self.OPS_V2 = OPS_V2
        self.components = []
        self.frequency = print_frequency
        self.fp_precision = fp_precision
        self.filename = 'reductions.h'
        self.NaN_check = NaN_check
        if hasattr(SimulationDataType.dtype, 'opsc'):
            self.dtype = SimulationDataType.opsc()
        else:
            self.dtype = 'double'  # Default to double
        return

    def initial_setup(self, blocks, arrays, probe_locations, output_file):
        # Check whether monitoring an array or a single value from a reduction already performed
        if isinstance(blocks, SimulationBlock):
            self.blocks = [blocks]
            arrays, probe_locations = [arrays], [probe_locations]
        elif isinstance(blocks, MultiBlock):
            self.blocks = [b for b in blocks.blocks]
        else:
            raise ValueError("SimulationMonitor input: pass either a SimulationBlock or MultiBlock class.")
        print(arrays)
        print(probe_locations)
        # exit()
        self.ndim = self.blocks[0].ndim
        self.nblocks = len(self.blocks)
        self.output_files = ['block%d_' % b.blocknumber + output_file for b in self.blocks]
        self.array_monitors = []
        self.scalar_monitors = []
        # Loop over all the blocks in the problem
        for block_id, b in enumerate(self.blocks):
            # if not any(x is None for x in probe_locations):
            self.array_monitors += [Monitor(b, var, loc, index) for index, (var, loc) in enumerate(zip(arrays[block_id], probe_locations[block_id])) if isinstance(loc, tuple)]
            self.scalar_monitors += [ScalarMonitor(b, var, output=loc) for index, (var, loc) in enumerate(zip(arrays[block_id], probe_locations[block_id])) if not isinstance(loc, tuple)]
            # Check if the scalar monitors have to be scaled before printing
            for SM in self.scalar_monitors:
                if SM.output == 'residual':
                    grid_factor = '*'.join([str(b.ranges[i][1]) for i in range(b.ndim)])
                    SM.output = 'sqrt(%s_out/(%s))' % (SM.scalar, grid_factor)
        return

    def add_components(self, components):
        """ Adds the given components to the simulation monitoring
        :param components: the components to be added to the timers, this can be a list or an individual component
        :return: None. """
        if isinstance(components, list):
            self.components += components
        else:
            self.components += [components]
        return

    def write_latex(self, latex):
        latex.write_string("Simulation monitoring start\\\\\n")
        for c in self.components:
            c.write_latex(latex)
        latex.write_string("Simulation monitoring end\\\\\n")
        return

    @property
    def opsc_code(self):
        # Loop over all the blocks in the problem
        code = []
        for b in self.blocks:
            code += self.opsc_start(b)
            code += self.opsc_middle(b)
            code += self.opsc_end(b)
        return code

    def initial_declarations(self, M):
        """ Initial declarations for the OPS data access."""
        name, number = str(M.flow_var), M.probe_no
        declarations = ["// Monitoring of %s" % name]
        declarations += ["ops_reduction reduce_%d_%s = ops_decl_reduction_handle(sizeof(%s), \"%s\", \"reduction_%d_%s\");" % (number, name, self.dtype, self.dtype, number, name)]
        declarations += ["%s %s_%d_output = 0.0;" % (self.dtype, name, number)]
        return declarations

    def reduction_range(self, M, idx):
        """ Sets up the grid locations and iteration range."""
        name, number = str(M.flow_var), M.probe_no
        bn = M.block.blocknumber
        if self.ndim == 1:
            i = M.probe_loc
            return ["int i%d%d = %s;" % (bn, idx, str(i))] + ["int monitor_range_%d_%s[] = {i%d%d, i%d%d+1};" % (number, name, bn, idx, bn, idx)]
        elif self.ndim == 2:
            i, j = M.probe_loc
            return ["int i%d%d = %s, j%d%d = %s;" % (bn, idx, str(i), bn, idx, str(j))] + ["int monitor_range_%d_%s[] = {i%d%d, i%d%d+1, j%d%d, j%d%d+1};" % (number, name, bn, idx, bn, idx, bn, idx, bn, idx)]
        elif self.ndim == 3:
            i, j, k = M.probe_loc
            return ["int i%d%d = %s, j%d%d = %s, k%d%d = %s;" % (bn, idx, str(i), bn, idx, str(j), bn, idx, str(k))] + ["int monitor_range_%d_%s[] = {i%d%d, i%d%d+1, j%d%d, j%d%d+1, k%d%d, k%d%d+1};" % (number, name, bn, idx, bn, idx, bn, idx, bn, idx, bn, idx, bn, idx)]

    def par_loop_declaration(self, M):
        """ Defines the parallel loop templates for the reductions."""
        name, number = str(M.flow_var), M.probe_no
        if self.ndim == 1:
            stencil_name = 'stencil_%d_00_1' % M.block.blocknumber
        elif self.ndim == 2:
            stencil_name = 'stencil_%d_00_00_2' % M.block.blocknumber
        else:
            stencil_name = 'stencil_%d_00_00_00_3' % M.block.blocknumber
        output_code = ["ops_par_loop(monitor_%d_%s, \"Reduction %s_%d\", %s, %d, monitor_range_%d_%s," % (number, name, name, number, M.block.blockname, self.ndim, number, name)]
        output_code += ["ops_arg_dat(%s, 1, %s, \"%s\", OPS_READ)," % (name, stencil_name, self.dtype)]
        output_code += ["ops_arg_reduce(reduce_%d_%s, 1, \"%s\", OPS_INC));" % (number, name, self.dtype)]
        return output_code

    def generate_result(self, M):
        """ Obtains the result of the OPS reduction."""
        name, number = str(M.flow_var), M.probe_no
        return ["ops_reduction_result(reduce_%d_%s, &%s_%d_output);\n" % (number, name, name, number)]

    def generate_kernel_code(self, M):
        """ Generates the kernel definitions to be added to the reductions header file."""
        name, number = str(M.flow_var), M.probe_no
        indices = ','.join(['0' for _ in range(self.ndim)])
        if self.OPS_V2:
            code = ["void monitor_%d_%s(const ACC<%s> &%s, %s *reduce_%d_%s){\n" % (number, name, self.dtype, name, self.dtype, number, name)]
            code += ["*reduce_%d_%s = %s(%s);\n}" % (number, name, name, indices)] + ["\n\n"]
        else:
            code = ["void monitor_%d_%s(const %s *%s, %s *reduce_%d_%s){\n" % (number, name, self.dtype, name, self.dtype, number, name)]
            code += ["*reduce_%d_%s = %s[OPS_ACC0(%s)];\n}" % (number, name, name, indices)] + ["\n\n"]
        return code

    def add_NaN_check(self, block):
        if self.NaN_check == None:
            return ['ops_NaNcheck(%s);' % str(block.block_datasets['rho_B%d' % block.blocknumber])] # default to density if no array given
        else:
            if "_B%d" % block.blocknumber in self.NaN_check:
                name = self.NaN_check
            else:
                name = '%s_B%d' % (self.NaN_check, block.blocknumber)
            return ['ops_NaNcheck(%s);' % (name)]

    @property
    def write_reductions_file(self):
        """ Creates a new header file to define the reduction kernels."""
        f = open(self.filename, 'w')
        f.write("#ifndef REDUCTIONS_H\n")
        f.write("#define REDUCTIONS_H\n")
        for M in self.array_monitors:
            f.write(''.join(self.generate_kernel_code(M)))
        f.write('#endif\n')
        f.close()
        return

    def format_output(self, block):
        """ Controls the printing format for the output."""
        AM, SM = [A for A in self.array_monitors if A.block.blocknumber == block.blocknumber], [S for S in self.scalar_monitors if S.block.blocknumber == block.blocknumber]
        placeholders = ', '.join(["%d"] + ["%%.%de" % self.fp_precision for _ in range(len(AM)+1)] + ["%%.%de" % self.fp_precision for _ in range(len(SM))])
        iterations = ['iter+1', '%s' % 'simulation_time']
        variables = ["%s_%d_output" % (str(M.flow_var), M.probe_no) for M in AM]
        # Scalar variables
        variables += ["%s" % str(M.output) for M in SM]
        # Normalise mean quantities by the number of iterations
        for i, var in enumerate(variables):
            if 'mean' in var:
                variables[i] = '%s/%s' % (var, '(iter+1)')
        variables = ', '.join(iterations + variables)
        # Check if the output should be written directly to a log file
        if len(self.output_files) > 0:
            # output_print = ["ops_fprintf(f%d, \"%s\\n\", %s);" % (block.blocknumber, placeholders, variables)] + ["fflush(f%d);" % block.blocknumber]
            output_print = ["ops_fprintf(f%d, \"%s\\n\", %s);" % (block.blocknumber, placeholders, variables)]
        else:
            if self.nblocks > 1:
                raise ValueError("For multi-block problems, please specify a output file name for the simulation monitor output.")
            # output_print = ["ops_printf(\"%s\\n\", %s);" % (placeholders, variables)] + ["fflush(stdout);"]
            output_print = ["ops_printf(\"%s\\n\", %s);" % (placeholders, variables)]
        return ["// Write the output values"] + output_print

    def generate_reduction_loops(self, block):
        """ Creates a block of code for each flow variable being monitored."""
        output_code = []
        AM = [A for A in self.array_monitors if A.block.blocknumber == block.blocknumber]
        for index, M in enumerate(AM):
            output_code += self.initial_declarations(M)
            output_code += self.reduction_range(M, index)
            output_code += self.par_loop_declaration(M)
            output_code += self.generate_result(M)
        return output_code

    def initial_print(self, block):
        """ Prints the headers at the top of the output once at the start."""
        headers = ['Iteration', 'Time'] + ['%s%s' % (str(M.flow_var), str(M.probe_loc)) for M in self.array_monitors if len(self.array_monitors) > 0 and M.block.blocknumber == block.blocknumber]
        headers += ['%s' % (str(M.scalar)) for M in self.scalar_monitors if len(self.scalar_monitors) > 0 and M.block.blocknumber == block.blocknumber]
        headers = ', '.join(headers)
        if len(self.output_files) > 0:
            output_code = ["if (iter == 0){\nops_fprintf(f%d, \"%s\\n\");}" % (block.blocknumber, headers)]
        else:
            output_code = ["if (iter == 0){\nops_printf(\"%s\\n\");}" % (headers)]
        return output_code

    def opsc_start(self, block):
        starting_code = ["// Data access for simulation monitoring"]
        if self.NaN_check:
            starting_code += self.add_NaN_check(block)
        starting_code += self.initial_print(block)
        return starting_code

    def opsc_middle(self, block):
        middle_code = []
        middle_code += self.generate_reduction_loops(block)
        return middle_code

    def opsc_end(self, block):
        end_code = self.format_output(block)
        self.write_reductions_file
        return end_code
