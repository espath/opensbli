"""@brief Algorithm generation
   @author Satya Pramod Jammy
   @contributors Satya Pramod Jammy and David J Lusher
   @details Implements the tree-based structure with the attribute (components) controls the
   depth of a node."""

# from sympy.core.compatibility import is_sequence
from sympy.printing.precedence import precedence
from sympy.utilities.iterables import is_sequence
from sympy.printing.ccode import C99CodePrinter
# from sympy.printing.c import C99CodePrinter
from sympy.core.relational import Equality
from opensbli.core.opensbliobjects import ConstantObject, ConstantIndexed, Constant, DataSetBase, GroupedPiecewise, ReductionVariable, DataObject, DataSet, WhileLoop, ForLoop
from sympy import Symbol, flatten, Rational, nsimplify
from opensbli.core.grid import GridVariable
from opensbli.core.datatypes import SimulationDataType
from opensbli.core.datatypes import FloatC, Double
from sympy import Pow, Idx, pprint, count_ops, Piecewise
import os
import logging
from collections import OrderedDict
LOG = logging.getLogger(__name__)
BUILD_DIR = os.getcwd()


class RationalCounter():

    # Counter for the kernels
    def __init__(self):
        self.name = 'rc%d'
        self.rational_counter = 0
        self.existing = OrderedDict()

    @property
    def increase_rational_counter(self):
        self.rational_counter = self.rational_counter + 1
        return

    def get_next_rational_constant(self, numerical_value):
        from opensbli.core.kernel import ConstantsToDeclare
        # name = self.name % self.rational_counter
        name = self.name
        self.increase_rational_counter
        ret = ConstantObject(name, rational=True) # Don't write rational constants to the HDF files
        ret.value = numerical_value
        self.existing[numerical_value] = ret
        ConstantsToDeclare.add_constant(ret)
        return ret


rc = RationalCounter()


class OPSCCodePrinter(C99CodePrinter):

    """ Prints OPSC code. """
    dataset_accs_dictionary = {}
    settings_opsc = {'rational': True, 'kernel': False, 'order': 'none'}

    def __init__(self, settings={}):
        """ Initialise the code printer. """
        self.settings_opsc = settings
        # Mixed precision settings
        if 'arrays_to_cast' in settings.keys():
            if len(settings['arrays_to_cast']) > 0:
                self.cast_precision = True
            else:
                self.cast_precision = False
        else:
            self.cast_precision = False
        C99CodePrinter.__init__(self, settings={'order':'none'})

    def return_args(self, expr):
        """ Retrieves the arguments of an input expression. Used for modifying the function call with custom code printers."""
        settings = {'kernel': True}
        for const in expr.atoms(Idx).union(expr.atoms(ConstantObject)):
            if hasattr(const, "main_file"):
                settings = {'kernel': False}
                break
        args = map(lambda x: ccode(x, settings), expr.args)
        args = [x for x in args]
        result = ','.join(args)
        return result

    def _print_ReductionVariable(self, expr):
        return '*%s' % str(expr)

    def _print_Rational(self, expr):
        """ Settings: if rational is True then rational numbers are printed as they are.
        Otherwise optimisations will be performed for rational constants that are evaluated
        at the start of the program to reduce divisions."""
        expr = nsimplify(expr)
        p, q = int(expr.p), int(expr.q)
        if isinstance(SimulationDataType.dtype(), FloatC):
            return '(%d.0f/%d.0f)' % (p, q)
        else:
            return '(%d.0/%d.0)' % (p, q)

    def _print_Mod(self, expr):
        """ All modulus functions are expressed as fmod currently and no integer values."""
        result = 'fmod(%s)' % self.return_args(expr)
        return result

    def _print_Float(self, expr):
        if isinstance(SimulationDataType.dtype(), FloatC):
            return str(float(expr)) + 'f'
        else:
            return super()._print_Float(expr)

    def _print_sin(self, expr):
        if isinstance(SimulationDataType.dtype(), FloatC):
            return 'sinf(%s)' % self.return_args(expr)
        else:
            return 'sin(%s)' % self.return_args(expr)

    def _print_cos(self, expr):
        if isinstance(SimulationDataType.dtype(), FloatC):
            return 'cosf(%s)' % self.return_args(expr)
        else:
            return 'cos(%s)' % self.return_args(expr)

    def _print_tan(self, expr):
        if isinstance(SimulationDataType.dtype(), FloatC):
            return 'tanf(%s)' % self.return_args(expr)
        else:
            return 'tan(%s)' % self.return_args(expr)

    def _print_sinh(self, expr):
        if isinstance(SimulationDataType.dtype(), FloatC):
            return 'sinhf(%s)' % self.return_args(expr)
        else:
            return 'sinh(%s)' % self.return_args(expr)

    def _print_cosh(self, expr):
        if isinstance(SimulationDataType.dtype(), FloatC):
            return 'coshf(%s)' % self.return_args(expr)
        else:
            return 'cosh(%s)' % self.return_args(expr)

    def _print_tanh(self, expr):
        if isinstance(SimulationDataType.dtype(), FloatC):
            return 'tanhf(%s)' % self.return_args(expr)
        else:
            return 'tanh(%s)' % self.return_args(expr)

    def _print_GridVariable(self, expr):
        """Prints the grid variable"""
        return str(expr)


    def _print_Abs(self, expr):
        if isinstance(SimulationDataType.dtype(), FloatC):
            return 'fabsf(%s)' % self.return_args(expr)
        else:
            return 'fabs(%s)' % self.return_args(expr)

    def _print_Max(self, expr):
        """MAXIMUM of the arguments, can handle any number of arguments:
        Max(a,b,c,d) is written as max(a, max(max(b,c),d))"""
        if isinstance(SimulationDataType.dtype(), FloatC):
            func_call = 'fmaxf(%s, %s)'
        else:
            func_call = 'fmax(%s, %s)'

        nargs = len(expr.args)
        args_code = [self._print(a) for a in expr.args]
        for i in range(nargs-1):
            # Max of the last 2 arguments in the array
            template = func_call % (args_code[-2], args_code[-1])
            # Remove the last 2 entries and append the max of the last 2
            del args_code[-2:]
            args_code.append(template)
        return str(args_code[0])

    def _print_Min(self, expr):
        """MINIUM of the arguments, can handle any number of arguments:
        Min(a,b,c,d) is written as min(a, min(min(b,c),d))"""
        if isinstance(SimulationDataType.dtype(), FloatC):
            func_call = 'fminf(%s, %s)'
        else:
            func_call = 'fmin(%s, %s)'

        nargs = len(expr.args)
        args_code = [self._print(a) for a in expr.args]
        for i in range(nargs-1):
            # Max of the last 2 arguments in the array
            template = func_call % (args_code[-2], args_code[-1])
            # Remove the last 2 entries and append the max of the last 2
            del args_code[-2:]
            args_code.append(template)
        return str(args_code[0])

    def _print_sign(self, expr):
        """signum function using ternary operators"""

        args = map(ccode, expr.args)
        args = [x for x in args]
        result = ','.join(args)
        result = '(%s > 0) ? 1 : ((%s < 0) ? -1 : 0)' % (result, result)
        return result

    def _print_DataObject(self, expr):
        """Raise error if a DataObject is found in the equation"""
        raise TypeError("Data object found in code generation, convert it to a dataset first, %s" % expr)

    def _print_Globalvariable(self, expr):
        # This should be handled in a different way if writing a kernel, if it is uses in the main cpp file it should be handled differently, only Global variables are supported and not arrays
        is_kernel = self.settings_opsc.get('kernel', False)
        if is_kernel:
            return ' *%s' % (str(expr))
        else:
            return "%s" % (str(expr))

    def _print_DataSetBase(self, expr):
        return str(expr)

    def _print_Pow(self, expr):
        """ Replace pow function calls with direct multiplication."""
        if isinstance(SimulationDataType.dtype(), FloatC):
            sqrt, pow_func, one = 'sqrtf(', 'powf(', '1.0f'
        else:
            sqrt, pow_func, one = 'sqrt(', 'pow(', '1.0'
        PREC = precedence(expr)
        if expr.exp in range(2, 7):
            return '(' + '*'.join([self.parenthesize(expr.base, PREC)] * int(expr.exp)) + ')'
        elif expr.exp in range(-6, 0):
            return '%s/(' % one + ('*'.join([self.parenthesize(expr.base, PREC)] * int(-expr.exp))) + ')'
        elif expr.exp == Rational(3,2):
            return '*'.join([self.parenthesize(expr.base, PREC)] + [sqrt + self.parenthesize(expr.base, PREC) + ')'])
        elif expr.exp == Rational(1,2):
            return '*'.join([sqrt + self.parenthesize(expr.base, PREC) + ')'])
        else:
            if isinstance(SimulationDataType.dtype(), FloatC):
                return '*'.join([pow_func + self.parenthesize(expr.base, PREC) + ', ' + str(expr.exp) + ')'])
            else:
                return super()._print_Pow(expr)

    def _print_Equality(self, expr):
        from opensbli.equation_types.opensbliequations import OpenSBLIEquation
        if isinstance(expr, OpenSBLIEquation):
            return "%s = %s" % (self._print(expr.lhs), self._print(expr.rhs))
        else:
            return "%s == %s" % (self._print(expr.lhs), self._print(expr.rhs))

    def _print_DataSet(self, expr):
        """ Prints the OpenSBLI dataset in the OPS format with the access numbers provided.
        Access numbers are updated for each kernel, see writing kernel in the OPSC class. """
        base = expr.base
        indices = expr.get_grid_indices
        out = ""
        # print(expr.__dict__)
        if self.cast_precision: # Change precision on the RHS of the equation only for certain quantites
            if hasattr(expr, "cast_precision") and not expr.cast_precision:
                pass # No casting on left hand side of equations
            elif hasattr(expr, "cast_precision") and expr.cast_precision: # Force the computation to be performed in lower precision on the RHS of the equations
                if base in self.settings_opsc['arrays_to_cast']:
                    # Change precision to that of the global simulation datatype
                    out += "(%s)" % SimulationDataType.dtype().opsc()
            else:
                print("Quantity: {} does not have this attribute.".format(expr))

        # Write the C code for this DataSet
        if self.settings_opsc.get('OPS_V2', True):
            out += "%s(%s)" % (self._print(base), ','.join([self._print(i) for i in indices]))
            return out
        else:
            if self.dataset_accs_dictionary[base]:
                out += "%s[%s(%s)]" % (self._print(base), self.dataset_accs_dictionary[base].name, ','.join([self._print(i) for i in indices]))
            else:
                raise ValueError("Did not find the OPS Access for DataSet %s " % expr.base)

#MBCHANGE
    def _print_IndexedBase(self, expr):
        return str(expr)
#MBCHANGE

    def _print_Indexed(self, expr):
        """ Print out an Indexed object.

        :arg expr: The Indexed expression.
        :returns: The indexed expression, as OPSC code.
        :rtype: str
        """
        indices = [ind for ind in expr.indices]
        for number, index in enumerate(indices):
            for sym in index.atoms(Symbol):
                indices[number] = indices[number].subs({sym: 0})
        out = "%s[%s]" % (self._print(expr.base.label), ','.join([self._print(index) for index in indices]))
        return out

    def _print_Grididx(self, expr):
        out = "%s[%s]" % (self._print(expr.base.label), ','.join([self._print(expr.number)]))
        return out

    def _print_UserFunction(self, fn):
        """ Prints the user defined funtion, we do not check if the function exists in the header files.
        This allows for users to write their own OPS functions and use them in the code."""
        return self._print(fn.args[0]) + '(%s)' % (', '.join([self._print(arg) for arg in fn.args[1:]]))


def pow_to_constant(expr):
    """Finds the negative powers of constant objects and evaluates them to a new constant object to reduce divisions."""
    from sympy.core.function import _coeff_isneg
    # Only negative powers i.e they correspond to division and they are stored into constant arrays
    inverse_terms = {}
    # Change the name of Rational counter to rcinv
    orig_name = rc.name

    for at in expr.atoms(Pow):
        # Remove common 1 / (gama - 1) factors
        if at is (1 / (ConstantObject('gama') - 1)):
            # print(at)
            rc.name = 'inv_' + 'gamma_m1'
            if at in rc.existing:
                inverse_terms[at] = rc.existing[at]
            else:
                inverse_terms[at] = rc.get_next_rational_constant(at)
        if _coeff_isneg(at.exp) and isinstance(at.base, ConstantObject):
            if (at.exp < -1):
                name = 'inv' + str(int(abs(at.exp))) + str(at.base)
                rc.name = name
            else:
                name = 'inv' + str(at.base)
                rc.name = name
            if at in rc.existing:
                inverse_terms[at] = rc.existing[at]
            else:
                inverse_terms[at] = rc.get_next_rational_constant(at)
    expr = expr.subs(inverse_terms)
    rc.name = orig_name  # change it back to the original name of Rational counter
    return expr


def ccode(expr, settings={}):
    """ Create an OPSC code printer object and write out the expression as an OPSC code string.

    :arg expr: The expression to translate into OPSC code.
    :arg Indexed_accs: Indexed OPS_ACC accesses.
    :arg constants: Constants that should be defined at the top of the OPSC code.
    :returns: The expression in OPSC code.
    :rtype: str."""
    if isinstance(expr, Equality):
        if 'boolean_equality' in settings.keys():
            if settings['boolean_equality']:
                equals = ' == '
            else:
                equals = ' = '
        else:
            equals = ' = '
        if 'rational' in settings.keys():
            pass
        else:
            expr = pow_to_constant(expr)
        code_print = OPSCCodePrinter(settings)
        code = code_print.doprint(expr.lhs) \
            + equals + OPSCCodePrinter(settings).doprint(expr.rhs)
        if isinstance(expr.lhs, GridVariable):
            code = code
        return code
    else:
        return OPSCCodePrinter(settings).doprint(expr)


class WriteString(object):
    def __init__(self, string):
        if isinstance(string, list):
            self.components = string
        elif isinstance(string, str):
            self.components = [string]
        else:
            raise ValueError("")
        return

    def __str__(self):
        return '\n'.join(self.components)

    def _print(self):
        return str(self)

    @property
    def opsc_code(self):
        return ['\n'.join(self.components)]


class OPSAccess(object):
    def __init__(self, no):
        self.name = "OPS_ACC%d" % no
        return


def indent_code(code_lines):
    """ Indent the code.

    :arg code_lines: The string or list of strings of lines of code to indent.
    :returns: A list of the indented line(s) of code.
    :rtype: list """

    p = C99CodePrinter()
    return p.indent_code(code_lines)


class OPSC(object):
    def __init__(self, algorithm, operation_count=False, OPS_diagnostics=1, OPS_V2=True, mixed_precision_config=None):
        """ Generating an OPSC code from the algorithm class.
        :arg object algorithm: An OpenSBLI algorithm class.
        :arg bool operation_count: If True, prints the number of arithmetic operations per kernel.
        :arg int OPS_diagnostics: OPS performance diagnostics. The default of 1 provides no kernel-based timing output
        A value of 5 gives a kernel breakdown of computational kernel and MPI exchange time.
        :arg bool OPS_V2: If True, uses the concise syntax in OPS for datasets (no explicit OPS_ARG number in kernels)."""
        if OPS_V2:
            self.ops_headers = {'input': "const ACC<%s> &%s", 'output': 'ACC<%s> &%s', 'inout': 'ACC<%s> &%s'}
        else:
            self.ops_headers = {'input': "const %s *%s", 'output': '%s *%s', 'inout': '%s *%s'}
        self.OPS_V2 = OPS_V2
        # if not algorithm.MultiBlock:
        self.operation_count = operation_count
        self.OPS_diagnostics = OPS_diagnostics
        self.MultiBlock = False
        self.nblocks = len(algorithm.block_descriptions)
        self.const_fname = 'constants.h'
        # Check if the simulation monitoring should be written to an output log file
        if algorithm.simulation_monitor:
            if len(algorithm.simulation_monitor.output_files) > 0:
                self.monitoring_output_file = True
            else:
                self.monitoring_output_file = False
        else:
            self.monitoring_output_file = False
        # Process any mixed precision customisations
        self.mixed_precision_config = mixed_precision_config
        self.arrays_to_cast = []
        if self.mixed_precision_config is not None:
            if 'casting' in self.mixed_precision_config:
                if self.mixed_precision_config['casting']:
                    self.cast_precision = True
                else:
                    self.cast_precision = False
                del self.mixed_precision_config['casting']
            else:
                self.cast_precision = False
            self.modify_dataset_precision(algorithm)
        else:
            self.cast_precision = False
        # First write the kernels, with this we will have the Rational constants to declare
        self.write_kernels(algorithm)
        def_decs = self.opsc_def_decs(algorithm)
        end = self.ops_exit(algorithm)
        algorithm.prg.components = def_decs + algorithm.prg.components + end
        code = algorithm.prg.opsc_code
        code = self.before_main(algorithm) + code
        f = open('opensbli.cpp', 'w')
        f.write('\n'.join(code))
        f.close()
        print("Successfully generated the OPS C code.")
        return

    def modify_dataset_precision(self, algorithm):
        """ Apply mixed precision options - change precision of certain quantities relative to the global simulation precision."""
        simulation_dsets = []
        # Get all the datasets defined in the simulation
        for d in algorithm.definitions_and_declarations.components:
            if isinstance(d, DataSetBase):
                simulation_dsets.append(d)
        # Add any missing ones
        for b in algorithm.blocks:
            for k, v, in b.block_datasets.items():
                simulation_dsets.append(v)
        # Remove duplicates
        simulation_dsets = list(set(simulation_dsets))
        # Change casting behaviour - force the quantities to lower precision in all RHS calculations
        if self.mixed_precision_config is not None:
            # Process the different input strategies to perform the precision changes
            for strategy, inputs in self.mixed_precision_config.items():
                # Get the inputs
                store_dsets = []
                arrays, modified_precision = [x.base for x in flatten(inputs[0])], inputs[1]
                # Different preset strategies
                # Time advance arrays (rho, rhou, rhov, rhow, rhoE)
                if strategy == 'q_vector':
                    lhs = [x.base for x in algorithm.time_advance_arrays]
                    for d in simulation_dsets:
                        if d in lhs:
                            d.datatype = modified_precision
                            store_dsets.append(d)
                # Work arrays used for temporary derivative calculations (StoreSome, and others)
                elif strategy == 'wk_arrays':
                    for d in simulation_dsets:
                        if 'wk' in str(d):
                            store_dsets.append(d)
                            d.datatype = modified_precision
                # Residual arrays used for time-advancement
                elif strategy == 'residuals':
                    for d in simulation_dsets:
                        if 'Residual' in str(d):
                            store_dsets.append(d)
                            d.datatype = modified_precision
                # Intermediate arrays used for time-stepping, filters
                elif strategy == 'RK_arrays':
                    RK_arrays = []
                    for b in flatten(algorithm.blocks):
                        for label, sc in b.discretisation_schemes.items():
                            if sc.schemetype == 'Temporal':
                                RK_arrays.append(sc.temp_RK_arrays)
                    RK_arrays = [x.base for x in flatten(RK_arrays)]
                    for d in simulation_dsets:
                            if d in RK_arrays:
                                store_dsets.append(d)
                                d.datatype = modified_precision
                # Custom input, user specified arrays
                elif strategy == 'custom':
                    for d in simulation_dsets:
                        if d in arrays:
                            store_dsets.append(d)
                            # print("Before:", d.datatype.opsc())
                            d.datatype = modified_precision
                            # print("After:", d.datatype.opsc())
                else:
                    raise ValueError("Unknown mixed precision preset: {}. Please choose from: q_vector, wk_arrays, residuals, RK_arrays, or custom.".format(strategy))
                store_dsets = sorted(store_dsets, key=lambda x: str(x))
                # Casting behaviour for the mixed-precision strategies
                if self.cast_precision:
                    self.arrays_to_cast += store_dsets
                print("Performed mixed precision on: {} - Modified precision of: {} from {} to {}.".format(strategy, store_dsets, SimulationDataType.dtype().opsc(), modified_precision.opsc()))

                # For preset values (non-custom) - update the list of arrays that had their precision modified
                self.mixed_precision_config[strategy] = (store_dsets, modified_precision)
        else: # No mixed precision
            pass    
        return

    def wrap_long_lines(self, code_lines):
        """ """
        limit = 120
        space = ' '
        formatted_code = []
        for code in code_lines:
            if len(code) > limit:
                # count the leading spaces so that we can use it at the end
                Leading_spaces = len(code) - len(code.lstrip())
                codelstrip = code.lstrip()
                length_ofstring = Leading_spaces
                split_lines = []
                # start the string with the leading spaces
                string = Leading_spaces*space
                for s in codelstrip.split(" "):
                    length_ofstring = len(string) + len(s)
                    if length_ofstring >= limit:
                        split_lines += [string]
                        # Increase the indentation by 4 spaces
                        string = (Leading_spaces + Leading_spaces)*space + s
                    else:
                        string = string + " " + s
                formatted_code += split_lines + [string]
                # Check the correctness of the code
                s1 = ''.join(split_lines + [string])
                s1 = "".join(s1.split())
                s2 = "".join(code.split())
                if s1 == s2:
                    pass
                else:
                    print(code)
                    print('\n'.join(split_lines + [string]))
                    raise ValueError("The code and the formatted line are not same")
            else:
                formatted_code += [code]
        return formatted_code

    def check_failed_central(self, code, kernel):
        """ Warns the user and exits code-generation if one of the CentralDerivatives failed to discretize and a code was not produced."""
        for line in code:
            if "CentralDerivative" in line:
                print('\33[91m' + "WARNING: Code generation failed to discretize derivative: {} in kernel: '{}', this code will not compile.".format(line, kernel.computation_name) + '\033[0m')
                exit()
            elif "Not supported in C" in line:
                print('\33[91m' + "WARNING: Code generation failed to discretize derivative: {} in kernel: '{}', this code will not compile.".format(line, kernel.computation_name) + '\033[0m')
                exit()
        return

    def kernel_header(self, tuple_list, idx_constants):
        code = []
        ins, outs, inouts = [x for x in tuple_list if x[1] == 'input'], [x for x in tuple_list if x[1] == 'output'], [x for x in tuple_list if x[1] == 'inout']
        ins, outs, inouts = sorted(ins, key=lambda x: str(x[0])), sorted(outs, key=lambda x: str(x[0])), sorted(inouts, key=lambda x: str(x[0]))
        tuple_list = ins + outs + inouts + idx_constants
        for key, val in (tuple_list):
            if str(key) == 'rkA' or str(key) == 'rkB' or str(key) == 'rkold' or str(key) == 'rknew': # RK coefficients in the kernel header
                if hasattr(key, "datatype") and key.datatype:
                    code += ['const %s *%s' % (key.datatype.opsc(), key)]
                else:
                    code += ['const %s *%s' % (SimulationDataType.opsc(), key)]
            elif str(key) == 'iter': # current iteration counter
                code += ['const int *%s' % key]
            elif isinstance(key, ReductionVariable):
                if key.reduction_type == 'OPS_INC': # summation reduction variables
                    code += ['%s *%s' % (key.datatype.opsc(), key)]
                elif val == 'input':
                    code += ['const %s *%s' % (key.datatype.opsc(), key)]
                else:
                    code += ['%s *%s' % (key.datatype.opsc(), key)]
            else:
                # Argument is a DataSet
                # Re-apply mixed precision if needed - some DataSets revert back to the simulation datatype - why?
                if hasattr(key, "datatype") and key.datatype:
                    if self.mixed_precision_config is not None:
                        for k, v in self.mixed_precision_config.items():
                            to_modify = [str(x) for x in flatten(v[0])]
                            if str(key) in to_modify:
                                key.datatype = v[1]
                    code += [self.ops_headers[val] % (key.datatype.opsc(), key)]
                else:
                    raise ValueError("Dataset: {} is missing its datatype.".format(key))
        code = ', '.join(code)
        return code

    def add_casting_switch(self, input_eqn):
        RHS_args = input_eqn.rhs.args
        for argument in RHS_args:
            for dset in argument.atoms(DataSet):
                if str(dset) in [str(x) for x in self.arrays_to_cast]:
                    dset.cast_precision = True
                else:
                    dset.cast_precision = False
        return input_eqn

    def kernel_computation_opsc(self, kernel):
        """ Function to write the out the contents of each computational kernel."""
        ins = kernel.rhs_datasetbases
        outs = kernel.lhs_datasetbases
        inouts = ins.intersection(outs)
        ins = ins.difference(inouts)
        outs = outs.difference(inouts)
        # eqs = kernel.equations
        all_dataset_inps = list(ins) + list(outs) + list(inouts)
        all_dataset_types = ['input' for i in ins] + ['output' for o in outs] + ['inout' for io in inouts]
        # add the global variables to the inputs and outputs
        global_ins, global_outs = kernel.global_variables
        if global_ins.intersection(global_outs):
            raise NotImplementedError("Input output of global variables is not implemented")
        all_dataset_inps += list(global_ins) + list(global_outs)
        all_dataset_types += ['input' for i in global_ins] + ['output' for o in global_outs]
        # Add any reduction variables present in the kernel, to generate the kernel headers
        reduction_ins, reduction_outs = kernel.reduction_variables
        all_dataset_inps += list(reduction_ins) + list(reduction_outs)
        all_dataset_types += ['input' for i in reduction_ins] + ['output' for o in reduction_outs]
        # Use list of tuples as dictionary messes the order
        header_dictionary = list(zip(all_dataset_inps, all_dataset_types))
        idx_constants = []
        if kernel.IndexedConstants:
            for i in sorted(kernel.IndexedConstants, key=lambda x: str(x)):
                idx_constants += [tuple([(i.base), 'input'])]
        other_inputs = ""
        # Local i, j, k index object (ignores MPI)
        if kernel.grid_indices_used:
            other_inputs += ", const int *idx"  # WARNING hard coded here
        else:
            other_inputs = ''
        code = ["void %s(" % kernel.kernelname + self.kernel_header(header_dictionary, idx_constants) + other_inputs + ')' + '\n{']
        ops_accs = [OPSAccess(no) for no in range(len(all_dataset_inps))]
        OPSCCodePrinter.dataset_accs_dictionary = dict(zip(all_dataset_inps, ops_accs))
        # Find all the grid variables and declare them at the top
        gridvariables = set()
        out = []
        for eq in kernel.equations:
            default_kernel_settings = {'kernel': True, 'OPS_V2': self.OPS_V2, 'arrays_to_cast' : self.arrays_to_cast}
            bool_settings = {'kernel': True, 'OPS_V2': self.OPS_V2, 'arrays_to_cast' : self.arrays_to_cast, 'boolean_equality' : True}
            # Note which DataSets are used on the LHS of equations
            if self.cast_precision:
                if isinstance(eq, GroupedPiecewise):
                    for i, (expr, condition) in enumerate(eq.args):
                        # Process all grouped equations within this condition
                        for single_eqn in expr:
                            # Don't cast LHS
                            LHS_of_equation = single_eqn.lhs
                            if isinstance(LHS_of_equation, DataSet):
                                single_eqn.lhs.cast_precision = False
                            # Loop over all arguments of this equation
                            single_eqn = self.add_casting_switch(single_eqn)
                            # Make sure LHS is not cast
                            # if not isinstance(LHS_of_equation, GridVariable):
                            #     single_eqn.lhs.cast_precision = False
                elif isinstance(eq, WhileLoop):
                    raise ValueError("Mixed precision not implemented yet for WhileLoop.")
                elif isinstance(eq, ForLoop):
                    raise ValueError("Mixed precision not implemented yet for ForLoop.")
                else: # Regular equations
                    # Never cast precision of left-hand side assignments
                    LHS_of_equation = eq.lhs
                    if isinstance(eq.lhs, DataSet):
                        eq.lhs.cast_precision = False
                    # # Check for Piecewise conditions
                    if isinstance(eq.rhs, Piecewise):
                        eq.lhs.cast_precision = False
                    else: # LHS of the equation is a DataSet
                        # # Quantity does not appear on the right hand side of the equation also
                        # if LHS_of_equation not in eq.rhs.atoms(DataSet):
                        #     eq.lhs.cast_precision = True
                        if isinstance(eq.rhs, Piecewise):
                            for pairs in eq.rhs.args:
                                pw_expr = pairs[0]
                                for dset in pw_expr.atoms(DataSet):
                                    dset.cast_precision = True
                            eq.lhs.cast_precision = False
                        else: # Regular equation
                            eq = self.add_casting_switch(eq)
                            # Make sure LHS is not cast
                            if not isinstance(LHS_of_equation, GridVariable):
                                eq.lhs.cast_precision = False

            # Get the grid variables
            gridvariables = gridvariables.union(eq.atoms(GridVariable))
            # Get the reduction variables and detect whether they are input or output
            if isinstance(eq, Equality):
                for rv in eq.lhs.atoms(ReductionVariable):
                    rv.usage = 'lhs'
                for rv in eq.rhs.atoms(ReductionVariable):
                    rv.usage = 'rhs'

            if isinstance(eq, Equality):
                out += [ccode(eq, settings=default_kernel_settings) + ';\n']
            elif isinstance(eq, GroupedPiecewise):
                for i, (expr, condition) in enumerate(eq.args):
                    if i == 0:
                        out += ['if (%s)' % ccode(condition, settings=bool_settings) + '{\n']
                        if is_sequence(expr):
                            for eqn in expr:
                                out += [ccode(eqn, settings=default_kernel_settings) + ';\n']
                        else:
                            out += [ccode(expr, settings=default_kernel_settings) + ';\n']
                        out += ['}\n']
                    elif condition != True:
                        out += ['else if (%s)' % ccode(condition, settings=bool_settings) + '{\n']
                        if is_sequence(expr):
                            for eqn in expr:
                                out += [ccode(eqn, settings=default_kernel_settings) + ';\n']
                        else:
                            out += [ccode(expr, settings=default_kernel_settings) + ';\n']
                        out += ['}\n']
                    else:
                        out += ['else{\n']
                        if is_sequence(expr):
                            for eqn in expr:
                                out += [ccode(eqn, settings=default_kernel_settings) + ';\n']
                        else:
                            out += [ccode(expr, settings=default_kernel_settings) + ';\n']
                        out += ['}\n']
            elif isinstance(eq, WhileLoop):
                for i, (expr, condition) in enumerate(eq.args):
                    if i == 0:
                        out += ['while (%s)' % ccode(condition, settings=bool_settings) + '{\n']
                        if is_sequence(expr):
                            for eqn in expr:
                                out += [ccode(eqn, settings=default_kernel_settings) + ';\n']
                        else:
                            out += [ccode(expr, settings=default_kernel_settings) + ';\n']
                        out += ['}\n']
                    elif i == 1:
                        pass
                    else:
                        raise ValueError("While Loop should only have two conditions.")
            elif isinstance(eq, ForLoop):
                for i, (expr, condition) in enumerate(eq.args):
                    if i == 0:
                        start = condition
                        evaluate = expr
                    elif i == 1:
                        end = expr
                    else:
                        raise ValueError("While Loop should only have two conditions.")
                # Create the C code for the for loop and populate it with equations
                iteration_index = start.lhs      
                out += ['for (int %s; %s; %s++)' % ((ccode(start, settings=default_kernel_settings), ccode(end, settings=bool_settings), ccode(iteration_index, settings=bool_settings))) + '{\n']
                # Add equations inside the for loop
                if is_sequence(evaluate):
                    for eqn in evaluate:
                        out += [ccode(eqn, settings=default_kernel_settings) + ';\n']
                else:
                    out += [ccode(evaluate, settings=default_kernel_settings) + ';\n']
                out += ['}\n']
            else:
                pprint(eq)
                raise TypeError("Unclassified type of equation.")
        # Sort the gridvariables to fix the order
        gridvariables = sorted(list(gridvariables), key=lambda x: str(x))
        for gv in gridvariables:
            code += ["%s %s = 0.0;" % (SimulationDataType.opsc(), str(gv))]
        # if '\n' in out[-1]:
        #     out[-1] = out[-1].replace('\n', '', out[-1].count(' \n ')-1)
        code += out + ['}']  # close Kernel
        OPSCCodePrinter.dataset_accs_dictionary = {}
        return code

    def write_kernels(self, algorithm):
        """ A function to write out the kernels header file definining all of the computations to be performed."""
        from opensbli.core.kernel import Kernel
        kernels = self.loop_alg(algorithm, Kernel)
        # Count the number of operations per kernel
        if self.operation_count:
            total = 0
            for k in kernels:
                n_operations = count_ops(k.equations)
                total += n_operations
                print([k.kernel_no, k.computation_name, 'operations: %d' % n_operations])
            print("Total operation count is: %d" % total)

        files = [open('%s_kernels.h' % b.block_name, 'w') for b in algorithm.block_descriptions]
        for i, f in enumerate(files):
            name = ('%s_kernel_H' % algorithm.block_descriptions[i].block_name).upper()
            f.write('#ifndef %s\n' % name)
            f.write('#define %s\n' % name)
        for k in kernels:
            out = self.kernel_computation_opsc(k) + ['\n']
            out = indent_code(out)
            out = self.wrap_long_lines(out)
            files[k.block_number].write('\n'.join(out))
            # Check for any failed discretisation
            self.check_failed_central(out, k)
        for f in files:
            f.write("#endif\n")
        files = [f.close() for f in files]
        return

    def ops_exit(self, algorithm):
        """ Exits the OPS program with optional kernel-based timing output."""
        output = []
        if self.OPS_diagnostics > 1:
            output += [WriteString("ops_timing_output(std::cout);")]
        if self.monitoring_output_file:
            for i in range(len(algorithm.simulation_monitor.output_files)):
                index = algorithm.simulation_monitor.blocks[i].blocknumber
                output += [WriteString("fclose(f%d);" % index)]
        output += [WriteString("ops_exit();")]
        return output

    def before_main(self, algorithm):
        """ Adds the required preamble to the main opensbli.cpp file and declares the simulation constants."""
        out = ['#include <stdlib.h> \n#include <string.h> \n#include <math.h> \n#include "constants.h"' ]
        from opensbli.core.kernel import ConstantsToDeclare
        # Declare a restart flag and loop variables globally
        constant_declarations = ["%s %s;" % ('int', 'restart')]
        constant_declarations += ["%s %s;" % ('int', 'iter')]
        constant_declarations += ["%s %s;" % ('int', 'stage')]
        constant_declarations += ["%s %s;" % (SimulationDataType.opsc(), 'tstart')]

        # Write the constants to a separate file instead
        for d in sorted(ConstantsToDeclare.constants, key=lambda x: str(x)):
            if isinstance(d, ConstantObject):
                constant_declarations += ["%s %s;" % (d.datatype.opsc(), d)]
            elif isinstance(d, ConstantIndexed):
                if not d.inline_array:
                    indices = ''
                    for s in d.shape:
                        indices = indices + '[%d]' % s
                    constant_declarations += ["%s %s%s;" % (d.datatype.opsc(), d.base.label, indices)]
        # Write the constant declarations to a separate file
        const_file = open(self.const_fname, 'w')
        const_file.write('\n'.join(['// Declaration of global constants'] + flatten(constant_declarations)))
        const_file.close()
        # Declare the simulation blocks
        out += ['#define OPS_%dD' % algorithm.block_descriptions[0].ndim]
        if self.OPS_V2:
            out += ['#define OPS_API 2']
        out += ['#include \"ops_seq.h\"']
        for b in algorithm.block_descriptions:
            out += ['#include \"%s_kernels.h\"' % b.block_name]
        # IO and constant functions
        out += ['#include "io.h"']
        # Include optional simulation monitoring reductions file
        if algorithm.simulation_monitor:
            self.opened = []
            out += ['#include \"%s\"' % algorithm.simulation_monitor.filename]
            if len(algorithm.simulation_monitor.output_files) > 0:
                for i in range(len(algorithm.simulation_monitor.output_files)):
                    index = algorithm.simulation_monitor.blocks[i].blocknumber
                    out += ['FILE *f%d = fopen(\"%s\", \"a\");' % (index, str(algorithm.simulation_monitor.output_files[i]))]
        return out

    def opsc_def_decs(self, algorithm):
        """ Declares the datasets and stencils required by the program."""
        from opensbli.core.kernel import StencilObject, ConstantsToDeclare
        from opensbli.core.boundary_conditions.exchange import Exchange
        output = []
        # Add OPS_init to the declarations
        output += self.ops_init()
        # Sort the constants to a consistent ordering
        ConstantsToDeclare.sort_constants()
        # First process all the constants in the definitions
        output += self.set_constant_values(ConstantsToDeclare.constants)
        # OPS declaration of the constants
        output += self.declare_ops_constants(ConstantsToDeclare.constants)
        # Once the constants are done define and declare OPS dats
        # Define and declare blocks
        for b in algorithm.block_descriptions:
            output += self.declare_block(b)
        # Define and declare datasets on each block
        f = open('defdec_data_set.h', 'w')
        datasets_dec = []
        output += [WriteString("#include \"defdec_data_set.h\"")]
        # Sort the declarations alphabetically before writing out
        store_stencils, store_dsets, store_reductions = [], [], []
        for d in algorithm.definitions_and_declarations.components:
            if isinstance(d, DataSetBase):
                store_dsets.append(d)
            elif isinstance(d, StencilObject):
                store_stencils.append(d)
            elif isinstance(d, ReductionVariable):
                store_reductions.append(d)
            else:
                print(d)
                print(type(d))
                raise TypeError("Quantity: {} is not defined within the simulation and cannot be declared.".format(d))
        dsets_to_declare = dict([(str(x), x) for x in store_dsets])

        for name in sorted(dsets_to_declare, key=str.lower):
            d = dsets_to_declare[name]
            datasets_dec += self.declare_dataset(d, algorithm)
        f.write('\n'.join(flatten([dset.opsc_code for dset in datasets_dec])))
        f.write('\n')
        f.close()
        # Declare stencils
        output += [WriteString("// Define and declare stencils")]
        f = open('stencils.h', 'w')
        output += [WriteString("#include \"stencils.h\"")]
        stencil_declarations = []
        for d in store_stencils:
            stencil_declarations += self.ops_stencils_declare(d)
        f.write('\n'.join(flatten([x.opsc_code for x in stencil_declarations])))
        f.write('\n')
        f.close()

        # Define reduction operation handles (global min, max reductions ...)
        if len(store_reductions) > 0:
            output += [WriteString("// Define and declare OPS reduction handles")]
            for rv in store_reductions:
                output += self.declare_reduction(rv)

        # Loop through algorithm components to include any halo exchanges
        exchange_list = self.loop_alg(algorithm, Exchange)
        if exchange_list:
            f = open('bc_exchanges.h', 'w')  # write BC_exchange code to a separate file
            exchange_code = []
            for e in exchange_list:
                call, code = self.bc_exchange_call_code(e)
                exchange_code += [code]
            f.write('\n'.join(flatten(exchange_code)))
            f.write('\n')
            f.close()
            output += [WriteString("#include \"bc_exchanges.h\"")]  # Include statement in the code
        # Write HDF5 I/O calls to a separate file
        io_file = open('io.h', 'w')
        io_file.close()
        output += self.ops_partition()
        # Notify whether the simulation is being restarted or not
        output += self.restart_notification()
        # This MUST be done after the partition command to avoid MPI HDF5 errors
        output += self.restart_simulation()
        return output

    def restart_simulation(self):
        """ Initialises the simulation time and iteration number from file if restarting the simulation."""
        out = [WriteString('// Constants from HDF5 restart file')]
        out += [WriteString('if (restart == 1){')]
        for c in self.restarted_constants:
            if c.name == 'start_iter':
                out += [WriteString('ops_get_const_hdf5(\"%s\", 1, \"%s\", (char*)&%s, "restart.h5");' % ('iter', c.datatype.opsc(), c.name))]
            else:
                out += [WriteString('ops_get_const_hdf5(\"%s\", 1, \"%s\", (char*)&%s, "restart.h5");' % (c.name, c.datatype.opsc(), c.name))]
        out += [WriteString('}')]
        out += [WriteString('else {')]
        for c in self.restarted_constants:
            out += [WriteString("%s = %s;" % (str(c), c.value))]
        out += [WriteString('}')]
        out += [WriteString('tstart = simulation_time;\n')]
        return out

    def ops_stencils_declare(self, s):
        out = []
        dtype = s.datatype.opsc()
        name = s.name + 'temp'
        sorted_stencil = s.sort_stencil_indices()
        out = [self.declare_inline_array(dtype, name, [st for st in flatten(sorted_stencil) if not isinstance(st, Idx)])]
        out += [WriteString('ops_stencil %s = ops_decl_stencil(%d,%d,%s,\"%s\");' % (s.name, s.ndim, len(s.stencil), name, name))]
        return out

    def ops_partition(self):
        """ Initialise an OPS partition for the purpose of multi-block and/or MPI partitioning.

        :returns: The partitioning code in OPSC format. Each line is a separate list element.
        :rtype: list"""
        output = [WriteString('// Init OPS partition')]
        # Add timers to MPI partition time
        output += [WriteString('double partition_start0, elapsed_partition_start0, partition_end0, elapsed_partition_end0;')]
        output += [WriteString('ops_timers(&partition_start0, &elapsed_partition_start0);')]
        output += [WriteString('ops_partition(\"\");')]
        output += [WriteString('ops_timers(&partition_end0, &elapsed_partition_end0);')]
        output += [WriteString('ops_printf("-----------------------------------------\\n MPI partition and reading input file time: %lf\\n -----------------------------------------\\n", elapsed_partition_end0-elapsed_partition_start0);')]
        # output += [WriteString('fflush(stdout);\n')]
        return output

    def restart_notification(self):
        """ Notifies the user whether the simulation is being restarted from file or not."""
        out = [WriteString('// Restart procedure')]
        out += [WriteString('ops_printf("\\033[1;32m\");')]
        out += [WriteString('if (restart == 1){')]
        ### Add the simulation time afterwards ###
        out += [WriteString('ops_printf("OpenSBLI is restarting from the input file: restart.h5\\n");')]
        out += [WriteString("}")]
        # Else clause
        out += [WriteString('else {')]
        out += [WriteString('ops_printf("OpenSBLI is starting from the initial condition.\\n");')]
        out += [WriteString("}")]
        out += [WriteString('ops_printf("\\033[0m");')]
        return out

    def ops_init(self):
        """ The default diagnostics level is 1, which offers no diagnostic information and should be used for production runs.
        Refer to OPS user manual for more information.
        :returns: The call to ops_init.
        :rtype: list """
        out = [WriteString('// Initializing OPS ')]
        return out + [WriteString('ops_init(argc,argv,%d);' % (self.OPS_diagnostics))]

    def bc_exchange_call_code(self, instance):
        """ Generates the code for OPS exchanges. Used for example in the periodic boundary condition."""
        off = 0
        halo = 'halo'
        # instance.transfer_size = instance.transfer_from
        # Name of the halo exchange
        name = instance.name
        # self.halo_exchange_number = self.halo_exchange_number + 1
        code = ['// Boundary condition exchange code on %s direction %s %s' % (instance.block_name, instance.direction, instance.side)]
        code += ['ops_halo_group %s %s' % (name, ";")]
        code += ["{"]
        code += ['int halo_iter[] = {%s}%s' % (', '.join([str(s) for s in instance.transfer_size]), ";")]
        code += ['int from_base[] = {%s}%s' % (', '.join([str(s) for s in instance.transfer_from]), ";")]
        code += ['int to_base[] = {%s}%s' % (', '.join([str(s) for s in instance.transfer_to]), ";")]
        # dir in OPSC. WARNING: Not sure what it is, but 1 to ndim works.
        from_dir = [ind+1 for ind in range(len(instance.transfer_to))]
        to_dir = [ind+1 for ind in range(len(instance.transfer_to))]
        # MBCHANGE - flip the direction if True
        if instance.flip[-1]:
            to_dir[instance.flip[1]] = -to_dir[instance.flip[1]]
        # MBCHANGE
        code += ['int from_dir[] = {%s}%s' % (', '.join([str(ind) for ind in from_dir]), ";")]
        code += ['int to_dir[] = {%s}%s' % (', '.join([str(ind) for ind in to_dir]), ";")]
        # Process the arrays
        for no, arr in enumerate(instance.transfer_arrays):
            from_array = instance.from_arrays[no]
            to_array = instance.to_arrays[no]
            code += ['ops_halo %s%d = ops_decl_halo(%s, %s, halo_iter, from_base, to_base, from_dir, to_dir)%s'
                     % (halo, off, from_array.base, to_array.base, ";")]
            off = off+1
        code += ['ops_halo grp[] = {%s}%s' % (','.join([str('%s%s' % (halo, of)) for of in range(off)]), ";")]
        code += ['%s = ops_decl_halo_group(%d,grp)%s' % (name, off, ";")]
        code += ["}"]
        # Finished OPS halo exchange, now get the call
        instance.call_name = 'ops_halo_transfer(%s)%s' % (name, ";")
        call = ['// Boundary condition exchange calls', 'ops_halo_transfer(%s)%s' % (name, ";")]
        for no, c in enumerate(code):
            code[no] = WriteString(c).opsc_code
        return call, code

    def loop_alg(self, algorithm, type_of_component):
        type_list = []

        def _generate(components, type_list):
            for component1 in components:
                if hasattr(component1, 'components'):
                    _generate(component1.components, type_list)
                elif isinstance(component1, type_of_component):
                    if component1 in type_list:
                        pass
                    else:
                        type_list += [component1]
        for c in algorithm.prg.components:
            _generate([c], type_list)
        return type_list

    def define_block(self, b):
        if not self.MultiBlock:
            return [WriteString("ops_block %s;" % b.block_name)]
        else:
            raise NotImplementedError("")

    def declare_block(self, b):
        if not self.MultiBlock:
            out = [WriteString("// Define and Declare OPS Block")]
            out += [WriteString('ops_block %s = ops_decl_block(%d, \"%s\");' % (b.block_name, b.ndim, b.block_name))]
            return out
        else:
            raise NotImplementedError("")

    def set_constant_values(self, constants):
        """ Declares all of the constants required by the simulation at the start of the program."""
        # First restart any constants from HDF5 if required
        out = []
        out += [WriteString('// Set restart to 1 to restart the simulation from HDF5 file')]
        restart = [c for c in constants if str(c) == 'restart'][0]
        out += [WriteString("%s = %s;" % (str(restart), restart.value))]
        constants.remove(restart)
        # Find which constants to restart
        self.restarted_constants = [x for x in constants if x.restart]
        init_constants = [c for c in constants if not c.restart]
        out += [WriteString('// User defined constant values')]
        # Give priority to lengths used to calculate grid spacings
        lengths = ['Lx', 'Ly', 'Lz', 'Lx0', 'Lx1', 'Lx2']
        for c in sorted(init_constants, key=lambda x: str(x))[::-1]:
            if str(c) in lengths:
                init_constants.remove(c)
                init_constants.insert(0, c)
        # Write the rest of the constants
        for c in init_constants:
            if isinstance(c, ConstantObject):
                if not isinstance(c.value, str):
                    out += [WriteString("%s = %s;" % (str(c), ccode(c.value, settings={'rational': True})))]
                else:
                    out += [WriteString("%s=%s;" % (str(c), c.value))]

            elif isinstance(c, ConstantIndexed):
                if c.value:
                    if len(c.shape) == 1:
                        if c.inline_array:
                            values = [ccode(c.value[i], settings={'rational': True}) for i in range(c.shape[0])]
                            out += [WriteString("%s %s[] = {%s};" % (c.datatype.opsc(), c.base.label, ', '.join(values)))]
                        else:
                            indices = ''
                            for s in c.shape:
                                indices = indices + '[%d]' % s
                            out += [WriteString("%s %s%s;" % (c.base.label, indices))]
                            for i in range(c.shape[0]):
                                out += [WriteString("%s[%d] = %s;" % (str(c.base.label), i, ccode(c.value[i], settings={'rational': True})))]
        return out

    def declare_ops_constants(self, input_constants):
        """ Calls the OPS declare constant function for all of the defined constants."""
        OPS_constant_declarations = []
        for c in sorted(input_constants, key=lambda x: str(x)):
            if isinstance(c, ConstantObject):
                OPS_constant_declarations += [WriteString("ops_decl_const(\"%s\" , 1, \"%s\", &%s);" % (str(c), c.datatype.opsc(), str(c)))]
        return OPS_constant_declarations

    def declare_inline_array(self, dtype, name, values):
        return WriteString('%s %s[] = {%s};' % (dtype, name, ', '.join([str(s) for s in values])))

    def update_inline_array(self, name, values):
        out = []
        for no, v in enumerate(values):
            out += [WriteString("%s[%d] = %s;" % (name, no, v))]
        return out

    def define_dataset(self, dset):
        if not self.MultiBlock:
            return [WriteString("ops_dat %s;" % (dset))]

    def get_max_halos(self, halos):
        halo_m = []
        halo_p = []
        for direction in range(len(halos)):
            if halos[direction][0]:
                hal = [d.get_halos(0) for d in halos[direction][0]]
                halo_m += [int(min(hal))]
            else:
                halo_m += [0]
            if halos[direction][1]:
                hal = [d.get_halos(1) for d in halos[direction][1]]
                halo_p += [int(max(hal))]
            else:
                halo_p += [0]
        return halo_m, halo_p

    def declare_reduction(self, rv):
        """ Declare a reduction variable in the code with the necessary handles."""
        dtype = SimulationDataType.dtype()
        variable_declaration = WriteString("%s %s = 0.0;" % (dtype.opsc(), str(rv.value)))
        handle_declaration = WriteString('ops_reduction %s = ops_decl_reduction_handle(sizeof(%s), \"%s\", \"reduction_%s\");' % (str(rv), dtype.opsc(), dtype.opsc(), str(rv)))
        out = [variable_declaration, handle_declaration]
        return out

    def initialize_dataset(self, dset, dtype):
        """ Initialize a dataset to zeros, not from a restart file."""
        # Residual and time-advance arrays do not require halos unless using shock filter
        if ('Residual' in str(dset) or 'tempRK' in str(dset) or 'RKold' in str(dset)):
            hm, hp = [-5 for _ in range(len(dset.size))], [5 for _ in range(len(dset.size))]
        else:
            hm, hp = self.get_max_halos(dset.halo_ranges)
        halo_p = self.declare_inline_array("int", "halo_p", hp)
        halo_m = self.declare_inline_array("int", "halo_m", hm)
        sizes = self.declare_inline_array("int", "size", [str(s) for s in (dset.size)])
        base = self.declare_inline_array("int", "base", [0 for i in range(len(dset.size))])
        value = WriteString("%s* value = NULL;" % dtype.opsc())
        temp = '%s = ops_decl_dat(%s, 1, size, base, halo_m, halo_p, value, \"%s\", \"%s\");' % (dset,
                                                                                                 dset.block_name, dtype.opsc(), dset)
        return [halo_p, halo_m, sizes, base, value, WriteString(temp)]

    def restart_dataset(self, dset, dtype, fname):
        """ Initialize a dataset from a restart HDF5 file."""
        temp = '%s = ops_decl_dat_hdf5(%s, 1, \"%s\", \"%s\", \"%s\");' % (dset,
                                                                           dset.block_name, dtype.opsc(), dset, fname)
        return [WriteString(temp)]

    def declare_dataset(self, dset, algorithm):
        """ Allocates memory for the storage arrays used by the simulation."""
        # Left hand side of the equations, quantities advanced in time
        time_advance_arrays = algorithm.time_advance_arrays
        # Coordinates evaluated in the code (not read from a grid file), which should also be restarted
        coordinates_to_restart = []
        for b in algorithm.block_descriptions:
            coordinates_to_restart += [str(x) for x in b.coordinate_arrays_to_restart]
        # Set the datatype of the array to declare
        if dset.datatype:
            dtype = dset.datatype
        else:
            dtype = SimulationDataType.dtype()
        # Create the code segment
        declaration = WriteString("ops_dat %s;" % dset)
        out = [declaration, WriteString("{")]
        # Add a restart flag to make it easier to restart the time advance arrays
        time_advance_arrays = [str(x) for x in time_advance_arrays]
        # Restart flag to the time advance arrays
        if str(dset) in time_advance_arrays:
            if not dset.write_to_hdf5:
                print('\33[91m' + "WARNING: Time-advance array: {} has not been set to be written to the restart file.".format(str(dset)) + '\033[0m')
            out += [WriteString('if (restart == 1){')]
            out += self.restart_dataset(dset, dtype, 'restart.h5')
            out += [WriteString("}")]
            # Else clause
            out += [WriteString('else {')]
            out += self.initialize_dataset(dset, dtype)
            out += [WriteString("}")]
        # Coordinates evaluated inside the code, not via an external grid file
        elif str(dset) in coordinates_to_restart:
            if not dset.write_to_hdf5:
                print('\33[91m' + "WARNING: Coordinate array: {} is computed during the initialisation kernel but has not been set to be written to the restart file.".format(str(dset)) + '\033[0m')
            out += [WriteString('if (restart == 1){')]
            out += self.restart_dataset(dset, dtype, 'restart.h5')
            out += [WriteString("}")]
            # Else clause
            out += [WriteString('else {')]
            out += self.initialize_dataset(dset, dtype)
            out += [WriteString("}")]
        # All other arrays
        else:
            # Externally provided grid file or restart file
            if dset.read_from_hdf5:
                out += self.restart_dataset(dset, dtype, dset.input_file_name)
            else:
                out += self.initialize_dataset(dset, dtype)
        out += [WriteString("}")]
        return out
