"""@brief
   @authors Satya Pramod Jammy, David J Lusher
   @contributors 
   @details
"""
from sympy import flatten, simplify, symbols, factor, count_ops, pprint, Piecewise, Equality, simplify, Pow, Mul, collect, Rational, factor_terms, gcd_terms, Tuple
from sympy.functions.elementary.piecewise import ExprCondPair
from opensbli.core.opensbliobjects import ConstantObject, CoordinateObject, DataObject, DataSet, GroupedPiecewise, GroupedCondition, ConstantIndexed, Grididx
from opensbli.core.grid import GridVariable
from opensbli.core.opensblifunctions import CentralDerivative
from opensbli.core.kernel import Kernel
from opensbli.schemes.spatial import Central
from opensbli.equation_types.opensbliequations import OpenSBLIEq, SimulationEquations, ConstituentRelations
from opensbli.equation_types.metric import MetricsEquation
from collections import OrderedDict
import re


class StoreSome(Central):
    """ Low-storage algorithms to reduce memory intensity and the number of global storage arrays.
        S.P. Jammy et al. Journal of Computational Science. Vol 36, September 2019 10.015."""

    def __init__(self, order, der_fns_to_store, merged=True, group_stored=False, level=1):
        """ Set up the scheme.
        :arg int order: The order of accuracy of the scheme."""
        Central.__init__(self, order)
        self.fns = der_fns_to_store
        self.merged = merged
        self.group_stored = group_stored
        self.factor = False
        return

    def generate_derivatives_to_store(self, coordinates, block):
        """ Creates CentralDerivative objects of the derivatives to be stored."""
        data_objects = flatten([symbols('%s' % self.fns, **{'cls': DataObject})])
        data_sets = [block.location_dataset(str(d)) for d in data_objects]
        # coords = [c for c in coordinates if not c.get_coordinate_type()]
        # Fix for Neil's issue with multiple passive scalars
        coords = [c for c in coordinates if not str(c) == 't']
        coords = sorted(coords, key=lambda x: x.direction)
        self.derivatives_to_store = dict(zip([i for i in range(block.ndim)], [[] for _ in range(block.ndim)]))
        for a in data_sets:
            for i, b in enumerate(coords):
                self.derivatives_to_store[i] += [CentralDerivative(a, b)]
        return

    def discretise(self, type_of_eq, block):
        """ This is the main calling function from opensbli equations.spatial_discretisation, which is called from block.discretise."""
        self.set_halos(block)
        if isinstance(type_of_eq, SimulationEquations):
            # Simulation equations are always solved as sbli_rhs_discretisation currently
            self.sbli_rhs_discretisation(type_of_eq, block)
            return self.required_constituent_relations
        else:
            discretised_eq = self.SS(type_of_eq, block, None, group=False)
            if discretised_eq:
                # discretised_eq = self.merge_conditionals(discretised_eq, block)
                discretisation_kernel = Kernel(block, computation_name="%s evaluation" % type_of_eq.__class__.__name__)
                discretisation_kernel.set_grid_range(block)
                for eq in discretised_eq:
                    discretisation_kernel.add_equation(eq)
                discretisation_kernel.update_block_datasets(block)
                type_of_eq.Kernels += [discretisation_kernel]
                return self.required_constituent_relations
            else:
                pass
            return self.required_constituent_relations

    def check_missing_constituent_relations(self, block, list_of_eq):
        """ Checks that there are no missing constituent relations kernels. These are CRs
        that don't include derivatives."""
        arrays = []
        for eq in flatten(list_of_eq):
            arrays += list(eq.atoms(DataSet))
        arrays = set(arrays)
        undefined = arrays.difference(self.required_constituent_relations.keys())
        # Find the user's constituent relations
        CR = [eqn_class for eqn_class in block.list_of_equation_classes if isinstance(eqn_class, ConstituentRelations)]
        CR_LHS = set([str(eqn.lhs) for eqn in flatten(CR[0].equations)])

        for dset in undefined:
            if str(dset) in CR_LHS:
                self.required_constituent_relations[dset] = Kernel(block, computation_name="CR%s" % dset)
                self.required_constituent_relations[dset].set_grid_range(block)
        return

    def generate_first_derivatives(self, equations, block):

        self.local_kernels = OrderedDict()
        self.required_constituent_relations = {}
        work_arrays = {}
        derivatives = []

        for dire in range(block.ndim):
            derivatives += self.derivatives_to_store[dire]

        for no, der in enumerate(derivatives):
            der.update_work(block)
            work_arrays[der] = der.work
            ker = Kernel(block, kernel_name='FD_%d' % no)
            ker.set_computation_name("Derivative evaluation %s " % (der))
            self.update_range_of_constituent_relations(der, block)
            v = der
            expr = OpenSBLIEq(v.work, v._discretise_derivative(self, block))
            if self.merged: # combine the branch conditions
                expr = self.merge_conditionals([expr], block)
            ker.add_equation(expr)
            ker.set_grid_range(block)
            self.local_kernels[v] = ker
            derivatives[no] = v

        # Get all the derivatives and traverse them to update the range of already evaluated derivatives
        for d in self.get_local_function(equations):
            if d.atoms(CentralDerivative).intersection(derivatives):
                expr, self.local_kernels = self.traverse(d, self.local_kernels, block)

        # Combine the kernels per direction?
        # if self.group_stored:
        #     store_local = self.local_kernels
        #     self.local_kernels = {}
        #     for dire in range(block.ndim):
        #         # Get the kernels in this direction
        #         kernels = [store_local[d] for d in self.derivatives_to_store[dire]]
        #         equations = flatten([ker.equations for ker in kernels])
        #         # Group conditionals if needed
        #         equations = flatten(self.merge_conditionals(equations, block))
        #         # Create new merged kernel
        #         ker = Kernel(block)
        #         ker.set_computation_name("StoreSome evaluations direction %d" % (dire))
        #         # ker.update_kernel_name('Testing')
        #         ker.add_equation(equations)
        #         # self.update_range_of_constituent_relations(der, block)
        #         ker.set_grid_range(block)
        #         self.local_kernels[dire] = ker    
        return work_arrays


    def sbli_rhs_discretisation(self, type_of_eq, block):
        """ Function that performs the discretisation to the convective and viscous terms for the Navier-Stokes equations."""
        block.reset_work_index
        equations = flatten(type_of_eq.equations)
        residual_arrays = [eq.residual for eq in equations]
        equations = [e._sanitise_equation for e in equations]
        coordinates = set()
        for e in equations:
            coordinates = coordinates.union(e.atoms(CoordinateObject))
        # User-specified first derivatives to store to global 2D/3D arrays
        self.generate_derivatives_to_store(coordinates, block)

        # Compute the first derivatives to store
        work_arrays = self.generate_first_derivatives(equations, block)

        # Update the main equations to solve using the local evaluation grid variables
        equations = [e.subs(work_arrays) for e in equations]
        rhs_eq = [e.rhs for e in equations]
        discretised_eq = [OpenSBLIEq(x, y) for x, y in zip(residual_arrays, rhs_eq)]

        # Distinguish between convective and viscous terms using Reynolds number parameter. Is there a better way?
        classify_parameter = ConstantObject("Re")
        viscous, convective = self.classify_equations_on_parameter(discretised_eq, classify_parameter)
        # Apply to the convective terms
        convective = [OpenSBLIEq(x, y) for x, y in zip(residual_arrays, convective)]
        convective_equations = self.SS(convective, block, 'Convective')
        # Apply to the viscous terms
        viscous = [OpenSBLIEq(x, x+y) for x, y in zip(residual_arrays, viscous)]
        viscous_equations = self.SS(viscous, block, 'Viscous')
        # Remove any non equations
        if viscous_equations is not None:
            viscous_equations = [x for x in viscous_equations if isinstance(x, OpenSBLIEq)]
        # Group conditionals for vectorisation into a grouped piecewise object instead
        if self.merged:
            if convective_equations is not None:
                convective_equations = self.merge_conditionals(convective_equations, block)
            if viscous_equations is not None:
                viscous_equations = self.merge_conditionals(viscous_equations, block)

        if convective_equations or viscous_equations:
            for der, ker in self.local_kernels.items():
                type_of_eq.Kernels += [ker]

        if convective_equations:
            # Factor common constants
            if self.factor:
                print("StoreSome: Factoring convective equations.")
                saving, total = 0, 0
                for i, eqn in enumerate(convective_equations):
                    lhs, rhs = eqn.lhs, eqn.rhs
                    before = count_ops(rhs)
                    total += before
                    rhs = factor_terms(eqn.rhs, clear=True)
                    after = count_ops(rhs)
                    # print("Operations - Before: {} After: {} Saving: {}".format(before, after, ))
                    if after < before:
                        saving += before - after
                        convective_equations[i] = OpenSBLIEq(lhs, rhs, evaluate=False)
                print("Reduced number of operations by: {:.2f}%.".format(saving/total * 100))
            convective_kernel = Kernel(block, computation_name="Convective terms")
            convective_kernel.set_grid_range(block)
            for eq in convective_equations:
                convective_kernel.add_equation(eq)
            convective_kernel.update_block_datasets(block)
            type_of_eq.Kernels += [convective_kernel]
        if viscous_equations:
            # Factor common constants
            if self.factor:
                print("StoreSome: Factoring viscous equations.")
                saving, total = 0, 0 
                for i, eqn in enumerate(viscous_equations):
                    lhs, rhs = eqn.lhs, eqn.rhs
                    before = count_ops(rhs)
                    total += before
                    rhs = factor_terms(eqn.rhs, clear=True)
                    after = count_ops(rhs)
                    # print("Operations - Before: {} After: {} Saving: {}".format(before, after, before - after))
                    if after < before:
                        saving += before - after
                        viscous_equations[i] = OpenSBLIEq(lhs, rhs, evaluate=False)
                print("Reduced number of operations by: {:.2f}%.".format(saving/total * 100))
            viscous_kernel = Kernel(block, computation_name="Viscous terms")
            viscous_kernel.set_grid_range(block)
            for eq in viscous_equations:
                viscous_kernel.add_equation(eq)
            viscous_kernel.update_block_datasets(block)
            type_of_eq.Kernels += [viscous_kernel]
        block.reset_work_index
        # Check for missing constituent relation kernels
        self.check_missing_constituent_relations(block, equations)
        return self.required_constituent_relations

    def merge_conditionals(self, input_equations, block):
        """ Optimisation to enable vectorisation by grouping the conditional expressions per direction (x, y, z)."""
        output_equations = []

        start_order = 0
        for direction in range(block.ndim):
            no_condition = dict() # no branches
            conditionals_d1 = dict() # first derivatives
            conditionals_d2 = dict() # second derivatives
            factor_dict = dict()
            for order, eqn in enumerate(input_equations):
                # Find conditional expressions and group them together based on their if condition
                if len(eqn.rhs.atoms(Piecewise)) > 0:
                    # Find the direction
                    idx = list(eqn.rhs.atoms(Grididx))
                    if len(idx) == 1:
                        pw_dire = idx[0].args[-1]
                    else:
                        # Compatibility path: some SymPy versions don't preserve Grididx
                        # in Piecewise conditions; infer direction from condition symbols.
                        pw_dire = None
                        for pw in eqn.rhs.atoms(Piecewise):
                            for pair in pw.args:
                                cond = pair.args[1]
                                for sym in getattr(cond, "free_symbols", set()):
                                    match = re.search(r"\[(\d+)\]$", str(sym))
                                    if match:
                                        pw_dire = int(match.group(1))
                                        break
                                if pw_dire is not None:
                                    break
                            if pw_dire is not None:
                                break
                        if pw_dire is None:
                            continue
                    if pw_dire == direction:
                        for const in eqn.rhs.atoms(ConstantObject):
                            if not isinstance(const, ConstantIndexed):
                                if const.rational:
                                    factor = const
                        lhs = eqn.lhs
                        pw = list(eqn.rhs.atoms(Piecewise))[0]
                        factor_dict[lhs] = factor
                        if 'd1_' in str(eqn.lhs):
                            for pair in pw.args:
                                expr, cond = pair.args[0], pair.args[1]
                                if cond in conditionals_d1:
                                    conditionals_d1[cond].append((lhs, expr, order+start_order))
                                else:
                                    conditionals_d1[cond] = [(lhs, expr, order+start_order)]
                        else:
                            for pair in pw.args:
                                expr, cond = pair.args[0], pair.args[1]
                                if cond in conditionals_d2:
                                    conditionals_d2[cond].append((lhs, expr, order+start_order))
                                else:
                                    conditionals_d2[cond] = [(lhs, expr, order+start_order)]
                else:
                    if direction == block.ndim - 1:
                        no_condition[order+start_order] = eqn

            # Processs conditional equations if they exist for first derivatives
            if len(conditionals_d1) > 0:
                # Construct the grouped equations
                conditions_to_evaluate = []
                for cond, exprs in conditionals_d1.items():
                    evaluations = []
                    if str(cond) != 'True': # avoid true condition until the end
                        for triple_value in exprs:
                            lhs, rhs, order = triple_value
                            factor = factor_dict[lhs]
                            evaluations.append(OpenSBLIEq(lhs, factor*rhs))
                        conditions_to_evaluate.append(GroupedCondition(Tuple(*evaluations), Equality(cond.lhs, cond.rhs)))
                # Add the default condition
                evaluations = []
                for triple_value in conditionals_d1[True]:
                    lhs, rhs, order = triple_value
                    factor = factor_dict[lhs]
                    evaluations.append(OpenSBLIEq(lhs, factor*rhs))
                conditions_to_evaluate.append(GroupedCondition(Tuple(*evaluations), True))
                # Create a GroupedPiecewise evalation object
                grouped = GroupedPiecewise(*conditions_to_evaluate)
                output_equations += [grouped]

            if len(conditionals_d2) > 0:
                # Construct the grouped equations
                conditions_to_evaluate = []
                for cond, exprs in conditionals_d2.items():
                    evaluations = []
                    if str(cond) != 'True': # avoid true condition until the end
                        for triple_value in exprs:
                            lhs, rhs, order = triple_value
                            factor = factor_dict[lhs]
                            evaluations.append(OpenSBLIEq(lhs, factor*rhs))
                        conditions_to_evaluate.append(GroupedCondition(Tuple(*evaluations), Equality(cond.lhs, cond.rhs)))
                # Add the default condition
                evaluations = []
                for triple_value in conditionals_d2[True]:
                    lhs, rhs, order = triple_value
                    factor = factor_dict[lhs]
                    evaluations.append(OpenSBLIEq(lhs, factor*rhs))
                conditions_to_evaluate.append(GroupedCondition(Tuple(*evaluations), True))
                # Create a GroupedPiecewise evalation object
                grouped = GroupedPiecewise(*conditions_to_evaluate)
                output_equations += [grouped]

            # Add the equations which have no branching conditions
            for key, val in no_condition.items():
                output_equations.append(val)
            start_order += order
        return output_equations

    def SS(self, type_of_eq, block, equation_type, group=True, level=1):
        """ Generates the GridVariables to store the local derivatives. Also sorts the
        derivatives to re-access the same arrays consecutively."""
        directions = ['x', 'y', 'z']
        if not isinstance(type_of_eq, list):
            equations = type_of_eq.equations
        else:
            equations = type_of_eq[:]

        discrete_equations = flatten(equations)[:]
        cds = self.get_local_function(flatten(equations))
        grid_variable_evaluations = []
        names = []
        if cds:
            # Sort by grouping variables
            if group:
                # Sort derivative evaluations by direction (x, y, z) and then by name
                if equation_type == 'Convective':
                    cds = sorted(cds, key=lambda x: (x.args[1].direction, str(x.args[0])))
                elif equation_type == 'Viscous':
                    cds = sorted(cds, key=lambda x: (x.args[1].direction, str(x.args[0]), x.order))
            for i, der in enumerate(cds):
                self.update_range_of_constituent_relations(der, block)
                if level == 1:
                    var_name = self.generate_name(der, block, i)
                    if len(der.args) == 2:
                        gv = GridVariable('d1_%s_d%s' % (var_name, directions[der.args[1].direction]))
                    elif len(der.args) == 3:
                        gv = GridVariable('d2_%s_d%s' % (var_name, directions[der.args[1].direction]))
                    else:
                        raise ValueError("Only first and second derivatives are supported in StoreSome.")
                    # Evaluate the expression and assign to the local grid variable
                    grid_variable_evaluations += [OpenSBLIEq(gv, der._discretise_derivative(self, block, type_of_eq=type_of_eq))]
                    if str(gv) in names:
                        print("WARNING: Duplicated derivative {} in StoreSome due to using old Skew operator.".format(str(gv)))
                    # assert str(gv) not in names # no repeated grid variable names
                    names.append(str(gv))
                    for no, c in enumerate(discrete_equations):
                        discrete_equations[no] = discrete_equations[no].subs(der, gv)
            # Check the each input derivative received a local grid variable
            assert len(cds) == len(names)
            return grid_variable_evaluations+discrete_equations
        else:
            return None

    def generate_name(self, der, block, identity):
        # Make a name for the local derivative evaluation
        input_vars = []
        input_args = der.args[0]
        dsets = sorted(list(input_args.atoms(DataSet)), key=lambda x: str(x))

        # Check for repeated variables as power
        if len(input_args.atoms(Pow)) == 0:
            var_name = ''.join([str(x).split('_B%d' % block.blocknumber)[0] for x in dsets])
        if len(input_args.atoms(Pow)) > 0:
            npow = list(input_args.atoms(Pow))[0].args[-1]
            if npow > 0:
                repeated = [(npow-1)*str(list(input_args.atoms(Pow))[0]).split('_B%d' % block.blocknumber)[0]]
            else:
                repeated = ['inv_']
            var_name = ''.join(repeated + [str(x).split('_B%d' % block.blocknumber)[0] for x in dsets])
        # var_name = str(identity) + '_' + var_name
        return var_name
