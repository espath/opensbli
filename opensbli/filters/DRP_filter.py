""" David J. Lusher 09/21. Dispersion Relation Preserving (DRP) explicit filters."""

from opensbli import *
from sympy import pprint, Piecewise, factor, Or, simplify, And
from opensbli.core.opensbliobjects import DataObject, ConstantObject, GroupedPiecewise
from opensbli.equation_types.opensbliequations import OpenSBLIEquation
from opensbli.postprocess.post_process_eq import *
from opensbli.code_generation.algorithm.common import *
from opensbli.utilities.user_defined_kernels import UserDefinedEquations
from opensbli.core.kernel import ConstantsToDeclare as CTD


class ExplicitFilter(object):
    """ Selective filtering from Bogey & Bailly, A family of low dispersive and low dissipative explicit
    schemes for flow and noise computations, JoCP (2004) 194-214."""
    def __init__(self, block, filter_directions, filter_type='DRP', q_vector = None, width=11, frequency=25, optimized=False, sigma=0.3333333, airfoil=False, multi_block=False):
        self.width, self.optimized = width, optimized
        directions = ['x', 'y', 'z']
        print("Using a %s filter with stencil width %d for block %d, in directions: %s." % (filter_type, self.width, block.blocknumber, [directions[x] for x in filter_directions]))
        self.depth = int(width/2.0)
        self.ndim = block.ndim
        self.block = block
        self.filter_directions = filter_directions
        self.wall_boundaries = [[False, False] for _ in range(self.ndim)]
        self.airfoil = airfoil
        if multi_block:
            self.nblocks = multi_block.nblocks
        else:
            self.nblocks = 1
        for x in filter_directions:
            assert isinstance(x, int)
        self.filter_type = filter_type
        # Arrays to filter
        # Conservative variables
        if q_vector == None:
            if block.conservative:
                if self.ndim == 2:
                    q = ['rho', 'rhou0', 'rhou1', 'rhoE']
                elif self.ndim == 3:
                    q = ['rho', 'rhou0', 'rhou1', 'rhou2', 'rhoE']
                self.q_vector = [block.location_dataset(x) for x in flatten(q)]
            else:
                if self.ndim == 2:
                    q = ['rho', 'u0', 'u1', 'Et']
                elif self.ndim == 3:
                    q = ['rho', 'u0', 'u1', 'u2', 'Et']
                self.q_vector = [block.location_dataset('rho')] + [block.location_dataset('rho')*block.location_dataset('%s' % x) for x in q[1:]]
                self.lhs = [block.location_dataset(x) for x in q]
        else: # User specifies which time-advance quantities need to be filtered
            q = flatten(q_vector)
            self.q_vector = [block.location_dataset(x) for x in q]
            print(self.q_vector)
            # exit()
        self.temp_arrays = [block.location_dataset('%s_RKold' % x) for x in q]
        self.freq = ConstantObject('filter_frequency')
        self.freq.value = frequency
        self.freq.datatype = Int()
        CTD.add_constant(self.freq)
        self.sigma = ConstantObject('DRP_filt')
        self.sigma.value = sigma
        # Generate the filter offset grid locations
        self.locations = [i for i in range(-int(self.width/2.0), int(self.width/2.0)+1)]
        # Generate the coefficients
        if self.filter_type == 'DRP':
            self.generate_DRP_weights()
        else:
            self.generate_Visbal_weights()
        # Create the filter equations
        self.create_filter(block)
        return

    def detect_wall_boundaries(self):
        """ The shock-filter is turned off in the near-wall region. This function detects which directions, if any, have
        wall boundary conditions."""
        try:
            for direction in range(self.ndim):
                for side in [0,1]:
                    if isinstance(self.block.boundary_types[direction][side], WallBC):
                        self.wall_boundaries[direction][side] = True
        except:
            raise ValueError("Please set boundary conditions on the block before calling the shock filter.")
        return

    def detect_interface_boundaries(self):
        """ The shock-filter is turned off close to block interfaces. This function detects which directions, if any, have
        interface boundary conditions."""
        self.interface_boundaries = [[False, False] for _ in range(self.ndim)]
        try:
            for direction in range(self.ndim):
                for side in [0,1]:
                    if isinstance(self.block.boundary_types[direction][side], InterfaceBC) or isinstance(self.block.boundary_types[direction][side], SharedInterfaceBC):
                        self.interface_boundaries[direction][side] = True
        except:
            raise ValueError("Please set boundary conditions on the block before calling the shock filter.")
        return

    def boundary_check(self, block):
        """ Checks if there are any non-periodic boundaries requirinig narrowing of the filter stencil."""
        self.non_periodic = [[False, False] for _ in range(block.ndim)]
        self.modify_directions = [False for _ in range(block.ndim)]
        for i, bc in enumerate(block.boundary_types):
            for side in [0,1]:
                if isinstance(bc[side], ModifyCentralDerivative):
                    self.non_periodic[i][side] = True
                    self.modify_directions[i] = True
        return

    def generate_DRP_weights(self):
        """ Weights are symmetric about the central point."""
        if self.width == 9:
            if self.optimized:
                self.weights = [0.008228661760,-0.045211119360,0.120007591680,-0.204788880640]
                self.weights += [0.243527493120] + self.weights[::-1]
            else:
                self.weights = [Rational(1,256),Rational(-1,32),Rational(7,64),Rational(-7,32)]
                self.weights += [Rational(35,128)] + self.weights[::-1]
        elif self.width == 11:
            if self.optimized:
                self.weights = [-0.002999540835,0.018721609157,-0.059227575576,0.123755948787,-0.187772883589]
                self.weights += [0.215044884112] + self.weights[::-1]
            else:
                self.weights = [Rational(-1,1024),Rational(5,512),Rational(-45,1024),Rational(15,128),Rational(-105,512)]
                self.weights += [Rational(63,256)] + self.weights[::-1]
        elif self.width == 13:
            raise ValueError("Six halos are required to use this filter width, update in Central Def_decs halos.")
            if self.optimized:
                self.weights = [0.001254597714,-0.008520738659,0.029662754736,-0.069975429105,0.123632891797,-0.171503832236]
                self.weights += [0.190899511506] + self.weights[::-1]
            else:
                self.weights = [Rational(1,4096),Rational(-3,1024),Rational(33,2048),Rational(-55,1024),Rational(495,4096),Rational(-99,512)]
                self.weights += [Rational(231,1024)] + self.weights[::-1]        
        return

    def generate_Visbal_weights(self):
        """ Weights are symmetric about the central point. Taken from M. Visbal, D. Gaitonde,
        On the use of higher-order finite-difference schemes on curvilinear and deforming meshes. JoCP 181, 155-185 (2002)."""
        if self.width == 3:
            self.weights = [Rational(1,2)] + [Rational(1,2)]
        elif self.width == 5:
            self.weights = [Rational(5,8)] + [Rational(-1,8), Rational(1,2)][::-1]
        elif self.width == 7:
            self.weights = [Rational(11,16)] + [Rational(1,32), Rational(-3,16), Rational(15,32)][::-1]
        elif self.width == 9:
            self.weights = [Rational(93,128)] + [Rational(-1,128), Rational(1,16), Rational(-7,32), Rational(7,16)][::-1]
        elif self.width == 11:
            self.weights = [Rational(193,256)] + [Rational(1,512), Rational(-5,256), Rational(45,512), Rational(-15,64), Rational(105,256)][::-1]
        return


    def create_stencil(self, direction):
        """ Indexes the datasets based on the width of the filter stencil."""
        output = []
        # Create conditions for adaptive filtering strength
        if self.airfoil: # Targeted adaptive filter for airfoil buffet cases
        # Spatially control the strength to reduce filter in the boundary-layer
            st = GridVariable('strength')
            error_indicator = ConstantObject('error_indicator')
            error_indicator.value = 0.5
            CTD.add_constant(error_indicator)
            locations = self.locations
            fmax = increment_dataset(self.Ren, direction, locations[0])
            for loc in locations[1:]:
                fmax = Max(increment_dataset(self.Ren, direction, loc), fmax)
            # Targeted filter to turn the filter off in the boundary-layer
            output += [OpenSBLIEq(gv('Ren'), fmax)]
            inner1 = ExprCondPair(3.0, gv('Ren') >= error_indicator)
            inner2 = ExprCondPair(0.01, True)
            inner_pw = OpenSBLIEq(st, Piecewise(*[inner1, inner2]))
            output += [GroupedPiecewise(ExprCondPair(inner_pw, Equality(self.filter_mask, 0)), ExprCondPair(OpenSBLIEq(st, 1.0), True))]
        else:
            # Regular application
            st = 1

        if self.filter_type == 'DRP':
            for dset_id, dset in enumerate(self.q_vector):
                stencil = []
                for i, location in enumerate(self.locations):
                    stencil.append(self.weights[i]*increment_dataset(dset, direction, location))
                output += [OpenSBLIEq(self.temp_arrays[dset_id], st*simplify(sum(stencil)))]
        elif self.filter_type == 'Visbal':
            for dset_id, dset in enumerate(self.q_vector):
                stencil = []
                total = 0
                for i in range(len(self.weights)):
                    total += Rational(1,2)*self.weights[i]*(increment_dataset(dset, direction, +i) + increment_dataset(dset, direction, -i))
                output += [OpenSBLIEq(self.temp_arrays[dset_id], st*factor(total))]
        else:
            raise ValueError("Wrong type of explicit filter specified.")
        # Restrict the filter if close to the wall, in the wall normal direction
        eqns = []
        buffer = self.depth
        if self.wall_boundaries[direction][0] or self.wall_boundaries[direction][1]:
            if self.wall_boundaries[direction][0]:
                check = self.block.grid_indexes[direction] < buffer
            if self.wall_boundaries[direction][1]:
                check = self.block.grid_indexes[direction] > self.block.ranges[direction][1] - (buffer+1)
            if self.wall_boundaries[direction][0] and self.wall_boundaries[direction][1]:
                check = Or(self.block.grid_indexes[direction] < buffer, self.block.grid_indexes[direction] > self.block.ranges[direction][1] - (buffer+1))
            zeroing = [OpenSBLIEq(self.temp_arrays[dset_id], 0.0) for dset_id, dset in enumerate(self.q_vector)]
            cond1 = ExprCondPair(zeroing, check)
            cond2 = ExprCondPair(output, True)
            eqns = [GroupedPiecewise(cond1, cond2)]
        else:
            eqns = output[:]
        return eqns

    def zero_temp_arrays(self):
        """ Ensure the temp arrays are zero everywhere."""
        zeroed = []
        for dset_id, dset in enumerate(self.temp_arrays):
            zeroed.append(OpenSBLIEq(dset, 0.0))
        return zeroed

    def create_equations(self, block, direction):
        # Create the indexed equations to calculate the filter
        application = self.create_stencil(direction)
        direction += 1
        # Update the q vector
        update = []
        if block.conservative:
            for dset_id, dset in enumerate(self.q_vector):
                if self.filter_type == 'DRP':
                    update += [OpenSBLIEq(dset, dset - self.sigma*self.temp_arrays[dset_id])]
                elif self.filter_type == 'Visbal':
                    update += [OpenSBLIEq(dset, dset - self.sigma*(dset - self.temp_arrays[dset_id]))]
                else:
                    raise ValueError("Wrong type of explicit filter specified.")
        else:
            if self.filter_type == 'DRP':
                update += [OpenSBLIEq(self.lhs[0], self.lhs[0] - self.sigma*self.temp_arrays[0])]
                inv_rho = GridVariable('inv_rho')
                update += [OpenSBLIEq(inv_rho, 1.0/self.lhs[0])]
                for dset_id, dset in enumerate(self.lhs[1:]):
                    update += [OpenSBLIEq(dset, dset - self.sigma*self.temp_arrays[dset_id+1]*inv_rho)]
            else:
                update += [OpenSBLIEq(self.lhs[0], self.lhs[0] - self.sigma*(self.lhs[0] - self.temp_arrays[0]))]
                inv_rho = GridVariable('inv_rho')
                update += [OpenSBLIEq(inv_rho, 1.0/self.lhs[0])]
                for dset_id, dset in enumerate(self.lhs[1:]):
                    update += [OpenSBLIEq(dset, dset - self.sigma*(dset - self.temp_arrays[dset_id+1]*inv_rho))]
        return application, update


    def pressure_correction(self, block, order):
        # Pressure gradient sensor for targeted filtering
        SS = ShockSensor()
        output_eqns, kappa = SS.Ren_sensor(block, name='kappa')
        Ren_output = output_eqns[-1].rhs
        self.Ren = block.location_dataset('Ren_sensor')
        output_eqns = [OpenSBLIEq(self.Ren, Ren_output)]
        UDF = UserDefinedEquations()
        UDF.algorithm_place = InTheSimulation(frequency=self.freq)
        UDF.computation_name = 'Block %d: pressure sensor evaluation' % (block.blocknumber)
        UDF.order = order
        UDF.add_equations(output_eqns)
        return UDF

    def create_UDF(self, block, equations, direction, order, UDF_type):
        UDF = UserDefinedEquations()
        UDF.algorithm_place = InTheSimulation(frequency=self.freq)
        if order == 0:
            # Mark as an explicit filter, to be used for full halo swaps
            UDF.full_swap = True
        if UDF_type == 'Zeroing':
            UDF.computation_name = 'Block %d: Zero the filter array' % block.blocknumber
        elif UDF_type == 'Calculation':
            UDF.computation_name = 'Block %d: %s filter calculation direction %s' % (block.blocknumber, self.filter_type, block.direction_labels[direction])
        elif UDF_type == 'Update':
            UDF.computation_name = 'Block %d: %s filter update direction %s' % (block.blocknumber, self.filter_type, block.direction_labels[direction])
        else:
            raise ValueError("The UDF should be one of the above actions.")
        # Ordering of the filter operations to fix the order of the kernel calls in the code
        UDF.order = order
        UDF.add_equations(equations)
        # Turn off filtering in the boundary-layer if required for airfoil cases
        if self.airfoil and self.wall_boundaries[1][0]:
            self.reduce_grid_range(block, UDF)
        return UDF

    def reduce_grid_range(self, block, filt_class):
        original = copy.deepcopy(block.ranges)
        direction = 1
        start = 1
        # original_start = original[direction][0]
        # Edit the start_index, currently assumes the filter should be applied to the end of the iteration range in that direction
        original[direction][0] = start
        filt_class.custom_grid_range = original
        return

    def airfoil_mask(self, block):
        """ Mask the spatial areas to apply filtering to."""
        UDF = UserDefinedEquations()
        UDF.algorithm_place = BeforeSimulationStarts()
        UDF.computation_name = 'Airfoil_filter_mask'
        # Filter regions
        x, y = block.location_dataset('x0'), block.location_dataset('x1')
        filter_condition = [ExprCondPair(0, And(Abs(y) < 0.5, x > 0.1))]
        filter_condition += [ExprCondPair(1, True)]
        self.filter_mask = block.location_dataset('filter_mask')
        output_eqns = [OpenSBLIEq(self.filter_mask, Piecewise(*filter_condition))]
        UDF.add_equations(output_eqns)
        UDF.order = 10000
        return UDF

    def create_filter(self, block):
        self.equation_classes = []
        # Zero the arrays
        zeroed = self.zero_temp_arrays() #### TODO: zero before each direction one by one
        self.equation_classes += [self.create_UDF(block, zeroed, 0, 0+block.blocknumber, 'Zeroing')]

        # Create a mask if airfoil problem
        if self.airfoil:
            mask_UDF = self.airfoil_mask(block)
            # Error indicator
            error_UDF = self.pressure_correction(block, order=block.blocknumber*10 + self.nblocks - 1) # evaluate before the filters
            self.equation_classes += [mask_UDF, error_UDF]

        # Check for non-periodic boundaries
        self.detect_wall_boundaries()
        # Create a kernel at the end of the time loop, every iteration (no frequency)
        start_number = block.blocknumber*10 + self.nblocks
        for direction in self.filter_directions:
            # Create the equations
            application, update = self.create_equations(block, direction)
            filter1 = self.create_UDF(block, application, direction, start_number, 'Calculation')
            start_number += 1
            filter2 = self.create_UDF(block, update, direction, start_number, 'Update')
            start_number += 1
            self.equation_classes += [filter1, filter2]
        return
