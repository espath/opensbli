from opensbli.multiblock.blockcollection import MultiBlock as MB
from sympy import flatten, pprint, Idx, Equality, Or
from opensbli.code_generation.latex import LatexWriter
from opensbli.core.kernel import ConstantsToDeclare as CTD
from opensbli.equation_types.opensbliequations import SimulationEquations, NonSimulationEquations, ConstituentRelations
from opensbli.equation_types.metric import MetricsEquation
from opensbli.core.opensbliobjects import Constant, DataSetBase, ConstantObject
from opensbli.core.datatypes import Int
from opensbli.code_generation.algorithm.common import BeforeSimulationStarts, AfterSimulationEnds, InTheSimulation
from opensbli.code_generation.algorithm import MainPrg, Condition, Condition, MainPrg, Loop, DoLoop, Timers, DefDecs
import copy
from opensbli.core.datatypes import SimulationDataType


class BlockDescription(object):
    def __init__(self, block):
        block.copy_block_attributes(self)
        return

class TraditionalAlgorithmRKMB(object):
    """ It is where the algorithm is generated, This is a seperate layer
    which gives user control to do any modifications for extra functionality that
    is to be performed like, doing some post processing for every time loop or
    sub rk loop. """

    def __init__(self, blocks, simulation_monitor=None):
        self.block_descriptions = []
        self.ntimers = 0
        self.MultiBlock = True
        self.simulation_monitor = simulation_monitor
        # For restart flag in the definitions to restart the time-advance arrays
        self.time_advance_arrays = []
        from opensbli.equation_types.opensbliequations import SimulationEquations
        for b in flatten([blocks.blocks]):
            for eqn_class in b.list_of_equation_classes:
                if isinstance(eqn_class, SimulationEquations):
                    self.time_advance_arrays += flatten(eqn_class.time_advance_arrays)
        self.datatype = SimulationDataType.dtype()
        self.check_temporal_scheme(blocks)
        self.prg = MainPrg()
        self.add_block_names(blocks)
        defdecs = self.get_definitions_declarations(blocks)
        self.definitions_and_declarations = defdecs
        self.spatial_solution(blocks)
        return

    def add_block_names(self, blocks):
        if self.MultiBlock:
            for b in range(blocks.nblocks):
                self.block_descriptions += [BlockDescription(blocks.get_block(b))]
        return

    def add_def_decs(self, defdecs):
        self.prg.add_components(defdecs)
        return

    def get_definitions_declarations(self, blocks):
        defdecs = DefDecs()
        if self.MultiBlock:
            for block_number in range(blocks.nblocks):
                b = blocks.get_block(block_number)
                # defdecs.add_components(list(b.constants.values()))
                defdecs.add_components(list(b.Rational_constants.values()))
                defdecs.add_components(list(b.block_datasets.values()))
                defdecs.add_components(list(b.block_stencils.values()))
                defdecs.add_components(list(b.block_reductions.values()))
        return defdecs

    def comapre_no_sims(self, s1, s2):
        return cmp(s1.order, s2.order)
    
    def spatial_solution(self, blocks):
        """ Add the spatial kernels to the temporal solution i.e temporalscheme.solution. """ 
        print("Writing algorithm")
        fname = 'algorithm.tex'
        latex = LatexWriter()
        latex.open(fname, "Algorithm for the equations")
        if self.MultiBlock:
            bc_kernels, inner_temporal_advance_kernels, temporal_start, temporal_end, spatial_kernels = [], [], [], [], []
            before_time, after_time, in_time, non_simulation_eqs = [], [], [], []
            metrics, tloop_blocks, inner_loop_blocks = [], [], []
            # Loop over the multiple blocks in turn
            for block_number in range(blocks.nblocks):
                b = blocks.get_block(block_number)
                for sc in b.get_temporal_schemes:
                    # Iteration counter for any conditional expressions
                    temporal_iteration = sc.temporal_iteration
                    temporal_iteration.main_file = True
                    inner_loop_blocks += [sc.stage]
                    tloop_blocks += [temporal_iteration]
                    for key, value in iter(sc.solution.items()):
                        #print(key, value)
                        if isinstance(key, SimulationEquations):
                            # Solution advancement kernels
                            temporal_start += sc.solution[key].start_kernels
                            temporal_end += sc.solution[key].end_kernels
                            inner_temporal_advance_kernels += sc.solution[key].kernels
                            bc_kernels += key.boundary_kernels
                            spatial_kernels += key.all_spatial_kernels(b)
                        elif isinstance(key, MetricsEquation):
                            metrics += [key]
                        elif isinstance(key, NonSimulationEquations):  # Add all other types of equations
                            non_simulation_eqs += [key]
                        else:
                            if not isinstance(key, ConstituentRelations):
                                print("NOT classified", type(key))
                                raise ValueError("Equations class can not be classified: %s" % key)
            # Place any non-simulation equation classes (statistics, filters, metric evaluations, ...)
            for key in sorted(non_simulation_eqs, key=lambda x: x.order):
                for place in key.algorithm_place:
                    if isinstance(place, BeforeSimulationStarts):
                        if place.start_condition is not None:
                            cond = Condition(place.start_condition)
                            cond.add_components(key.Kernels)
                            before_time += [cond]
                        else:
                            before_time += key.Kernels
                    elif isinstance(place, AfterSimulationEnds):
                        after_time += key.Kernels
                    else:
                        if place.frequency: # iteration frequency condition
                            t = Equality((temporal_iteration + 1) % key._place[0].frequency, 0)
                            cond = Condition(t)
                            cond.add_components(key.Kernels)
                            in_time += [cond]
                        elif place.execution_condition is not None:  # boolean condition
                            cond = Condition(place.execution_condition)
                            cond.add_components(key.Kernels)
                            in_time += [cond]
                        else: # no condition, always evaluate
                            in_time += key.Kernels

            # Add optional simulation monitors
            if self.simulation_monitor is not None:
                t = (Or(Equality((temporal_iteration + 1) % self.simulation_monitor.frequency, 0), Equality(temporal_iteration, 0)))
                cond = Condition(t)
                cond.add_components(self.simulation_monitor)
                in_time += [cond]

            # Process the metrics we will control here it self later we will move this to multi block
            # The first derivatives this includes bc application
            for m in metrics:
                before_time += m.fd_kernels
            # Now the second derivatives
            for m in metrics:
                before_time += m.sd_kernels
            #tloop_blocks = set(tloop_blocks)
            #inner_loop_blocks = set(inner_loop_blocks)
            if len(set(tloop_blocks)) != 1 or len(set(inner_loop_blocks)) != 1:
                raise ValueError("")
            # RK time stepping loop
            innerloop = DoLoop(inner_loop_blocks[0])
            innerloop.add_components(spatial_kernels + inner_temporal_advance_kernels + bc_kernels)

            print("Found %d kernels." % len(spatial_kernels + inner_temporal_advance_kernels + bc_kernels + in_time))
            ## Add BC kernels to temporal start
            temporal_start = bc_kernels + temporal_start
            #temporal_iteration = sc.temporal_iteration
            tloop = DoLoop(tloop_blocks[0])
            self.get_io_mb(blocks, before_time, in_time, after_time, tloop_blocks[0])
            tloop.add_components(temporal_start)
            tloop.add_components(innerloop)
            tloop.add_components(in_time)
            tloop.add_components(temporal_end)
            # Create a time loop and add the components to it
            timed_tloop = self.add_timers(tloop)
            self.prg.add_components(before_time)
            self.prg.add_components(timed_tloop)
            self.prg.add_components(after_time)
        latex.close()
        return
    
    def get_io_mb(self, blocks, before_time, in_time, after_time, temporal_iteration):
        for block_number in range(blocks.nblocks):
            b = blocks.get_block(block_number)
            for io in b.InputOutput:
                io.block_number = block_number # MBCHANGE
                # Check all of the datasets added to the IO classes have been defined elsewhere in the simulation
                io.check_datasets(b)
                for place in io.algorithm_place:
                    if isinstance(place, BeforeSimulationStarts):
                        io_copy = copy.deepcopy(io)
                        io_copy.dynamic_fname = False
                        before_time += [io_copy]
                    elif isinstance(place, AfterSimulationEnds):
                        io_copy = copy.deepcopy(io)
                        io_copy.dynamic_fname = False
                        after_time += [io_copy]
                    elif isinstance(place, InTheSimulation):
                        if place.initial_condition:
                            t = Or(Equality((temporal_iteration + 1) % place.frequency, 0), Equality(temporal_iteration, 0))
                        else:
                            t = (Equality((temporal_iteration + 1) % place.frequency, 0, evaluate=False))
                        cond = Condition(t)
                        io_copy = copy.deepcopy(io)
                        io_copy.dynamic_fname = True
                        io_copy.control_parameter = temporal_iteration + 1
                        cond.add_components(io_copy)
                        in_time += [cond]
                    else:
                        raise NotImplementedError("In Nonsimulation equations")
        return 

    def add_timers(self, components):
        timer = Timers(self.ntimers)
        self.ntimers += 1
        timer.add_components(components)
        return timer

    def check_temporal_scheme(self, blocks):
        """ If Multi-block this checks the temporal scheme is the same for all the blocks."""
        if self.MultiBlock:
            sc = set()
            for b in range(blocks.nblocks):
                for s in blocks.get_block(b).get_temporal_schemes:
                    sc.add(s.__class__.__name__)
            if len(sc) != 1:
                raise ValueError("More than one temporal scheme")
        else:
            if len(blocks[0].get_temporal_schemes) > 1:
                raise ValueError("More than one temporal scheme for a block")
        return
