"""Reference paper: Richez, Leguille and Marquez, "Selective frequency damping method for steady RANS solutions of
turbulent separated flows around an airfoil at stall", Computer and Fluids, Vol. 132, pp. 51-61 (2016)

Also from: Akervik (2006)

chifilt : damping coefficient (normally 0.1)
omegafilt : cut frequency
dt: time step
vfilt: filtered solution
v : conservative variables

Ported to OpenSBLI (D.Lusher 05/2020). Original version by Dr. Andrea Sansica.
Updated with periodic reset function (D. Lusher 08/2023).
--------------------------------------------------------------------------------------------
Some notes:
    initialize vfilt = v at the first time iteration
    write out vfilt for restart with sfd."""

from opensbli import *
from sympy import symbols, exp, pprint, Mod
from opensbli.core.opensbliobjects import DataObject, ConstantObject
from opensbli.equation_types.opensbliequations import OpenSBLIEquation
from opensbli.postprocess.post_process_eq import *
from opensbli.core.kernel import ConstantsToDeclare as CTD
from opensbli.code_generation.algorithm.common import *
from opensbli.utilities.user_defined_kernels import UserDefinedEquations

class SFD(object):
    def __init__(self, block, chifilt=0.1, omegafilt=1.0/0.75, formulation='standard'):
        self.block = block
        self.formulation = formulation
        print("Using Selective Frequency Damping (SFD) on block {}".format(self.block.blocknumber))
        # Arrays for the filtered solution
        self.create_arrays()
        self.equation_classes = []
        # Create a kernel to initialize vfilt = v
        self.generate_initial_condition()
        # Create the filter equations
        self.create_filter(chifilt, omegafilt)
        return

    def create_arrays(self):
        ndim = self.block.ndim
        # Conservative variables
        if ndim == 2:
            cons_vars = ['rho', 'rhou0', 'rhou1', 'rhoE']
        elif ndim == 3:
            cons_vars = ['rho', 'rhou0', 'rhou1', 'rhou2', 'rhoE']
        self.cons_arrays = [DataObject('%s' % var) for var in cons_vars]
        # Filter variables
        if self.formulation == 'standard':
            self.f_old = [GridVariable('%s_SFDold' % var) for var in cons_vars]
        else:
            self.f_old = [DataObject('%s_SFDold' % var) for var in cons_vars]
        self.f_new = [DataObject('%s_SFD' % var) for var in cons_vars]
        return

    def generate_initial_condition(self):
        """ Initializes the conservative variables with the freestream conditions."""
        initial_class = UserDefinedEquations()
        initial_class.computation_name = 'Initialize the SFD temporal filter'
        initial_class.algorithm_place = BeforeSimulationStarts()
        # Ensure that the evaluation comes after the initial condition
        initial_class.order = 10000000000
        # Create the equations from the conservative variables
        initial_equations = [OpenSBLIEq(left, right) for (left, right) in zip(self.f_new, self.cons_arrays)]
        initial_class.add_equations(initial_equations)
        self.equation_classes.append(initial_class)
        return

    def create_filter(self, chifilt, omegafilt):
        # Create a kernel at the end of the time loop, every iteration (no frequency)
        filter_class = UserDefinedEquations()
        apply_SFD = ConstantObject('apply_SFD', integer=True)
        apply_SFD.value = 1
        CTD.add_constant(apply_SFD)
        filter_class.algorithm_place = InTheSimulation(execution_condition=Equality(apply_SFD, 1))
        filter_class.computation_name = 'SFD application'
        filter_class.order = 10000000000000
        # Create the filtered equation for all of the conservative variables
        cons_vars, f_old, f_new = self.cons_arrays, self.f_old, self.f_new
        # Copy the value from the previous time-step
        equations = [OpenSBLIEq(old, new) for (old, new) in zip(f_old, f_new)]
        # Filter coefficients as stored constants, so they can be modified in the C code
        chi, omega, dt = symbols('chi_filt omega_filt dt', **{'cls': ConstantObject})
        chi.value, omega.value = chifilt, omegafilt
        CTD.add_constant([chi, omega])
        # Apply the filter
        equations += [OpenSBLIEq(filt, (v*(1.0-exp(-(chi+omega)*dt))+filt*(chi/omega+exp(-(chi+omega)*dt)))/(chi/omega+1.0)) for (filt, v) in zip(f_new, cons_vars)]
        # Second part
        equations += [OpenSBLIEq(v, (v*(exp(-(chi+omega)*dt)*chi/omega+1)+old_filt*chi/omega*(1.0-exp(-(chi+omega)*dt)))/(chi/omega+1.0)) for (old_filt, v) in zip(f_old, cons_vars)]
        # Apply periodic reset feature or not?
        if self.formulation != 'standard':
            equations = self.SFD_reset(equations)
        filter_class.add_equations(equations)
        self.equation_classes.append(filter_class)
        return

    def reset_q(self):
        """ Set the SFD filter array equal to the current flow field (q-vector)."""
        equations = []
        for i, term in enumerate(self.f_old):
            equations += [OpenSBLIEq(term, self.cons_arrays[i])]
        return equations

    def reset_f(self):
        """ Set the SFD filter array equal to the filter state at the previous periodic reset."""
        equations = []
        # Save the current filter state
        equations += [OpenSBLIEq(GridVariable('%s_reset_save' % str(x)), x) for x in self.f_new] # rhoSFD_new_reset_save
        # Reset the filter to its state at the previous periodic reset
        equations += [OpenSBLIEq(x, y) for (x,y) in zip(self.f_new, self.f_old)]
        # Set the long term storage to be used at next reset periodic reset
        equations += [OpenSBLIEq(x, GridVariable('%s_reset_save' % str(y))) for (x,y) in zip(self.f_old, self.f_new)] # rhoSFD_new_reset_save
        return equations

    def SFD_reset(self, input_equations):
        # Reset frequency parameter
        reset_control = ConstantObject('reset_SFD_frequency', integer=True)
        reset_control.value = 1
        current_iter = Globalvariable("iter", integer=True)
        CTD.add_constant(reset_control)
        check = Equality(Mod(current_iter+1, reset_control), 0)
        # Which reset procedure to apply?
        if self.formulation == 'reset_f':
            reset_equations = self.reset_f()
        elif self.formulation == 'reset_q':
            reset_equations = self.reset_q()
        # Create the conditional expression based on the reset frequency
        cond1 = ExprCondPair(reset_equations, check)
        cond2 = ExprCondPair(OpenSBLIEq(gv('temp'), 0.0), True)
        input_equations += [GroupedPiecewise(cond1, cond2)]
        return input_equations
