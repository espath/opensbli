"""@brief This contains the Physics Object for Navier-Stokes equations.
   @authors Luis Espath
   @contributors
   @details
    Dimentional dual-velocity:
      ∂t ρ + Div(ρ vl) = Div(ρ jl)
      ∂t(ρ vl) + Div(ρ vl⊗vl + p I) = Div(S + 1/2 ρ (vl⊗jl + jl⊗vl)) + b
      ∂t(ρ e) + Div((ρ e + p) vl)
        = Div(S^T vl - q + 1/2 (ρ vl ∧ jl) vl + (ρ e + p) jl)
        + 1/2 (ρ vl ∧ jl) : Wl + r + b·vl

      q = -kappa(ρ,θ) ∇θ,  p = R ρ θ
      jl = iota(ρ,θ) ∇ ln p,
      S = 2 μ(ρ,θ) D + ζ(ρ,θ) (tr D) I,  ζ = -2/3 μ
      Wl = 1/2 (∇vl - (∇vl)^T)
      D  = 1/2 (∇vl + (∇vl)^T)

    Dimensionless dual-velocity:
      ∂t ρ + Div(ρ vl) = Div(ρ jl)
      ∂t(ρ vl) + Div(ρ vl⊗vl + p I) = Div(S + 1/2 ρ (vl⊗jl + jl⊗vl))
      ∂t(ρ e) + Div((ρ e + p) vl)
        = Div(S^T vl - q + 1/2 (ρ vl ∧ jl) vl + (ρ e + p) jl)
        + 1/2 (ρ vl ∧ jl) : Wl

      θ = p / ρ
      μ = μ0 θ^s
      kappa(θ) = γ/((γ-1) Pr) * μ(θ)
      q = -kappa(θ) ∇θ
      jl = Ml ∇ ln p
      S = 2 μ(θ) D + ζ(θ) (tr D) I,  ζ = -2/3 μ
"""
from opensbli.core.grid import GridVariable
from opensbli.core.opensbliobjects import DataObject, ConstantObject, DataSetBase
from sympy import S, Rational, pprint, Eq
from opensbli.utilities.helperfunctions import dot
from opensbli.core.parsing import EinsteinEquation
import itertools


class Physics(object):
    """A base class defining physical model's that are to be used for
    codegeneration"""


class PhysicsVariable(object):
    blocknumber = None
    shape = None

    def __init__(self, name):
        self._variable = self.block.location_dataset(name)
        self._relation = None
        self._conservative_relation = None
        return

    @property
    def relation(self):
        """Conver the dataobjects to datasets"""
        return self._relation

    @relation.setter
    def relation(self, formula):
        self._relation = formula
        return

    @property
    def conservative_relation(self):
        return self._conservative_relation

    @conservative_relation.setter
    def conservative_relation(self, formula):
        self._conservative_relation = formula

    @property
    def variable(self):
        return self._variable

    @property
    def datasetbase(self):
        return DataSetBase(self._variable, self.shape, self.blocknumber)


class DualVelocityPhysics(Physics):
    eqobject = EinsteinEquation()

    def __init__(self, block, **settings):
        self.ndim = block.ndim
        PhysicsVariable.blocknumber = block.blocknumber
        PhysicsVariable.shape = block.shape
        PhysicsVariable.block = block
        self.set_names(**settings)
        return

    def set_names(self, **settings):
        self._density = PhysicsVariable('rho')
        self._density.relation = self._density.variable
        m = self.eqobject.expand('rhou_i', self.ndim, "x", substitutions=[], constants=[])
        self._momentum = [PhysicsVariable(m1) for m1 in m]  # Expand rhou_i to number of dimensions
        self._totalenergy = PhysicsVariable('rhoE')
        self._pressure = PhysicsVariable('p')
        self._temperature = PhysicsVariable('T')
        self._speed_of_sound = PhysicsVariable('a')
        m = self.eqobject.expand('vl_i', self.ndim, "x", substitutions=[], constants=[])
        self._velocity = [PhysicsVariable(m1) for m1 in m]  # Expand vl_i to number of dimensions
        m = self.eqobject.expand('jl_i', self.ndim, "x", substitutions=[], constants=[])
        self._diffusive_mass_flux = [PhysicsVariable(m1) for m1 in m]  # Expand jl_i to number of dimensions
        m = self.eqobject.expand('gradlnp_i', self.ndim, "x", substitutions=[], constants=[])
        self._log_pressure_gradient = [PhysicsVariable(m1) for m1 in m]  # Expand gradlnp_i to number of dimensions
        self._viscosity = PhysicsVariable('mu')
        self._thermal_conductivity = PhysicsVariable('kappa')

        # Set the relations for the lagrangian velocity
        for d in range(self.ndim):
            self._velocity[d].relation = self.momentum()[d]/self.density()

        # jl = Ml * gradlnp
        for d in range(self.ndim):
            self._diffusive_mass_flux[d].relation = self.mass_mobility()*self.log_pressure_gradient()[d]

        # Set the relations for the momentum
        for d in range(self.ndim):
            self._momentum[d].relation = self.velocity()[d]*self.density()

        self._Reynoldsnumber = ConstantObject("Re")
        self._CpbyCv = ConstantObject("gama")  # Ratio of specific heats
        self._Prandtlnumber = ConstantObject("Pr")
        self._Mach = ConstantObject("Minf")
        self._MassMobility = ConstantObject("Ml")
        self._ViscosityRef = ConstantObject("mu0")
        self._ViscosityExponent = ConstantObject("s")

        self._speed_of_sound.relation = (self.specific_heat_ratio()*self.pressure()/self.density())

        self._pressure.relation = (self.specific_heat_ratio() - S.One)*(self.total_energy() - Rational(1, 2)*self.density()*(dot(self.velocity(), self.velocity())))
        self._pressure.conservative_relation = (self.specific_heat_ratio() - S.One)*(self.total_energy() - Rational(1, 2)*(dot(self.momentum(), self.momentum()))/self.density())

        # Dimensionless dual-velocity thermodynamic closure: theta = p / rho
        self._temperature.relation = self.pressure(relation=True)/self.density()
        self._temperature.conservative_relation = self.pressure(relation=True, conservative=True)/self.density()

        # Transport model from the dual-velocity header:
        #   mu = mu0 * theta^s
        #   kappa = gamma/((gamma-1)Pr) * mu
        self._viscosity.relation = self.viscosity_reference() * self.temperature(relation=True)**self.viscosity_exponent()
        self._thermal_conductivity.relation = (
            self.specific_heat_ratio()/((self.specific_heat_ratio() - S.One)*self.prandtl_number())
        ) * self.viscosity(relation=True)
        return

    def specific_heat_ratio(self):
        return self._CpbyCv

    def mach_number(self):
        return self._Mach

    def prandtl_number(self):
        return self._Prandtlnumber

    def mass_mobility(self):
        return self._MassMobility

    def viscosity_reference(self):
        return self._ViscosityRef

    def viscosity_exponent(self):
        return self._ViscosityExponent

    def density(self, relation=False):
        if not relation:
            return self._density.variable
        else:
            return self._density.relation

    def momentum(self, relation=False):
        if not relation:
            return [m.variable for m in self._momentum]
        else:
            return [m.relation for m in self._momentum]

    def total_energy(self, relation=False):
        if not relation:
            return self._totalenergy.variable
        else:
            raise NotImplementedError("")

    def pressure(self, relation=False, conservative=False):
        if not relation:
            return self._pressure.variable
        else:
            if conservative:
                return self._pressure.conservative_relation
            else:
                return self._pressure.relation

    def speed_of_sound(self, relation=False):
        if not relation:
            return self._speed_of_sound.variable
        else:
            return self._speed_of_sound.relation

    def velocity(self, relation=False):
        if not relation:
            return [m.variable for m in self._velocity]
        else:
            return [m.relation for m in self._velocity]

    def lagrangian_velocity(self, relation=False):
        return self.velocity(relation=relation)

    def diffusive_mass_flux(self, relation=False):
        if not relation:
            return [m.variable for m in self._diffusive_mass_flux]
        else:
            return [m.relation for m in self._diffusive_mass_flux]

    def log_pressure_gradient(self, relation=False):
        if not relation:
            return [m.variable for m in self._log_pressure_gradient]
        else:
            return [m.relation for m in self._log_pressure_gradient]

    def temperature(self, relation=False, conservative=False):
        if not relation:
            return self._temperature.variable
        else:
            if conservative:
                return self._temperature.conservative_relation
            else:
                return self._temperature.relation

    def viscosity(self, relation=False):
        if not relation:
            return self._viscosity.variable
        else:
            return self._viscosity.relation

    def thermal_conductivity(self, relation=False):
        if not relation:
            return self._thermal_conductivity.variable
        else:
            return self._thermal_conductivity.relation


# Backward-compatible alias for callers expecting NSphysics-style naming.
NSphysics = DualVelocityPhysics


class DualVelocitySplit(object):
    """Dual-velocity PDE builder.

    This mirrors NS_Split usage style and returns expanded mass, momentum,
    and energy equations including jl terms from the dual-velocity theory.
    """

    def __init__(
        self,
        ndim,
        constants,
        coordinate_symbol="x",
        conservative=True,
        viscosity="dynamic",
        include_skew_work=False,
        use_reynolds=True,
    ):
        self.ndim = ndim
        self.constants = constants
        self.coordinate_symbol = coordinate_symbol
        self.conservative = conservative
        self.viscosity = viscosity
        self.include_skew_work = include_skew_work
        self.use_reynolds = use_reynolds
        self.EE = EinsteinEquation()

        if self.conservative:
            self.mom_lhs = "rhou"
            self.energy_lhs = "rhoE"
        else:
            self.mom_lhs = "rho*vl"
            self.energy_lhs = "Et"

        self.substitutions = self.diffusive_terms()
        self.mass = self.mass_eq()
        self.momentum = self.momentum_eq()
        self.energy = self.energy_eq()

    def mass_eq(self):
        eq = "Eq(Der(rho, t), -Conservative(rho*vl_j, x_j) + Conservative(rho*jl_j, x_j))"
        return self.EE.expand(eq, self.ndim, self.coordinate_symbol, self.substitutions, self.constants)

    def momentum_eq(self):
        eq = (
            "Eq(Der(%s_i, t), "
            "-Conservative(%s_i*vl_j + KD(_i,_j)*p, x_j) "
            "+Conservative((1/2)*rho*(vl_i*jl_j + jl_i*vl_j), x_j) "
            "+Der(tau_i_j, x_j))"
        ) % (self.mom_lhs, self.mom_lhs)
        return self.EE.expand(eq, self.ndim, self.coordinate_symbol, self.substitutions, self.constants)

    def energy_eq(self):
        base = (
            "Eq(Der(%s, t), "
            "-Conservative((%s + p)*vl_j, x_j) "
            "+Conservative((%s + p)*jl_j, x_j) "
            "+Der(vl_i*tau_i_j, x_j) + Der(q_j, x_j))"
        ) % (self.energy_lhs, self.energy_lhs, self.energy_lhs)

        # Optional source term for (1/2) (rho vl ^ jl) : Wl when a model-specific
        # scalar source is provided externally as skew_work.
        if self.include_skew_work:
            base = base[:-1] + " + skew_work)"

        return self.EE.expand(base, self.ndim, self.coordinate_symbol, self.substitutions, self.constants)

    def diffusive_terms(self):
        # Avoid Der(log()) here: OpenSBLI sanitisation expects special
        # derivative function objects and can fail on raw SymPy log atoms.
        substitutions = ["Eq(jl_i, Ml*Der(p, x_i)/p)"]
        if self.viscosity == "inviscid":
            return substitutions

        if self.viscosity == "constant":
            if self.use_reynolds:
                stress = "Eq(tau_i_j, (1.0/Re)*(Der(vl_i,x_j)+Der(vl_j,x_i)-(2/3)*KD(_i,_j)*Der(vl_k,x_k)))"
                heat = "Eq(q_j, (1.0/Re)*Der(T,x_j))"
            else:
                stress = "Eq(tau_i_j, (Der(vl_i,x_j)+Der(vl_j,x_i)-(2/3)*KD(_i,_j)*Der(vl_k,x_k)))"
                heat = "Eq(q_j, Der(T,x_j))"
        else:
            # Dynamic dual-velocity transport closures.
            mu_eq = "Eq(mu, mu0*T**s)"
            kappa_eq = "Eq(kappa, (gama/((gama-1)*Pr))*mu)"
            if self.use_reynolds:
                stress = "Eq(tau_i_j, (mu/Re)*(Der(vl_i,x_j)+Der(vl_j,x_i)-(2/3)*KD(_i,_j)*Der(vl_k,x_k)))"
                heat = "Eq(q_j, (kappa/Re)*Der(T,x_j))"
            else:
                stress = "Eq(tau_i_j, mu*(Der(vl_i,x_j)+Der(vl_j,x_i)-(2/3)*KD(_i,_j)*Der(vl_k,x_k)))"
                heat = "Eq(q_j, kappa*Der(T,x_j))"
            substitutions += [mu_eq, kappa_eq]

        substitutions += [stress, heat]
        return substitutions


# class StatsVariable(object):
#     def __init__(self, var):
#         if not isinstance(var, DataObject):
#             raise ValueError("")
#         self._variable = var
#         self._relation = None
#         self._conservative_relation = None
#         return

#     @property
#     def relation(self):
#         """Conver the dataobjects to datasets"""
#         return self._relation

#     @relation.setter
#     def relation(self, formula):
#         self._relation = formula
#         return

#     @property
#     def conservative_relation(self):
#         return self._conservative_relation

#     @conservative_relation.setter
#     def conservative_relation(self, formula):
#         self._conservative_relation = formula

#     @property
#     def variable(self):
#         # returns the dataset of the current variable at the location
#         dsetbase = self.datasetbase
#         return dsetbase[dsetbase.location()]

#     @property
#     def datasetbase(self):
#         return DataSetBase(self._variable)


# class NSPhysics_Stats(NSphysics, StatsVariable):
#     def __init__(self, ndim, **settings):
#         NSphysics.__init__(self, ndim)
#         self.ndim = ndim
#         # Create the components required for stats
#         self.create_stat_components(**settings)
#         self.niter = ConstantObject('niter', integer=True)
#         from sympy import Int
#         self.niter.datatype = Int()
#         # Equations to initialise stat arrays to zeros
#         self.init = self.init_stats()
#         # During time loop stat collection equations
#         self.collection = self.stat_collection()
#         # Reynolds or Favre averaging after the simulation
#         self.average = self.favre_average()
#         # self.average = self.reynolds_average()
#         return

#     def create_stat_components(self, **settings):
#         # rho mean equation
#         self._rho_mean = StatsVariable(DataObject('rhomean'))
#         self._rho_mean.relation = self._rho_mean.variable + self._density.relation
#         # momentum component means
#         self._momentum_means = [StatsVariable(DataObject('rhou%dmean' % i)) for i in range(self.ndim)]
#         for no, component in enumerate(self._momentum_means):
#             component.relation = component.variable + self._momentum[no].variable
#         # Create mean Reynolds stresses
#         indices = [i for i in range(self.ndim)]
#         stress_components = sorted(set(tuple(sorted(t)) for t in set(itertools.product(indices, indices))), key=lambda element: (element[0], element[1]))
#         self._reynolds_stress_means = [StatsVariable(DataObject('rhou%d%dmean' % (i, j))) for (i, j) in stress_components]
#         for no, indices in enumerate(stress_components):
#             self._reynolds_stress_means[no].relation = self._reynolds_stress_means[no].variable + \
#                                                         (self._momentum[indices[0]].variable*self._momentum[indices[1]].variable)/self._density.variable
#         return

#     def init_stats(self):
#         equations = [Eq(self._rho_mean.variable, 0.0)]
#         equations += [Eq(mean.variable, 0.0) for mean in self._momentum_means]
#         equations += [Eq(mean.variable, 0.0) for mean in self._reynolds_stress_means]
#         for eqn in equations:
#             pprint(eqn)
#         return equations

#     def stat_collection(self):
#         equations = [Eq(self._rho_mean.variable, self._rho_mean.relation)]
#         equations += [Eq(mean.variable, mean.relation) for mean in self._momentum_means]
#         equations += [Eq(mean.variable, mean.relation) for mean in self._reynolds_stress_means]
#         for eqn in equations:
#             pprint(eqn)
#         return

#     def reynolds_average(self):
#         equations = [Eq(self._rho_mean.variable, self._rho_mean.variable/self.niter)]
#         equations += [Eq(mean.variable, mean.variable/(self.niter)) for mean in self._momentum_means]
#         equations += [Eq(mean.variable, mean.variable/(self.niter)) for mean in self._reynolds_stress_means]

#         return

#     def favre_average(self):
#         equations = [Eq(self._rho_mean.variable, self._rho_mean.variable/self.niter)]
#         equations += [Eq(GridVariable('rmean'), self._rho_mean.variable)]
#         equations += [Eq(mean.variable, mean.variable/(self.niter*GridVariable('rmean'))) for mean in self._momentum_means]
#         equations += [Eq(mean.variable, mean.variable/(self.niter*GridVariable('rmean'))) for mean in self._reynolds_stress_means]
#         for eqn in equations:
#             pprint(eqn)
#         return equations
