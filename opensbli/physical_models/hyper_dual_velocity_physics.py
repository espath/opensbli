"""@brief This contains the Physics Object for the dimensionless
   dual-velocity compressible Navier-Stokes equations.

   @authors Luis Espath
   @contributors
   @details

   Dimensionless dual-velocity:

     ∂t ρ + Div(ρ vl) − Div(ρ jl) = 0,

     ∂t(ρ vl)
       + Div(ρ vl⊗vl + p I)
       − Div( 1/Re S + 1/2 ρ (vl⊗jl + jl⊗vl) )
       − b = 0,

     ∂t(ρ e)
       + Div((ρ e + p) vl)
       − Div( 1/Re ( S^T vl − 1/((γ−1)Ma^2 Pr) q )
              + 1/2 (ρ vl ∧ jl) vl
              + (ρ e + p) jl )
       − 1/2 (ρ vl ∧ jl) : Wl
       − r
       − b·vl = 0.

   Constitutive relations:

     μ(θ) = θ^(3/2) * (1 + SuthT/RefT) / (θ + SuthT/RefT)
     q = −kappa(ρ, θ) ∇θ,
     p = 1/(γ Ma^2) ρ θ,
     jl = (Ml/(γ Ma^2)) iota(ρ, θ) ∇p,
     S = 2 μ(ρ, θ) Dl + ζ(ρ, θ) (tr D) I,
     ζ(ρ, θ) = −2/3 μ(ρ, θ),
     Wl = 1/2 (∇vl − (∇vl)^T),
     Dl  = 1/2 (∇vl + (∇vl)^T).
"""
from opensbli.core.grid import GridVariable
from opensbli.core.opensbliobjects import DataObject, ConstantObject, DataSetBase
from sympy import S, Rational
from opensbli.utilities.helperfunctions import dot
from opensbli.core.parsing import EinsteinEquation
from opensbli.physical_models.split_forms import NS_Split
from opensbli.equation_types.opensbliequations import OpenSBLIEq


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


class HyperDualVelocityPhysics(Physics):
    eqobject = EinsteinEquation()

    def __init__(self, block, **settings):
        self.ndim = block.ndim
        PhysicsVariable.blocknumber = block.blocknumber
        PhysicsVariable.shape = block.shape
        PhysicsVariable.block = block
        self.set_names(**settings)
        return

    def set_names(self, **settings):
        self._density = PhysicsVariable("rho")
        self._density.relation = self._density.variable
        m = self.eqobject.expand("rhou_i", self.ndim, "x", substitutions=[], constants=[])
        self._momentum = [PhysicsVariable(m1) for m1 in m]
        self._totalenergy = PhysicsVariable("rhoE")
        self._pressure = PhysicsVariable("p")
        self._temperature = PhysicsVariable("T")
        self._speed_of_sound = PhysicsVariable("a")
        m = self.eqobject.expand("vl_i", self.ndim, "x", substitutions=[], constants=[])
        self._velocity = [PhysicsVariable(m1) for m1 in m]
        m = self.eqobject.expand("jl_i", self.ndim, "x", substitutions=[], constants=[])
        self._diffusive_mass_flux = [PhysicsVariable(m1) for m1 in m]
        m = self.eqobject.expand("gradp_i", self.ndim, "x", substitutions=[], constants=[])
        self._pressure_gradient = [PhysicsVariable(m1) for m1 in m]
        self._iota = PhysicsVariable("iota")
        self._viscosity = PhysicsVariable("mu")
        self._thermal_conductivity = PhysicsVariable("kappa")

        for d in range(self.ndim):
            self._velocity[d].relation = self.momentum()[d] / self.density()
            self._momentum[d].relation = self.velocity()[d] * self.density()

        self._CpbyCv = ConstantObject("gama")
        self._Prandtlnumber = ConstantObject("Pr")
        self._Mach = ConstantObject("Minf")
        self._MassMobility = ConstantObject("Ml")
        self._ViscosityRef = ConstantObject("mu0")
        self._ViscosityExponent = ConstantObject("s")
        self._SutherlandTemperature = ConstantObject("SuthT")
        self._ReferenceTemperature = ConstantObject("RefT")

        self._speed_of_sound.relation = self.specific_heat_ratio() * self.pressure() / self.density()

        self._pressure.relation = (self.specific_heat_ratio() - S.One) * (
            self.total_energy() - Rational(1, 2) * self.density() * (dot(self.velocity(), self.velocity()))
        )
        self._pressure.conservative_relation = (self.specific_heat_ratio() - S.One) * (
            self.total_energy() - Rational(1, 2) * (dot(self.momentum(), self.momentum())) / self.density()
        )

        # p = rho*T/(gama*Minf^2) -> T = gama*Minf^2*(p/rho)
        self._temperature.relation = (
            self.specific_heat_ratio() * self.mach_number() ** 2
        ) * self.pressure(relation=True) / self.density()
        self._temperature.conservative_relation = (
            self.specific_heat_ratio() * self.mach_number() ** 2
        ) * self.pressure(relation=True, conservative=True) / self.density()

        self._iota.relation = S.One
        for d in range(self.ndim):
            self._diffusive_mass_flux[d].relation = (
                self.mass_mobility()
                / (self.specific_heat_ratio() * self.mach_number() ** 2)
                * self.iota(relation=True)
                * self.pressure_gradient()[d]
            )

        # Sutherland viscosity law in dimensionless temperature theta (= T variable here).
        self._viscosity.relation = (
            self.temperature(relation=True) ** Rational(3, 2)
            * (S.One + self.sutherland_temperature() / self.reference_temperature())
            / (self.temperature(relation=True) + self.sutherland_temperature() / self.reference_temperature())
        )
        # Kept as a standalone variable to preserve q = -kappa grad(theta) form in substitutions.
        self._thermal_conductivity.relation = self.viscosity(relation=True)
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

    def sutherland_temperature(self):
        return self._SutherlandTemperature

    def reference_temperature(self):
        return self._ReferenceTemperature

    def density(self, relation=False):
        return self._density.relation if relation else self._density.variable

    def momentum(self, relation=False):
        return [m.relation for m in self._momentum] if relation else [m.variable for m in self._momentum]

    def total_energy(self, relation=False):
        if relation:
            raise NotImplementedError("")
        return self._totalenergy.variable

    def pressure(self, relation=False, conservative=False):
        if not relation:
            return self._pressure.variable
        if conservative:
            return self._pressure.conservative_relation
        return self._pressure.relation

    def speed_of_sound(self, relation=False):
        return self._speed_of_sound.relation if relation else self._speed_of_sound.variable

    def velocity(self, relation=False):
        return [m.relation for m in self._velocity] if relation else [m.variable for m in self._velocity]

    def lagrangian_velocity(self, relation=False):
        return self.velocity(relation=relation)

    def pressure_gradient(self, relation=False):
        return [m.relation for m in self._pressure_gradient] if relation else [m.variable for m in self._pressure_gradient]

    def iota(self, relation=False):
        return self._iota.relation if relation else self._iota.variable

    def diffusive_mass_flux(self, relation=False):
        return [m.relation for m in self._diffusive_mass_flux] if relation else [m.variable for m in self._diffusive_mass_flux]

    def temperature(self, relation=False, conservative=False):
        if not relation:
            return self._temperature.variable
        if conservative:
            return self._temperature.conservative_relation
        return self._temperature.relation

    def viscosity(self, relation=False):
        return self._viscosity.relation if relation else self._viscosity.variable

    def thermal_conductivity(self, relation=False):
        return self._thermal_conductivity.relation if relation else self._thermal_conductivity.variable


class HyperDualVelocitySplit(object):
    """Dimensionless hyper-dual-velocity PDE builder.

    Uses NS_Split for the baseline convective/diffusive operator and adds
    dual-velocity correction fluxes. This guarantees Ml=0 is on the same
    discrete NS path as the selected split form.
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
        split_type="KGP",
        energy_formulation="enthalpy",
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
            self.mom_lhs = "u"
            self.energy_lhs = "Et"

        # Baseline NS operator (same discrete convective core as NS case).
        if self.use_reynolds:
            ns_viscosity = self.viscosity
        else:
            # NS_Split always uses Re-scaled diffusion. Keep compatibility by
            # routing no-Re mode through constant form; this path is currently
            # unused by production apps.
            ns_viscosity = "constant"
        ns = NS_Split(
            split_type,
            self.ndim,
            self.constants,
            coordinate_symbol=self.coordinate_symbol,
            conservative=self.conservative,
            viscosity=ns_viscosity,
            energy_formulation=energy_formulation,
            debug=False,
        )
        self.mass = ns.mass
        self.momentum = ns.momentum
        self.energy = ns.energy

        self.substitutions = self.dual_substitutions()
        self.add_dual_corrections()

    def dual_substitutions(self):
        return [
            "Eq(iota, 1.0)",
            "Eq(gradp_i, Der(p, x_i))",
            "Eq(jl_i, (Ml/(gama*Minf*Minf))*iota*gradp_i)",
        ]

    def add_dual_corrections(self):
        mass_corr = self.EE.expand(
            "Eq(Der(rho, t), Conservative(rho*jl_j, x_j))",
            self.ndim,
            self.coordinate_symbol,
            self.substitutions,
            self.constants,
        )
        mom_corr = self.EE.expand(
            "Eq(Der(%s_i, t), Conservative((1/2)*rho*(u_i*jl_j + jl_i*u_j), x_j))" % self.mom_lhs,
            self.ndim,
            self.coordinate_symbol,
            self.substitutions,
            self.constants,
        )
        ene_corr = self.EE.expand(
            "Eq(Der(%s, t), Conservative((%s + p)*jl_j, x_j))" % (self.energy_lhs, self.energy_lhs),
            self.ndim,
            self.coordinate_symbol,
            self.substitutions,
            self.constants,
        )

        mass_base = self._as_list(self.mass)
        mass_add = self._as_list(mass_corr)
        mom_base = self._as_list(self.momentum)
        mom_add = self._as_list(mom_corr)
        ene_base = self._as_list(self.energy)
        ene_add = self._as_list(ene_corr)

        mass_out = [OpenSBLIEq(m.lhs, m.rhs + c.rhs) for m, c in zip(mass_base, mass_add)]
        mom_out = [OpenSBLIEq(m.lhs, m.rhs + c.rhs) for m, c in zip(mom_base, mom_add)]
        ene_out = [OpenSBLIEq(m.lhs, m.rhs + c.rhs) for m, c in zip(ene_base, ene_add)]

        self.mass = self._unpack(mass_out)
        self.momentum = self._unpack(mom_out)
        self.energy = self._unpack(ene_out)

        if self.include_skew_work:
            ene_src = self.EE.expand(
                "Eq(Der(%s, t), skew_work)" % self.energy_lhs,
                self.ndim,
                self.coordinate_symbol,
                self.substitutions,
                self.constants,
            )
            ene_base = self._as_list(self.energy)
            ene_add = self._as_list(ene_src)
            ene_out = [OpenSBLIEq(e.lhs, e.rhs + s.rhs) for e, s in zip(ene_base, ene_add)]
            self.energy = self._unpack(ene_out)

    @staticmethod
    def _as_list(eq):
        return eq if isinstance(eq, list) else [eq]

    @staticmethod
    def _unpack(eq_list):
        return eq_list[0] if len(eq_list) == 1 else eq_list
