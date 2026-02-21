#!/usr/bin/env python
from opensbli import *
import copy
import os
import re
from opensbli.utilities.helperfunctions import substitute_simulation_parameters

# Dual-velocity 1D shock structure case inspired by:
#   ms-dg-dual-continuum/code/test_model.py + test_shock1d.py
ndim = 1
coordinate_symbol = "x"
eq = EinsteinEquation()
def _float_env(name, default):
    val = os.getenv(name)
    return default if val is None else float(val)

def _int_env(name, default):
    val = os.getenv(name)
    return default if val is None else int(val)

def _bool_env(name, default):
    val = os.getenv(name)
    if val is None:
        return default
    v = val.strip().lower()
    if v in ("1", "true", "yes", "on"):
        return True
    if v in ("0", "false", "no", "off"):
        return False
    raise ValueError(f"Invalid boolean for {name}: {val}")

def _tag(v):
    return f"{v:g}".replace(".", "p")


def _patch_tvd_defaults_cpp(delta_value, eps_value, simulation_name="opensbli"):
    """Override TVD defaults that are hardcoded in generated C++."""
    cpp_path = f"./{simulation_name}.cpp"
    with open(cpp_path) as f:
        txt = f.read()
    txt, _ = re.subn(r"(delta_TVD\s*=\s*)([^;]+);", rf"\g<1>{delta_value};", txt, count=1)
    txt, _ = re.subn(r"(eps_TVD\s*=\s*)([^;]+);", rf"\g<1>{eps_value};", txt, count=1)
    with open(cpp_path, "w") as f:
        f.write(txt)

MA_VALUE = _float_env("MA_VALUE", 1.55)
ML_VALUE = _float_env("ML_VALUE", 0.4)
TVD_KAPPA_VALUE = _float_env("TVD_KAPPA_VALUE", 1.0)
TVD_DELTA_VALUE = _float_env("TVD_DELTA_VALUE", 0.5)
TVD_EPS_VALUE = _float_env("TVD_EPS_VALUE", 1.0e-8)
DT_VALUE = _float_env("DT_VALUE", 2.0e-4)
NITER_VALUE = _int_env("NITER_VALUE", 500000)
NP0_VALUE = _int_env("NP0_VALUE", 201)
M1_VALUE = MA_VALUE
DELTA_IC_VALUE = 1.0
FILTER_METHOD = os.getenv("SHOCK_FILTER", "").strip().lower()
if FILTER_METHOD:
    if FILTER_METHOD == "tvd":
        USE_TVD_FILTER = True
        USE_WENO_FILTER = False
    elif FILTER_METHOD == "weno":
        USE_TVD_FILTER = False
        USE_WENO_FILTER = True
    elif FILTER_METHOD in ("none", "off"):
        USE_TVD_FILTER = False
        USE_WENO_FILTER = False
    else:
        raise ValueError(f"Invalid SHOCK_FILTER='{FILTER_METHOD}'. Use tvd|weno|none.")
else:
    # Backward-compatible env controls.
    USE_WENO_FILTER = _bool_env("USE_WENO_FILTER", False)
    USE_TVD_FILTER = _bool_env("USE_TVD_FILTER", True)

CASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = os.path.join(CASE_DIR, "outputs")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Constants used by the dual-velocity model.
constants = ["gama", "Pr", "Ml", "mu0", "s"]

# Governing equations.
DV = DualVelocitySplit(
    ndim,
    constants,
    coordinate_symbol=coordinate_symbol,
    conservative=True,
    viscosity="dynamic",
    include_skew_work=False,
    use_reynolds=False,
)

simulation_eq = SimulationEquations()
simulation_eq.add_equations(DV.mass)
simulation_eq.add_equations(DV.momentum)
simulation_eq.add_equations(DV.energy)

# Constituent relations.
constituent = ConstituentRelations()
constituent.add_equations(eq.expand("Eq(vl_i, rhou_i/rho)", ndim, coordinate_symbol, [], constants))
constituent.add_equations(
    eq.expand(
        "Eq(p, (gama-1)*(rhoE - (1/2)*rho*(KD(_i,_j)*vl_i*vl_j)))",
        ndim,
        coordinate_symbol,
        [],
        constants,
    )
)
constituent.add_equations(eq.expand("Eq(T, p/rho)", ndim, coordinate_symbol, [], constants))
constituent.add_equations(eq.expand("Eq(jl_i, Ml*Der(p, x_i)/p)", ndim, coordinate_symbol, [], constants))
constituent.add_equations(eq.expand("Eq(mu, mu0*T**s)", ndim, coordinate_symbol, [], constants))
constituent.add_equations(
    eq.expand("Eq(kappa, (gama/((gama-1)*Pr))*mu)", ndim, coordinate_symbol, [], constants)
)
constituent.add_equations(eq.expand("Eq(a, (gama*p/rho)**0.5)", ndim, coordinate_symbol, [], constants))

block = SimulationBlock(ndim, block_number=0)

# Initial condition: smoothed RH jump over x in [-10, 10].
initial = GridBasedInitialisation()
x0 = "Eq(DataObject(x0), -10.0 + block.deltas[0]*block.grid_indexes[0])"
rhoL = "Eq(GridVariable(rhoL), gama)"
pL = "Eq(GridVariable(pL), 1.0)"
uL = "Eq(GridVariable(uL), %.16f)" % M1_VALUE

pR = "Eq(GridVariable(pR), pL*(1.0 + 2.0*gama/(gama+1.0)*((%.16f)**2.0 - 1.0)))" % M1_VALUE
rhoR = "Eq(GridVariable(rhoR), rhoL*((gama+1.0)*(%.16f)**2.0)/((gama-1.0)*(%.16f)**2.0 + 2.0))" % (M1_VALUE, M1_VALUE)
uR = "Eq(GridVariable(uR), %.16f*rhoL/rhoR)" % M1_VALUE

tau0 = "Eq(GridVariable(tau0), 0.5*(1.0 + tanh(DataObject(x0)/%.16f)))" % DELTA_IC_VALUE
d = "Eq(GridVariable(d), rhoL + (rhoR-rhoL)*tau0)"
p = "Eq(GridVariable(p), pL + (pR-pL)*tau0)"
u0 = "Eq(GridVariable(vl0), uL + (uR-uL)*tau0)"

rho = "Eq(DataObject(rho), d)"
rhou0 = "Eq(DataObject(rhou0), d*vl0)"
rhoE = "Eq(DataObject(rhoE), p/(gama-1.0) + 0.5*d*(vl0**2.0))"

eqns = [x0, rhoL, pL, uL, pR, rhoR, uR, tau0, d, p, u0, rho, rhou0, rhoE]
local_dict = {"block": block, "GridVariable": GridVariable, "DataObject": DataObject}
initial_equations = [parse_expr(e, local_dict=local_dict) for e in eqns]
initial.add_equations(initial_equations)

# Dirichlet far-field states from RH relations.
dl, ul, pl = symbols("dl ul pl", **{"cls": GridVariable})
left_eval = [
    OpenSBLIEq(dl, Symbol("gama")),
    OpenSBLIEq(ul, M1_VALUE),
    OpenSBLIEq(pl, 1.0),
    OpenSBLIEq(DataObject("rho"), dl),
    OpenSBLIEq(DataObject("rhou0"), dl * ul),
    OpenSBLIEq(DataObject("rhoE"), pl / (Symbol("gama") - 1.0) + 0.5 * dl * ul**2),
]

dr, ur, pr = symbols("dr ur pr", **{"cls": GridVariable})
right_eval = [
    OpenSBLIEq(dr, Symbol("gama") * ((Symbol("gama") + 1.0) * M1_VALUE ** 2.0) / (((Symbol("gama") - 1.0) * M1_VALUE ** 2.0) + 2.0)),
    OpenSBLIEq(ur, M1_VALUE * Symbol("gama") / dr),
    OpenSBLIEq(pr, 1.0 * (1.0 + 2.0 * Symbol("gama") / (Symbol("gama") + 1.0) * (M1_VALUE ** 2.0 - 1.0))),
    OpenSBLIEq(DataObject("rho"), dr),
    OpenSBLIEq(DataObject("rhou0"), dr * ur),
    OpenSBLIEq(DataObject("rhoE"), pr / (Symbol("gama") - 1.0) + 0.5 * dr * ur**2),
]

boundaries = [DirichletBC(0, 0, left_eval), DirichletBC(0, 1, right_eval)]

schemes = {}
cent = StoreSome(4, "vl0")
schemes[cent.name] = cent
rk = RungeKuttaLS(3, formulation="SSP")
schemes[rk.name] = rk

block.set_block_boundaries(boundaries)
output_name = os.path.join(OUTPUT_DIR, f"opensbli_output_ma{_tag(MA_VALUE)}_ml{_tag(ML_VALUE)}.h5")
kwargs = {"iotype": "Write", "name": output_name}
h5 = iohdf5(**kwargs)
h5.add_arrays(simulation_eq.time_advance_arrays)
h5.add_arrays([DataObject("x0"), DataObject("p"), DataObject("vl0"), DataObject("T"), DataObject("jl0")])
block.setio(copy.deepcopy(h5))

if USE_WENO_FILTER and USE_TVD_FILTER:
    raise ValueError("Enable only one shock filter at a time: WENO or TVD.")

if USE_TVD_FILTER:
    TVD_filter = TVDFilter(block, airfoil=False)
    block.set_equations(TVD_filter.equation_classes)
elif USE_WENO_FILTER:
    WF = WENOFilter(
        block,
        order=5,
        formulation="Z",
        flux_type="LLF",
        airfoil=False,
        optimize=False,
    )
    block.set_equations(WF.equation_classes)

if USE_TVD_FILTER:
    filter_name = "TVD"
elif USE_WENO_FILTER:
    filter_name = "WENO"
else:
    filter_name = "NONE"
print(f"Shock filter selection: {filter_name}")

block.set_equations([copy.deepcopy(constituent), copy.deepcopy(simulation_eq), initial])
block.set_discretisation_schemes(schemes)
block.discretise()

alg = TraditionalAlgorithmRK(block)
SimulationDataType.set_datatype(Double)
OPSC(alg)

constants = [
    "gama",
    "Pr",
    "Ml",
    "mu0",
    "s",
    "dt",
    "niter",
    "block0np0",
    "Delta0block0",
    "eps",
    "inv_rfact0_block0",
    "delta_TVD",
    "eps_TVD",
    "kappa_TVD",
]
values = [
    "1.6666666666666667",
    "0.6666666666666666",
    f"{ML_VALUE}",
    "1.0",
    "0.75",
    f"{DT_VALUE}",
    f"{NITER_VALUE}",
    f"{NP0_VALUE}",
    "20.0/(block0np0-1)",
    "1.0e-16",
    "1.0/Delta0block0",
    f"{TVD_DELTA_VALUE}",
    f"{TVD_EPS_VALUE}",
    f"{TVD_KAPPA_VALUE}",
]
substitute_simulation_parameters(constants, values)
_patch_tvd_defaults_cpp(TVD_DELTA_VALUE, TVD_EPS_VALUE)
print(f"TVD parameters: kappa_TVD={TVD_KAPPA_VALUE}, delta_TVD={TVD_DELTA_VALUE}, eps_TVD={TVD_EPS_VALUE}")
print_iteration_ops(NaN_check="rho")
