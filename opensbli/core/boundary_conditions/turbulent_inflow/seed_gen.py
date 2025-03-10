""" Seed generation for turbulent inflow RNG
    - Seeds are generated randomly for each fresh code generation
    - For replicatibility, users may wish to set constant values - this is done by passing fixed_seed=True
    - Generates 6 seed arrays with size ny*nz
    - Seed values are uniformly distributed between 0 and 1E7 """

from sympy import pprint
import numpy as np


def get_seeds(ny, nz, niter, fixed_seed=False):
    if fixed_seed:
        np.random.seed(seed=123456789)
    seeds = 1.E07*np.random.rand(ny,nz,6)
    return seeds

