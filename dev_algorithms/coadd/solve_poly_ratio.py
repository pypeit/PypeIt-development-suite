

import numpy as np
import scipy
import matplotlib.pyplot as plt
import os

def solve_poly_fn(theta, xvector, polyfunc, nback = None):

    if nback is None:
        acoeff = theta
    else:
        acoeff = theta[0:-nback]
        bcoeff = theta[-nback:]

    ymult = utils.func_val(acoeff, xvector, polyfunc, minx=wave_min, maxx=wave_max)
    if nback is not None:
        ymult = utils.func_val(acoeff, xvector, polyfunc, minx=wave_min, maxx=wave_max)
