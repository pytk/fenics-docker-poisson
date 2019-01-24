from __future__ import print_function
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import scipy
import bias
import matplotlib.animation as animation
import mpmath

# Warning: from fenics import * will import both `sym` and
# `q` from FEniCS. We therefore import FEniCS first and then
# overwrite these objects.
from dolfin import *        

def subbandContribution(fermienergy, eigenvalue, cons, device):
    """
    to calculate number of electron in i-th subband dominated

    return 2d numpy array of electron number
    """
    hbar=1.05457266 * 10**-34
    constant_value = hbar*hbar*cons.PI
    effective_mass = device.material["electron_effective_mass"] * cons.M
    dos = effective_mass / constant_value
    state = dos*cons.KB*cons.T*math.log(math.exp((fermienergy + eigenvalue)/(cons.KB*cons.T)) + 1)
    return state

def electronDensity(eigenvector, state):
	"""
	print the electron density function(x)
	"""
	electron_density = []
	for i in range(len(eigenvector)):
		kth_term = [(j * j * state[i]) for j in eigenvector[i]]
		electron_density.append(kth_term)

	n = np.sum(electron_density, axis=0)
	return n