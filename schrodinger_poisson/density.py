import mpmath
import scipy
import numpy as np
import math
import constant as cons

# constant value
# dirac consant value (J)
hbar=1.05457266 * 10**-34
# germanium effective mass
mass = 0.041 * cons.M
# m* / (pi * hbar^2)
temp1 = mass / (cons.PI * hbar**2)
# Kb * T
kb = 8.6173303 * 10**-5
temp2 = kb * cons.T

Q = cons.Q


# function to return the electron occuypation state nk
# conduction electron concentrairion per unit box / m^-3
def electronOccupationState(eigenvalues, fermiEnergy):
    """
    return electron occupation state for each subband meaning just number of suband value sets
    """
    nk = []
    for eigenvalue in eigenvalues:
        result = temp1*mpmath.quad(lambda x: 1/(1 + mpmath.exp((x-fermiEnergy)/temp2)), [eigenvalue/Q, float('inf')])
        nk.append(result)
    return nk

# function  to calculate the electron density
def electronDensityFunction(eigenvectors, nk, eigenvalues):
    """
    return electron density numpy n(x) list not each subband

    args: 
        - eigenvectors: 1d wavefunction array of each slice
        - nk: electron occupation state
    """
    result = []
    for i in range(len(eigenvalues)):
        kth_term = [(wavefunction * wavefunction * nk[i]) for wavefunction in eigenvectors]
        result.append(kth_term)

    # sigma with each eigen value
    n = np.sum(result, axis=0)
    return n

