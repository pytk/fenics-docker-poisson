import mpmath
import scipy
import numpy as np
import math

# function to return the electron occuypation state nk
# conduction electron concentrairion per unit box / m^-3
def electronOccupationState(eigenvalues, fermi_energy, cons):
    """
    return electron occupation state for each subband meaning just number of suband value sets
    """
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

    nk = []
    for eigenvalue in eigenvalues:
        result = temp1*mpmath.quad(lambda x: 1/(1 + mpmath.exp((x-fermi_energy)/temp2)), [eigenvalue, 2*eigenvalue + fermi_energy])
        nk.append(result)
    nk = np.array(nk)
    return nk

# function  to calculate the electron density
def electronDensityFunction(each_x_eigenvectors, nk, eigenvalues):
    """
    return electron density numpy n(x) list not each subband

    args: 
        - eigenvectors: Dict{subband::String, Array}
        - nk: electron occupation state (Float)

    Note: we can get electron density distribution, hence you can figure out
    the final percentage of doner n[x] = Nd[x] 
    """
    each_x_density = []
    for subband in range(len(eigenvalues)-1):
        kth_term = [(each_z_eigenvector * each_z_eigenvector * nk[subband]) for each_z_eigenvector in each_x_eigenvectors[str(subband+1)]]
        # result :: Array[Array, Array, Array]
        each_x_density.append(kth_term)

    # sigma with each eigen value
    each_x_density = np.array(each_x_density)

    # sum electron density for each subband
    each_x_density = np.sum(each_x_density, axis=0)
    # return number of electron in a cell (1 nm)
    return each_x_density

