from __future__ import print_function
from scipy.integrate import simps
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import sys

# Warning: from fenics import * will import both `sym` and
# `q` from FEniCS. We therefore import FEniCS first and then
# overwrite these objects.
from dolfin import *


def potential(x):
    return (1/1000) * x

def makeHamiltonian(ny, dy, potential):
    hamiltonian = np.zeros((ny, ny))
    for i in range(0, ny):
        for dx in range(-1, 1):
            j = i + dx
            v = 0.0
            if dx == 0:
                v = (2/dy**2 + potential[i])
            elif dx == 1:
                v = -1/dy**2
            elif dx == -1:
                v = -1/dy**2

            if j >= 1 and j <= ny+1:
                hamiltonian[i, j] = v
    print(hamiltonian)
    return hamiltonian

def schrodinger():
    """
    Schrodinger equation solver with FEniCS
    if given potential calculated by Poisson equation,
    it can compute wave function and eigen vector meaning
    fixed potential.

    Args
    ---------------
    potential: (1D numpy array) alculated by poisson solver
    """

    N = 1000 
    L = 10
    dx = 0.01
    x, dy = np.linspace(0, L, N), L / N
    # plot with matplotlib instead of paraview

    # crate potential array
    vector = np.array([potential(i) for i in range(1000)])

    H = makeHamiltonian(N, dx, vector)
    w, v = np.linalg.eigh(H)
    print("----------")
    print(w[:3])
    wavefunction = v.T[1]
    print("----------")
    print(v)

    plt.plot(x, wavefunction)
    plt.plot(x, vector)
    plt.savefig("wavefunction.png")

    print("finish!!!!!!!!!!!!!!")

if __name__ == "__main__":
    print("start")
    schrodinger()
    print("finish!!")