from __future__ import print_function
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import math
from scipy.integrate import simps

def makeHamiltonian(ny, dy, potential):
    hamiltonian = np.zeros((ny+1, ny+1))
    for i in range(1, ny+1):
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
    return hamiltonian

# test function
def potential(ny):
    vector = np.zeros(ny+1)
    for i in range(1, ny+1):
        vector[i] = (1/ny)*i
    return vector



# main function
if __name__ == "__main__":
    ny = 100
    L = 1
    y, dy = np.linspace(-L, L, ny), 2*L/ny
    vector = potential(ny)
    H = makeHamiltonian(ny, dy, vector)
    w, v = np.linalg.eigh(H)

    plt.plot(v[0:ny, 0:3])
    plt.show()
    plt.savefig("schrodinger.png")
