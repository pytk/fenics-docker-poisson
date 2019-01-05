from __future__ import print_function
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as 
import math

def makeHamiltonian(mesh, potential):
    hamiltonian = np.zeros(math.sqrt(mesh.num_vertices()))
    for i in range(mesh.num_vertices()):
        for 
