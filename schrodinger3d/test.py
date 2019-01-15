from __future__ import print_function
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

# Warning: from fenics import * will import both `sym` and
# `q` from FEniCS. We therefore import FEniCS first and then
# overwrite these objects.
from fenics import *

if __name__ == "__main__":
    mesh = UnitCubeMesh(10, 10, 10)

    V = FunctionSpace(mesh, 'Lagrange', 1)
    u = TrialFunction(V)  # Note: not TrialFunction!
    v = TestFunction(V)

    a = dot(grad(u), grad(v))*dx

    A = PETScMatrix()
    assemble(a, tensor=A)

    # Compute solution
    # Create eigensolver
    eigensolver = SLEPcEigenSolver(A)
    #eigensolver.parameters["problem_type"] = "gen_hermitian"
    eigensolver.parameters['spectrum'] = 'smallest magnitude'
    eigensolver.parameters['solver']   = 'lapack'
    eigensolver.parameters['tolerance'] = 1.e-15

    #solve for eigenvalues
    eigensolver.solve()

    u = Function(V)
    r, c, rx, cx = eigensolver.get_eigenpair(0)
    print('eigenvalue: '+ str(r))

    #assign eigenvector to function
    u.vector()[:] = rx

    plot(u)
    plt.savefig("schrodinger-fenics.png")

    # Save solution in VTK format
    file = File("wavefunction.pvd")
    file << u

