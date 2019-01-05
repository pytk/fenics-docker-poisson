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

# Use SymPy to compute f from the manufactured solution u
def schrodinger(yin, yfi, ny, potential):
    """
    Schrodinger equation solver with FEniCS
    if given potential calculated by Poisson equation,
    it can compute wave function and eigen vector meaning
    fixed potential.

    Args
    ---------------
    potential: (1D numpy array) alculated by poisson solver
    """
    # Create mesh and define function space
    mesh = IntervalMesh(ny, yin, yfi)
    V = FunctionSpace(mesh, 'P', 3)

    def boundary(x, on_boundary):
        return on_boundary

    bc = DirichletBC(V, 0, boundary)

    # Define variational problem
    u = TrialFunction(V)  # Note: not TrialFunction!
    v = TestFunction(V)

    # potential calculated by poisson solver
    p = function(V)

    p.vector()[:] = np.array([i for i in potential.vector()])

    a = (inner(grad(u), grad(v)) + p*u*v)*dx
    m = u*v*dx

    A = PETScMatrix()
    assemble(a, tensor=A)
    M = PETScMatrix()
    assemble(m, tensor=M)
    bc.apply(A)          # apply the boundary conditions
    bc.apply(M)

    # Compute solution
    #create eigensolver
    eigensolver = SLEPcEigenSolver(A,M)
    eigensolver.parameters['spectrum'] = 'smallest magnitude'
    eigensolver.parameters['solver']   = 'lapack'
    eigensolver.parameters['tolerance'] = 1.e-15

    #solve for eigenvalues
    eigensolver.solve()

    u = Function(V)
    for i in range(0,3):
        #extract next eigenpair
        r, c, rx, cx = eigensolver.get_eigenpair(i)
        print('eigenvalue: '+ str(r))

        #assign eigenvector to function
        u.vector()[:] = rx

        #plot eigenfunction
        plot(u)
        plt.savefig("schrodinger-fenics.png")


