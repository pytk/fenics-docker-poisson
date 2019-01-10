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
def schrodinger_2d(mesh, potential):
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
    # mesh = IntervalMesh(ny, yin, yfi)
    print("Starting Schrodinger Wave Function")
    V = FunctionSpace(mesh, 'P', 3)

    def boundary(x, on_boundary):
        return on_boundary

    bc = DirichletBC(V, 0, boundary)

    # Define variational problem
    u = TrialFunction(V)  # Note: not TrialFunction!
    v = TestFunction(V)

    # potential calculated by poisson solver
    p = Function(V)

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
    # Create eigensolver
    eigensolver = SLEPcEigenSolver(A,M)
    eigensolver.parameters['spectrum'] = 'smallest magnitude'
    eigensolver.parameters['solver']   = 'lapack'
    eigensolver.parameters['tolerance'] = 1.e-15

    #solve for eigenvalues
    eigensolver.solve()

    u = Function(V)
    for i in range(0,2):
        #extract next eigenpair
        r, c, rx, cx = eigensolver.get_eigenpair(i)
        print('eigenvalue: '+ str(r))

        #assign eigenvector to function
        u.vector()[:] = rx

        #plot eigenfunction
        plot(u)
        plt.savefig("schrodinger-fenics.png")


def schrodinger(mesh, potentials):
    """
    Schrodinger equation solver with FEniCS
    if given potential calculated by Poisson equation,
    it can compute wave function and eigen vector meaning
    fixed potential.

    Args
    ---------------
    potential: (1D numpy array) alculated by poisson solver
    """
    subband_number = 3

    potentials = potentials.T
    
    # numpy array for storing potential
    eigen_vector = np.empty((subband_number, potentials.shape[0], potentials.shape[1]))

    eigen_value = np.empty((subband_number, potentials.shape[0], potentials.shape[1]))

    for (index, potential) in enumerate(potentials):
        # Create mesh and define function space
        # mesh = IntervalMesh(ny, yin, yfi)
        V = FunctionSpace(mesh, 'P', 1)

        def boundary(x, on_boundary):
            return on_boundary

        bc = DirichletBC(V, 0, boundary)

        # Define variational problem
        u = TrialFunction(V)  # Note: not TrialFunction!
        v = TestFunction(V)

        # potential calculated by poisson solver
        p = Function(V)

        p.vector()[:] = np.array([i for i in potential])

        a = (inner(grad(u), grad(v)) + p*u*v)*dx
        m = u*v*dx

        A = PETScMatrix()
        assemble(a, tensor=A)
        M = PETScMatrix()
        assemble(m, tensor=M)
        bc.apply(A)          # apply the boundary conditions
        bc.apply(M)

        # Compute solution
        # Create eigensolver
        eigensolver = SLEPcEigenSolver(A,M)
        eigensolver.parameters['spectrum'] = 'smallest magnitude'
        eigensolver.parameters['solver']   = 'lapack'
        eigensolver.parameters['tolerance'] = 1.e-15

        #solve for eigenvalues
        eigensolver.solve()

        u = Function(V)
        eigenvalue = [ [0] ] * subband_number
        eigenvector = [ [0] ] * subband_number
        for i in range(0,2):
            #extract next eigenpair
            r, c, rx, cx = eigensolver.get_eigenpair(i)
            print('eigenvalue: '+ str(r))
            eigenvalue[i] = r
            eigenvector[i] = rx

            #assign eigenvector to function
            u.vector()[:] = rx

            #plot eigenfunction
            plot(u)
            plt.savefig("schrodinger-fenics.png")

        for i in range(0,subband_number):
            eigen_value[i][index][:] = eigenvalue[i]
            eigen_vector[i][index][:] = eigenvector[i]
        """
        save electric field along with z-axis
        after putting the eigen value into numpy array,
        it will be transposed for z-axis.
        """

    electric_field = np.empty((subband_number, potentials.shape[0] ,potentials.shape[1]))

    for i in range(0,subband_number):
        arr = np.array(eigen_vector[i][:][:])
        arr = arr.T
        eigen_vector[i][:][:] = arr

        for (index, energy) in enumerate(eigen_value[i][:][:]):
            electric_field[i][index][:] = np.gradient(eigen_value[i][index][:])

        arr = np.array(electric_field[i][:][:])
        arr = arr.T
        electric_field[i][:][:] = arr

    return electric_field, eigen_vector

        


    


