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
from fenics import *

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


def schrodinger_test(mesh, potential, device, cons):

    V = FunctionSpace(mesh, 'CG', 1)

    # plot with matplotlib instead of paraview
    n = V.dim()
    d = mesh.geometry().dim()

    dof_coordinates = V.tabulate_dof_coordinates()
    dof_coordinates.resize((n, d))
    dof_x = dof_coordinates[:, 0]
    dof_y = dof_coordinates[:, 1]

    N = device.ny 
    L = device.yfi
    x, dx = np.linspace(0, L, N), L / N

    potential = np.reshape(potential, (device.ny+1,device.nx+1))

    #potential = potential.flatten()
    y, dy = np.linspace(0, L, N+1), L / N
    wavefunction = np.empty((device.ny+1, device.nx+1))

    for index in range(device.nx + 1):
        vector = potential[:, index]
        H = makeHamiltonian(N, dx, vector)
        np.savetxt("hamiltonian.csv", H)
        sys.exit()
        w, v = np.linalg.eigh(H)
        wavefunction[:, index] = v[:,0]
        print(v[:,0])

    """
    #plt.plot(x, v.T[0] / simps(v.T[0]**2, x)**0.5, label="gradn state")
    wavefunction = wavefunction.T
    """

    wavefunction = wavefunction.flatten()

    fig = plt.figure()
    ax = fig.gca(projection="3d")
    ax.plot_trisurf(dof_x, dof_y, wavefunction, linewidth=0.2, antialiased=True, cmap=plt.cm.CMRmap)
    ax.view_init(10, -220)
    plt.savefig("wavefunction.png")

    print("finish!!!!!!!!!!!!!!")



def schrodinger(mesh, potential, device, cons):
    """
    Schrodinger equation solver with FEniCS
    if given potential calculated by Poisson equation,
    it can compute wave function and eigen vector meaning
    fixed potential.

    Args
    ---------------
    potential: (1D numpy array) alculated by poisson solver
    """
    wavefunction = np.empty((device.ny+1, device.nx+1))
    potential = np.reshape(potential, (device.ny+1,device.nx+1))

    # plot with matplotlib instead of paraview
    V = FunctionSpace(mesh, 'Lagrange', 1)
    n = V.dim()
    d = mesh.geometry().dim()

    dof_coordinates = V.tabulate_dof_coordinates()
    dof_coordinates.resize((n, d))
    dof_x = dof_coordinates[:, 0]
    dof_y = dof_coordinates[:, 1]

    for index in range(device.nx + 1):
        mesh = IntervalMesh(device.ny, 0, device.yfi)
        vector = potential[:, index]
        V = FunctionSpace(mesh, 'Lagrange', 1)

        def boudary(x, on_boundary):
            return on_boundary

        bc = DirichletBC(V, 0, boudary)

        # difine function
        u = TrialFunction(V)
        v = TestFunction(V)
        Potenial = Function(V)
        Potenial.vector().array()[:] = vector

        # difine problem
        temp = -1 * pow(cons.HBAR, 2) / (2*device.material["electron_effective_mass"]*cons.M)
        a = (inner(temp * grad(u), grad(v)) + Potenial*u*v)*dx
        m = u*v*dx

        # assemble stiffness matrix
        A = PETScMatrix()
        assemble(a, tensor=A)
        M = PETScMatrix()
        assemble(m, tensor=M)
        bc.apply(A)
        bc.apply(M)

        # create eigensolver
        eigensolver = SLEPcEigenSolver(A, M)
        eigensolver.parameters['spectrum'] = 'smallest magnitude'
        eigensolver.parameters['solver'] = 'lapack'
        eigensolver.parameters['tolerance'] = 1e-15

        #solve for eigenvalues
        eigensolver.solve()

        r, c, rx, cx = eigensolver.get_eigenpair(0)

        print("eigen value : " + str(r))

        #assign eigenvector to function
        wavefunction[:, index] = rx


    wavefunction = wavefunction.flatten()
    fig = plt.figure()
    ax = fig.gca(projection="3d")
    ax.plot_trisurf(dof_x, dof_y, wavefunction, linewidth=0.2, antialiased=True, cmap=plt.cm.CMRmap)
    ax.view_init(10, -220)
    plt.savefig("wavefunction.png")

    print("finish!!!!!!!!!!!!!!")
    

"""
        for i in range(0,subband_number):
            eigen_value[i][index][:] = eigenvalue[i]
            eigen_vector[i][index][:] = eigenvector[i]

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
"""

        


    


