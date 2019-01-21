from __future__ import print_function
from scipy.integrate import simps
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

# Warning: from fenics import * will import both `sym` and
# `q` from FEniCS. We therefore import FEniCS first and then
# overwrite these objects.
from dolfin import *


# make Hamiltonian for eigen value problem
def makeHamiltonian(ny, dy, potential):
    """
    return Hamiltonian for solving eigen value problem

    Args
        - ny : number of vector in Interval Mesh
        - dy : mesh particion in Interval Mesh
        - potential : 1d potential in each slice of 2d electro static potential
    """
    hamiltonian = np.zeros((ny+1, ny+1))
    for i in range(0, ny+1):
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


def schrodinger(mesh, potential, device, cons):
    """
    return normalized hamiltonian with each eigen value (n = 1, 2, 3) and wave function in rectangler mesh

    Args
        - mesh : Rectangler Mesh (Dolfin Class)
        - potential : 2d electro static potential calculated in Poisson Equation 
        - device : Original Class of device structure
        - cons : Original Class of constant value
    """
    subbands = 4
    for subband in range(1, subbands):
        # Function space of rectangle mesh
        V = FunctionSpace(mesh, 'CG', 1)

        # plot with matplotlib instead of paraview
        n = V.dim()
        d = mesh.geometry().dim()

        # Excerpt of x-axis and y-axis
        dof_coordinates = V.tabulate_dof_coordinates()
        dof_coordinates.resize((n, d))
        dof_x = dof_coordinates[:, 0]
        dof_y = dof_coordinates[:, 1]

        # dimention of rectangle mesh
        N = device.ny 
        L = device.yfi
        x, dx = np.linspace(0, L, N), L / N

        # reshape potential from 2d rectangle shape to 1d array
        potential = np.reshape(potential, (device.ny+1,device.nx+1))

        y, dy = np.linspace(0, L, N+1), L / N
        wavefunction = np.empty((device.ny+1, device.nx+1))

        # calculate eigen vector and eigen value for each slice of rectangle mesh
        for index in range(device.nx + 1):
            vector = potential[:, index]
            H = makeHamiltonian(N, dx, vector)
            w, v = np.linalg.eigh(H)
            wavefunction[:, index] = v[:,subband]

        # reshape 2d wavefunction array to 1d array
        
        X = np.linspace(device.xin, device.xfi, device.nx+1)
        Y = np.linspace(device.yin, device.yfi, device.ny+1)
        X, Y = np.meshgrid(X, Y)
        print(X.shape)
        print(wavefunction.shape)

        # plot wavefunction
        fig = plt.figure()
        ax = fig.gca(projection="3d")
        ax.plot_surface(X, Y, wavefunction, linewidth=0.2, antialiased=True, cmap=plt.cm.coolwarm)
        ax.view_init(30, -120)
        plt.savefig("img/wavefunction_" + str(subband) + ".png")


# schrodinger equation to solve eigen value in triangle quantum well by using 
# fenics eigen solver tool

def schrodinger_fenics(mesh, potential, device, cons):
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

        r, c, rx, cx = eigensolver.get_eigenpair(2)

        print("eigen value : " + str(r))

        #assign eigenvector to function
        wavefunction[:, index] = rx


    wavefunction = wavefunction.flatten()
    fig = plt.figure()
    ax = fig.gca(projection="3d")
    ax.plot_trisurf(dof_x, dof_y, wavefunction, linewidth=0.2, antialiased=True, cmap=plt.cm.CMRmap)
    ax.view_init(10, -120)
    plt.savefig("wavefunction.png")

    print("finish!!!!!!!!!!!!!!")