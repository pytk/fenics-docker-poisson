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
    N = device.ny 
    L = device.yfi
    x, dx = np.linspace(0, L, N), L / N

    potential = np.array([-1 * i for i in potential.vector()[:]])
    potential = np.reshape(potential, (device.ny+1,device.nx+1))

    np.savetxt("potential.csv", potential)

    #potential = potential.flatten()

    #wavefunction = np.empty((potential.shape[0], potential.shape[1]))

    """
    for (index, p) in enumerate(potential):
        vector = p
        H = makeHamiltonian(ny, dy, vector)
        w, v = np.linalg.eigh(H)
        wavefunction[index][:] = v[:,0]
    """

    #plt.plot(x, v.T[0] / simps(v.T[0]**2, x)**0.5, label="gradn state")
    plt.savefig("wavefunction.png")

    """
    np.savetxt("array.csv", wavefunction)
    wavefunction = wavefunction.T
    wavefunction = wavefunction.flatten()

    rectangle_mesh = device.RectangleMeshCreate()
    V = FunctionSpace(rectangle_mesh, 'CG', 1)

    u = Function(V)

    u.vector()[:] = np.array([i for i in wavefunction])
    plot(u)
    plt.savefig("wavefunction.png")

    file = File("wavefunction.pvd")
    file << u

    return u
    """





# Use SymPy to compute f from the manufactured solution u
def schrodinger_2d(mesh, potential, device, cons):
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
    V = FunctionSpace(mesh, 'CG', 1)

    def boundary(x, on_boundary):
        return on_boundary

    # Initialize mesh function for interior domains
    domains = CellFunction("size_t", mesh)
    domains.set_all(0)

    bc = DirichletBC(V, 0, boundary)

    dx = Measure("dx", subdomain_data=domains)

    # Define variational problem
    u = TrialFunction(V)  # Note: not TrialFunction!
    v = TestFunction(V)
    Vpot = Function(V)

    Vpot.vector()[:] = np.array([-1 * i for i in potential.vector()])

    # p.vector()[:] = np.array([i for i in potential.vector()])
    temp = -1 * pow(cons.HBAR, 2) / (2*device.material["electron_effective_mass"]*cons.M)

    a = (inner(temp * grad(u), grad(v)) + Vpot*u*v)*dx(0)
    m = u*v*dx(0)

    A = PETScMatrix()
    M = PETScMatrix()
    _ = PETScVector()
    L = Constant(0.0)*v*dx(0)

    assemble_system(a, L, A_tensor=A, b_tensor=_)
    assemble_system(m, L, A_tensor=M, b_tensor=_)

    # Compute solution
    # Create eigensolver
    eigensolver = SLEPcEigenSolver(A,M)
    #eigensolver.parameters["problem_type"] = "gen_hermitian"
    eigensolver.parameters['spectrum'] = 'smallest magnitude'
    eigensolver.parameters['solver']   = 'lapack'
    eigensolver.parameters['tolerance'] = 1.e-15

    #solve for eigenvalues
    eigensolver.solve()

    u = Function(V)
    eigenvectors = []
    eigenvalues = []
    rx_list = []

    # extract first eigenpair
    r, c, rx, cx = eigensolver.get_eigenpair(0)

    #assign eigenvector to function
    u.vector()[:] = rx

    plot(u)
    plt.savefig("schrodinger-fenics.png")

    # Save solution in VTK format
    file = File("wavefunction.pvd")
    file << u

    i = 0
    while (i < eigensolver.get_number_converged() and r < 1.8):#fermiEnergy):
        #extract next eigenpair
        r, c, rx, cx = eigensolver.get_eigenpair(i)
        #assign eigenvector to function
        #rx = np.array([q*xsi_factor for q in rx])
        #if(r > 0.24):
        #  break
        
        u.vector()[:] = rx
        eigenvalues.append(r)
        eigenvectors.append([j for j in u.vector()])
        
        rx_list.append(rx) 
        
        #plot eigenfunction
        print('eigenvalue: ' + str(r))
        #plot(u).update(u)
        #interactive()
    
        #increment i
        i = i+1


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
    subband_number = 3

    potential = np.array([i for i in potential.vector()[:]])
    potential = np.reshape(potential, (device.ny+1,device.nx+1))

    potential = potential.T

    #potential = potential.flatten()

    wavefunction = np.empty((potential.shape[0], potential.shape[1]))

    for (index, p) in enumerate(potential):
        # Create mesh and define function space
        # mesh = IntervalMesh(ny, yin, yfi)

        print("Starting Schrodinger Wave Function")
        V = FunctionSpace(mesh, 'CG', 1)


        # Initialize mesh function for interior domains
        domains = CellFunction("size_t", mesh)
        domains.set_all(0)

        dx = Measure("dx", subdomain_data=domains)

        # Define variational problem
        u = TrialFunction(V)  # Note: not TrialFunction!
        v = TestFunction(V)
        Vpot = Function(V)

        Vpot.vector()[:] = np.array([i for i in p])

        # p.vector()[:] = np.array([i for i in potential.vector()])
        temp = -1 * pow(cons.HBAR, 2) / (2*device.material["electron_effective_mass"]*cons.M)

        a = (inner(temp * grad(u), grad(v)) + Vpot*u*v)*dx(0)
        m = u*v*dx(0)

        """
        A = PETScMatrix()
        assemble(a, tensor=A)
        M = PETScMatrix()
        assemble(m, tensor=M)
        bc.apply(A)          # apply the boundary conditions
        bc.apply(M)
        """

        A = PETScMatrix()
        M = PETScMatrix()
        _ = PETScVector()
        L = Constant(0.0)*v*dx(0)

        assemble_system(a, L, A_tensor=A, b_tensor=_)
        assemble_system(m, L, A_tensor=M, b_tensor=_)

        # Compute solution
        # Create eigensolver
        eigensolver = SLEPcEigenSolver(A,M)
        eigensolver.parameters['spectrum'] = 'smallest magnitude'
        eigensolver.parameters['solver']   = 'lapack'
        eigensolver.parameters['tolerance'] = 1.e-15

        #solve for eigenvalues
        eigensolver.solve()

        r, c, rx, cx = eigensolver.get_eigenpair(0)

        print("eigen value : " + str(r))

        #assign eigenvector to function
        wavefunction[index][:] = rx

    wavefunction = wavefunction.T
    print(wavefunction)
    wavefunction = wavefunction.flatten()

    rectangle_mesh = device.RectangleMeshCreate()
    V = FunctionSpace(rectangle_mesh, 'CG', 1)

    u = Function(V)

    u.vector()[:] = np.array([i for i in wavefunction])
    plot(u)
    plt.savefig("wavefunction.png")

    file = File("wavefunction.pvd")
    file << u
        
    return u
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

        


    


