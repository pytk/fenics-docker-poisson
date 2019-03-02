from __future__ import print_function
from scipy.integrate import simps
import numpy as np
import math

# Warning: from fenics import * will import both `sym` and
# `q` from FEniCS. We therefore import FEniCS first and then
# overwrite these objects.
from dolfin import *

# make Hamiltonian for eigen value problem
def makeHamiltonian(ny, dy, potential, device, cons):
    """
    return Hamiltonian for solving eigen value problem

    Args
        - ny : number of vector in Interval Mesh
        - dy : mesh particion in Interval Mesh
        - potential : 1d potential in each slice of 2d electro static potential
    """

    hbar=1.05457266 * 10**-34
    hamiltonian = np.zeros((ny+1, ny+1))
    constant_value = (hbar*hbar)/(2.0*dy*dy)
    effective_mass = device.material["electron_effective_mass"] * cons.M
    for i in range(0, ny+1):
        for dx in range(-1, 1):
            j = i + dx
            v = 0.0
            if dx == 0:
                v = (2*constant_value/effective_mass + potential[i]*cons.Q)
            elif dx == 1:
                v = -1*constant_value/effective_mass
            elif dx == -1:
                v = -1*constant_value/effective_mass
            else:
                v = 0
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

    Return
        - wavefunction_dict: Dict{"subband_number": [nx, nz], "subband_number": [nx, nz], ......}
    """
    subbands = device.subband_number+1
    # reshape potential from 2d rectangle shape to 1d array
    #potential = np.array([i for i in potential])
    #potential = np.reshape(potential, (device.ny+1,device.nx+1))

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
    N = device.nz 
    L = device.zfi
    x, dx = np.linspace(0, L, N), L / N

    y, dy = np.linspace(0, L, N+1), L / N
    wavefunction = np.zeros((device.nz+1, device.nx+1))
    eigenvalue = np.zeros((subbands, device.nx+1))

    wavefunction_dict = {}
    eigenvalue_dict = {}

    for subband in range(0, subbands):
        # calculate eigen vector and eigen value for each slice of rectangle mesh
        for index in range(device.nx + 1):
            #if(device.gate_ini < index*device.dx and index*device.dx < device.gate_fin):
            vector = potential[:, index]
            H = makeHamiltonian(N, dx, vector, device, cons)
            w, v = np.linalg.eigh(H)
            temp = v[:, subband]

            # eigen vector is already normalized!!!
            wavefunction[:, index] = -1*temp
            eigenvalue[subband, index] = w[subband] / cons.Q

        wavefunction_dict[subband+1] = wavefunction
        eigenvalue_dict[subband+1] = eigenvalue[subband][:]

        """
        # reshape 2d wavefunction array to 1d array
        if("schrodinger" in device.flag):
            X = np.linspace(device.xin, device.xfi, device.nx+1)
            Y = np.linspace(device.yin, device.yfi, device.ny+1)
            X, Y = np.meshgrid(X, Y)

            # plot wavefunction
            fig = plt.figure()
            ax = fig.gca(projection="3d")
            ax.plot_surface(X, Y, wavefunction, linewidth=0.2, antialiased=True, cmap=plt.cm.coolwarm)
            ax.view_init(10, -120)
            plt.savefig("img/wave/wavefunction_" + str(iterate) + "-" + str(subband) + ".png")
        """

    print("Schrodinger Equation got finished!")

    return wavefunction_dict, eigenvalue, eigenvalue_dict
