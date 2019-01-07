from __future__ import print_function
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D


def plot_potential_distribution(device, mesh, u):
    # output potential calculated with FEniCS
    """
    output = open("out.txt", "w+")
    cordinate = mesh.coordinates()
    u_array = u.vector()[:]
    for i in range(mesh.num_vertices()-1):
        output.write('x:%8g, y:%8g, potential %g\n' % (cordinate[i][0], cordinate[i][1], u_array[i]))
    output.close()
    """

    # plot 3d figure in x-y electro static potential with matplot lib
    x = np.array([])
    y = np.array([])
    z = np.array([])
    u_array = u.vector()[:]
    grid = mesh.coordinates()
    for i in range(mesh.num_vertices()):
        x = np.append(x, grid[i][0])
        y = np.append(y, grid[i][1])
        z = np.append(z, u_array[i])
    X = np.reshape(x, (device.nx+1,device.ny+1)) 
    Y = np.reshape(y, (device.nx+1,device.ny+1))
    Z = np.reshape(z, (device.nx+1,device.ny+1))
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    surf = ax.plot_surface(X ,Y, Z, cmap='bwr', linewidth=0)
    fig.colorbar(surf)
    ax.set_title("Electro Static Potential by FEniCS(PDEs)")
    fig.show()
    plt.savefig("non-liner-poisson3d.png")


def plot_wave_function(device, mesh, wavefunction):
    x = np.array([])
    y = np.array([])
    grid = mesh.coordinates()
    for i in range(mesh.num_vertices()):
        x = np.append(x, grid[i][0])
        y = np.append(y, grid[i][1])
    X = np.reshape(x, (device.nx+1,device.ny+1)) 
    Y = np.reshape(y, (device.nx+1,device.ny+1))

    for i in range(0,3):
        Z = wavefunction[i][:][:]
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        surf = ax.plot_surface(X ,Y, Z, cmap='bwr', linewidth=0)
        fig.colorbar(surf)
        ax.set_title("Electro Static Potential by FEniCS(PDEs)")
        fig.show()
        plt.savefig("wave-function3d" + "-" + str(i) + ".png")

def plot_electric_field(device, mesh, electric_field):
    x = np.array([])
    y = np.array([])
    grid = mesh.coordinates()
    for i in range(mesh.num_vertices()):
        x = np.append(x, grid[i][0])
        y = np.append(y, grid[i][1])
    X = np.reshape(x, (device.nx+1,device.ny+1)) 
    Y = np.reshape(y, (device.nx+1,device.ny+1))

    for i in range(0,3):
        Z = electric_field[i][:][:]
        wfig = plt.figure()
        wax = wfig.add_subplot(111, projection='3d')
        wsurf = wax.plot_surface(X ,Y, Z, cmap='bwr', linewidth=0)
        wfig.colorbar(wsurf)
        wax.set_title("wave function by FEniCS(PDEs)")
        wfig.show()
        plt.savefig("wave-function3d" + "-" + str(i) + ".png")