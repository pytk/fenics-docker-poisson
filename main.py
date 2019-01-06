from __future__ import print_function
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from fenics import *

import schrodinger_fenics
import poisson
import plotter

class Device(object):
    def __init__(self, xin, yin, xfi, yfi, nx, ny):
        self.xin = xin
        self.yin = yin
        self.xfi = xfi
        self.yfi = yfi
        self.nx = nx
        self.ny = ny

    def IntervalMeshCreate(self):
        mesh = IntervalMesh(self.ny, self.yin, self.yfi)
        return mesh

    def RectangleMeshCreate(self):
        mesh = RectangleMesh(Point(self.xin, self.yin), Point(self.xfi, self.yfi), self.nx, self.ny)
        return mesh


if __name__ == "__main__":

    # create mesh
    xin = 0
    yin = 0
    xfi = 30
    yfi = 30
    nx = 30
    ny = 30

    device = Device(xin, yin, xfi, yfi, nx, ny)
    rectangle_mesh = device.RectangleMeshCreate()
    potential = poisson.poissonSolver(rectangle_mesh)
    plotter.plot_potential_distribution(device, rectangle_mesh, potential)

    u = potential.vector()[:]
    temp = np.array([])
    for i in range(rectangle_mesh.num_vertices()):
        temp = np.append(temp, u[i])
    potential = np.reshape(temp, (device.nx+1,device.ny+1))

    interval_mesh = device.IntervalMeshCreate()
    electric_field, wavefunction = schrodinger_fenics.schrodinger(interval_mesh, potential)
    plotter.plot_wave_function(device, rectangle_mesh, wavefunction)