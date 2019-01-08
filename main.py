from __future__ import print_function
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from fenics import *

import schrodinger_fenics
import poisson_test
import plotter

class Device(object):
    def __init__(self, xin, yin, xfi, yfi, nx, ny, gate_ini, gate_fin):
        self.xin = xin
        self.yin = yin
        self.xfi = xfi
        self.yfi = yfi
        self.nx = nx
        self.ny = ny
        self.dx = xfi / nx
        self.dy = yfi / ny
        self.gate_ini = gate_ini
        self.gate_fin = gate_fin

    def IntervalMeshCreate(self):
        mesh = IntervalMesh(self.ny, self.yin, self.yfi)
        return mesh

    def RectangleMeshCreate(self):
        mesh = RectangleMesh(Point(self.xin, self.yin), Point(self.xfi, self.yfi), self.nx, self.ny)
        return mesh

    def craeteChargeDistribution(self, doner, constant):
        """
        Args
        ------
        doner - Array
        doner[  {xi: ..., xf: ..., yi: ..., yf: ..., nplus: ...} ]
        """
        temp = (self.nx + 1) * (self.ny + 1)
        arr = np.arange(temp).reshape(self.nx+1, self.ny+1)

        doner_density = 5.0 * 10**1

        surfacepotential = 1.0 * 10**3
        backcharge = 1.0 * 10**3
        
        arr = np.full_like(arr, doner_density)

        for d in doner:
            d["xi"] = int(d["xi"] / self.dx)
            d["xf"] = int(d["xf"] / self.dx)
            d["yi"] = int(d["yi"] / self.dy)
            d["yf"] = int(d["yf"] / self.dy)
            arr[d["xi"] : d["xf"] , d["yi"] : d["yf"] ] = float(d["nplus"])

        gate_ini = int(self.gate_ini / self.dx)
        gate_fin = int(self.gate_fin / self.dx)
        xfi = int(self.xfi / self.dx)
        arr[0: gate_ini, 0] = surfacepotential
        arr[gate_fin: xfi, 0] = surfacepotential
        arr[0: xfi, 0] = backcharge

        print(arr)
        result = arr.flatten()
        result[:]

        return result





class Constant(object):
    def __init__(self):
        self.Q = 1.60217662e-19
        self.HBAR = 5.682119514e-16
        self.H = 4.135667662e-15
        self.PI = 3.14159265359
        self.EPS = 8.854187817 * 10**-12

if __name__ == "__main__":

    # create mesh
    xin = 0
    yin = 0
    xfi = 30
    yfi = 30
    nx = 30
    ny = 30
    gate_ini = 10
    gate_fin = 10

    # doner density about n++ layer
    doner = [
        {
            "xi": 0.0,
            "xf": 3.0 * 10**-9,
            "yi": 0.0,
            "yf": 5.0 * 10**-9,
            "nplus": 1.0 * 10**21
        },
        {
            "xi": 27 * 10**-9,
            "xf": 30 * 10**-9,
            "yi": 0.0,
            "yf": 5.0 * 10**-9,
            "nplus": 1.0 * 10 **21
        }
    ]

    constant = Constant()

    device = Device(xin, yin, xfi, yfi, nx, ny, gate_ini, gate_fin)

    rectangle_mesh = device.RectangleMeshCreate()

    dopant = device.craeteChargeDistribution(doner, constant)

    potential = poisson_test.poissonSolver(rectangle_mesh, dopant, device)
    plotter.plot_potential_distribution(device, rectangle_mesh, potential)

    u = potential.vector()[:]
    temp = np.array([])
    for i in range(rectangle_mesh.num_vertices()):
        temp = np.append(temp, u[i])
    potential = np.reshape(temp, (device.nx+1,device.ny+1))

    interval_mesh = device.IntervalMeshCreate()
    electric_field, wavefunction = schrodinger_fenics.schrodinger(interval_mesh, potential)
    plotter.plot_wave_function(device, rectangle_mesh, wavefunction)