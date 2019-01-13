from __future__ import print_function
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from fenics import *

import schrodinger_fenics
import poisson_bcs
import plotter
import sys

class Device(object):
    def __init__(self, xin, yin, xfi, yfi, nx, ny, gate_ini, gate_fin, src, drain, material, oxyde, acceptor):
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
        self.src = src
        self.drain = drain
        self.material = material
        self.oxyde = oxyde
        self.acceptor = acceptor

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
        arr = np.arange(temp).reshape(self.ny+1, self.nx+1)

        doner_density = 1.0 * 10**3

        #surfacepotential = 1.0 * 10**19
        backcharge = -2.0 * 10**3
        
        arr = np.full_like(arr, doner_density)


        for d in doner:
            d["xi"] = int(d["xi"] / self.dx)
            d["xf"] = int(d["xf"] / self.dx)
            d["yi"] = int(d["yi"] / self.dy)
            d["yf"] = int(d["yf"] / self.dy)
            arr[d["yi"] : d["yf"] , d["xi"] : d["xf"] ] = float(d["nplus"])

        gate_ini = int(self.gate_ini / self.dx)
        gate_fin = int(self.gate_fin / self.dx)
        xfi = int(self.xfi / self.dx)
        
        #arr[0: gate_ini, 0] = surfacepotential
        #arr[gate_fin: xfi, 0] = surfacepotential
        arr[0, :] = backcharge

        arr = arr / constant.EPS / 10**-9
        arr = arr * constant.Q

        np.savetxt("array.csv", arr)
        
        result = arr.flatten()
        
        return result

class Constant(object):
    def __init__(self):
        self.Q = 1.60217662e-19
        self.HBAR = 6.582119514 * 10**-16
        self.H = 6.626070040 * 10**-34
        self.PI = 3.14159265359
        self.EPS = 8.854187817 * 10**-12
        self.T = 300
        self.M = 9.10938356 * 10**-31
        self.KB = 1.38064852* 10**-23

if __name__ == "__main__":
    if (len(sys.argv) > 1) and (sys.argv[1] == "debug"):
        import ptvsd
        print("waiting...")
        ptvsd.enable_attach("my_secret", address=('0.0.0.0', 53005))
        ptvsd.wait_for_attach()

    # create mesh
    nm = 1 * 10**-9
    xin = 0
    yin = 0
    xfi = 300 * nm
    yfi = 50 * nm
    nx = 100
    ny = 50
    gate_ini = 125 * nm
    gate_fin = 175 * nm
    src = 100 * nm
    drain = 200 * nm
    oxyde = 3 * nm
    acceptor = 1 * 10**21

    # doner density about n++ layer
    doner = [
        {
            "xi": 0,
            "xf": 100 * nm,
            "yi": 0 * nm,
            "yf": 50 * nm,
            "nplus": -1.0 * 10**4
        },
        {
            "xi": 200 * nm,
            "xf": 300 * nm,
            "yi": 0 * nm,
            "yf": 50 * nm,
            "nplus": -1.0 * 10 **4
        }
    ]

    material = {
        # Germanium
        "electron_effective_mass": 0.041,
        "hevy_hole_effective_mass": 0.28,
        "band_gap": 0.66,
        "workf_function": 4.05,
        "conduction_band_min": 0.66,
        "valence_band_min": 0.0,
        "electron_affinity": 1.2,
        "intrinsic_carrier_concentration": 2.0 * 10**19,
        "effective_conduction_band_density": 1.0 * 10**25,
        "effective_valence_band_density": 5.0 * 10**24,
        # metal (TaN)
        "metal_work_function": 5.43
    }

    constant = Constant()

    device = Device(xin, yin, xfi, yfi, nx, ny, gate_ini, gate_fin, src, drain, material, oxyde, acceptor)

    rectangle_mesh = device.RectangleMeshCreate()

    dopant = device.craeteChargeDistribution(doner, constant)

    potential = poisson_bcs.poissonSolver(rectangle_mesh, dopant, device, constant)
    plotter.plot_potential_distribution(device, rectangle_mesh, potential)

    V = FunctionSpace(rectangle_mesh, 'CG', 1)
    new_potential = Function(V)
    new_potential.vector()[:] = np.array([-constant.Q*j for j in potential.vector()[:]])
    """
    u = potential.vector()[:]
    temp = np.array([])
    for i in range(rectangle_mesh.num_vertices()):
        temp = np.append(temp, u[i])
    potential = np.reshape(temp, (device.nx+1,device.ny+1))
    
    potential = np.reshape(new_potential, (device.nx+1,device.ny+1))

    interval_mesh = device.IntervalMeshCreate()
    electric_field, wavefunction = schrodinger_fenics.schrodinger(interval_mesh, potential)
    plotter.plot_wave_function(device, rectangle_mesh, wavefunction)
    """