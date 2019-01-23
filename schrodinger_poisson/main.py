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
    def __init__(self, structure, material):
        self.xin = structure["xin"]
        self.yin = structure["yin"]
        self.xfi = structure["xfi"]
        self.yfi = structure["yfi"]
        self.nx = structure["nx"]
        self.ny = structure["ny"]
        self.dx = self.xfi / self.nx
        self.dy = self.yfi / self.ny
        self.gate_ini = structure["gate_ini"]
        self.gate_fin = structure["gate_fin"]
        self.src = structure["src"]
        self.drain = structure["drain"]
        self.material = material
        self.doner = structure["doner"]

    def IntervalMeshCreate(self):
        mesh = IntervalMesh(self.ny, self.yin, self.yfi)
        return mesh

    def RectangleMeshCreate(self):
        mesh = RectangleMesh(Point(self.xin, self.yin), Point(self.xfi, self.yfi), self.nx, self.ny)
        return mesh

    def ChangeDolfinVectorToNumpy(self, array1, array2, array3):
        """
        change fenics PETSc vector into numpy array. The original PETSc vector itself is really
        unuseful to excerpt or get slice from 2d original vector.

        Args:
            dof_coordinates = V.tabulate_dof_coordinates()
            dof_coordinates.resize((n, d))
            dof_x = dof_coordinates[:, 0]
            dof_y = dof_coordinates[:, 1]

            - array1: dof_x
            - array2: dof_y
            - array3: potential

        Return: 
            array: (2D numpy array include electro static potetial)

            you can use returned array without any change
        """
        size = np.array(array3).shape[0]
        # check if all array has same dimention
        assert (np.array(array1).shape[0] == np.array(array3).shape[0] and np.array(array1).shape[0] == np.array(array3).shape[0]), "Not mutch array size!"
        array = [[a, b, c] for a, b, c in zip(array1, array2, array3)]

        # sort array with composite key value (x, y)
        sorted_array = sorted(array, key=lambda x: (x[0], x[1]))
        array = np.array(sorted_array)

        # reshape like [x, y, potential]
        array = np.reshape(array, (size, 3))

        # excerpt only potential
        array = array[:, 2]

        array = np.reshape(array, (self.nx+1, self.ny+1))
        array = array.T

        return array

    def craeteChargeDistribution(self, doner, constant):
        """
        Args
        ------
        doner - Array
        doner[  {xi: ..., xf: ..., yi: ..., yf: ..., nplus: ...} ]
        """
        temp = (self.nx + 1) * (self.ny + 1)
        arr = np.arange(temp).reshape(self.ny+1, self.nx+1)

        doner_density = self.doner
        
        arr = np.full_like(arr, doner_density)


        for d in doner:
            d["xi"] = int(d["xi"] / self.dx)
            d["xf"] = int(d["xf"] / self.dx)
            d["yi"] = int(d["yi"] / self.dy)
            d["yf"] = int(d["yf"] / self.dy)
            arr[d["yi"] : d["yf"] , d["xi"] : d["xf"] ] = float(d["nplus"])
        
        result = arr.flatten()
        
        return result

class ConstantValue(object):
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
        address = ('127.0.0.1', 8000)
        ptvsd.enable_attach(address)
        ptvsd.wait_for_attach()
        print("attached")

    nm = 1 * 10**-9

    # create mesh
    structure = {
        "nm": 1 * 10**-9,
        "xin": 0,
        "yin": 0,
        "xfi": 300 * nm,
        "yfi": 20 * nm,
        "nx": 300,
        "ny": 200,
        "gate_ini": 50 * nm,
        "gate_fin": 250 * nm,
        "src": 40 * nm,
        "drain": 260 * nm,
        "doner": 1.0 * 10**21 
    }

    # doner density about n++ layer
    doner = [
        {
            "xi": 0,
            "xf": 50 * nm,
            "yi": 0 * nm,
            "yf": 12 * nm,
            "nplus": 1.0 * 10**21
        },
        {
            "xi": 250 * nm,
            "xf": 300 * nm,
            "yi": 0 * nm,
            "yf": 12 * nm,
            "nplus": 1.0 * 10 **21
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

    constant = ConstantValue()

    device = Device(structure, material)

    rectangle_mesh = device.RectangleMeshCreate()

    dopant = device.craeteChargeDistribution(doner, constant)

    #plotter.electron_density(device, rectangle_mesh, dopant)

    potential = poisson_bcs.poissonSolverTest(rectangle_mesh, dopant, device, constant)
    #plotter.plot_potential_distribution(device, rectangle_mesh, potential)

    """
    V = FunctionSpace(rectangle_mesh, 'CG', 1)
    new_potential = Function(V)
    new_potential.vector()[:] = np.array([-constant.Q*j for j in potential.vector()[:]])
    u = potential.vector()[:]
    temp = np.array([])
    for i in range(rectangle_mesh.num_vertices()):
        temp = np.append(temp, u[i])
    potential = np.reshape(temp, (device.nx+1,device.ny+1))
    new_potential = np.reshape(new_potential, (device.nx+1,device.ny+1))
    """
    
    schrodinger_fenics.schrodinger(rectangle_mesh, potential, device, constant)
    