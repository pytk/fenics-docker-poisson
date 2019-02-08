from __future__ import print_function
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from fenics import *

import schrodinger_fenics
import density
import poisson_bcs
import plotter
import sys
import time

class Device(object):
    def __init__(self, structure, material, flag):
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
        self.nplus = structure["nplus"]
        self.flag = flag
        self.subband_number = structure["subband_number"]

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
        self.Number_of_Particle = 100000

if __name__ == "__main__":
    start = time.time()
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
        "doner": 1.0 * 10**21,
        "nplus": [50 * nm, 250 * nm],
        "subband_number": 4
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
        "non-parabolicity": 0.3,
        # metal (TaN)
        "metal_work_function": 5.43,
        "oxyde_affinity": 0.75
    }

    # assume initial fermi energy is 4.3 eV
    fermi_energy = 3.5

    constant = ConstantValue()

    plot_flag = []
    for arg in sys.argv:
        plot_flag.append(arg)

    device = Device(structure, material, plot_flag)

    rectangle_mesh = device.RectangleMeshCreate()

    dopant = device.craeteChargeDistribution(doner, constant)

    V = FunctionSpace(rectangle_mesh, 'CG', 1)

    electron_distribution = Constant(0.0)

    iterate = 0

    while iterate < 5:
        print("Loop No." + str(iterate) + " get started")
        print("----------------------------------------")

        potential = poisson_bcs.poissonSolverTest(rectangle_mesh, dopant, device, constant, electron_distribution)
        #plotter.plot_potential_distribution(device, rectangle_mesh, potential)

        if ("poisson" in device.flag):
            plotter.electrostaticPotential(device, potential, iterate)
        
        wavefunction, eigenvalue = schrodinger_fenics.schrodinger(rectangle_mesh, potential, device, constant, iterate)

        # reinitialize electron distribution
        electron_density = np.zeros((device.ny+1, device.nx+1))

        wavefunction_dict = {}
        for index in range(device.nx + 1):
            # calculate occupation state of each subband
            # return a numpy array of occuaption state in each subband
            nk = density.electronOccupationState(eigenvalue[:, index], fermi_energy)
            for subband in range(device.subband_number):
                wavefunction_dict[str(subband)] = wavefunction[subband][:, index]
            n = density.electronDensityFunction(wavefunction_dict, nk, eigenvalue[:, index])
            electron_density[:, index] = n

        # plot electron distribution
        if ("density" in device.flag):
            plotter.electronDistribution(device, electron_density, iterate)

        # check if electron density get convergence or not
        density_error = []
        electron_density = electron_density.flatten()
        if (iterate > 0):
            for el in range(len(electron_density)):
                density_error.append(electron_distribution.vector()[el] - electron_density[el])
            density_error = abs(sum(density_error)/len(density_error))

            print("Electron Density Error is: " + str(density_error[0]))

        electron_distribution = Function(V)
        # add electron to source term of Poisson equation
        electron_distribution.vector()[:] = np.array([i for i in electron_density])  

        iterate  = iterate + 1
        

    end = time.time()
    print("elpsed_time:{}".format(end - start) + "[sec]")
