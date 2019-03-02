from dolfin import *
import numpy as np

class Device(object):
    def __init__(self, structure, material, doner, fermi_energy):
        self.xin = structure["xin"]
        self.zin = structure["zin"]
        self.xfi = structure["xfi"]
        self.zfi = structure["zfi"]
        self.nx = structure["nx"] - 1
        self.nz = structure["nz"] - 1
        self.dx = self.xfi / self.nx
        self.dz = self.zfi / self.nz
        self.gate_ini = structure["gate_ini"]
        self.gate_fin = structure["gate_fin"]
        self.src = structure["src"]
        self.drain = structure["drain"]
        self.material = material
        self.doner = doner
        self.fermi_energy = fermi_energy
        self.subband_number = structure["subband_number"]
        self.carrier_per_superparticle = structure["carrier_per_superparticle"]
        self.scat_list = structure["scat_list"]

    def IntervalMeshCreate(self):
        mesh = IntervalMesh(self.ny, self.zin, self.zfi)
        return mesh

    def RectangleMeshCreate(self):
        mesh = RectangleMesh(Point(self.xin, self.zin), Point(self.xfi, self.zfi), self.nx, self.nz)
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

        array = np.reshape(array, (self.nx+1, self.nz+1))
        array = array.T

        return array

    def craeteChargeDistribution(self, doner, constant):
        """
        Args
        ------
        doner - Array
        doner[  {xi: ..., xf: ..., yi: ..., yf: ..., nplus: ...} ]
        """
        temp = (self.nx + 1) * (self.nz + 1)
        arr = np.arange(temp).reshape(self.nz+1, self.nx+1)

        # flatter doner density with some resonable value which is much less than nplus density in N plus layer
        doner_density = self.doner[0]["nplus"] * self.dx * self.dz
        arr = np.full_like(arr, doner_density * (10**-2))


        for d in doner:
            d["xi"] = int(d["xi"] / self.dx)
            d["xf"] = int(d["xf"] / self.dx)
            d["zi"] = int(d["zi"] / self.dz)
            d["zf"] = int(d["zf"] / self.dz)
            arr[d["zi"] : d["zf"] , d["xi"] : d["xf"] ] = float(d["nplus"])
        
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