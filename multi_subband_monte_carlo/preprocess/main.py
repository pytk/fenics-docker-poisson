from __future__ import print_function
import numpy as np
from fenics import *
from multi_subband_monte_carlo.preprocess import schrodinger
from multi_subband_monte_carlo.preprocess import density
from multi_subband_monte_carlo.preprocess import poisson
from multi_subband_monte_carlo.plot import plotter

import sys
import time

def main(constant, device, mesh, dopant):

    V = FunctionSpace(mesh, 'CG', 1)

    electron_distribution = Constant(10000)

    iterate = 0
    start = time.time()
    # preprocess get started
    #================================================================#
    while iterate < 1:
        print("Preprocess Loop No." + str(iterate) + " get started")
        print("----------------------------------------")

        # caluculate poisson equation
        potential = poisson.poissonSolver(mesh, dopant, device, constant, electron_distribution)

        """
        if ("poisson" in device.flag):
            plotter.electrostaticPotential(device, potential, iterate)
        """
        
        # calculate schrodinger equation
        wavefunction, eigenvalue, eigenvalue_dict = schrodinger.schrodinger(mesh, potential, device, constant)

        # reinitialize electron distribution
        electron_density = np.zeros((device.nz+1, device.nx+1))

        # calculate electron density for each subband 
        wavefunction_dict = {}
        for index in range(device.nx + 1):
            # calculate occupation state of each subband
            # return a numpy array of occuaption state in each subband
            nk = density.electronOccupationState(eigenvalue[:, index], device.fermi_energy, constant)
            for subband in range(device.subband_number):
                wavefunction_dict[str(subband+1)] = wavefunction[subband+1][:, index]
            n = density.electronDensityFunction(wavefunction_dict, nk, eigenvalue[:, index])
            electron_density[:, index] = n

        # plot electron distribution
        plotter.electronDistribution(device, electron_density, iterate)

        # check if electron density get convergence or not
        density_error = []
        electron_density_flat = electron_density.flatten()
        if (iterate > 0):
            for el in range(len(electron_density_flat)):
                density_error.append(electron_distribution.vector()[el] - electron_density_flat[el])
            density_error = abs(sum(density_error)/len(density_error))

            print("Electron Density Error is: " + str(density_error[0]))


        electron_distribution = Function(V)
        # add electron to source term of Poisson equation
        electron_distribution.vector()[:] = np.array([i for i in electron_density_flat])  

        iterate  = iterate + 1
        
    #================================================================#
    # pre process got done!
    end = time.time()
    print("elpsed_time:{}".format(end - start) + "[sec]")

    return potential, wavefunction, eigenvalue_dict, electron_density

