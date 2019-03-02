from __future__ import print_function
import numpy as np
from fenics import *
import time

import os
import sys

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

from multi_subband_monte_carlo import simulation_condition
# preprocess
from multi_subband_monte_carlo.preprocess import main as preprocess

# main process
from multi_subband_monte_carlo.mainprocess import main as mainprocess

def main():
    # create device structure
    nm = 1 * 10**-9
    structure = {
        "xin": 0,
        "zin": 0,
        "xfi": 300 * nm,
        "zfi": 20 * nm,
        "nx": 300,
        "nz": 20,
        "gate_ini": 50 * nm,
        "gate_fin": 250 * nm,
        "src": 40 * nm,
        "drain": 260 * nm,
        "subband_number": 3,
        "carrier_per_superparticle": 100000,
        "scat_list": ["acoustic", "non-polor-optics"]
    }

    # doner density about n++ layer
    doner = [
        {
            "xi": 0,
            "xf": 50 * nm,
            "zi": 0 * nm,
            "zf": 12 * nm,
            "nplus": 1.0 * 10**21
        },
        {
            "xi": 250 * nm,
            "xf": 300 * nm,
            "zi": 0 * nm,
            "zf": 12 * nm,
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
        "elasctic_constant": 12.60,
        "deformation_potential": 3.0*10**8,
        "phonon_energy": 0.03704,
        # metal (TaN)
        "metal_work_function": 5.43,
        "oxyde_affinity": 0.75
    }

    # assume initial fermi energy is 4.3 eV
    fermi_energy = 3.5

    constant = simulation_condition.ConstantValue()

    device = simulation_condition.Device(structure, material, doner, fermi_energy)

    mesh = device.RectangleMeshCreate()

    dopant = device.craeteChargeDistribution(doner, constant)

    # preprocess get started, which is a self consistent loop for schrodinger-poisson
    electron_potential, wave_functions, eigen_values, electron_density = preprocess.main(constant, device, mesh, dopant)

    # main process get started, which is a self consistent loop for schrodinger-poisson, boltzman transport equation and scattering process based on felmi's golden rule
    mainprocess.multiSubBand(electron_potential, wave_functions, eigen_values, electron_density, mesh, device, constant)
    print("All of the processes has done!")
    
# whole process get started
if __name__ == "__main__":
    if (len(sys.argv) > 1) and (sys.argv[1] == "debug"):
        import ptvsd
        print("waiting for being attached from local machine...")
        address = ('0.0.0.0', 3000)
        ptvsd.enable_attach(address)
        ptvsd.wait_for_attach()
        print("attached!!")

        plot_flag = []
        for arg in sys.argv:
            plot_flag.append(arg)

    main()
