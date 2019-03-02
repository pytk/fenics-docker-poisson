from multi_subband_monte_carlo.mainprocess.schrodinger_poisson import schrodinger
from multi_subband_monte_carlo.mainprocess.schrodinger_poisson import poisson
import time

def multiSubBand(electron_potential, wave_functions, eigen_values, electron_density, mesh, device, constant):
    """
    this function is a hub file to execute all main loop including schrodinger and poisson but we will execute juliat function in hub.jl

    Args:
        - electro static potential (Numpy Array)
        - electron density (Numpy Array)
        - doner distribution (Numpy Array)
        - eigen vector (Numpy Array)
        - eigen value (Numpy Array)
    """
    from julia import Julia
    julia = Julia()
    julia.eval("@eval Main import Base.MainInclude: include")
    from julia import Main
    Main.include("./mainprocess/initialization/main.jl")

    start = time.time()

    # convert from python object to python dictionary to calculate in julia
    device = device.__dict__

    particles, scattering_rate, Gm = Main.initialize(device, wave_functions, eigen_values, electron_density, constant)

    Main.include("./mainprocess/transport/main.jl")
    Main.include("./mainprocess/transport/scattering_rate.jl")

    print("start transport loop for each particle")

    # loop for all time step
    while(time < max_time):

        # ensemble monte carlo including drift and scat, finally it return electron density charge
        electron_density, each_subband_density, particles = Main.ensembleMonteCarlo(particles, device, Gm, electron_density, scattering_rate)

        # poisson equation
        potential = poisson_bcs.poissonSolverTest(mesh, device, constant, electron_density)

        # schrodinger equation  
        wavefunction, eigen_values = schrodinger_fenics.schrodinger(mesh, potential, device)

        # calculate scattering rate based on updated wave function
        scattering_rate, Gm = Main.getScatteringRate(wavefunction, eigen_values, device)

    elapsed_time = time.time() - start
    print("Elapsed Time(Python): %.2f [s]" % elapsed_time)

        