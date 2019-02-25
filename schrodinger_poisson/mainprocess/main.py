import schrodinger
import poisson
import time

def multiSubBand(device, wave_function, eigen_values, electron_density, mesh):
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
    Main.include("./../initialization/main.jl")

    start = time.time()

    particles, scattering_rate, Gm = Main.initialize(device, wave_function, eigen_values, electron_density)

    Main.include("./transport/main.jl")
    Main.include("./transport/scattering_rate.jl")

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

        