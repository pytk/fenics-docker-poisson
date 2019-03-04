from multi_subband_monte_carlo.mainprocess.schrodinger_poisson import schrodinger
from multi_subband_monte_carlo.mainprocess.schrodinger_poisson import poisson
import time

def multiSubBand(electric_field, wave_functions, eigen_values, electron_density, mesh, device, constant):
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

    time_step = device.time_step
    final_time = device.final_time
    current_time = 0.0

    # convert from python object to python dictionary to calculate in julia
    device_dictionary = device.__dict__
    
    particles, Gm = Main.initialize(device_dictionary, wave_functions, eigen_values, electron_density, constant)

    # make scattering rate global scope to avoid heavy trafic in a process of passing data as args

    Main.include("./mainprocess/transport/ensemble_monte_carlo.jl")
    Main.include("./mainprocess/transport/scattering_rate.jl")

    print("start transport loop for each particle")

    # loop for all time step
    while(current_time < final_time):
        # ensemble monte carlo including drift and scat, finally it return electron density charge
        #electron_density, each_subband_density, particles = Main.ensembleMonteCarlo(particles, device_dictionary, Gm, electron_density, scattering_rate, current_time)
        electron_density, each_subband_density, particles = Main.ensembleMonteCarlo(particles, device_dictionary, electron_density, current_time, electric_field, Gm)

        print("poisson equation for main loop")
        # poisson equation
        potential = poisson.poissonSolver(mesh, device, constant, electron_density)

        print("Schrodinger equation for main loop")
        # schrodinger equation  
        wavefunction, eigen_values = schrodinger.schrodinger(mesh, potential, device, constant)

        print("Caluculation for scattering rate in main loop")
        # calculate scattering rate based on updated wave function
        Main.getScatteringRate(wavefunction, eigen_values, device_dictionary, constant)

        if(current_time + time_step >= final_time):
            time_step = final_time - current_time
        current_time = current_time + time_step

    elapsed_time = time.time() - start
    print("Elapsed Time(Python): %.2f [s]" % elapsed_time)

        