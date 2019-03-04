include("./initialize_particle.jl")
include("./scattering_rate.jl")
function initialize(device, wave_function, eigen_values, electron_density, constant)
    #=
    Initialize particle information, firstly we have to calculate scattering rate so that we need to set electron free flight time to each particle

    Args: 
        - device (object)
        - wave_function (numpy array) from fenics project preprocess
        - eigen_value (numpy array) from fenics project preprocess
        - electron_density (numpy array) from fenics project preprocess

    Return:
        - particle (julia dictionaly)
    =#
    # scattering_rate has too much size so we have to make it global variable and then access it without passing data as args (scattering_rate means each slice scattering_rate)
    each_slice_Gm = getScatteringRate(wave_function, eigen_values, device, constant)

    particles = initializeParticle(electron_density, device, each_slice_Gm, eigen_values, constant)

    return particles, each_slice_Gm
end