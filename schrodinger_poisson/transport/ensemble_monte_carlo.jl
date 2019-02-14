function ensembleMonteCarlo(particle, device, Gm)
    #=
    Args:
        - particle: Dict{Int32, Dict{String, Float64}}
        particle[particle_number][property]
            - property:
                -subband: Int32 the subband the particle belong to
                -kx: Float64 wave vector along x axis in wave space
                -ky: Float64 wave vctor along y axis in wave space
                -kz: Float64 wave vctor along z axis in wave space
                -x: Float64 
                -z: Float64
                -ts: Float64 time which couse scattering rate

        - device: Object that devecice structure and simulation condition

        - Gm: Dict{Int32, Float64}

    Return:
        - particle: Dict{Int32, Dict{String, Float64}} -> new particle information
        particle[particle_number][property]
            - property:
                -kx: Float64 wave vector along x axis in wave space
                -ky: Float64 wave vctor along y axis in wave space
                -kz: Float64 wave vctor along z axis in wave space
                -x: Float64 
                -z: Float64
                -ts: Float64 time which couse scattering rate
    =#

    # dx , dz
    dx = device["dx"]
    dz = device["dz"]

    nx = device["nx"]
    nz = device["nz"]

    start_time = 0.0
    time_step =  start_time + device["time_step"]

    particle_number = 1
    time = start_time
    scattering_time = particle[particle_number]["ts"]

    temp_particle = particle

    while true
        τ = scattering_time - time

        # drift process from the point to other point that particle can move during the free flight time and update only particle position
        temp_particle = drift(τ, temp_particle[particle_number])

        # just take the x position and z position after drift process
        i = trunc(Int, temp_particle["x"]/dx)+1
        j = trunc(Int, temp_particle["z"]/dz)+1

        if i <= 1 i=1 end
        if j <= 1 j=1 end
        if i >= nx i=nx end
        if j >= nz j=nz end

        # scat process where a scat will get energy of particle changed based on the scattering rate on the particle's position and update only particle energy
        temp_particle = scat(temp_particle)

        # set next happen scat time from random value and subband the particle belong to
        scattering_time = time - log(rand())/Gm[i][temp_particle["subband"]]

        # if the time scat process occuerd later over time step of simulation, it must be forced to move on to next particle process and next time step
        if time_step <= scattering_time break end

    end

    τ = scattering_time - time

    # TODO: function to calculate drift process which update particle information
    temp_particle = drift(τ, temp_particle[particle_number])

    # check if a particle is going out from the right edge of the device







    