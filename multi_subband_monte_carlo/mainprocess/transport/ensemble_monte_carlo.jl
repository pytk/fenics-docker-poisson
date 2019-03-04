include("./create_particle.jl")
include("./drift.jl")
include("./scat.jl")
include("./electron_charge.jl")

function ensembleMonteCarlo(particles, device, electron_density, current_time, electric_field, Gm)
    #=
    Args:
        - particle: Dict{Int32, Dict{String, Float64}}
        particles[particle_number][property]
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

        - density: Array{Float64, 2 x 2}

        - doner: the number of particle in n++ region

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

    # time means tdt in Archimedes
    time = current_time + device["time_step"]

    # initial particle number for accessing Julia dictionaly
    particle_number = 1

    while true
        # copy of all particle information
        println("current particle number is ", particle_number)
        temp_particle = particles[particle_number]
        temp_time = current_time
        while temp_particle["ts"] <= time
            τ = temp_particle["ts"] - temp_time

            # drift process from the point to other point that particle can move during the free flight time and update only particle position
            temp_particle = drift(τ, temp_particle, device, electric_field)

            # just take the x position and z position after drift process
            i = trunc(Int, temp_particle["x"]/dx)+1
            j = trunc(Int, temp_particle["z"]/dz)+1

            if i <= 1 i=1 end
            if j <= 1 j=1 end
            if i >= nx i=nx end
            if j >= nz j=nz end

            # scat process where a scat will get energy of particle changed based on the scattering rate on the particle's position and update only particle energy
            temp_particle = scat(temp_particle, device)
            temp_time = temp_particle["ts"]

            # set next happen scat time from random value and subband the particle belong to
            temp_particle["ts"] = temp_time - log(rand())/Gm[i][temp_particle["subband"]]

            # if the time scat process occuerd later over time step of simulation, it must be forced to move on to next particle process and next time step
        end

        τ = time - temp_time

        # TODO: add electric_field
        temp_particle = drift(τ, temp_particle, device, electric_field)

        if temp_particle["subband"] != 9
            particles[particle_number] = temp_particle
            particle_number += 1
        elseif temp_particle["subband"] == 9
            x = particle["x"]
            particles[particle_number] = particles[end]
            deleteat!(particle, length(particle))
            # the point where disappared particle is
            density_for_creation = doner - density[x][2]
            for j in 1:density_for_creation
                created_particle = createParticle(device, Gm, eigen_values, x)
                particles[particle_number] = created_particle
                particle_number += 1
            end
        end

        if particle_number < length(particles)
            break
        end
    end
    electron_density, each_subband_density = electronCharge(particles, device)
    return electron_density, each_subband_density, particles
end