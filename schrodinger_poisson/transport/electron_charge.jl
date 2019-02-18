function electronCharge(particles, device)
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

    Return:
        - electron_density[subband][x][z] = Float64
    =#

    dx = device["dx"]
    dz = device["dz"]

    nx = device["nx"]
    nz = device["nz"]

    electron_density = Dict{Int32,Dict{Int32, Dict{Int32, Float64}}}

    for particle in particles
        x = particle["x"]
        z = particle["z"]
        i = trunc(Int, x/dx)+1
        j = trunc(Int, z/dz)+1
        subband = particle["subband"]

        cloud_density = (1.0-(x-(i-1.0)))*(1.0-(y-(j-1.0)))

        # clould in cell method
        electron_density[i][j][subband] += cloud_density
        if i <= nx
            electron_density[i+1][j][subband] += cloud_density
        end
        if i <= ny
            electron_density[i][j+1][subband] += cloud_density
        end
        if i <= nx && j <= ny
            electron_density[i+1][j+1][subband] += cloud_density
        end
    end

    for x in 1:nx
        for z in 1:nz
            electron_density[x][z] *= number_of_particle_in_superparticle / (dx*dz)
        end
    end
end


