function electronCharge(particles, device)
    #=
    To return electron charge based on cell cloud method

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
        - electron_density[nx, nz] = Float64
        - each_subband_density[subband_number, nx, nz] = Float64
    =#

    dx = device["dx"]
    dz = device["dz"]

    nx = device["nx"]
    nz = device["nz"]

    subband_number = device["subband_number"]
    
    each_subband_density = Array{Float64}(undef, subband_number, nx, ny)
    each_subband_density = zero(each_subband_density)

    electron_density = Array{Float64}(undef, nx, ny)
    electron_density = zero(electron_density)

    carrier_per_superparticle = device["carrier_per_superparticle"]

    for particle in particles
        x = particle["x"]
        z = particle["z"]
        i = trunc(Int, x/dx)+1
        j = trunc(Int, z/dz)+1
        subband = particle["subband"]

        cloud_density = (1.0-(x-(i-1.0)))*(1.0-(y-(j-1.0)))

        # clould in cell method
        if 0 < subband && subband < subband_number
            each_subband_density[subband, i, j] += cloud_density
        end
        electron_density[i, j] += cloud_density

        if i <= nx
            if 0 < subband && subband < subband_number
                each_subband_density[subband, i+1, j] += cloud_density
            end
            electron_density[i+1, j] += cloud_density
        end
        if i <= ny
            if 0 < subband && subband < subband_number
                each_subband_density[subband, i, j+1] += cloud_density
            end
            electron_density[i, j+1] += cloud_density
        end
        if i <= nx && j <= ny
            if 0 < subband && subband < subband_number
                each_subband_density[subband, i+1, j+1] += cloud_density
            end
            electron_density[i+1, j+1] += cloud_density
        end
    end

    temp = map(e -> e*(carrier_per_superparticle/(dx*dz)), electron_density)
    electron_density = temp
    temp = map(e -> e*(carrier_per_superparticle/(dx*dz)), each_subband_density)
    each_subband_density = temp
    
    return electron_density, each_subband_density
end