function checkIfParticleOut(particle, device)
    #=
    Args:
        - particle: Dict{String, Float64}
        particle[property]
            - property:
                -subband: Int32 the subband the particle belong to
                -kx: Float64 wave vector along x axis in wave space
                -ky: Float64 wave vctor along y axis in wave space
                -kz: Float64 wave vctor along z axis in wave space
                -x: Float64 
                -z: Float64
                -ts: Float64 time which couse scattering rate

            - scattering_rate: Dict{Int32, Dict{Int32, Dict{Int32, Array{Float64, scattering_number}}}}
                - scattering_rate[nx][subband][energy]
    =#

    # dx , dz
    dx = device["dx"]
    dz = device["dz"]

    nx = device["nx"]
    nz = device["nz"]

    kx = particle["kx"]
    kz = particle["kz"]

    xfi = device["xfi"]
    zfi = device["zfi"]

    subband = particle["subband"]

    src = device["src"]
    drain = device["drain"]

    gate_ini = device["gate_ini"]
    gate_fin = device["gate_fin"]

    x = particle["x"]
    z = particle["z"]

    # check if a particle is going out from the right or left edge of the device
    if x < 0
        x = -x
        kx = -kx
    elseif x > xfi
        x = xfi - (x - xfi)
        kx = -kx
    end

    if z < 0
        if x < src || x > doner
            subband = 9
        elseif x > gate_ini && x < gate_fin
            subband = 9
        else
            z = -z
            kz = -kz
        end
    elseif z > nz
        z = zfi - (z-zfi)
        kz = -kz
    end

    particle["x"] = x
    particle["z"] = z
    particle["kx"] = kx
    particle["kz"] = kz
    particle["subband"] = subband

    return particle
end


