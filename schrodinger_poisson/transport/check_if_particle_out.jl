function checkIfParticle(particle, device)
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

    src = device["src"]
    drain = device["drain"]

    gate_ini = device["gate_ini"]
    gate_fin = device["gate_fin"]

    x = trunc(Int, particle["x"]/dx + 1.5)
    z = trunc(Int, particle["z"]/dz + 1.5)

    # check if a particle is going out from the right or left edge of the device
    if x < 0
        x = -x
        kx = -kx
    else if x > nx+1
        x = nx - (x - nx)
        kx = -kx
    end

    if j < 0
        if x < src || x > doner
            subband = 9
        else if x > gate_ini && x < gate_fin
            subband = 9
        else
            z = -z
            kz = -kz
        end
    else if z > nz
        z = nz - (z-nz)
        kz = -kz
    end

    particle["x"] = x
    particle["z"] = z
    particle["kx"] = kx
    particle["kz"] = kz
    particle["subband"] = subband

    return particle
end


