function drift(τ, particle, device, electron_field)
    #=
    Note: This doesn't include non-prabolicityt

    Args:
        - τ: Float64 free flight time
        - particle: Object
        - device: Object
        - electric_filed: python dictionaly Dict{Int32, Array{Float64, 1}}
    =#
    # effective mass in L band
    m = device["electron_effective_mass"]

    # Dirac's constant (Js)
    ħ = 1.054571800*10^-34

    e = 1.6021766208*10^-19

    F = electric_filed

    # dx , dz
    dx = device["dx"]
    dz = device["dz"]

    nx = device["nx"]
    nz = device["nz"]

    i = trunc(Int, particle["x"]/dx)+1
    j = trunc(Int, particle["z"]/dz)+1

    if i <= 1 i=1 end
    if j <= 1 j=1 end
    if i >= nx i=nx end
    if j >= nz j=nz end

    # electron drift process
    dkx =  e*F[j]*τ/ħ

    kx = particle["kx"]
    ky = particle["ky"]
    kz = particle["kz"]

    x = particle["x"]
    z = particle["z"]

    # calculate group velocity along with x-axis. parhaps, It must be considered z-axis
    particle["x"] = x + ħ/m*τ*(kx+0.5*dkx)
    particle["z"] = z + ħ*kz/m*τ

    particle["kx"] = kx+dkx

    return particle
end


