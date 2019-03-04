include("./check_if_particle_out.jl")

function drift(τ, particle, device, electric_field)
    #=
    Note: This doesn't include non-prabolicityt

    Args:
        - τ: Float64 free flight time
        - particle: Object
        - device: Object
        - electric_filed: python dictionaly Dict{Int32, Array{Float64, 1}}
    =#
    println("drift process get started")

    # effective mass in L band
    m = device["material"]["electron_effective_mass"]

    # Dirac's constant (Js)
    ħ = 6.582119514*10^-16

    e = 1.6021766208*10^-19

    F = electric_field

    # dx , dz
    dx = device["dx"]
    dz = device["dz"]

    nx = device["nx"]
    nz = device["nz"]

    kx = particle["kx"]
    ky = particle["ky"]
    kz = particle["kz"]

    x = particle["x"]
    z = particle["z"]

    i = trunc(Int, x/dx)+1
    j = trunc(Int, z/dz)+1

    if i <= 1 i=1 end
    if j <= 1 j=1 end
    if i >= nx i=nx end
    if j >= nz j=nz end
    # electron drift process
    dkx =  e*F[j][i]*τ/ħ

    # calculate group velocity along with x-axis. parhaps, It must be considered z-axis
    particle["x"] = x + ħ/m*τ*(kx+0.5*dkx)
    particle["z"] = z + ħ*kz/m*τ

    particle["kx"] = kx+dkx
    # check if particle goes out from device range
    particle = checkIfParticleOut(particle, device)
    println("drift process done")
    return particle
end


