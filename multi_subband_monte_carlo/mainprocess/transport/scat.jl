include("./../initialization/scattering_rate.jl")
function scat(particle, device)
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

            - scattering_rate: Dict{Int32, Dict{Int32, Dict{Int32, Array{Float64, scattering_number}}}}
                - scattering_rate[nx][subband][energy]
    =#
    println("scat process get started")

    kx = particle["kx"]
    ky = particle["ky"]
    kz = particle["kz"]

    dx = device["dx"]
    dz = device["dz"]

    x = particle["x"]
    x = trunc(Int, x/dx)+1

    subband = particle["subband"]

    # energy step to calculate scattering rate from 0.0 eV to 1.0 eV per 0.001 eV
    # WARNING: if you wanna change this value you have to change same stuff also in scattering_rate process
    de = 0.001
    ne = 1000

    # effective mass in L band
    m = device["material"]["electron_effective_mass"]

    phonon_energy = device["material"]["phonon_energy"]

    # Dirac's constant (Js)
    ħ = 1.054571800 * 10^-15

    k = kx^2 + ky^2 + kz^2
    superparticle_energy = (ħ^2)*k/(2m)

    if superparticle_energy <= 0.0 return particle end
    println(superparticle_energy)
    ie = trunc(Int, superparticle_energy/de)+1

    if ie > ne ie = ne end

    # select scattering process
    #==================================#
    r = rand()

    # non-polor optical phonon scattering
    # acoustic
    if r <= scattering_rate[x][subband][ie][1]
        final_energy = superparticle_energy
        if final_energy <= 0.0 return end
    
    # absourb
    elseif r <= scattering_rate[x][subband][ie][2]
        final_energy = superparticle_energy + phonon_energy
    
    # emission
    elseif r <= scattering_rate[x][subband][ie][3]
        final_energy = superparticle_energy - phonon_energy
    else
        println("x: ", x)
        println("subband: ", subband)
        println("ie: ", ie)
        println("r: ", r)
        println(scattering_rate[x][subband][ie][3])
    end
    #==================================#

    # determine final wave vector based on the calculated energy which got changed with scattering process
    Ek = final_energy
    # wave bector (E-k) without non-parabolicity
    k = √(2m*Ek)/ħ

    ϕ = 2π*rand()

    cosθ = 1 - 2rand()
    sinθ = √(1 - cosθ^2)
    cosϕ = cos(ϕ)
    sinϕ = sin(ϕ)

    kz = k*sinθ*sinϕ
    kx = k*cosϕ
    ky = k*sinθ*cosϕ

    particle["kx"] = kx
    particle["ky"] = ky
    particle["kz"] = kz
    println("scattering process done")
    return particle

end




    
    
