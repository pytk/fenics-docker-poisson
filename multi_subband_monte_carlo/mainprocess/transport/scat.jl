function scat(particle, device, scattering_rate)
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

    kx = particle["kx"]
    ky = particle["ky"]
    kz = particle["kz"]

    x = trunc(Int, particle["x"]/dx)+1
    subband = particle["subband"]

    de = device["energy_step"]
    ne = device["number_of_energy"]

    number_of_scat = device["number_of_scat"]

    # effective mass in L band
    m = device["electron_effective_mass"]

    phonon_energy = device["material"]["phonon_energy"]

    # Dirac's constant (Js)
    ħ = 1.054571800*10^-34

    k = kx^2 + ky^2 + kz^2
    superparticle_energy = (ħ^2)*k/(2m)

    if superparticle_energy <= 0.0 return end
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

    return particle

end




    
    
