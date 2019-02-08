function initializeParticle(electron_density, device, scattering_rate):
    #=
    Return particle information like wave number (kx, ky)
    position, subband number, etc ... for each particle.
    as an args, we can get potential and eigen value from schrodinger-poisson loop process, so we have to assign each particle to specific subband firstly
    when it comes to kenetic energy, we assign random value to each particle

    Args:
        - electron_density[x, z] electron density distribution from self-consistent schrodinger-poisson loop

        - p[particle number][i]
            - i = 1: subband number
            - i = 2: kx
            - i = 3: ky
            - i = 4: ts (free flight time)
            - i = 5: x (x position)
            - i = 6: z (z position)

    Note: we are only taking consider 2d particle gas not 3d particle,
          hence, we don't need to consider kz because the energy along the z-axis is discretized
    =#
    
    # constant value

    # Boltzman constant (JK^-1)
    kb = 1.38064852*10^-23

    # device temperture
    T = device["T"]

    # effective mass in L band
    m = device["electron_effective_mass"]

    # Dirac's constant (Js)
    ħ = 1.054571800*10^-34

    particle_number = 0

    # particle array 
    particle = Any[]

    # dx , dz
    dx = device["dx"]
    dz = device["dz"]

    for i in 1:device["nx"]
        for j in 1:device["nz"]
            number_of_particle = electron_density[i, j]
            for np in 1:number_of_particle
                particle_number += 1

                # we assume that the particles are initially at near thermal equilibrium
                Ek = -kb*T*rand()

                # wave bector (E-k) without non-parabolicity
                k = sqrt(2m*Ek)/ħ

                ϕ = 2π*rand()

                cosθ = 1 - 2rand()
                sinθ = sqrt(1 - cosθ^2)
                cosϕ = cos(ϕ)
                sinϕ = sin(ϕ)

                # build the new particle information as julia dictionaly
                dict = Dict("subband" => 1, "kx" => k*sinθ*cosϕ, "ky" => k*sinθ*sinϕ, "ts" => "TODO", "x" => dx*(rand()+i-1.5), "z" => dz*(rand()+j-1.5))

                # append No.n particle dict into particle array
                push!(particle, dict, nothing)
            end
        end
    end

    @printf "Initial number of electron super-particles is %d" particle_number
    return particle

end
                
            
