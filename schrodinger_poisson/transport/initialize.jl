function initializeParticle(electron_density, device, Gm, eigen_values)
    #=
    Return particle information like wave number (kx, ky)
    position, subband number, etc ... for each particle.
    as an args, we can get potential and eigen value from schrodinger-poisson loop process, so we have to assign each particle to specific subband firstly
    when it comes to kenetic energy, we assign random value to each particle

    Args:
        - electron_density[x, z] electron density distribution from self-consistent schrodinger-poisson loop

        - eigen_values [subband, nx] eigen value from schrodinger equation for each slice along x-axsis

        - p[particle number][i]
            - i = 1: subband 1 ~ 3 but 8 means 3d particle and 9 means the particle in ohmic or schottky layer so it would be deleted
            - i = 2: kx
            - i = 3: ky
            - i = 4: ts (free flight time)
            - i = 5: x (x position)
            - i = 6: z (z position)
        - scattering_rate

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
    particle = Dict{Int32, Dict{String, Float64}}()

    # eigen values
    eigen_value = eigen_values[1]

    # dx , dz
    dx = device["dx"]
    dz = device["dz"]

    nx = device["nx"]
    nz = device["nz"]

    gate_ini = device["gate_ini"]
    gate_fin = device["gate_fin"]

    number_of_carrier_per_superparticle = device["carrier_per_superparticle"]


    particle_count = 0
    for i in 1:nx
        for j in 1:nz
            number_of_particle = electron_density[i][j]*dx*dz/number_of_carrier_per_superparticle
            for np in 1:number_of_particle
                particle_number += 1

                # we assume that the particles are initially at near thermal equilibrium
                Ek = -kb*T*rand()

                # wave bector (E-k) without non-parabolicity
                k = √(2m*Ek)/ħ

                ϕ = 2π*rand()

                cosθ = 1 - 2rand()
                sinθ = √(1 - cosθ^2)
                cosϕ = cos(ϕ)
                sinϕ = sin(ϕ)

                # kz from eigen value by schrodinger equation
                kx = k*sinθ*cosϕ
                ky = k*sinθ*sinϕ
                z = dz*(rand()+j-1.5)
                ts = -log(rand())/Gm[trunc(Int, x/dx)]
                x = dx*(rand()+i-1.5)

                # When it comes to z-axis in wave space, it depends on the position of particle. if a particle is within under the gate, the particle would be 2d electron
                if gate_ini <= x && x <= gate_fin
                    kz = √(2m*eigen_value[i])/ħ
                    subband = 1
                else
                    kz = k*cosθ
                    subband = 8

                # build the new particle information as julia dictionaly
                dict = Dict("subband" => subband, "kx" => kx, "ky" => ky, "kz" => kz, "ts" => ts, "x" => x, "z" => z)

                # insert No.n particle dict into particle array
                particle_count += 1
                particle[particle_count] = dict
            end
        end
    end

    @printf "Initial number of electron super-particles is %d" particle_count
    return particle

end