function getScatteringRate(wavefunction, eigen_values, device)
    #=
    Return scattering rate about 2D electron gass by calculating wave function and eigen values with Schrodinger equation

    Args: 
        - wavefunction[subband, nx, nz] wave function from schrodinger equation for each cell

        - eigen_values[subband, nx] eigen valus from schrodinger equation

        - device device instance

    Return:
        - scattering_rate_first[{}]
        return scattering rate for each enrgy step.
        you can access specific value like scattering_rate[energy]['scat_name']['subband_number']
    =#

    nx = device["nx"]
    subband_number = device["subband_number"]
    # energy step to calculate scattering rate from 0.0 eV to 1.0 eV per 0.001 eV
    de = 0.001
    ne = 1000

    # Boltzman constant (JK^-1)
    kb = 1.38064852*10^-23

    # device temperture
    T = device["T"]

    # effective mass in L band
    m = device["electron_effective_mass"]

    # Dirac's constant (Js)
    ħ = 1.054571800*10^-34

    # phonon energy
    ħω = device["phonon_energy"]

    # deformation energy
    D = device["deformation_potential"]

    cL = device["elasctic_constant"]

    scattering = ["acoustic", "non-polor-optics"]

    function density(energy, eigen_values)
        #= 
        Args:
            - energy: eneryg step to compare with subband energy
            - eigen_values: [E2, E1] E2 means higher subband value, E1 means smaller subband energy
        =#
        if enrgy >= eigen_values[1]
            Nn = (2m^1.5)/(4(π^2)(ħ^3)) * (eigen_values[subband+1, x] -eigen_values[subband, x])
        else
            Nn = 0
        end
        return Nn
    end
    
    # scattering rate array
    scattering_rate_first = Array{Float64, 2}[undef, nx, ne]
    scattering_rate_second = Array{Float64, 2}[undef, nx, ne]
    scattering_rate_third = Array{Float64, 2}[undef, nx, ne]

    # calculate Gnm and φ which is required to calculate scattering rate for each kind of scattering that's why I compute these constant for each of first step of calculating scattering rate for each step
    for x in 1:nx+1
        # calculate Gn,m which is an integration of wave function along z axsis
        Gnm = Dict{Int32, Array{Float64, 1}}
        φ = Array{Float64, 1}
        for subband in 1:subband_number
            qzmax = eigen_values[subband+1, x] - eigen_values[subband, x]
            
            sum = 0
            for qz in -qzmax:de:qzmax
                for z in 1:nz+1
                    sum = sum + (wavefunction[subband][:, z]*wavefunction[subband+1][:, z]) * exp.((z*qz)im)
                end
                push!(Gnm, sum)
            end
            #Gnm = map(Gnm) do x x^2 end
            φ = map(Gnm) do t t^2 end

            # φ has number of subband element
            # φ include all pairs (1, 2), (2, 3), (3, 4)
            if Gnm in φ
                continue
            else
            # φ depends on only energy itself so that we don't add existing value
                push!(φ, Gnm)
            end

            # compute each type of scattering rate
            acoustics = []
            non_polor_absorbs = []
            non_polor_emitions = []
            scattring_rate_temp = 0
            for e in 0:ne
                
                # acoustic scattering 
                Nn = density(e)
                acoustic = 2*π*(D^2)*kb*T/(ħ*cL)*φ*Nn
                push!(acoustics, acoustic)

                scattring_rate_temp += acoustic

                # non-polor optical phonon absorb
                Nn = density(e+ħω)
                non_polor_absorb = 0
                push!(non_polor_absorbs, non_polor_absorb)

                scattring_rate_temp += non_polor_absorb

                # non-polor optical phonon emition
                Nn =density(e-ħω)
                non_polor_emition = 0
                push!(non_polor_emitions, non_polor_emition)

                scattring_rate_temp += non_polor_emition

                

            # now that we have only acoustic inttra valley phonon scattering
            # you have to set scattering that you set finally!!!!!
            max_value = maximum(acoustics)

            # normalize scattering rate
            temp = map(scattring_rate_temp) do t t/max_value end

            # add normalized scattering rate for each slic 
            if subband == 1
                push!(scattering_rate_first, temp)
            else if subband == 2
                push!(scattering_rate_second, temp)
            else
                push!(scattering_rate_third, temp)
            end
        end
    end

    return scattering_rate_first, scattering_rate_second, scattering_rate_third
end             
                    