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
    scattering_rate = Array{Float64, 4}[0.0, nx, subband_number, ne, scattering_number]

    # calculate Gnm and φ which is required to calculate scattering rate for each kind of scattering that's why I compute these constant for each of first step of calculating scattering rate for each step

    # we have to make dictionary for each loop and get integrated finally
    each_slice = Dict{Int32, Dict{Int32, Dict{Int32, Array{Float64, scattering_number}}}}()
    for x in 1:nx+1
        each_subband = Dict{Int32, Dict{Int32, Array{Float64, scattering_number}}}
        for subband in 1:subband_number

            # this process is independent but other scattering rate is based on calculated value in this area
            #######==========================================#########
            Gnm = Dict{Int32, Array{Float64, 1}}
            φ = Array{Float64, 1}
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
            #=========================================#

            # compute each type of scattering rate
            each_energy = Dict{Int32, Array{Float64, scattering_number}}

            # Array to normalize scattering rate by calculatin sum of final number of scat
            final_number_of_scat = Array{Float64, 1}

            for energy in 1:energies
                # Inititalize scattering rate
                # Array to store each type of scattering rate
                each_scat = Array{Float64, scattering_number}
                scattring_rate_temp = 0

                # TODO: each kind of scattering rate should be calculated in each component as separete file
                # acoustic scattering 
                Nn = density(e)
                acoustic = 2*π*(D^2)*kb*T/(ħ*cL)*φ*Nn
                scattring_rate_temp += acoustic
                push!(each_scat, scattring_rate_temp)

                # Non-polor optical phonon absorb
                Nn = density(e+ħω)
                non_polor_absorb = 0
                scattring_rate_temp += non_polor_absorb
                push!(each_scat, scattring_rate_temp)

                # Non-polor optical phonon emition
                Nn =density(e-ħω)
                non_polor_emition = 0
                scattring_rate_temp += non_polor_emition
                push!(each_scat, scattring_rate_temp)

                # this is speacial only for final number of scat
                push!(final_number_of_scat, scattring_rate_temp)

                #===========================================#

                each_energy[energy] = each_scat
            end

            # normalize scattering rate for each energy and each type of scat
            max_scattering_rate = maximum(final_number_of_scat)
            for energy in 1:energies
                temp = map(each_energy[energy]) do el el/max_scattering_rate end
                each_energy[energy] = temp
            #=================================================================#

            each_subband[subband] = each_energy
        end

        each_slice[x] = each_subband
    end

    return each_slice
end


