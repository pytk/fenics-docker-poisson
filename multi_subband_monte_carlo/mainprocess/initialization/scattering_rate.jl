function getScatteringRate(wavefunction, eigen_values, device, constant)
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
    nz = device["nz"]
    subband_number = device["subband_number"]
    # energy step to calculate scattering rate from 0.0 eV to 1.0 eV per 0.001 eV
    de = 0.001
    ne = 1000

    # max energies
    energies = 1.00

    # Boltzman constant (JK^-1)
    kb = 8.6173303*10^-5

    # device temperture
    T = constant["T"]

    # effective mass in L band
    m = device["material"]["electron_effective_mass"]

    # Dirac's constant (Js)
    ħ = 6.582119514*10^-16

    # phonon energy
    ħω = device["material"]["phonon_energy"]

    # deformation energy
    D = device["material"]["deformation_potential"]

    cL = device["material"]["elasctic_constant"]

    scattering = ["acoustic", "non-polor-optics"]

    function density(energy, eigen_values, subband, x)
        #= 
        Args:
            - energy: eneryg step to compare with subband energy
            - eigen_values: [E2, E1] E2 means higher subband value, E1 means smaller subband energy
        =#
        if energy >= (eigen_values[subband+1][x] - eigen_values[subband][x])
            Nn = (2m^1.5)/(4*(π^2)*(ħ^3)) * (eigen_values[subband+1][x] -eigen_values[subband][x])
        else
            Nn = 0.0
        end
        return Nn
    end

    # calculate Gnm and φ which is required to calculate scattering rate for each kind of scattering that's why I compute these constant for each of first step of calculating scattering rate for each step
    # we have to make dictionary for each loop and get integrated finally
    #scattering_rate = Dict{Int32, Dict{Int32, Dict{Float64, Array{Float64, 1}}}}()
    # store.jl is a module to store global variable
    each_slice_scattering_rate = Dict{Int32, Dict{Int32, Dict{Int32, Array{Float64, 1}}}}()
    each_slice_Gm = Dict{Int32, Dict{Int32, Float64}}()
    for x in 1:nx+1
        each_subband = Dict{Int32, Dict{Int32, Array{Float64, 1}}}()
        each_subband_Gm = Dict{Int32, Float64}()
        each_subband_Gnm = Dict{Int32, Float64}()
        for subband in 1:subband_number

            # this process is independent but other scattering rate is based on calculated value in this area
            #######==========================================#########
            each_qz_Gnm = Array{Float64, 1}(undef, 0)
            φ = Array{Float64, 1}(undef, 0)
            qzmax = eigen_values[subband+1][x] - eigen_values[subband][x]
            temp = 0.00
            for qz in -qzmax:de:qzmax
                for z in 1:nz+1
                    temp = temp + (wavefunction[subband][z, x]*wavefunction[subband+1][z, x]) * exp(z*qz*10^-9)
                end
                push!(each_qz_Gnm, temp)
            end
            φ = sum(map(each_qz_Gnm) do t t^-2 end)

            each_subband_Gnm[subband] = φ

            # compute each type of scattering rate
            each_energy = Dict{Int32, Array{Float64, 1}}()

            # Array to normalize scattering rate by calculatin sum of final number of scat
            final_number_of_scat = Array{Float64, 1}(undef, 0)
            for energy in 0:de:energies
                # Inititalize scattering rate
                # Array to store each type of scattering rate
                each_scat = Array{Float64, 1}(undef, 0)
                scattring_rate_temp = 0.00

                # TODO: each kind of scattering rate should be calculated in each component as separete file
                # acoustic scattering 
                Nn = density(energy, eigen_values, subband, x)
                acoustic = 2*π*(D^2)*kb*T/(ħ*cL)*φ*Nn
                scattring_rate_temp = scattring_rate_temp + acoustic
                push!(each_scat, scattring_rate_temp)

                # Non-polor optical phonon absorb
                Nn = density(energy+ħω, eigen_values, subband, x)
                non_polor_absorb = 0
                scattring_rate_temp = scattring_rate_temp + non_polor_absorb
                push!(each_scat, scattring_rate_temp)

                # Non-polor optical phonon emition
                Nn =density(energy-ħω, eigen_values, subband, x)
                non_polor_emition = 0
                scattring_rate_temp = scattring_rate_temp + non_polor_emition
                push!(each_scat, scattring_rate_temp)

                # this is speacial only for final number of scat
                push!(final_number_of_scat, scattring_rate_temp)
                energy = trunc(Int, energy/de)+1
                each_energy[energy] = each_scat
            end

            # normalize scattering rate for each energy and each type of scat
            max_scattering_rate = maximum(final_number_of_scat)
            # this is for initialize each particle info in initialize.jl         
            each_subband_Gm[subband] = max_scattering_rate
            for energy in 0:de:energies
                energy = trunc(Int, energy/de)+1
                each_energy[energy] = map(each_energy[energy]) do el el/max_scattering_rate end
            end
            each_subband[subband] = each_energy
        end
        # scattering rate is a global variable
        each_slice_scattering_rate[x] = each_subband
        each_slice_Gm[x] = each_subband_Gm
    end
    global scattering_rate
    scattering_rate = each_slice_scattering_rate
    println("scattering rate done!!!!!!!!!!!!!!!!!!")
    return each_slice_Gm
end