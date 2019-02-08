function getScatteringRate(wavefunction, eigen_values, device)
    #=
    Return scattering rate about 2D electron gass by calculating wave function and eigen values with Schrodinger equation

    Args: 
        - wavefunction[subband, nx, nz] wave function from schrodinger equation for each cell

        - eigen_values[subband, nx] eigen valus from schrodinger equation

        - device device instance

    Return:
        - scattering_rate[ {'scat name': {'grand': 4398, 'second': 4387, 'third': 4387343}},  {'scat name': {'grand': 4398, 'second': 4387, 'third': 4387343}}, .....]
        return scattering rate for each enrgy step.
        you can access specific value like scattering_rate[energy]['scat_name']['subband_number']
    =#

    nx = device["nx"]
    subband_number = device["subband_number"]
    # energy step to calculate scattering rate from 0.0 eV to 1.0 eV per 0.001 eV
    de = 0.001
    ne = 1000
    
    # scattering rate array
    scattering_rate = Array{Float64, 2}[undef, nx, ne]

    # calculate Gnm and φ which is required to calculate scattering rate for each kind of scattering that's why I compute these constant for each of first step of calculating scattering rate for each step
    for e in 1:ne+1
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
                    Gnm = map(Gnm) do x x^2 end
                end
                # φ has number of subband element
                push!(φ, Gnm)
            end

            ##### here to write next step for each slice and store the calculated scattering rate into array



        end






    end

                
                     
                
                
                


    


        


