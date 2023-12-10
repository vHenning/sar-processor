function deconvolve(smallSignals, chirpFFT, pulseSamples)
    shape = size(smallSignals)

    # add zero padding at the beginning of each pulse echo (each column is an echo)
    cimg = vcat(zeros(Complex{Float32},(pulseSamples,shape[2])),
                Complex{Float32}.(smallSignals));

    fft!(cimg,(1)); # perform an FFT on each column (each pulse echo)

    cimg =  cimg .* conj.(chirpFFT) ; # convolution with chirp signal performed in frequency domain

    ifft!(cimg,(1)); # perform an inverse FFT on each column (each pulse echo)
    cimg = Complex{Float16}.(cimg');  #transpose so that each echo is a horizontal line
    smallSignals = [] # free up memory of smallSignals

    return cimg;
end

function rangeMigration(cimg, c, sampleRate, PRF, wavelength, R0)
    # run an fft on each column of cimg (echos are rows here)
    cimg = Complex{Float16}.(fft(Complex{Float32}.(cimg),(1)));

    #now we want to shift each frequency in range space, so make each frequency bin a column for speed
    cimg = cimg'
    shape = size(cimg)

    slantRes = 1/2*c/sampleRate

    for i = 1:shape[2]
        n = i
        #the "highest frequencies" are actually the negative frequencies aliased up!
        if n>shape[2]/2
            n = shape[2]-i
        end
        fn = (n-1)/shape[2]*PRF    #check this but pretty sure
        
        #range migration distance in meters
        ΔR = wavelength^2*R0*fn^2/(8*Vr^2)
        cellshift = Integer(round(ΔR/slantRes))
        
        #interpolation
        #NEAREST NEIGHBOR - bad!
        cimg[:,i] = vcat(cimg[cellshift+1:shape[1],i],zeros(Complex{Float64},cellshift))
    end
    return cimg'
end