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
