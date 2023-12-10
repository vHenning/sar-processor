using FFTW

function parseMetadata(pathname, metadataName, rangeCells)
    file = open("$pathname/$metadataName.0__A")
    rec = parseFile(file,datasetSummaryRecordScheme,720)
    close(file)
    
    PRF = parse(Float64,rec.fields[rec.key["PRF"]])/1000
    
    sampleRate = parse(Float64,rec.fields[rec.key["samplingRate"]])*1e6     #Hz
    pulseSamples = let
        pulseLength = parse(Float64,rec.fields[rec.key["pulseLength"]])*1e-6   #s
        Integer(floor(pulseLength*sampleRate))
    end
    
    #process header
    chirpFFT = let
        fdot = -parse(Float64,split(rec.fields[rec.key["coeffs"]])[2])
        sampleRate = parse(Float64,rec.fields[rec.key["samplingRate"]])*1e6     #Hz
        pulseLength = parse(Float64,rec.fields[rec.key["pulseLength"]])*1e-6   #s
        pulseSamples = Integer(floor(pulseLength*sampleRate))
    
        Sif(t) = exp(pi*im*fdot*t^2)
        t = 1/sampleRate* (range(1, stop = pulseSamples) |> collect)
        t = t.-maximum(t)/2
        sig = Sif.(t)/sqrt(pulseSamples)
        chirp = vcat(sig,zeros(Complex{Float32}, rangeCells))
    
        fft(chirp)
    end
    
    wavelength = parse(Float32,rec.fields[rec.key[  "wavelength" ]]) #m 

    return (chirpFFT, pulseSamples, sampleRate, PRF, wavelength);
end

function readData(pathname, imagename, rangeCells)
    file = open("$pathname/$imagename.0__A")

    imageDescriptor = parseFile(file,imageFileDescrictorScheme)
    sarDataBytes    = parse(Int64,imageDescriptor.fields[imageDescriptor.key["sarDataBytes"]])
    numSignals      = parse(Int64,imageDescriptor.fields[imageDescriptor.key["numSignals"]])
    
    signalRecords = [ ]
    
    for i = 1:(numSignals-1)
        parsedRecord = parseFile(file,signalDataRecordScheme,720+sarDataBytes*i)
        signal = parsedRecord.fields[parsedRecord.key["signalData"]]
        
        #TODO  - is 15 or 16 better?? compute mean of signals?
        #TODO  - DEAL WITH VALUES > 0x1f!!
        signal = map(x->min(0x1f,x), signal)
        
        signal = Int8.(signal)-Int8.(15*ones(size(signal)))
        I = signal[1:2:length(signal)]
        Q = signal[2:2:length(signal)]
    
        parsedRecord.key["I"] = length(parsedRecord.fields)+1;
        parsedRecord.key["Q"] = length(parsedRecord.fields)+2;
        push!(parsedRecord.fields, I)
        push!(parsedRecord.fields, Q)
        
        parsedRecord.fields[parsedRecord.key["signalData"]] = []  #remove original signalData from record!
        
        push!(signalRecords,parsedRecord)
    end
    
    close(file)
    
    #sub image formation:
    sampleNum = rangeCells #how many samples of each echo to keep - keep all range cells by default
    echoNum = 35000  #how many echos to keep
    echostart = 1

    # smallSignals is matrix containing a subset (defined by sampleNum and echoNum)
    # of the echo signals in signalRecords
    smallSignals = zeros(Complex{Float16},sampleNum,echoNum)

    #This is soooooo much faster than hcatS!!!!
    #https://stackoverflow.com/questions/38308053/julia-how-to-fill-a-matrix-row-by-row-in-julia
    for i = echostart+1:echostart+echoNum
        line = signalRecords[i]
        I = Float16.(line.fields[line.key["I"]])
        Q = Float16.(line.fields[line.key["Q"]])
        
        smallSignals[:,i-echostart] = Complex.(I[1:sampleNum],Q[1:sampleNum])
    end

    signalRecords = []  #free up signalRecords - now using smallSignals
    GC.gc();

    return smallSignals;
end