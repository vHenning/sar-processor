using ImageView
using Serialization
using FFTW

pathname = "../data/"    # relative path of folder containing L1.0 data
imagename = "IMG-HH-ALPSRP273530900-H1"  # which image polarization to use - IMG-HH or IMG-HV for PALSAR

function bytesToString(file, start, stop)
    len = stop-start+1
    seek(file,start-1)
    String(read(file,len))
end

function bytesToUInt32(file, start, stop)
    #TODO check if range matches size
    seek(file, start-1)
    num = try
        read(file,UInt32)
    catch
        0
    end
    Int32(ntoh(num))
end

function bytesToUInt16(file, start, stop)
    #TODO check if range matches size
    seek(file, start-1)
    Int16(ntoh(read(file,UInt16)))
end

function bytesToArray(file, start, stop)
    len = stop-start+1
    seek(file, start-1)
    read(file,len)
end

#Image File Descriptor
imageFileDescrictorScheme = [
        ("length", (9, 12), bytesToUInt32),
        ("numSignals", (181, 186), bytesToString),
        ("sarDataBytes", (187, 192), bytesToString),
        ("bitsPerSample", (217, 220), bytesToString),
        ("samplesPerPixel", (221, 224), bytesToString), #could be 1 or 2 (IQ or A)
        #245-272 - borders and interleaving
        ("physicalRecordsPerLine", (273, 274), bytesToString),
        #281-400 - unexplored
        ("sarDataFormat", (401, 428), bytesToString),
        ("sarDataFormatTypeCode", (429,432), bytesToString)
        
        ]

signalDataRecordScheme = [
        ("recordSequenceNumber", (1,4), bytesToUInt32),
        ("lengthOfRecord", (9,12), bytesToUInt32),
        ("sarImageDatalineNumber", (13, 16), bytesToUInt32),
        ("pixelCount", (25, 28), bytesToUInt32),
        #33-56 - instrument settings
        #57-84 - chirp characteristics
        ("PRF (mHz)", (57,60), bytesToUInt32),
        ("Chirp type", (67,68), bytesToUInt16),    #0 for linear
        
        #93-388 - platform info - important!!!
        ("Valid", (97,100),bytesToUInt32),
        ("Electronic Squint 1",(101,104),bytesToUInt32),
        ("Electronic Squint 2",(109,112),bytesToUInt32),
        ("Mechanical Squint",(113,116),bytesToUInt32),
        ("FirstSampleRange (m)",(117,120),bytesToUInt32),
        ("Platform Altitude (m)",(141,144),bytesToUInt32),
        ("Platform Ground Speed (cm/s)",(145,148),bytesToUInt32),
        ("signalData", (413, 10800), bytesToArray)
        ]

datasetSummaryRecordScheme = [
    ("wavelength", (501,516), bytesToString),
    ("chirpType", (519,534), bytesToString),
    ("coeffs",(534,694),bytesToString),
    ("samplingRate",(711,726), bytesToString),
    ("pulseLength",(743,758),bytesToString),      #us
    ("PRF", (935,950), bytesToString),            #mHz
    ("doppler", (1415,1654), bytesToString)
]

function parseFile(file,scheme,offset=0)
    keyMap = Dict([])
    for i=1:length(scheme)
        keyMap[scheme[i][1]] = i  #make dictionary mapping field names to indices
    end
    results = []
    for (name, (start, stop), f) in scheme
        len  = stop-start+1
        push!(results, f(file,start+offset,stop+offset))
    end
    return (key=keyMap, fields=results)
end

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
    if(i%10000==0)
        print(i)
        print("   ")
        print(Base.summarysize(parsedRecord))
        print("B    ")
        print(Base.summarysize(signalRecords)/1e6)
        println(" MB")
    end
end
print(Base.summarysize(signalRecords)/1e6)
println(" MB")

close(file)

# save signalRecords for later use
Serialization.serialize(open("$pathname/$imagename-signal-records.ser","w"),signalRecords)
signalRecords = [];  # free up the memory of the signalRecords object - it'll be read back in later

file = open("$pathname/LED-ALPSRP273530900-H1.0__A")
rec = parseFile(file,datasetSummaryRecordScheme,720)
close(file)

PRF = parse(Float64,rec.fields[rec.key["PRF"]])/1000
rangeCells = 5000


sampleRate = parse(Float64,rec.fields[rec.key["samplingRate"]])*1e6     #Hz
pulseSamples = let
    pulseLength = parse(Float64,rec.fields[rec.key["pulseLength"]])*1e-6   #s
    #PRF = parse(Float64,rec.fields[rec.key["PRF"]])/1000
    Integer(floor(pulseLength*sampleRate))
end

#process header
chirpFFT = let
    f = parse(Float64,split(rec.fields[rec.key["coeffs"]])[1])
    fdot = -parse(Float64,split(rec.fields[rec.key["coeffs"]])[2])
    sampleRate = parse(Float64,rec.fields[rec.key["samplingRate"]])*1e6     #Hz
    pulseLength = parse(Float64,rec.fields[rec.key["pulseLength"]])*1e-6   #s
    #PRF = parse(Float64,rec.fields[rec.key["PRF"]])/1000
    pulseSamples = Integer(floor(pulseLength*sampleRate))
    #pulseSamples = 400
    print("Max IF freq: ")
    println(fdot*pulseLength)

    Sif(t) = exp(pi*im*fdot*t^2)

    t = 1/sampleRate* (range(1, stop = pulseSamples) |> collect)

    t = t.-maximum(t)/2

    sig = Sif.(t)/sqrt(pulseSamples)

    chirp = vcat(sig,zeros(Complex{Float32}, rangeCells))

    fft(chirp)
end

print("Ready")

R0 = 848665     #m  
altitude = 628000 #m, nominal
c = 3e8

vorbital = let 
    a = (6371+628)*1000.0  #m            #https://www.eorc.jaxa.jp/ALOS-2/en/about/overview.htm
    G = 6.67408e-11
    Me= 5.97219e24       #kg
    P = sqrt(4*pi^2/(G*Me)*a^3)
    print("Period: ")
    print(round(P/60,digits=4))
    println(" m")
    2pi*a/P            #These are terrible assumptions!!!!!!!
end
vorbital = 7593
Vr = vorbital  # TODO: remove reference to redundant var Vr
println("Orbital Velocity: ",vorbital," m/s")

wavelength = parse(Float32,rec.fields[rec.key[  "wavelength" ]]) #m 
#antenna length
La = 8.9 #m

#beamwidth
bw = 0.886*wavelength/La

nadirAngle = acos(altitude/R0)
println(nadirAngle*180/pi)
RD = 1/2*c*rangeCells/16000000   # TODO where does 16,000,000 come from? 16 MHz sampling rate??
rangeCellLength = 1/rangeCells * (sqrt((R0+RD)^2-altitude^2)-sqrt(R0^2-altitude^2))

# compare the width and height of a pixel (cell)
# to get an approximate scaling for image formation
println("Azimuth:     ", round(vorbital/PRF, digits=3), " m")
println("Range Cell: ",round(rangeCellLength, digits=3), " m")

# approximate "y/x" scaling needed to render images without stretching
aspectRatio = rangeCellLength/(vorbital/PRF)
println("Aspect Ratio: ",round(aspectRatio,digits=1))

# load serialized version of signalRecords if not loaded from file above
signalRecords = Serialization.deserialize(open("$pathname/$imagename-signal-records.ser","r"))
print("Loaded signalRecords")

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

print("Done")

# TODO: serialize smallSignals?
#Serialization.serialize(open("$pathname/smallSignalsBigger.ser","w"),smallSignals)

shape = size(smallSignals)
rawMagnitude = abs.(view(smallSignals,1:10:shape[1],1:40:shape[2]));
rawMagnitude = reverse(rawMagnitude,dims=1)
imshow(rawMagnitude);
# TODO: might want to extend the width of the image by pulseSamples before downsizing,
# as was done in the original version. See snippet below:
# vcat(zeros(Complex{Float16},pulseSamples),smallSignals[:,i])

# deconvolution by the chirp signal

shape = size(smallSignals)

# add zero padding at the beginning of each pulse echo (each column is an echo)
cimg = vcat(zeros(Complex{Float32},(pulseSamples,shape[2])),
             Complex{Float32}.(smallSignals));

fft!(cimg,(1)); # perform an FFT on each column (each pulse echo)

cimg =  cimg .* conj.(chirpFFT) ; # convolution with chirp signal performed in frequency domain

ifft!(cimg,(1)); # perform an inverse FFT on each column (each pulse echo)
cimg = Complex{Float16}.(cimg');  #transpose so that each echo is a horizontal line
smallSignals = [] # free up memory of smallSignals

shape = size(cimg)
rangeCompressedMagnitude = abs.(view(cimg,1:40:shape[1],1:10:shape[2]));
rangeCompressedMagnitude = reverse(rangeCompressedMagnitude,dims=1)
imshow(rangeCompressedMagnitude);

Serialization.serialize(open("$pathname/$imagename.rcc","w"),cimg)
cimg = [];
GC.gc();

cimg = Serialization.deserialize(open("$pathname/$imagename.rcc","r"))
print("Loaded")

# run an fft on each column of cimg (echos are rows here)
cimg = Complex{Float16}.(fft(Complex{Float32}.(cimg),(1)));

# show fourier transformed cimg
# the curves that will be corrected
# by range cell migration should be visible
shape = size(cimg)
rccftpre = (abs.(view(cimg,
                2:100:shape[1],
                3600+54:1:3600+473)))
imshow(rccftpre);

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
    
    if n%10000==0
        print(n)
        print("  ")
        println(cellshift)
    end
end
cimg = cimg'
println("Done")

# show fourier transformed cimg
# the curves that will be corrected
# by range cell migration should be visible
shape = size(cimg)
rccftpost = (abs.(view(cimg,
                2:100:shape[1],
                3600+54:1:3600+473)))
imshow(rccftpost);

theta(s,R) = atan(vorbital*s/R)

#one way beam pattern:
p(a) = sinc(a*La/wavelength)
w(s,R) = p(theta(s,R))^2

Rc = R0+10000        # TODO where does 10000 come from? I think this is a focusing dist choice

R(s) = Rc - 1/2*vorbital^2/Rc*s^2

C(s) = exp(-4pi*im/wavelength*R(s))*w(s,Rc)

complexAzimuthFFT = let
    width = 200   # TODO where does this come from?
    s = 1/PRF*(range(-width, stop = width) |> collect)

    sig = C.(s)/sqrt(width)
    
    azimuth = vcat(sig, zeros(Complex{Float32},size(cimg)[1]-length(sig)))

    fft(azimuth)
end

print("Ready")

for i = 1:size(cimg)[2]
    line = Complex{Float32}.(cimg[:,i])
    
    ####### Azimuth Compression
    
    #lineFFT = fft(line)
    lineFFT = cimg[:,i]
    #lineFFT = fft(Complex{Float32}.(cimg[:,i]))
    crossCorrelated = AbstractFFTs.ifft(conj.(complexAzimuthFFT).*lineFFT)
    
    ####### End Azimuth Compression
    
    
    #complex = Complex.(I,Q)
    result = abs.( crossCorrelated )
    cimg[:,i] = result
    if i%1000 == 1
        print("#")
    end
end

shape = size(cimg)
azcompmag = abs.(view(cimg,1:16:shape[1],1:4:shape[2]));
azcompmag = reverse(azcompmag,dims=1)
imshow(azcompmag);
