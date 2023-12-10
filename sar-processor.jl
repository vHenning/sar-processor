using ImageView
using Serialization
using FFTW

include("tools.jl");
include("processingFunctions.jl");

pathname = "../data/"    # relative path of folder containing L1.0 data
imagename = "IMG-HH-ALPSRP273530900-H1"  # which image polarization to use - IMG-HH or IMG-HV for PALSAR

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

println("Done")

print("Starting Deconvolution... ");
cimg = deconvolve(smallSignals, chirpFFT, pulseSamples);
println("done");

GC.gc();

print("Starting Range Migration... ");
cimg = rangeMigration(cimg, c, sampleRate, PRF, wavelength, R0);
println("done");

# show fourier transformed cimg
# the curves that will be corrected
# by range cell migration should be visible
shape = size(cimg)
rccftpost = (abs.(view(cimg,
                2:100:shape[1],
                3600+54:1:3600+473)))

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

saveImg(azcompmag, "finalImage.png");