using ImageView
using Serialization
using FFTW

include("tools.jl");
include("processingFunctions.jl");
include("parser.jl")

pathname = "../data/"    # relative path of folder containing L1.0 data
imagename = "IMG-HH-ALPSRP273530900-H1"  # which image polarization to use - IMG-HH or IMG-HV for PALSAR

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

print("Reading Data... ");
smallSignals = readData(pathname, imagename, rangeCells);
println("done")

print("Starting Deconvolution... ");
cimg = deconvolve(smallSignals, chirpFFT, pulseSamples);
println("done");

GC.gc();

print("Starting Range Migration... ");
cimg = rangeMigration(cimg, c, sampleRate, PRF, wavelength, R0, Vr);
println("done");

print("Starting Azimuth Compression... ");
cimg = azimuthCompression(cimg, vorbital, La, wavelength, R0, PRF);
println("done");

shape = size(cimg)
azcompmag = abs.(view(cimg,1:16:shape[1],1:4:shape[2]));
azcompmag = reverse(azcompmag,dims=1)

saveImg(azcompmag, "finalImage.png");