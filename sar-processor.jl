using Distributed

cpus = 10;
newProcesses = cpus - size(workers(),1);
addprocs(newProcesses);

include("tools.jl");
include("processingFunctions.jl");
include("parser.jl")

pathname = "../data/"    # relative path of folder containing L1.0 data
imagename = "IMG-HH-ALPSRP273530900-H1"  # which image polarization to use - IMG-HH or IMG-HV for PALSAR
metadataName = "LED-ALPSRP273530900-H1"

R0 = 848665     #m
c = 3e8
Vr = 7593
La = 8.9 #m
rangeCells = 5000

print("Parsing Metadata... ");
(pulseSamples, sampleRate, PRF, wavelength, fdot, sampleRate, pulseLength) = parseMetadata(pathname, metadataName, rangeCells);
println("done");

print("Reading Data... ");
smallSignals = readData(pathname, imagename, rangeCells);
println("done")

print("Starting Deconvolution... ");
cimg = deconvolve(smallSignals, fdot, sampleRate, pulseLength, pulseSamples);
println("done");

GC.gc();

print("Starting Range Migration... ");
cimg = rangeMigration(cimg, c, sampleRate, PRF, wavelength, R0, Vr);
println("done");

print("Starting Azimuth Compression... ");
cimg = azimuthCompression(cimg, Vr, La, wavelength, R0, PRF);
println("done");

shape = size(cimg)
azcompmag = abs.(view(cimg,1:16:shape[1],1:4:shape[2]));
azcompmag = reverse(azcompmag,dims=1)

saveImg(azcompmag, "finalImage.png");