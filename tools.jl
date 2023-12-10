using Images

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

function saveImg(image, name)
    data_scaled = (image .- minimum(image)) / (maximum(image) - minimum(image))
    save("$name.png", colorview(Gray, data_scaled));
end