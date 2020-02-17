#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 11:49:34 2020

@author: florie

%SM4READER takes a file ID of an RHK .sm4 file, then reads

SEE => 

SM4 DATA FILE FORMAT V5.pdf

AND =>

http://unh2d.weebly.com/using-sm4-files-in-matlab.html
http://unh2d.weebly.com/sm4readerm-and-the-mat-file-format.html

"""

import numpy as np
from matplotlib.patches import Circle


def string(array):
    array = list(map(lambda x: chr(x), array))
    array = ("".join(array))
    return array


filepath='/home/florie/Documents/LABO/Python/SM4_File_Format/FCg01-2-010719_0061.sm4'
def readCitsSm4Bin(filepath):

    ObjectIDCode = ['Undefined', 'Page Index Header', 'Page Index Array',
                    'Page Header', 'Page Data', 'Image Drift Header',
                    'Image Drift', 'Spec Drift Header',
                    'Spec Drift Data (with X,Y coordinates)', 'Color Info',
                    'String data', 'Tip Track Header', 'Tip Track Data', 'PRM',
                    'Thumbnail', 'PRM Header', 'Thumbnail Header', 'API Info',
                    'History Info', 'Piezo Sensitivity',
                    'Frequency Sweep Data', 'Scan Processor Info', 'PLL Info',
                    'CH1 Drive Info', 'CH2 Drive Info', 'Lockin0 Info',
                    'Lockin1 Info', 'ZPI Info', 'KPI Info', 'Aux PI Info',
                    'Low-pass Filter0 Info', 'Low-pass Filter1 Info']
    f = open(filepath, "rb")

    # File Header
    Header_size = int(np.fromfile(f, dtype=np.uint16, count=1)[0])
    Signature = string(np.fromfile(f, dtype=np.uint16, count=18))
#    print(Signature)
    Total_Pagecount = int(np.fromfile(f, dtype=np.uint32, count=1)[0])
    ObjectListCount = int(np.fromfile(f, dtype=np.uint32, count=1)[0])
    ObjectFieldSize = int(np.fromfile(f, dtype=np.uint32, count=1)[0])
    Reserved = int(np.fromfile(f, dtype=np.uint32, count=2)[0])

    # iterate over the known objects from file header
    ObjectlistName = []
    ObjectlistOffset = []
    ObjectlistSize = []
    for i in range(ObjectListCount):
        ObjectlistName.append(ObjectIDCode[np.fromfile(f, dtype=np.uint32, count=1)[0]])
        ObjectlistOffset.append(np.fromfile(f, dtype=np.uint32, count=1))
        ObjectlistSize.append(np.fromfile(f, dtype=np.uint32, count=1))
    
    # Read and record the Page Index Header
    PageIndexHeader_PageCount = int(np.fromfile(f, dtype=np.uint32, count=1)[0]) # the number of pages in the Page Index Array
    PageIndexHeader_ObjectListCount = int(np.fromfile(f, dtype=np.uint32, count=1)[0])# %the count of objects stored after the Page Index Header
                                                                        #currently there is just one: Page Index Array    
    PageIndexHeader_Reserved = int(np.fromfile(f, dtype=np.uint32, count=2)[0]) #two fields reserved for future use 


    # Read and record the Page Index Array
    PageIndexHeader_ObjectID = int(np.fromfile(f, dtype=np.uint32, count=1)[0])  
    PageIndexHeader_Offset = int(np.fromfile(f, dtype=np.uint32, count=1)[0])   
    PageIndexHeader_Size = int(np.fromfile(f, dtype=np.uint32, count=1)[0])


    # Get info on each page
    PageIndex = []
    for j in range(PageIndexHeader_PageCount):
        PageIndex.append({'PageID': np.fromfile(f, dtype=np.uint16, count=8)[0],
                                'PageDataType': int(np.fromfile(f, dtype=np.uint32, count=1)[0]),
                                'PageSourceType': int(np.fromfile(f, dtype=np.uint32, count=1)[0]),
                                'ObjectListCount': int(np.fromfile(f, dtype=np.uint32, count=1)[0]),
                                'MinorVersion': int(np.fromfile(f, dtype=np.uint32, count=1)[0]),
                                'ObjectList': []
                                })
        for i in range(PageIndex[j]['ObjectListCount']):
            dict1 = {'ObjectID': ObjectIDCode[int(np.fromfile(f, dtype=np.uint32, count=1)[0])],
                     'Offset': int(np.fromfile(f, dtype=np.uint32, count=1)[0]),
                     'Size': int(np.fromfile(f, dtype=np.uint32, count=1)[0])}
            PageIndex[j]['ObjectList'].append(dict1)


    # Read and record each Pages Headers and Data
    PageHeader = []
    TextStrings = []
    Spatial = []
    Spectral = []
    SpectralInfo = []
    SpatialInfo = []
    topocount = 0
    currentcount = 0
    dIdV_Line_Speccount = 0
    for j in range(PageIndexHeader_PageCount):
        #f.seek(nbytes,0 = a partir du début) va au byte n + 1
        f.seek(PageIndex[j]['ObjectList'][0]['Offset'], 0)
    
    # Read Page Header    
        PageHeader.append({'FieldSize': int(np.fromfile(f, dtype=np.uint16, count=1)[0]),
                           'StringCount': int(np.fromfile(f, dtype=np.uint16, count=1)[0]),
                           'PageType_DataSource': int(np.fromfile(f, dtype=np.uint32, count=1)[0]),
                          'DataSubSource': int(np.fromfile(f, dtype=np.uint32, count=1)[0]),
                          'LineType': int(np.fromfile(f, dtype=np.uint32, count=1)[0]),
                          'XCorner': int(np.fromfile(f, dtype=np.uint32, count=1)[0]),
                          'YCorner': int(np.fromfile(f, dtype=np.uint32, count=1)[0]),
                          'Width': int(np.fromfile(f, dtype=np.uint32, count=1)[0]),
                          'Height': int(np.fromfile(f, dtype=np.uint32, count=1)[0]),
                          'ImageType': int(np.fromfile(f, dtype=np.uint32, count=1)[0]),
                          'ScanDirection': int(np.fromfile(f, dtype=np.uint32, count=1)[0]),
                          'GroupID': int(np.fromfile(f, dtype=np.uint32, count=1)[0]),
                          'PageDataSize': int(np.fromfile(f, dtype=np.uint32, count=1)[0]),
                          'MinZValue': int(np.fromfile(f, dtype=np.uint32, count=1)[0]),
                          'MaxZValue': int(np.fromfile(f, dtype=np.uint32, count=1)[0]),
                          'XScale': np.float(np.fromfile(f, dtype=np.float32, count=1)[0]),#single = 4 bytes Floating-point numbers
                          'YScale': np.float(np.fromfile(f, dtype=np.float32, count=1)[0]),
                          'ZScale': np.float(np.fromfile(f, dtype=np.float32, count=1)[0]),
                          'XYScale': np.float(np.fromfile(f, dtype=np.float32, count=1)[0]),
                          'XOffset': np.float(np.fromfile(f, dtype=np.float32, count=1)[0]),
                          'YOffset': np.float(np.fromfile(f, dtype=np.float32, count=1)[0]),
                          'ZOffset': np.float(np.fromfile(f, dtype=np.float32, count=1)[0]),
                          'Period': np.float(np.fromfile(f, dtype=np.float32, count=1)[0]),
                          'Bias': np.float(np.fromfile(f, dtype=np.float32, count=1)[0]),
                          'Current': np.float(np.fromfile(f, dtype=np.float32, count=1)[0]),
                          'Angle': np.float(np.fromfile(f, dtype=np.float32, count=1)[0]),
                          'ColorInfoListCount': int(np.fromfile(f, dtype=np.uint32, count=1)[0]),
                          'GridXSize': int(np.fromfile(f, dtype=np.uint32, count=1)[0]),
                          'GridYSize': int(np.fromfile(f, dtype=np.uint32, count=1)[0]),
                          'ObjectListCount': int(np.fromfile(f, dtype=np.uint32, count=1)[0]),
                          'ObjectList': []
                          })
        # Skip flags and reserved data
        f.seek(1+3+60,1);
        # Read the Object List
        for i in range(PageHeader[j]['ObjectListCount']):
            dict1 = {'ObjectID': int(np.fromfile(f, dtype=np.uint32, count=1)[0]),
                     'Offset': int(np.fromfile(f, dtype=np.uint32, count=1)[0]),
                     'Size': int(np.fromfile(f, dtype=np.uint32, count=1)[0])}
            PageHeader[j]['ObjectList'].append(dict1)

        # PageheaderObjectList
        c = int(np.fromfile(f, dtype=np.uint16, count=1)[0])
        d = {'strLabel': string(np.fromfile(f, dtype=np.uint16, count=c))} # eg "current image"
        count = int(np.fromfile(f, dtype=np.uint16, count=1)[0])
        d['strSystemText'] = string(np.fromfile(f, dtype=np.uint16, count=count))

        count = int(np.fromfile(f, dtype=np.uint16, count=1)[0])
        d['strSessionText'] = string(np.fromfile(f, dtype=np.uint16, count=count))

        count = int(np.fromfile(f, dtype=np.uint16, count=1)[0])
        d['strUserText'] = string(np.fromfile(f, dtype=np.uint16, count=count))

        count = int(np.fromfile(f, dtype=np.uint16, count=1)[0])
        d['strPath'] = string(np.fromfile(f, dtype=np.uint16, count=count))

        count = int(np.fromfile(f, dtype=np.uint16, count=1)[0])
        d['strDate'] = string(np.fromfile(f, dtype=np.uint16, count=count))

        count = int(np.fromfile(f, dtype=np.uint16, count=1)[0])  # DAQ time
        d['strTime'] = string(np.fromfile(f, dtype=np.uint16, count=count))

        count = int(np.fromfile(f, dtype=np.uint16, count=1)[0])  # physical units of x axis
        d['strXUnits'] = string(np.fromfile(f, dtype=np.uint16, count=count))

        count = int(np.fromfile(f, dtype=np.uint16, count=1)[0])
        d['strYUnits'] = string(np.fromfile(f, dtype=np.uint16, count=count))

        count = int(np.fromfile(f, dtype=np.uint16, count=1)[0])
        d['strZUnits'] = string(np.fromfile(f, dtype=np.uint16, count=count))

        count = int(np.fromfile(f, dtype=np.uint16, count=1)[0])
        d['strYLabel'] = string(np.fromfile(f, dtype=np.uint16, count=count))

        count = int(np.fromfile(f, dtype=np.uint16, count=1)[0])
        d['strStatusChannelText'] = string(np.fromfile(f, dtype=np.uint16, count=count))

        count = int(np.fromfile(f, dtype=np.uint16, count=1)[0])
        # contains last saved line count for an image data page
        d['strCompletedLineCount'] = string(np.fromfile(f, dtype=np.uint16, count=count))

        count = int(np.fromfile(f, dtype=np.uint16, count=1)[0])
        # Oversampling count for image data pages
        d['strOverSamplingCount'] = string(np.fromfile(f, dtype=np.uint16, count=count))

        count = int(np.fromfile(f, dtype=np.uint16, count=1)[0])
        # voltage at which the sliced image is created from the spectra page.  empty if not a sliced image
        d['strSlicedVoltage'] = string(np.fromfile(f, dtype=np.uint16, count=count))

        count = int(np.fromfile(f, dtype=np.uint16, count=1)[0])
        # PLLPro status text: blank, master or user
        d['strPLLProStatus'] = string(np.fromfile(f, dtype=np.uint16, count=count))

        count = int(np.fromfile(f, dtype=np.uint16, count=1)[0])
        # ZPI controller item's set-point unit
        d['strSetpointUnit'] = string(np.fromfile(f, dtype=np.uint16, count=count))

        count = int(np.fromfile(f, dtype=np.uint16, count=1)[0])
        # stores value of CH1 and CH2 if they are in hardware space
        d['strCHDriveValues'] = string(np.fromfile(f, dtype=np.uint16, count=count))
    
        TextStrings.append(d)
        
        # Page Header – Sequential Data Page - Get Datas
        Data = []
        f.seek(PageIndex[j]['ObjectList'][1]['Offset'],0)
        Data.append(np.fromfile(f, dtype=np.uint32, count=int(round(PageIndex[j]['ObjectList'][1]['Size']/4))))
        print(np.shape(Data))
        #/4 because total data size is divided by the number of bytes that use each 'long' data
        ScaleData = PageHeader[j]['ZOffset']+Data[0]*PageHeader[j]['ZScale']
        ScaleData = np.reshape(ScaleData,(PageHeader[j]['Width'],PageHeader[j]['Height']), order="F")
        #order Fortran = "F" to match readCITSsm4File function
        print(np.shape(ScaleData))

        ###################### Spatial Data
        # is it label topo ? # cf p.28 of SM4 DATA FILE FORMAT V5.pdf
        print(TextStrings[j]['strLabel'])
        if TextStrings[j]['strLabel'] == 'Topography' and PageHeader[j]['PageType_DataSource'] == 1:
            topocount += 1
            Spatialpagenumber = j  # records last spatial page; to read the rest of the spatial data later
# !!! rotate ? transpose ?
            Spatial.append({'TopoData': np.rot90(ScaleData,3)})
            SpatialInfo.append({'TopoUnit': TextStrings[j]['strZUnits']})
            if topocount == 1:
                FImage = Spatial[-1]['TopoData']  # Image forward
            elif topocount ==2:
                BImage = Spatial[-1]['TopoData']  # Image bacward
            else:
                print('there is more topo data than expected')
        # is it spatial Current data ?
        elif TextStrings[j]['strLabel'] == 'Current' and PageHeader[j]['PageType_DataSource'] == 2:
            currentcount += 1
            Spatial.append({'CurrentData': np.rot90(ScaleData,3)})  
            Spatial[-1]['CurrentUnit'] = TextStrings[j]['strZUnits']
            Spatialpagenumber = j
            if currentcount == 1:
                FImage_I = Spatial[-1]['CurrentData']  # Image forward
            elif currentcount ==2:
                BImage_I = Spatial[-1]['CurrentData']  # Image bacward
            else:
                print('there is more current data than expected')

        ###################### Spectral Data - can be Point or Line
        # Is it Spectral Point(38) ? Not taken in charge currently
        elif PageHeader[j]['PageType_DataSource'] == 38:
            print('You didnt load a CITS. Use sm4_reader to read this data')
            
        # CITS is recarded as a line of LIA
        # Is it Spectral Line(16) LIA ?
        elif TextStrings[j]['strLabel'] == 'LIA' and PageHeader[j]['PageType_DataSource'] == 16:
            dIdV_Line_Speccount += 1;
            Spectral.append({'dIdV_Line_Data': ScaleData})
            Linespectrapagenumber = j
            if dIdV_Line_Speccount == 1 :
                SpectralData_y = Spectral[-1]['dIdV_Line_Data']
                print(PageHeader[j]['XScale'])
                #Get spectra coordinates offset
                for objectNbre in range(PageHeader[j]['ObjectListCount']):
                    # 7 == Tip track Info Header
                    if PageHeader[j]['ObjectList'][objectNbre]['ObjectID'] == 7:
                        TipTrackInfo_offset = PageHeader[j]['ObjectList'][objectNbre]['Offset']
                        TipTrackInfo_Size = PageHeader[j]['ObjectList'][objectNbre]['Size']
                    # 8 == Tip track Data
                    elif PageHeader[j]['ObjectList'][objectNbre]['ObjectID'] == 8:
                        TipTrackData_offset = PageHeader[j]['ObjectList'][objectNbre]['Offset']
                        TipTrackData_Size = PageHeader[j]['ObjectList'][objectNbre]['Size']
            else:
                print('there is more spectral data than expected')

        # We could add 'Spaectral line Current', or 'PLL Amplitude' ;
        # 'PLL Phase' ; 'dF' ; 'PLL Drive'spectra, AFM, ... in the same way
        else:
            print('Data Type '+str(TextStrings[j]['strLabel'])
                  + ' or PageType_DataSource '
                  + str(PageHeader[j]['PageType_DataSource'])
                  + ' not found. Check spelling or add it.')

    ###################### Get spectra coordinates
    # Nbre of measurements taken along a line :
    nbreScans = np.shape(SpectralData_y)[0]
    print(nbreScans)
    xCoord = np.zeros(nbreScans)
    print(np.shape(xCoord))
    yCoord = np.zeros(nbreScans)
    f.seek(TipTrackData_offset,0)  # Go to beggining of header
    for i in range(nbreScans):
        f.seek(4, 1)  # skip start time
        a = np.float(np.fromfile(f, dtype=np.float32, count=1)[0])
#        print(a)
        xCoord[i] = a
        yCoord[i] = np.float(np.fromfile(f, dtype=np.float32, count=1)[0])
        f.seek(16, 1)  # skip dx dy xcumul ycumul fields
    print(np.shape(xCoord))
    print(xCoord[0])
    print((PageHeader[Linespectrapagenumber]['Width']))
    SpectralData_x = PageHeader[Linespectrapagenumber]['XOffset'] + PageHeader[Linespectrapagenumber]['XScale'] * np.array(list(range(0,PageHeader[Linespectrapagenumber]['Width']))) * 1000.0#mV
    #each measurement is taken at a coordinate,
    #but several measurements (repetitions) are taken on the same spot.
    #Eg if repetitions = 4, 2 measurements, one forward, one backward (saved in right direction)
    repetitions = 0
    pointdiff = 0
    repindex = 0
    # step through data points until there is a difference between the current (x,y) point and the next (x1,y1) point
    while pointdiff == 0: 
        pointdiff = (xCoord[repindex+1]- xCoord[repindex]) + (xCoord[repindex+1]- xCoord[repindex]) #adds the x and y difference values together.  if this is anything other than 0, we have found the limit of the repetitions
        repetitions = repetitions+1
        repindex = repindex+1
    numberOfPlots = nbreScans/repetitions #this is the number of distinct plotting locations

    # Center coordinates and metric dimensions in nm
    xL = abs(PageHeader[Spatialpagenumber]['XScale']*PageHeader[Spatialpagenumber]['Width']) * (10 ** 9)
    yL = abs(PageHeader[Spatialpagenumber]['YScale']*PageHeader[Spatialpagenumber]['Height']) * (10 ** 9)
    xC = PageHeader[Spatialpagenumber]['XOffset'] * (10 ** 9)
    yC = PageHeader[Spatialpagenumber]['YOffset'] * (10 ** 9) * (10 ** 9)
    #repetition_index = 1 #0(forw),1(back),2(forw)... according to the number of repetitions.
    #size spatial data
    xPx = int(PageHeader[Spatialpagenumber]['Height'])
    yPx = int(PageHeader[Spatialpagenumber]['Width'])
    #size spectral data ( not necessarily the same in RHK )
    xSpec = int(np.sqrt(numberOfPlots))
    ySpec = int(np.sqrt(numberOfPlots))
    zPt = int(len(SpectralData_y[:, 0]))#nbre of points in spec data
#        x_m = np.linspace(0, xL*1e+9, xPx)#?
#        y_m = np.linspace(0, yL*1e+9, yPx)
    try:
        topo = np.zeros(shape=(xPx, yPx))
        m_data = np.zeros(shape=(repetitions, ySpec, xSpec, zPt))
        print(np.shape(m_data))
    except MemoryError:
        print("The data is too big ! Or the memory too small...")
        return False
        
    #in Spectraldata_y, dIdV info corresponds to the spec data saved from left to right and increasing y(downwards in RHK), with same nbre of repetions at each spot
    for r in range(repetitions):#even : forward, odd : backwards
        for y in range(ySpec):#len(SpectralData_y[:,0])):
            for x in range(xSpec):
                m_data[r][y][x] = SpectralData_y[:, (xSpec*y+x)*repetitions+r]
    
    patch = []
    for m in range(0, numberOfPlots, repetitions): #iterate over the number of locations
                        patch.append(Circle((-xC + xL/2 + xCoord[m]*1e+9,- yC + yL/2 + yCoord[m]*1e+9), yL/5000,facecolor='r',edgecolor='None'))

    channelList = ['Data {}'.format(i) for i in range(repetitions)]

    m_params = {"xPx": xSpec, "yPx": ySpec, "xL": xL, "yL": yL, "zPt": zPt,
                "vStart": SpectralData_x[0],"vEnd": SpectralData_x[-1], "dV": abs(SpectralData_x[-1] - SpectralData_x[0])/zPt,
                "Patch": patch}

    topo = FImage
    print(np.shape(topo))
    
    #!!! create average Data :
    average = True
    if average:
        average = np.zeros(shape=(ySpec, xSpec, zPt))
        for y in range(ySpec):#len(SpectralData_y[:,0])):
            for x in range(xSpec):
                for r in range(repetitions):
                    average[y][x] += SpectralData_y[:, (xSpec*y+x)*repetitions+r]/repetitions
        print('okk')
#        self.addChannel(average, "average")
    return True
        
readCitsSm4Bin(filepath)
