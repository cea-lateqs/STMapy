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
# -*- coding: utf-8 -*-
import struct
import sys
import codecs
import unittest
import numpy as np


def string(array):
    array = list(map(lambda x: chr(x), array))
    array = ("".join(array))
    return array


filepath='/home/florie/Documents/LABO/Python/SM4_File_Format/FCg01-2-010719_0061.sm4'
def loadCitsSm4Bin(filepath):

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


    # Read and record the Page Index Headers and Page Headers for each Page
    PageHeader = []
    TextStrings = []
    for j in range(PageIndexHeader_PageCount):
        #f.seek(nbytes,0 = a partir du début) va au byte n + 1
        f.seek(PageIndex[j]['ObjectList'][0]['Offset'], 0)
    
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

#    print(PageIndexHeader_PageCount)
#        print(PageIndex[6]['PageID'])
#        print(PageHeader[j]['ObjectListCount'])
#        print((PageHeader[j]['ObjectList']))
#        print(type(PageHeader[j]['ObjectList']))

        # Skip flags and reserved data
        f.seek(1+3+60,1);
        # Read the Object List
        for i in range(PageHeader[j]['ObjectListCount']):
            dict1 = {'ObjectID': int(np.fromfile(f, dtype=np.uint32, count=1)[0]),
                     'Offset': int(np.fromfile(f, dtype=np.uint32, count=1)[0]),
                     'Size': int(np.fromfile(f, dtype=np.uint32, count=1)[0])}
#            print(dict1)
            PageHeader[j]['ObjectList'].append(dict1)
#            print(np.shape(PageHeader[j]['ObjectList']))

#        print(PageHeader[j]['ObjectList'][20]['ObjectID'])
#        print(PageHeader[j]['ObjectList'][20]['Offset']-PageIndex[j]['ObjectList'][0]['Offset'])
# eg 10 : History info ; 5 : SectralData ; 0 : empty ?.

        # Read and record the Text Strings = (str lenghth ; str) for each Page

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

############################################################
    # Initialize variables that track how many pages of each type there are

    topocount = 0  # Topography
    currentcount = 0  # Spatial current
    liacurrentcount = 0  # Spatial LIA current
    IV_Point_Speccount = 0  # Spectral point current
    dIdV_Point_Speccount = 0  # Spectral point LIA (dIdV) current
    IV_Line_Speccount = 0  # Spectral line current
    dIdV_Line_Speccount = 0  # Spectral line LIA (dIdV) current
    Spatialpagenumber = 0  # Total number of spatial pages
    Pointspectrapagenumber = 0  # Total number of point line pages
    Linespectrapagenumber = 0  # Total number of line spectra pages
    Frequencyspectrapagenumber = 0  # Total number of frequency sweep pages


    # Page Header – Sequential Data Page - Get Datas
    Data = []
    ScaleData = []
    
    Spatial = []
    Spectral = []
    SpectralInfo = {}
    SpatialInfo = {}
    # There is only one type of data per page. But several page with same type of Data.
    # eg. Topo forward, topo backward
    for j in range(PageIndexHeader_PageCount):
        f.seek(PageIndex[j]['ObjectList'][1]['Offset'],0)
        Data.append(np.fromfile(f, dtype=np.uint32, count=int(round(PageIndex[j]['ObjectList'][1]['Size']/4))))
        print(np.shape(Data))
        #/4 because total data size is divided by the number of bytes that use each 'long' data
        a = PageHeader[j]['ZOffset']+Data[j]*PageHeader[j]['ZScale']
        print(np.shape(a))
        a = np.reshape(a,(PageHeader[j]['Width'],PageHeader[j]['Height']), order="F")
        #order Fortran = "F" to match readCITSsm4File function
        ScaleData.append(a)
        print(np.shape(ScaleData))


        ###################### Spatial Data
        # is it label topo ? # cf p.28 of SM4 DATA FILE FORMAT V5.pdf
        print(TextStrings[j]['strLabel'])
        if TextStrings[j]['strLabel'] == 'Topography' and PageHeader[j]['PageType_DataSource'] == 1:
            print('Topo!!')
            topocount += 1
# !!! rotate ? transpose ?
            Spatial.append({'TopoData': np.rot90(ScaleData[j],-1)})
            Spatial[-1]['TopoUnit'] = TextStrings[j]['strZUnits']
            Spatialpagenumber = j  # records last spatial page; to read the rest of the spatial data later
        
        # is it spatial Current data ?
        elif TextStrings[j]['strLabel'] == 'Current' and PageHeader[j]['PageType_DataSource'] == 2:
            currentcount += currentcount
            Spatial.append({'CurrentData': np.rot90(ScaleData[j],-1)})  
            Spatial[-1]['CurrentUnit'] = TextStrings[j]['strZUnits']
            Spatialpagenumber = j
        
        # is it spatial LIA current data ?
        elif TextStrings[j]['strLabel'] == 'LIA Current' and PageHeader[j]['PageType_DataSource'] == 2:
            # !!! Just 'LIA' ?
            # if the comparison returns a true, flattendebug is set to true
            liacurrentcount = liacurrentcount+1;
            Spatial.append({'LIACurrentData': np.rot90(ScaleData[j],-1)})
            Spatial[-1]['LIACurrentUnit'] = TextStrings[j]['strZUnits']
            Spatialpagenumber = j

        # 'PLL Amplitude' ; 'PLL Phase' ; 'dF' ;
        # 'PLL Drive' could be added in the same way

        ###################### Spectral Data - can be Point or Line
        # !!! Can there be several data types in one page ? elif = if ? No as ScaleData contains only one ? Auquel cas Spectral[-1] ou Spectral[new data index]
        # Is it Spectral Point(38) Current ?
        elif TextStrings[j]['strLabel'] == 'Current' and PageHeader[j]['PageType_DataSource'] == 38:
            IV_Point_Speccount = IV_Point_Speccount + 1;
            Spectral.append({'IV_Point_Data': ScaleData[j]})
            Spectral[-1]['IV_Point_DataUnit'] = 'A'
            spectralpagenumber = j
            Pointspectrapagenumber = j

        # Is it Spectral Point(38) LIA Current ?
        elif TextStrings[j]['strLabel'] == 'LIA' and PageHeader[j]['PageType_DataSource'] == 38:
            dIdV_Point_Speccount += 1
            Spectral.append({'dIdV_Point_Data': ScaleData[j]})
            Spectral[-1]['dIdV_Point_DataUnit'] = 'A'
            spectralpagenumber = j
            Pointspectrapagenumber = j

        # Is it Spectral Line(16) Current ?
        elif TextStrings[j]['strLabel'] == 'Current' and PageHeader[j]['PageType_DataSource'] == 16:
            IV_Line_Speccount = IV_Line_Speccount + 1;
            Spectral.append({'IV_Line_Data': ScaleData[j]})
            Spectral[-1]['IV_Line_DataUnit'] = 'A'
            spectralpagenumber = j
            Linespectrapagenumber = j

        # Is it Spectral Line(16) LIA Current ?
        elif TextStrings[j]['strLabel'] == 'LIA' and PageHeader[j]['PageType_DataSource'] == 16:
            dIdV_Line_Speccount = dIdV_Line_Speccount + 1;
            Spectral.append({'dIdV_Line_Data': ScaleData[j]})
            Spectral[-1]['dIdV_Line_DataUnit'] = 'A'
            spectralpagenumber = j
            Linespectrapagenumber = j

        # We could add AFM spectra in th same way

        else:
            print('Data Type '+str(TextStrings[j]['strLabel'])
                  + ' or PageType_DataSource '
                  + str(PageHeader[j]['PageType_DataSource'])
                  + ' not found. Check spelling or add it.')
    print(Spatialpagenumber)

    # If pages are empty :
    if (topocount + currentcount + liacurrentcount + IV_Point_Speccount +
       dIdV_Point_Speccount + IV_Line_Speccount + dIdV_Line_Speccount == 0):
        print('No data found in the file')
        return None

    # If there are any spatial pages ; Get Spatial information
    if topocount + currentcount + liacurrentcount > 0:
        SpatialInfo = PageHeader['Spatialpagenumber']
        SpatialInfo['widthheightUnit'] = 'm'
        SpatialInfo['biasUnit'] = 'V'
        SpatialInfo['currentUnit'] = 'A'
        SpatialInfo['xyoffsetUnit'] = 'm'
        SpatialInfo['width'] = abs(PageHeader['Spatialpagenumber']['XScale']*PageHeader['Spatialpagenumber']['Width'])
        SpatialInfo['height'] = abs(PageHeader['Spatialpagenumber']['YScale']*PageHeader['Spatialpagenumber']['Height'])

 

#    %There are three types of spatial data: frequency, point (IV,V,Z), line
#    %(IV).
  
    # If there are any spectral pages ; Get spectral information
    if (IV_Point_Speccount + dIdV_Point_Speccount
       + IV_Line_Speccount + dIdV_Line_Speccount) > 0:
        SpectralInfo = PageHeader['spectralpagenumber']
        SpectralInfo['biasUnit'] = 'V'
        SpectralInfo['currentUnit'] = 'A'
        SpectralInfo['xyoffsetUnit'] = 'm'
        SpectralInfo['xpoints'] =  PageHeader['spectralpagenumber']['Width']
        SpectralInfo['xscale'] = list(range(0,SpectralInfo['xpoints']))
        SpectralInfo['xdata'] =  PageHeader['spectralpagenumber']['XOffset'] + PageHeader['spectralpagenumber']['XScale']*SpectralInfo['xscale']
# !!! PH(p) ?? outfile.PageHeader(p).XOffset + outfile.PageHeader(p).XScale * xscale)
        #SpectralInfo['xdataUnit'] =  TextString(p).strXUnits;


    # If there are any point spectra pages ; Get each point spectra information
    if Pointspectrapagenumber > 0:
        for objectnumber in range(PageHeader['Pointspectrapagenumber']['ObjectListCount']):
            if PageHeader['Pointspectrapagenumber'][objectnumber]['ObjectID'] == 7:  # Not the good index ! should be 4 ?
                PageHeader['Pointspectrapagenumber'].SpecDriftHeader.Offset = PageHeader(Pointspectrapagenumber).Offset(objectnumber);
#                    outfile.PageHeader(Pointspectrapagenumber).SpecDriftHeader.Size = outfile.PageHeader(Pointspectrapagenumber).Size(objectnumber);


#        for objectnumber = 1:outfile.PageHeader(Pointspectrapagenumber).ObjectListCount
#            if outfile.PageHeader(Pointspectrapagenumber).ObjectID(objectnumber) == 8
#                    outfile.PageHeader(Pointspectrapagenumber).SpecDriftData.Offset = outfile.PageHeader(Pointspectrapagenumber).Offset(objectnumber);
#                    outfile.PageHeader(Pointspectrapagenumber).SpecDriftData.Size = outfile.PageHeader(Pointspectrapagenumber).Size(objectnumber);


 
#            %read the Spec data for tip start time and position
#            fseek(fileID,outfile.PageHeader(Pointspectrapagenumber).SpecDriftData.Offset,'bof');
#            formatoutfile.Spectral.startTime = fread(fileID,1,'single');
#            formatoutfile.Spectral.xCoord = fread(fileID,1,'single');
#            formatoutfile.Spectral.yCoord = fread(fileID,1,'single');
#            formatoutfile.Spectral.xyCoordUnit = 'm';
#            formatoutfile.Spectral.dx = fread(fileID,1,'single');
#            formatoutfile.Spectral.dy = fread(fileID,1,'single');
#            formatoutfile.Spectral.xCumulative = fread(fileID,1,'single');
#            formatoutfile.Spectral.yCumulative = fread(fileID,1,'single');
#            formatoutfile.Spectral.type = 'Point';
#   
#    
#            if Frequencyspectrapagenumber > 0
#                   if formatoutfile.Spectral.xdataUnit == 'Hz'
#                        formatoutfile.Spectral.type = 'Frequency';
#                   end
#
#                   if formatoutfile.Spectral.xdataUnit == 'V'
#                        formatoutfile.Spectral.type = 'Vspec';
#                   end
#
#                   if formatoutfile.Spectral.xdataUnit == 'm'
#                        formatoutfile.Spectral.type = 'Zspec';
#          
                
                
#    %If there are point spectra pages, add the specific point spectra
#    %information to formatoutfile.Spectral
#    if Linespectrapagenumber > 0
#    %need to find out which pages have the extra header data
#        for objectnumber = 1:outfile.PageHeader(Linespectrapagenumber).ObjectListCount
#            if outfile.PageHeader(Linespectrapagenumber).ObjectID(objectnumber) == 7
#                    outfile.PageHeader(Linespectrapagenumber).SpecDriftHeader.Offset = outfile.PageHeader(Linespectrapagenumber).Offset(objectnumber);
#                    outfile.PageHeader(Linespectrapagenumber).SpecDriftHeader.Size = outfile.PageHeader(Linespectrapagenumber).Size(objectnumber);



#        for objectnumber = 1:outfile.PageHeader(Linespectrapagenumber).ObjectListCount
#            if outfile.PageHeader(Linespectrapagenumber).ObjectID(objectnumber) == 8
#                    outfile.PageHeader(Linespectrapagenumber).SpecDriftData.Offset = outfile.PageHeader(Linespectrapagenumber).Offset(objectnumber);
#                    outfile.PageHeader(Linespectrapagenumber).SpecDriftData.Size = outfile.PageHeader(Linespectrapagenumber).Size(objectnumber);



#        [~, numberofscans] = size(outfile.ScaleData{Linespectrapagenumber}); %determine number of scans in the data set
#        formatoutfile.Spectral.xCoord = zeros(numberofscans,1);%initialize the arrays to hold the coordinates
#        formatoutfile.Spectral.yCoord = zeros(numberofscans,1);
#        formatoutfile.Spectral.xyCoordUnit = 'm';
#        
#        %read the Spec data for tip start time and position
#        fseek(fileID,outfile.PageHeader(Linespectrapagenumber).SpecDriftData.Offset,'bof'); %go to the beginning of the header
#        for k = 1:numberofscans %the header repeats a sequence of information, need to iterate over all points
#            fseek(fileID,4,'cof'); %skip the StartTime
#            formatoutfile.Spectral.xCoord(k) = fread(fileID,1,'single'); %Record x and y coordinates
#            formatoutfile.Spectral.yCoord(k) = fread(fileID,1,'single');
#            fseek(fileID,16,'cof'); %skip the dx dx xcumulative and ycumulative fields since we won't record them here
#        end
#           
#        fseek(fileID,outfile.PageHeader(Linespectrapagenumber).SpecDriftData.Offset,'bof');
#        formatoutfile.Spectral.startTime = fread(fileID,1,'single');
#        fseek(fileID,8,'cof'); %skip the x and y coord
#        formatoutfile.Spectral.dx = fread(fileID,1,'single');
#        formatoutfile.Spectral.dy = fread(fileID,1,'single');
#        formatoutfile.Spectral.xCumulative = fread(fileID,1,'single');
#        formatoutfile.Spectral.yCumulative = fread(fileID,1,'single');
#        formatoutfile.Spectral.type = 'Line';
#           
    
#    
#    %If there is PLL data present, create a PLL header and format to
#    %contain the relevant information.
#    if PLLTrue == 1 %our test for the presence of PLL data
#    
#        %need to find out which pages have the extra header data
#        for objectnumber = 1:outfile.PageHeader(PLLPagenumber).ObjectListCount
#            if outfile.PageHeader(PLLPagenumber).ObjectID(objectnumber) == 22
#                    outfile.SpatialPLL.Offset = outfile.PageHeader(PLLPagenumber).Offset(objectnumber);
#                    outfile.SpatialPLL.Size = outfile.PageHeader(PLLPagenumber).Size(objectnumber);
#            end
#        end
#        
#        fseek(fileID,outfile.SpatialPLL.Offset,'bof'); %go to the relevant offset
#        fseek(fileID,8,'cof'); %skip the beginning, not useful information
#        formatoutfile.PLL.DriveAmplitude = fread(fileID,1,'double'); %begin reading the data
#        formatoutfile.PLL.DriveAmplitudeUnit = 'V';
#        
#        formatoutfile.PLL.DriveRefFrequency = fread(fileID,1,'double');
#        formatoutfile.PLL.DriveRefFrequencyUnit = 'Hz';
#        
#        formatoutfile.PLL.LockinFrequencyOffset = fread(fileID,1,'double');
#        formatoutfile.PLL.LockinFrequencyOffsetUnit = 'Hz';
#        
#        formatoutfile.PLL.LockinHarmonicFactor = fread(fileID,1,'double');
#            
#        formatoutfile.PLL.LockinPhaseOffset = fread(fileID,1,'double');
#        formatoutfile.PLL.LockinPhaseUnit = 'deg';
#        
#        formatoutfile.PLL.PIGain = fread(fileID,1,'double');
#        formatoutfile.PLL.PIGainUnit = 'Hz/deg';
#        
#        formatoutfile.PLL.PIntCutOffFreq = fread(fileID,1,'double');
#        formatoutfile.PLL.PIntCutOffFreqUnit = 'Hz';
#        
#        formatoutfile.PLL.LowerBound = fread(fileID,1,'double');
#        formatoutfile.PLL.UpperBound = fread(fileID,1,'double');
#        formatoutfile.PLL.PIOutputUnit = 'Hz';
#        
#        formatoutfile.PLL.DissPIGain = fread(fileID,1,'double');
#        formatoutfile.PLL.DissPIGainUnit = 'V/V';
#        
#        formatoutfile.PLL.DissIntCutOffFreq = fread(fileID,1,'double');
#        formatoutfile.PLL.DissIntCutOffFreqUnit = 'Hz';
#        
#        formatoutfile.PLL.DissLowerBound = fread(fileID,1,'double');
#        formatoutfile.PLL.DissUpperBound = fread(fileID,1,'double');
#        formatoutfile.PLL.DissPIOutputUnit = 'V';
#        
        
    print('****')
    print(np.shape(Spatial))
    print(np.shape(Spectral))
    return Spatial, Spectral


loadCitsSm4Bin(filepath)

#        Spatial = loadedData['Spatial']
#        Spectral = loadedData['Spectral']
#        Image = Spatial['TopoData'][0, 0]
#        FImage = Image[0, 0]# Image forward
#        BImage = Image[0, 1]#Image bacward
#        spectraType = Spectral['type'][0, 0]#see scipy.org : scipy.io.loadmat
