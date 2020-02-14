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
        f.seek(PageIndex[j]['ObjectList'][0]['Offset'], 0)
#        fseek(fileID,outfile.PageIndex(j).Offset(1) ,'bof'); %use the offsets
#        %we read to earlier to find the beginning of the Page Header
    
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
#    print(PageIndex[6]['PageID'])
#    print(PageHeader[0]['XScale'])
#    
        # Skip flags and reserved data
        f.seek(1+3+60,1);

        # Read the Object List
        for i in range(PageIndex[j]['ObjectListCount']):
            dict1 = {'ObjectID': ObjectIDCode[int(np.fromfile(f, dtype=np.uint32, count=1)[0])],
                     'Offset': int(np.fromfile(f, dtype=np.uint32, count=1)[0]),
                     'Size': int(np.fromfile(f, dtype=np.uint32, count=1)[0])}
            PageHeader[j]['ObjectList'].append(dict1)

        # Read and record the Text Strings = (str lenghth ; str) for each Page

        # PageheaderObjectList
        count = int(np.fromfile(f, dtype=np.uint16, count=1)[0])
        d = {'strLabel': string(np.fromfile(f, dtype=np.uint16, count=count))} # eg "current image"

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


    # Page Header – Sequential Data Page
    Data = []
    ScaleData = []
    for j in range(PageIndexHeader_PageCount):
        f.seek(PageIndex[j]['ObjectList'][1]['Offset'],0)
        Data.append(np.fromfile(f, dtype=np.uint32, count=int(round(PageIndex[j]['ObjectList'][1]['Size']/4))))
        #/4 because total data size is divided by the number of bytes that use each 'long' data
        a = PageHeader[j]['ZOffset']+Data[j]*PageHeader[j]['ZScale']
        a = np.reshape(a,(PageHeader[j]['Width'],PageHeader[j]['Height']), order="F")
        #order Fortran = "F" to match readCITSsm4File function
        ScaleData.append(a)


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
    SpatialPLLpagenumber = 0  # Spatial pages
    PLLPageNumber = 0  # PLL Page
    PLLTrue = 0  # Logical for if there is PLL data, start at 0 (false)

    # number of PLL and AFM pages
    PLLAmpcount = 0
    PLLPhasecount = 0
    dFcount = 0
    dFspeccount = 0
    PLLDrivecount = 0
    PLLAmp_Speccount = 0
    PLLPhase_Speccount = 0
    PLLDrive_Speccount = 0

    Spatial = []
    for j in range(PageIndexHeader_PageCount):
        #Spatial Data : is label topo ?
        if TextStrings[j]['strLabel'] == 'Topography':
            print('Topo!!')
            topocount += 1



    Spatial = 1
    Spectral = 1
    return Spatial, Spectral

readCitsSm4Bin(filepath)
