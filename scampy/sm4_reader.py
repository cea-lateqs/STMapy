#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 11:49:34 2020

@author: florie

%SM4READER takes a file ID of an RHK .sm4 file, then reads

SEE => SM4 DATA FILE FORMAT V5.pdf

%Further information can be found at: 
%     http://unh2d.weebly.com/using-sm4-files-in-matlab.html
%     http://unh2d.weebly.com/sm4readerm-and-the-mat-file-format.html

"""
# -*- coding: utf-8 -*-
import struct
import sys
import codecs
import unittest
import numpy as np


filepath='/home/florie/Documents/LABO/Python/SM4_File_Format/FCg01-2-010719_0061.sm4'
def readCitsSm4Bin(filepath):
    nbrelines = 0
    error = 0
    byte = 0
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
    Signature = np.fromfile(f, dtype=np.uint16, count=18)
    Signature = list(map(lambda x: chr(x), Signature))
    Signature = ("".join(Signature))
    Total_Pagecount = int(np.fromfile(f, dtype=np.uint32, count=1)[0])
    ObjectListCount = int(np.fromfile(f, dtype=np.uint32, count=1)[0])
    ObjectFieldSize = int(np.fromfile(f, dtype=np.uint32, count=1)[0])
    Reserved = int(np.fromfile(f, dtype=np.uint32, count=1)[0])
    
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
    PageHeaderIndex = []
    for j in range(PageIndexHeader_PageCount):
        PageHeaderIndex.append({'PageID': np.fromfile(f, dtype=np.uint16, count=8)[0],
                                'PageDataType': int(np.fromfile(f, dtype=np.uint32, count=1)[0]),
                                'PageSourceType': int(np.fromfile(f, dtype=np.uint32, count=1)[0]),
                                'ObjectListCount': int(np.fromfile(f, dtype=np.uint32, count=1)[0]),
                                'MinorVersion': int(np.fromfile(f, dtype=np.uint32, count=1)[0]),
                                'ObjectList': []
                                })
        for i in range(PageHeaderIndex[j]['ObjectListCount']):
            dict1 = {'ObjectID': ObjectIDCode[int(np.fromfile(f, dtype=np.uint32, count=1)[0])],
                     'Offset': int(np.fromfile(f, dtype=np.uint32, count=1)[0]),
                     'Size': int(np.fromfile(f, dtype=np.uint32, count=1)[0])}
            PageHeaderIndex[j]['ObjectList'].append(dict1)


    # Read and record the Page Index Headers and Page Headers for each Page
    PageHeader = []
    for j in range(PageIndexHeader_PageCount):
        f.seek(PageHeaderIndex[j]['ObjectList'][0]['Offset'], 0)
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
                          'ObjectListCount': int(np.fromfile(f, dtype=np.uint32, count=1)[0])
                          })
    print(PageIndexHeader_PageCount)
    PageHeader[0]['XScale']
    
#       %Skip flags and reserved data
#       fseek(fileID,1+3+60,0);
#
#       %Read the Object List
#       for i=1:outfile.PageHeader(j).ObjectListCount
#        outfile.PageHeader(j).ObjectID(i) = fread(fileID,1,'uint32');
#        outfile.PageHeader(j).Offset(i) = fread(fileID,1,'uint32');     % 4 bytes
#        outfile.PageHeader(j).Size(i) = fread(fileID,1,'uint32');
#       end
#
#       
#       
#       %Read and record the Text Strings for each Page
#       %----------------------------------------------
#       %PageheaderObjectList
#        count=fread(fileID,1,'uint16'); %tells how long the next string is
#        outfile.TextString(j).strLabel = fread(fileID,count,'uint16=>char');%string that goes on top of plot window, such as "current image"
#        outfile.TextString(j).strLabel=transpose(outfile.TextString(j).strLabel);
#
#        count=fread(fileID,1,'uint16'); %tells how long the next string is
#        outfile.TextString(j).strSystemText = fread(fileID,count,'uint16=>char');%a comment describing the data
#        outfile.TextString(j).strSystemText=transpose(outfile.TextString(j).strSystemText);
#
#        count=fread(fileID,1,'uint16'); %tells how long the next string is
#        outfile.TextString(j).strSessionText = fread(fileID,count,'uint16=>char');%general session comments
#        outfile.TextString(j).strSessionText=transpose(outfile.TextString(j).strSessionText);
#
#        count=fread(fileID,1,'uint16'); %tells how long the next string is
#        outfile.TextString(j).strUserText = fread(fileID,count,'uint16=>char');%user comments
#        outfile.TextString(j).strUserText=transpose(outfile.TextString(j).strUserText);
#
#        count=fread(fileID,1,'uint16'); %tells how long the next string is
#        outfile.TextString(j).strPath = fread(fileID,count,'uint16=>char');%Path
#        outfile.TextString(j).strPath=transpose(outfile.TextString(j).strPath);
#
#        count=fread(fileID,1,'uint16'); %tells how long the next string is
#        outfile.TextString(j).strDate = fread(fileID,count,'uint16=>char');%DAQ date
#        outfile.TextString(j).strDate=transpose(outfile.TextString(j).strDate);
#
#        count=fread(fileID,1,'uint16'); %tells how long the next string is
#        outfile.TextString(j).strTime = fread(fileID,count,'uint16=>char');%DAQ time
#        outfile.TextString(j).strTime=transpose(outfile.TextString(j).strTime);
#
#        count=fread(fileID,1,'uint16'); %tells how long the next string is
#        outfile.TextString(j).strXUnits = fread(fileID,count,'uint16=>char');%physical units of x axis
#        outfile.TextString(j).strXUnits=transpose(outfile.TextString(j).strXUnits);
#
#        count=fread(fileID,1,'uint16'); %tells how long the next string is
#        outfile.TextString(j).strYUnits = fread(fileID,count,'uint16=>char');%physical units of Y axis
#        outfile.TextString(j).strYUnits=transpose(outfile.TextString(j).strYUnits);
#
#        count=fread(fileID,1,'uint16'); %tells how long the next string is
#        outfile.TextString(j).strZUnits = fread(fileID,count,'uint16=>char');%physical units of Z axis
#        outfile.TextString(j).strZUnits=transpose(outfile.TextString(j).strZUnits);
#
#        count=fread(fileID,1,'uint16'); %tells how long the next string is
#        outfile.TextString(j).strXLabel = fread(fileID,count,'uint16=>char');
#        outfile.TextString(j).strXLabel=transpose(outfile.TextString(j).strXLabel);
#
#        count=fread(fileID,1,'uint16'); %tells how long the next string is
#        outfile.TextString(j).strYLabel = fread(fileID,count,'uint16=>char');
#        outfile.TextString(j).strYLabel=transpose(outfile.TextString(j).strYLabel);
#
#        count=fread(fileID,1,'uint16'); %tells how long the next string is
#        outfile.TextString(j).strStatusChannelText = fread(fileID,count,'uint16=>char');%status channel text
#        outfile.TextString(j).strStatusChannelText=transpose(outfile.TextString(j).strStatusChannelText);
#
#        count=fread(fileID,1,'uint16'); %tells how long the next string is
#        outfile.TextString(j).strCompletedLineCount = fread(fileID,count,'uint16=>char');%contains last saved line count for an image data page
#        outfile.TextString(j).strCompletedLineCount=transpose(outfile.TextString(j).strCompletedLineCount);
#
#        count=fread(fileID,1,'uint16'); %tells how long the next string is
#        outfile.TextString(j).strOverSamplingCount = fread(fileID,count,'uint16=>char');%Oversampling count for image data pages
#        outfile.TextString(j).strOverSamplingCount=transpose(outfile.TextString(j).strOverSamplingCount);
#
#        count=fread(fileID,1,'uint16'); %tells how long the next string is
#        outfile.TextString(j).strSlicedVoltage = fread(fileID,count,'uint16=>char');%voltage at which the sliced image is created from the spectra page.  empty if not a sliced image
#        outfile.TextString(j).strSlicedVoltage=transpose(outfile.TextString(j).strSlicedVoltage);
#
#        count=fread(fileID,1,'uint16'); %tells how long the next string is
#        outfile.TextString(j).strPLLProStatus = fread(fileID,count,'uint16=>char');%PLLPro status text: blank, master or user
#        outfile.TextString(j).strPLLProStatus=transpose(outfile.TextString(j).strPLLProStatus);
#
#        count=fread(fileID,1,'uint16'); %tells how long the next string is
#        outfile.TextString(j).strSetpointUnit = fread(fileID,count,'uint16=>char');%ZPI controller item's set-point unit
#        outfile.TextString(j).strSetpointUnit=transpose(outfile.TextString(j).strSetpointUnit);
#
#        count=fread(fileID,1,'uint16'); %tells how long the next string is
#        outfile.TextString(j).strCHDriveValues = fread(fileID,count,'uint16=>char');%staores value of CH1 and CH2 if they are in hardware space
#        outfile.TextString(j).strCHDriveValues=transpose(outfile.TextString(j).strCHDriveValues);
#    end
    
    
#    file = f.read()
#    first_four_bytes = f.read(byte + 4)
#    print(first_four_bytes)
    #read Header 
#    Headersize = lines[1].decode('utf-8')
#    print(Headersize)
#    for line in file:
#        try:
#            L = line.decode("utf-8")
#            if 'Current' in L:
#                print ('youpi')
#            if 'LIA Current' in L:
#                print ('youpi')
#            nbrelines +=1
#        except UnicodeDecodeError : error +=1; pass
#    print(nbrelines, error)
#    print(lines[6])
#        # Header lines can be treated as regular strings
#        if header_end_not_found:
#            print("Problem while reading the file : could not find ':HEADER END:' in file")
#            f.close()
#            return False


readCitsSm4Bin(filepath)
