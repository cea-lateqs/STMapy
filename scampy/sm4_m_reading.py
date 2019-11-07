#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep  8 17:47:15 2019

@author: florie

This code is based on STSPlot.m matlab code that can be found at: 
%     http://unh2d.weebly.com/using-sm4-files-in-matlab.html
%     http://unh2d.weebly.com/sm4readerm-and-the-mat-file-format.html

It requires data which has been converted from .sm4 
to the .mat format by sm4reader.m, that will also be traduced to python later.
Meanwhile, one should run it through Octave interface. !'with '-v7' saving option!

See citswidget from scampy for more possibilities

#STSPLOT     
#     Plots STS data converted by sm4reader
#     function call: STSPLOT(yPlot,yOffset,repetitionNumber,varargin)
#
#The program is designed to aid in creating waterfall plots (with or without
#data averaging and/or smoothing), and also displays a flattened 
#topographic image (with the locations the STS data were taken 
#marked on the image, with colours coordinated with the STS spectra plots)=not yet.
#
#The STS files are either a set of point spectroscopy data (with possible 
#multiple repetitions of the measurement per file), or a single file which
#contains STS measurements at several points along a line (aka grid) on a sample
#surface (with possible repetitions).  (The associated topographical file is
#required to mark the physical locations, but can be turned off with an 
#optional flag.)=not yet
#
#yPlot is the types of plots you want to generate, in whichever order you 
#choose, such as {'dI/dV', 'Current (A)'}.
#
#yOffset is used to generate waterfall plots, the offset is used between 
#different point spectroscopy locations, or the series of locations in a 
#line for line spectroscopy.  The offsets are given as an 
#array in the same order of the plot types, e.g. [4e-12,1e-10] for two plot
#types.
#
#repetitionNumber controls which repetitions of each measurement you are
#plotting.  You can enter the following values:
#
#   an integer       --that repetition is plotted, e.g. 1, 2
#   'all'            --all repetitions per point are plotted simultaneously
#   'average'        --this plots the average of all the repetitions
#                      for each measurement location
#   'allwithaverage' --plots all repetitions, with a curve for the average 
#                      shown in black (note if the last measurement is in 
#                      black the average remains black as well)    
#
#After reptition number you can include optional arguments.  These include:
#
#*Data smoothing
#   Call by entering {'smooth',<smoothwidth>}, where <smoothwidth> is a
#   natural number which says what number neighbors you will use to average
#   at each point.  0 is no smooth, 1 is neareast neighbors only (3
#   points), 2 is two neighbors on each side (5 pts), etc..
#
#*Scan direction arrow = not yet
#   On the plotted topography data for line spectroscopy, you can include 
#   an arrow that shows what direction the plotting took place.  Include 
#   'arrow' as the argument to enable the arrow.
#
#*No spatial image
#   By including 'nospatial', the user will not be queried for a
#   topographic scan, and will not plot the topographic data.
#
#An example call using some optional fields is:
#STSplot({'Current (A)','dI/dV'},[0,0],'average', smoothwidth = 20, arrow = True )

#See also SM4READER, SM4TOMATLAB
"""
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from numbers import Number
from PyQt5 import QtCore, QtGui
import scipy.io
import tkinter as tk
from tkinter import filedialog
import Common.functions as fc
root = tk.Tk()
root.withdraw()
plt.style.use('/home/florie/Documents/LABO/Python/plot_kerelski/nature.mplstyle') #pour avoir les jolies ecritures
mpl.rcParams['contour.negative_linestyle'] = 'solid'


cmap = fc.getGwyCm()


##
#def loadwithsm4check(loadName,fileName):### NOT WORKING, WE ONLY IMPORT MAT FILES ###
        #check to see if the data is an .sm4 file, which will require reading
#        [~,~,fExtension] = fileparts(strcat(pwd,loadName));
#
#        if strcmp(fExtension,'.sm4') == 1;
#            fileID = fopen(fileName{1}); 
#            [loadedData,~] = sm4reader(fileID);
#            loadedData = load(loadName);
def Datatype (y,spectratype) :
       #Sets up a plot selection criteria, and determines the order
       datatypelabel = []
       for i in range (len(y)) :
            if y[i] == 'dI/dV':
                datatypelabel.append( 'dIdV_'+''.join(spectratype)+'_Data' );    
    
            else : 
                if y[i] == 'Current (A)' :
                    datatypelabel.append( 'IV_'+''.join(spectratype)+'_Data');
            
                else : print(str(y[i]) + ' Data type is not recognized. Remove it from STSplot arguments or check the .mat file loading/Saving')
       return datatypelabel

def STSplot(yPlot,yOffset,repetitionNumber, smoothwidth = False , arrow = False , nospatial = False ) :
    #must be initialized to 0 in order for the nospatial flag to work
    xcoords = 0;
    ycoords = 0;
    CITS = True

    #------------
    #Flags/checks
    #------------
    
    #check to see if RepetitionNumber is in fact a number
    repCheck = isinstance(repetitionNumber,Number)

    #check if repetition number is 'all'
    if (repetitionNumber == 'all') : allCheck , averageCheck ,  allwithaverageCheck = True , False , False ; 

    #check if repetition number is 'average'
    if (repetitionNumber =='average') : 
        allCheck , averageCheck ,  allwithaverageCheck = False, True, False
        if (smoothwidth == False) : print('you must choose smoothwith to get average done!!!')

    #check if repetition number is 'allwithaverage'
    if (repetitionNumber =='allwithaverage') : 
        allCheck , averageCheck ,  allwithaverageCheck = False, False, True
        if (smoothwidth == False) : print('you must choose smoothwith to get average done!!!')

    
    #------------
    #Choose files
    #------------
    
    #Prompt the user to select files
    fileSelectMessage = 'Select the file containing the STS data.';
    filename = filedialog.askopenfilename(title=fileSelectMessage)        
    
    
    #-------------------------------------
    #Determine what the user wants to plot
    #-------------------------------------
    
    #Determines how many plot types we want to do (e.g. dIdV or IV)
    numPlotTypes = len(yPlot);
    dataTypeLabel = [None]*10; #initializes a structure to hold the DataTypeLabels
        
    #Check to see what kind of data is present by loading the first file
    loadedData = scipy.io.loadmat(filename)
    Spectral = loadedData['Spectral']
    spectraType = Spectral['type'][0,0];#see scipy.org : scipy.io.loadmat    
    #I am making the assumption that the data is only point spectra or line
    #spectra, but not both
    dataTypeLabel = Datatype (yPlot,spectraType)
    
    if (nospatial == 0):
#Topographic Information
#-----------------------
        fileSelectMessage2 = 'Select the file containing the topographical data. May be the same as before.' #waits until you press okay before prompts the file select
        toposelectname = filedialog.askopenfilename(title=fileSelectMessage2);
        loadedTopo = scipy.io.loadmat(toposelectname)
        FImage,xcoords,ycoords,ImageWidth,ImageLength,fxscale,fyscale,centerx,centery = toposetup(loadedTopo); #generates the topographical information using a subfunction defined at the end of the main program

#-------------
#Point Spectra
#-------------
    if (spectraType =='Point') :       
        for PlotCount in range(numPlotTypes): #iterate over the plot types
            
            try :
                SpectralData_y = Spectral[dataTypeLabel[PlotCount]][0,0]
            except ValueError : print ( 'The file seems to contain no ' +str(dataTypeLabel[PlotCount]) + '. remove it from STSplot arguments or check the .mat file loading/Saving')
            SpectralData_x = Spectral['xdata'][0,0]*1000;#mV  
            
            repInFile = np.shape(SpectralData_y)[1];
#pull out how many repetitions are in this data by reading the matrix size (in the proper dimension) of the data folder

            #Creates a figure to plot in.
            fig,axp = plt.subplots()
            axp.set_xlabel('Voltage (mV)')
            axp.set_ylabel(yPlot[PlotCount])
                                
            
            yDataAverage = np.array(np.sum(SpectralData_y,axis=1)*1.0/repInFile);
            yDataAverage = linearsmooth(yDataAverage,smoothwidth);
            yData = linearsmooth(SpectralData_y,smoothwidth); #smoothing the data.  default is zero smoothing
       
            if repCheck == 1 : axp.plot(SpectralData_x,yData[:,repetitionNumber]);#this plots the curve
    
            
            if allCheck == 1 : #if the user entered 'all' for RepetitionNumber then we will plot every repetition for each file
                axp.plot(SpectralData_x,yData);#this plots all the curves for each data file
            
            if averageCheck == 1 : #if the user entered 'average' for RepetitionNumber then we plot the average of all plots in the file
                 axp.plot(SpectralData_x,yDataAverage);#this plots all the average of all the curves for the data file            
    
            
            if allwithaverageCheck == 1 : #if the user entered 'average' for RepetitionNumber then we plot the average of all plots in each file
                 axp.plot(SpectralData_x,yData);#this plots all the curves for each data file
                 axp.plot(SpectralData_x,yDataAverage);#this plots all the average of all the curves for the active data file
                 
                 
#------------
#Line Spectra
#------------
    #find out how many measurement locations were taken along a line :
    #each measurement is taken at a coordinate, but several measurements (repetitions) are taken on the same spot. Eg if repetitions = 4, 2 measurements, one forward, one backward (saved in right direction)
    numOfMeasurements = np.shape(Spectral['xCoord'][0,0])[0]; 
    repetitions = 0;
    pointdiff = 0;
    repindex = 0;
    while pointdiff == 0 :  #step through data points until there is a difference between the current (x,y) point and the next (x1,y1)point
        pointdiff = (Spectral['xCoord'][0,0][repindex+1]- Spectral['xCoord'][0,0][repindex]) + (Spectral['yCoord'][0,0][repindex+1]- Spectral['yCoord'][0,0][repindex]); #adds the x and y difference values together.  if this is anything other than 0, we have found the limit of the repetitions
        repetitions = repetitions+1;
        repindex = repindex+1;
    
    numberOfPlots = numOfMeasurements/repetitions; #this is the number of distinct plotting locations    
    
#We iterate over the number of plotting types (e.g. dI/dV)
    for PlotCount in range(numPlotTypes) : 
        #Spectral[dataTypeLabel[PlotCount]] = linearsmooth(Spectral[dataTypeLabel[PlotCount]],smoothwidth); #smooth the data, the default is no smoothing
        try :
            SpectralData_y = Spectral[dataTypeLabel[PlotCount]][0,0]
        except ValueError : print ( 'The file seems to contain no ' +str(dataTypeLabel[PlotCount]) + '. remove it from STSplot arguments or check the .mat file loading/Saving')
        SpectralData_x = Spectral['xdata'][0,0]*1000;#mV
        
        #Creates a figure to plot in.
        fig,axl = plt.subplots()
        axl.set_xlabel('Voltage (mV)')
        axl.set_ylabel(yPlot[PlotCount])
             
        
        #there are three possible plot types
        if repCheck == 1:    #if the user input a number as the RepetitionNumber field, we use this plot the user selected.  E.g. the second measurement at each location on the line if they had entered 2 (so repnum-1 in python).
            for m in range( (repetitionNumber-1),numOfMeasurements,repetitions):    #start at the repetition number the user wants, then step by the number of repetitions until you have plotted all the data.  E.g. There are three measurements at each of 5 points, and the user wants to plot the second.  The loop therefore is m = 2:3:15, or 2,5,8,11,14...the second measurement at each point.  [Point 1 has 1,2 and 3, Point 2 has 4, 5, and 6, etc.]
                
                axl.plot(SpectralData_x,SpectralData_y[:,m]+(int(m)/repetitions)*yOffset[PlotCount]);#this plots the curve
                    

#?                [xminindex(1,minindex), yminindex(1,minindex)] = pointsfinder(loadedData, dataTypeLabel,PlotCount,m,xcoords,ycoords);
        
        
        if allCheck == 1 :   #if the user input a number as the RepetitionNumber field, we use this plot, which plots only the repetition the user selected.  E.g. the second measurement at each location on the line if they had entered 2.
            for m in range (numOfMeasurements):
                axl.plot(SpectralData_x,SpectralData_y[:,m]+(int(m)/repetitions)*yOffset[PlotCount]);#this plots the curve
                    

#?                [xminindex(1,minindex), yminindex(1,minindex)] = pointsfinder(loadedData, dataTypeLabel,PlotCount,m,xcoords,ycoords);

            
        if averageCheck == 1  :  #if the user input a number as the RepetitionNumber field, we use this plot, which plots only the repetition the user selected.  E.g. the second measurement at each location on the line if they had entered 2.
            for m in range(0,numOfMeasurements,repetitions):    
                avedata =np.zeros(int(Spectral['points'][0,0])); #initialize the array which will store the averaged data

                for k in range( m,m+repetitions): #step through the measurements at a single location.  example of the indices.  There are 16 measurements with two repetitions at each point.  The outer loop will iterate m from 1 by 2 until 16: 1,3,5,7,9,11,13,15.  For the inner loop we see n goes from m:m+repetitions-1, in this case 1 to 1 + 2 -1, or 1 to 2.  m then goes forward to 3, and n goes from 3 to 4.
                    avedata = avedata + SpectralData_y[:,k];#as we iterate over the measurements at each point, avedata holds the sum of all the measurements
                #ok because these are arrays
                avedata = avedata*1.0/repetitions; #once we are done summing, avedata is divided by the number or repetitions, leading to the averaged data
               #see bottom avedata = linearsmooth(avedata,smoothwidth);
                axl.plot(SpectralData_x,avedata+(int(m)/repetitions)*yOffset[(PlotCount)]);#this plots the curve
                    


        if allwithaverageCheck == 1:    #if the user input a number as the RepetitionNumber field, we use this plot, which plots only the repetition the user selected.  E.g. the second measurement at each location on the line if they had entered 2.
            
            avedata =np.zeros(int(Spectral['points'][0,0])); #initialize the array which will store the averaged data

            for m in range(numOfMeasurements):    
                axl.plot(SpectralData_x,SpectralData_y[:,m]+(int(m)/repetitions)*yOffset[PlotCount]);#this plots the curve
                    
            for m in range(0,numOfMeasurements,repetitions):
                
                avedata =np.zeros(int(Spectral['points'][0,0]));
                for k in range( m,m+repetitions): #step through the measurements at a single location.  example of the indices.  There are 16 measurements with two repetitions at each point.  The outer loop will iterate m from 1 by 2 until 16: 1,3,5,7,9,11,13,15.  For the inner loop we see n goes from m:m+repetitions-1, in this case 1 to 1 + 2 -1, or 1 to 2.  m then goes forward to 3, and n goes from 3 to 4.
                    avedata = avedata + SpectralData_y[:,k];#as we iterate over the measurements at each point, avedata holds the sum of all the measurements
                
                avedata = avedata*1.0/repetitions; #once we are done summing, avedata is divided by the number or repetitions, leading to the averaged data
                #see bottom avedata = linearsmooth(avedata,smoothwidth);
                
                axl.plot(SpectralData_x,avedata+(int(m)/repetitions)*yOffset[PlotCount]);#this plots the curve

#PLOT CITS written to fit citswidget requirements (trying at least)
        fig2 = plt.figure()
        # fig.canvas.mplconnect()
        (ax2) = fig2.add_subplot(1, 1, 1, aspect='equal')
        if dataTypeLabel[PlotCount] == 'dIdV_Line_Data' and CITS == True : 
            Voltageindex = int(len(SpectralData_x)/2.) #test in the middle
            repetition_index = 1 #0,1,2... according to the number of repetitions.
            xPx = int(np.sqrt(numberOfPlots))
            yPx = int(np.sqrt(numberOfPlots))
            x_m = np.linspace(0, fxscale*1e+9, xPx)#?
            y_m = np.linspace(0, fyscale*1e+9, yPx)
            m_data = np.zeros( shape=( repetitions, yPx, xPx, int(len(SpectralData_y[:,0])) ) )
            mapData = np.ndarray(shape=(yPx, xPx))
#Spectral info is saved from right to left and up to down
            i=0
            for r in range (0,repetitions) :
                for y in range( 0 , yPx):#len(SpectralData_y[:,0])):
                    for x in range( 0, xPx):#int(numOfMeasurements/repetitions) ): #If you want bacward : columms+1 in spectralData_y
                        if i+r<4096 : m_data[r][y][x] = SpectralData_y[:,i+r]#check if int? #repetition_index
                        i+=1
            print('######')
            valMin = 10**100
            valMax = -10**(100)
            for y in range(0, yPx):
                for x in range(0, xPx):
                    val = m_data[repetition_index,y,x,Voltageindex] #repetition_index
                    mapData[y][x] = val
                    if val < valMin:
                        valMin = val
                    elif val > valMax:
                        valMax = val
                        
            print(valMin,valMax)
            print( np.shape(mapData) ) 
            XYmap = ax2.pcolormesh( mapData, cmap=cmap)#, vmin=valMin, vmax=valMax, clip_box
            ax2.axis([0, xPx, 0, yPx])
            ax2.set_title("V=" + str(SpectralData_x[Voltageindex])+'mV')
            # Image color scale is adjusted to the data:
            XYmap.set_clim(valMin,valMax)   
    
    
    if nospatial == 0 :
    #---------------
    #Plot Topography
    #---------------
    #Here we determine the optimal contrast to display the
            #topographical information.  We found a good value is to use RMS
            #+/- 5*STD for the contrast range.  These values are found here.

            imagevector = FImage[FImage!=0];#we create a vector with all the nonzero values in the image.  Actual zero values are unlikely unless there is empty data represented by zeros, so we drop them.
            imagerms = np.sqrt(np.sum(np.square(imagevector))/(1.0*len(imagevector))); #calculate the RMS
            imagestd=np.std(imagevector,ddof=1); #calculate the standard deviation. degree of freedom: N-1 instead of N


            #Plot the topography data
            #generate the image figure
            fig1,ax = plt.subplots()
            ax.set_xlabel('nm')
            ax.set_ylabel('nm')


            #show the topography with the contrast determined by the range of RMS +/- 5*STD
            ax.imshow(FImage,extent=[0,fxscale*1e+9,0,fyscale*1e+9], vmin = (imagerms - 5*imagestd) , vmax = (imagerms + 5*imagestd));#show the topography with the contrast determined by the range of RMS +/- 5*STD
            #mesh(imagex,imagey,FImage)
            
            
            #We plot our spec location points.
            locations = True
            if locations == True:
                patch = []
                for m in range (0,numOfMeasurements,repetitions): #iterate over the number of locations
                    print (m)
                    patch.append(Circle(((-centerx + fxscale/2 + Spectral['xCoord'][0,0][m])*1e+9,(- centery + fyscale/2 + Spectral['yCoord'][0,0][m])*1e+9),0.05,facecolor='r',edgecolor='None'));
                for n in range(len(patch)):
                    ax.add_patch(patch[n])

            #this is optional, to include use the flag 'arrow'; only works with
            #line spectra
#            if arrow == 1 and lineTrue == 1 :
#                #make an arrow showing direction, we want to do this first so that
#                #the location points sit on top of the arrow
#                arrowx = xminindex(1,1);    #starting x coordinate of the arrow is the first xminindex
#                arrowy = yminindex(1,1);    #starting y coordinate of the arrow is the first yminindex
#                arrowdx = xminindex(1,numberOfPlots)-xminindex(1,1);     #the dx is the difference between the first and last location points
#                arrowdy = yminindex(1,numberOfPlots)-yminindex(1,1);     #the dy is the difference between the first and last location points
#                quiver(arrowx, arrowy, arrowdx, arrowdy, 0, 'LineWidth',3,'MaxHeadSize',0.4,'Color',[0 0 0]); #this draws the arrow and, importantly, the field with 0 turns off the autoscaling so that the length stays just as we want it.
#
#                #make an arrow showing direction, we want to do this first so that
#                #the location points sit on top of the arrow
#                arrowx = xminindex(1,1);    #starting x coordinate of the arrow is the first xminindex
#                arrowy = yminindex(1,1);    #starting y coordinate of the arrow is the first yminindex
#                arrowdx = xminindex(1,numberOfPlots)-xminindex(1,1);     #the dx is the difference between the first and last location points
#                arrowdy = yminindex(1,numberOfPlots)-yminindex(1,1);     #the dy is the difference between the first and last location points
#                quiver(arrowx, arrowy, arrowdx, arrowdy, 0, 'LineWidth',2,'MaxHeadSize',0.3,'Color',[.4 .4 1]); #this draws the arrow and, importantly, the field with 0 turns off the autoscaling so that the length stays just as we want it.


##
def toposetup(filename):
    #flattening
    #Load the image file we want 
    Spatial = filename['Spatial']
    Image = Spatial['TopoData'][0,0];
    FImage = Image[0,0]# Image forward
    BImage = Image[0,1]#Image bacward 
    
    ImageWidth=Spatial['points'][0,0]; #gets the width of the image
    ImageLength=Spatial['lines'][0,0]; #length
    fxscale = Spatial['width'][0,0];
    fyscale = Spatial['height'][0,0];
    
    #Performs the x line-by-line flatten first
    for n in range(int(ImageLength)):
       #create an array of x-data
       x=np.arange(1,int(ImageWidth)+1)*1.0;
       #extract the y data from row n
       y=FImage[n];

       #fit this data with a linear fit
       p =np.polyfit(x,y,1);
       linecorrect =x*p[0]+p[1];
       FImage[n]=FImage[n]-linecorrect;
   
    #y correction is not line-by-line but averaged over y direction
    averageY = np.sum(FImage,axis=1);

    #I don't want unfinished portions of the scan with all zero values to
    #influence the scan, so I will remove them.
    averageY[averageY ==0] =np.nan; #replace zeros with NaN so I can use isnan

    #create an array of x-data
    x=np.arange(1,len(averageY)+1);#np.transpose?
    #extract the y data from row n
    y=averageY;
    include = ~(np.isnan(y)); #include indexes which points are NOT NaN
    p =np.polyfit(x[include],y[include],1);
    linecorrect =x*p[0]+p[1];
    for n in range(int(ImageWidth)) :
       FImage[:,n]=FImage[:,n]-linecorrect;

    #set the initially zero (empty) values back to zero
    FImage[~include,:] =0;
    print(np.shape(FImage))
    #calculate the x and y vectors for the image in the global coordinate
    #system of the STM
    centerx = Spatial['xoffset'][0,0];
    centery = Spatial['yoffset'][0,0];
    width = Spatial['width'][0,0];
    height = Spatial['height'][0,0];
    xpoints = Spatial['points'][0,0];
    ypoints = Spatial['lines'][0,0];
    
    #initialize arrays to hold spatial information
    xcoords = np.zeros(int(xpoints));
    ycoords = np.zeros(int(ypoints));
    
    
    #generate the spatial information for x
    for xpointnumber in range (int(xpoints)):
        xcoords[xpointnumber] = centerx - width/2 + width/(xpoints-1)*xpointnumber;
    
    #generate the spatial information for y, remember, as the y index goes
    #up, the number should move in the negative direction, this is the
    #reason for the inverted signs versus the x calculation.
    for ypointnumber in range(int(ypoints)):
        ycoords[ypointnumber] = centery + height/2 - height/(ypoints-1)*ypointnumber;

    return FImage,xcoords,ycoords,ImageWidth,ImageLength,fxscale,fyscale,centerx,centery


###
##this function finds the x and y coordinates of the spatial scan which most
##closely match the x and y coordinates of the spectra by subtracting the x
##and y coordinates from the x and y spatial coordinate vectors and finding
##the minimum difference to identify which is closest
#def pointsfinder(loadeddata, ~,~,m,xcoords,ycoords):
#    Spectral = loadeddata['Spectral']
#    xfinder = xcoords - Spectral['xCoord'][0,0][m];
#    [~,xmincoordinate]=min(abs(xfinder));
#
#    yfinder = ycoords - Spectral[.yCoord(m);
#    [~,ymincoordinate]=min(abs(yfinder));
#    return xmincoordinate, ymincoordinate

##
def linearsmooth(data,smoothwidth):
    smoothwidth = int(smoothwidth)
    #Part One: Gathering Information
    if smoothwidth == 0 : smootheddata=data;
    else :
    #find out how long a dataset is, and how many data sets there are
        numberofdatasets,datalength = np.shape(data);

    #Part Two: Create "padded" datasets, to average at the edges even though the dataset ends.

    #initialize the padded data array. E.g. if our smoothwidth is 2, we need 
    #two extra points before the first point (and two after the last) in order to average.
        paddeddata = np.pad(data, [(0,0),(smoothwidth,smoothwidth)], 'edge');

    #initalize the smoothed data array with same size as the original data sets
        smootheddata = np.zeros((numberofdatasets,datalength));


    #Part Three: Smoothing
        #step through the points within paddeddata, averaging over their neighbors, saving into smootheddata.
        counter = 0;
        for a in range ((smoothwidth),(smoothwidth+datalength)):
            smootheddata[counter] = np.sum(paddeddata[(a-smoothwidth):(a+smoothwidth)]) / (1+2*smoothwidth);
            counter = counter+1;
    
    return smootheddata

STSplot(['Current (A)','dI/dV'],[0,0],'all', smoothwidth = 2 , arrow = False , nospatial = False)