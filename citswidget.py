# -*- coding: utf-8 -*-
"""
Created on Mon Aug 03 12:04:19 2015

@author: LH242250
"""
from PyQt4.QtGui import *
from PyQt4.QtCore import *
from ui_citswidget import Ui_CitsWidget
from matplotlib.backends.backend_qt4 import NavigationToolbar2QT as NavigationToolbar
import numpy as np
import pylab
import os.path as osp
from os import listdir
import scipy as sp
import scipy.optimize
import scipy.signal
import matplotlib.patches as patches
import matplotlib.pyplot
import struct

class CitsWidget(QMainWindow, Ui_CitsWidget):
### Building methods
    def __init__(self,parent):
        """ Builds the widget with parent widget in arg """
        QMainWindow.__init__(self,parent)
        #Set up the user interface from Designer.
        self.setupUi(self)
        self.setAutoFillBackground(True)
        self.m_data=np.ndarray([])
        self.tot_data=np.ndarray([])
        self.nAvgSpectra=0
        self.m_params={}
        self.channelList=[]
        self.mapName=""
        #Set up figures
        self.toolbar_map = NavigationToolbar(self.m_mapWidget,self)
        self.toolbar_spec = NavigationToolbar(self.m_specWidget,self)
        self.map_layout.insertWidget(0,self.toolbar_map)
        self.spec_layout.insertWidget(0,self.toolbar_spec)
        self.fig_spec=self.m_specWidget.figure
        self.ax_spec=self.fig_spec.add_subplot(1,1,1)
        self.ax_spec.hold(True)
        self.fig_spec.subplots_adjust(left=0.125,right=0.95,bottom=0.15,top=0.92)
        #Variables linked to clicks on map
        self.pts_clicked=[[],[],[]]
        self.origin_x=0
        self.origin_y=0
        self.voltageLine=0
        #Colormaps
        self.m_colorBarBox.addItems(matplotlib.pyplot.colormaps())
        self.m_colorBarBox.setCurrentIndex(2)
        #Boolean that is True if a map is loaded
        self.dataLoaded=False
        self.wdir="E:\\PhD\\Experiments\\STM\\Lc12\\CITS\\"
        self.connect()
        #self.readTopo("H:\\Experiments\\STM\\Lc0 (Si-260)\\4K\\Spectro ASCII\\TestZ.txt")
        self.nSpectraDrawn=0
        self.spectrumColor=['b', 'g', 'r', 'c', 'm', 'y', 'k']
        #Set up layouts
        self.m_avgWidget.hide()
        self.m_cbarWidget.hide()
    
    def connect(self):
        """ Connects all the signals. Only called in the constructor """
        self.m_openButton.clicked.connect(self.loadCITS)
        self.m_channelBox.currentIndexChanged.connect(self.updateMap)
        self.m_colorBarBox.currentIndexChanged.connect(self.updateMap)
        self.m_voltageBox.valueChanged.connect(self.drawXYMap)
        self.m_mapWidget.mpl_connect('button_press_event',self.onPressOnMap)
        self.m_mapWidget.mpl_connect('button_release_event',self.onReleaseOnMap)
        self.m_specWidget.mpl_connect('button_press_event',self.updateToPointedX)
        self.m_clearSpec.clicked.connect(self.clearSpectrum)
        self.m_scaleVoltage.toggled.connect(self.clearSpectrum)
        #self.m_scaleMetric.toggled.connect(self.updateMap)
        self.m_averageCitsButton.clicked.connect(self.averageCITS)
        self.m_wholeLengthCutButton.clicked.connect(self.launchBigCut)
        self.m_magicButton.clicked.connect(self.magicFunction)
        self.m_avgSpec.clicked.connect(self.launchAvgSpectrum)
        self.m_fitSpec.clicked.connect(self.fitSpectrum)
        self.m_fitCustomCheckbox.stateChanged.connect(self.m_fitUpperBox.setEnabled)
        self.m_fitCustomCheckbox.stateChanged.connect(self.m_fitLowerBox.setEnabled)
        self.m_avgBox.toggled.connect(self.updateAvgVariables)
        self.m_vLineBox.toggled.connect(self.clearVoltageLine)
        self.m_avgVButton.clicked.connect(self.averageSpectrumWithValues)
        #Averaging with respect to values
        self.m_avgCheckBox.toggled.connect(self.m_avgWidget.setVisible)
        self.m_aboveBox.valueChanged.connect(self.updateAboveValue)
        self.m_belowBox.valueChanged.connect(self.updateBelowValue)
        #Cbar custom limits
        self.m_cbarCheckBox.toggled.connect(self.m_cbarWidget.setVisible)
        self.m_cbarCustomCheckbox.stateChanged.connect(self.m_cbarUpperBox.setEnabled)
        self.m_cbarCustomCheckbox.stateChanged.connect(self.m_cbarLowerBox.setEnabled)
        
        
### Reading and loading CITS methods
        
    def loadCITS(self):
        """ Slot that launches the reading of the CITS by asking the path of the file """
        filename=QFileDialog.getOpenFileName(self,"Choose a CITS file",self.wdir,"3D binary file (*.3ds);;Ascii file (*.asc);;Text file (*.txt)")
        extension=filename.split('.')[-1]        
        if(extension=="asc"):
            self.clearMap()
            self.dataLoaded=self.readCitsAscii(filename)
        elif(extension=="3ds"):
            self.clearMap()
            self.dataLoaded=self.readCitsBin(filename)
        elif(extension=="txt"):
            self.readTopo(filename)
        else:
            print("Extension not recognized")
            return
        if(self.dataLoaded):
            self.wdir=osp.dirname(filename)
            self.mapName=osp.basename(filename)
            print(self.mapName)
            self.updateWidgets()
        
    def averageCITS(self):
        """ Slot that averages the CITS chosen in the dialog box (they have to be of the same dimensions) """
        Cits_names=QFileDialog.getOpenFileNames(self,"Choose the CITS files to average","E:\\PhD\\Experiments\\STM\\Lc12\\CITS","3D binary file (*.3ds);;Ascii file (*.asc);;Text file (*.txt)")
        N_Cits=len(Cits_names)
        first=True
        for cits in Cits_names:
            extension=cits.split('.')[-1]  
            if(extension=="asc" or extension=="txt"):
                self.dataLoaded=self.readCitsAscii(cits)
            elif(extension=="3ds"):
                self.dataLoaded=self.readCitsBin(cits)
            else:
                print("Extension not recognized")
                return
            #Verify that the data was loaded. Returns otherwise
            if(not self.dataLoaded):
                return
            #Create mean_data array if first cits. Add to it if not
            if(first):
                mean_data=self.m_data/N_Cits
                first=False
                nPts=(cits.split('.')[-2]).split('_')[-2]
            else:
                mean_data+=self.m_data/N_Cits
        #If everthing went well, I clear the map and set the data to mean_data
        self.mapName="Average_"+nPts
        self.clearMap()
        self.dataLoaded=True
        self.m_data=mean_data
        self.updateWidgets()
        
        
    def readCitsAscii(self,filepath):
        """ Reads an Ascii CITS file (Omicron) and stores all the parameters"""
        f=open(filepath)
        divider=1
        while(True):
            #Read the parameters of the map until "Start of Data"
            line=f.readline()
            #Pixel dimensions in X
            if("x-pixels" in line):
                xPx=int(line.split()[-1])
            #Pixel dimensions in Y
            elif("y-pixels" in line):
                yPx=int(line.split()[-1])
            #Metric dimensions in X
            elif("x-length" in line):
                xL=float(line.split()[-1])
            #Metric dimensions in Y
            elif("y-length" in line):
                yL=float(line.split()[-1])
            #Number of points IN TOTAL therefore the division per 2 to have the number of points per channel
            elif("z-points" in line):
                zPt=int(line.split()[-1])
                #There is zPt/2 points for fwd and zPt/2 points for bwd
                zPt=zPt/2
            #Starting voltage of the spectro
            elif("Device_1_Start" in line):
                vStart=round(float(line.split()[-2]),6)
            #Final voltage of the spectro
            elif("Device_1_End" in line):
                vEnd=round(float(line.split()[-2]),6)
            #Any eventual divider
            elif("divider" in line):
                divider=int(line.split()[-1])
            elif("Start of Data" in line):
                break
        # Matlab convention : columns first then rows hence [y][x]
        # In Omicron CITS, there is only two channels : fwd and bwd so it is read as such
        self.channelList=("Data [Fwd]","Data [Bwd]")
        self.m_data=np.zeros(shape=(2,yPx,xPx,zPt))
        for y in range(0,yPx):
            for x in range(0,xPx):
                #The line read is an array containing the dI/dV (or I(V)) values indexed by the voltage index
                #Strip to remove the newline at the end and split to transform the string in a list
                data_list=f.readline().strip().split()
                #Forward data
                self.m_data[0][y][x]=data_list[0:zPt]
                #No need to reverse the backward data as it was from Vmin to Vmax in the file as the fwd data
                #Backward data
                self.m_data[1][y][x]=(data_list[zPt:2*zPt])
        f.close()
        #Store the parameters in a dictonnary to use them later
        self.m_params={"xPx":xPx,"yPx":yPx,"xL":xL,"yL":yL,"zPt":zPt,"vStart":vStart/divider,"vEnd":vEnd/divider,"dV":abs(vEnd-vStart)/(divider*zPt)}
        if(divider!=1):
            print "A divider of "+str(divider)+" was found and applied"
        return True
        
    def readCitsBin(self,filepath):
        """ Reads a binary CITS file (Nanonis) and stores all the parameters"""
        #The divider is already taken into account by Nanonis during the experiment so no need to process it again*
        f=open(filepath,"rb")
        zSpectro=False
        half=False
        while(True):
            #Read the header of the map until its end ("HEADER_END")
            line=f.readline()
            #Pixel dimensions
            if("Grid dim" in line):
                splitted_line=line.split('"')[1].split()
                xPx=int(splitted_line[0])
                yPx=int(splitted_line[-1])
            #Center coordinates and metric dimensions in nm (Grid settings also contains other data)
            elif("Grid settings" in line):
                xC=float(line.split(";")[0].split("=")[-1])*(10**9)
                yC=float(line.split(";")[1])*(10**9)
                xL=float(line.split(";")[-3])*(10**9)
                yL=float(line.split(";")[-2])*(10**9)
            elif("Sweep Signal" in line):
                if(line.split('"')[1]=='Z (m)'):
                    zSpectro=True
            #Number of points per channel
            elif("Points" in line):
                zPt=int(line.split('=')[-1])
            #Channels recorded
            elif("Channels" in line):
                self.channelList=line.split('"')[1].split(';')
                nChannels=len(self.channelList)
            #Experiment parameters. Not used for now, only the number is recorded to skip the corresponding bytes afterwards
            elif("Experiment parameters" in line):
                nbExpParams=len(line.split(';'))
            #End of the header
            elif(":HEADER_END:" in line):
                break
            
        #Reading vStart and vEnd (floats of 4 bytes each)
        try:
            reading=struct.unpack('>'+'f'*2,f.read(8))
        except struct.error:
            print("Problem while reading the file : number of bytes to read different than what was expected")
            return False
        #If it is a Z-Spectroscopy, put the Z boundaries in nm
        if(zSpectro):
            vStart=round(reading[0]*10**9,6)
            vEnd=round(reading[1]*10**9,6)
        else:
            vStart=round(reading[0],6)
            vEnd=round(reading[1],6)
        #Reading experiment parameters (nbExpParams*4 bytes)
        #f.read(nbExpParams*4)
        # Matlab convention : columns first then rows hence [y][x]
        try:
            self.topo=np.zeros(shape=(yPx,xPx))
            self.m_data=np.zeros(shape=(nChannels,yPx,xPx,zPt))
        except MemoryError:
            print("The data is too big ! Or the memory too small...\nI will take half of the channels...\n")
            half=True
        #If the first alloc didn't work, try to halve it
        if(half):
            try:
                self.topo=np.zeros(shape=(yPx,xPx))
                self.m_data=np.zeros(shape=(nChannels/2,yPx,xPx,zPt))
            except MemoryError:
                print("The data is REALLY too big ! Or the memory REALLY too small...\nI give up...\n")
                return False
        #Format string for unpacking zPt big-endians floats ('>f')
        fmtString='>'+'f'*zPt
        #zPt floats to read of 4 bytes each
        bytesToRead=4*zPt
        for y in range(0,yPx):
            for x in range(0,xPx):
                b=f.read(nbExpParams*4)
                try:
                    self.topo[y][x]=struct.unpack('>'+'f'*nbExpParams,b)[2]
                except struct.error:
                    print("Problem while reading the topo : number of bytes to read different than what was expected")
                    return False
                for chan in range(0,nChannels):
                # Each channel is written successively by sequences of 4*zPt bytes. I then read these bytes and unpack them as big-endians floats ('>f')
                    b=f.read(bytesToRead)
                    try:
                        if(not half or chan<nChannels/2):
                            self.m_data[chan][y][x]=struct.unpack(fmtString,b)
                    except struct.error:
                        print("Problem while reading the file : number of bytes to read different than what was expected")
                        return False
                #After each loop over channels, a new "experiment" begins so I need to skip the vStart, vEnd and experiments parameters floats that are written once again before the channels
                f.read(8)                
                #f.read(8+nbExpParams*4)   
        f.close()
        if(half):        
            self.channelList=self.channelList[0:nChannels/2]
        #Store the parameters in a dictonnary to use them later
        self.m_params={"xPx":xPx,"yPx":yPx,"xC":xC,"yC":yC,"xL":xL,"yL":yL,"zPt":zPt,"vStart":vStart,"vEnd":vEnd,"dV":(vEnd-vStart)/(zPt)}
        print(self.m_params)
        #Test to have currents in nA
        i=0
        for i in range(0,len(self.channelList)):
            chan=self.channelList[i]
            if("(A)" in chan):
                self.m_data[i]=np.abs(self.m_data[i])*10**9
                self.channelList[i]=chan.replace("(A)","(nA)")
        #Test
        if(zSpectro):
            self.extractSlope(0.01,0)
            #self.extractSlope(0.01,1)
        #Post-processing
        #self.extractOutOfPhase(1,2)
        #Topo : convert in nm and draw it
        self.topo=self.topo*10**9
        self.drawTopo()
        #pylab.close()
        return True
        
### Reading and loading topo images methods (not finished yet)
        
    def readTopo(self,filepath):
        """ Reads a topography file """
        f=open(filepath)
        topo_data=[]
        w=0
        h=0
        for line in f:
            #Treat the headers differently
            if(line[0]=='#'):
                if("Width" in line):
                    w=float(line.split()[-2])
                if("Height" in line):
                    h=float(line.split()[-2])
            else:
                topo_data.append(line.strip().split())
        f.close()
        #Set up figure
        xPx=len(topo_data[0])
        yPx=len(topo_data)
        fig=pylab.figure()
        self.ax_topo=fig.add_subplot(1,1,1,aspect=float(yPx)/xPx)
        fig.subplots_adjust(left=0.125,right=0.95,bottom=0.15,top=0.92)
        data=(np.asfarray(topo_data))*10**9
        if(w!=0 and h!=0):
            dx=w/xPx
            dy=h/yPx
            self.ax_topo.axis([0,w,0,h])
            self.ax_topo.invert_yaxis()
            # PLot data (y first (matlab convention)
            XYmap=self.ax_topo.pcolormesh(np.arange(0,w+dx,dx),np.arange(0,h+dy,dy),data,cmap=self.m_colorBarBox.currentText())
        else:
            self.ax_topo.axis([0,xPx,0,yPx])
            #self.ax_topo.invert_yaxis()
            XYmap=self.ax_topo.pcolormesh(data,cmap=self.m_colorBarBox.currentText())
        #Colorbar stuff
        cbar = fig.colorbar(XYmap, shrink=.9, pad=0.05, aspect=15)
        cbar.ax.yaxis.set_ticks_position('both')
        cbar.ax.tick_params(axis='y', direction='in')
        #self.addChannel(data,"Z (nm)")
        return False
        
    def drawTopo(self):
        """ Draws the topography read while opening the CITS """
        #If yPx==1, it is a Line Spectro so I need to call the specific method to draw the topo
        yPx=self.m_params["yPx"]
        if(yPx==1):             
            self.drawLineTopo()
            return
        lineFit=True
        w=self.m_params["xL"]
        h=self.m_params["yL"]
        #Set up figure
        xPx=self.m_params["xPx"]
        
        fig=pylab.figure()
        self.ax_topo=fig.add_subplot(1,1,1,aspect=float(yPx)/xPx)
        fig.subplots_adjust(left=0.125,right=0.95,bottom=0.15,top=0.92)
        #Numpy array to save the leveled topo
        topo_leveled=np.zeros(shape=(yPx,xPx))
        #Line fiting
        if(lineFit):
            fig.suptitle("Raw topo data")
            fitX=np.arange(0,xPx)
            fitF = lambda z,a,b: a*z+b
            #Figure to plot the leveled data
            fig2=pylab.figure()
            axes2=fig2.add_subplot(1,1,1,aspect=float(yPx)/xPx)
            fig2.suptitle("With line leveling")
            fig2.subplots_adjust(left=0.125,right=0.95,bottom=0.15,top=0.92)   
            for y in range(0,yPx):            
                fitY=self.topo[y]
                f = sp.interpolate.InterpolatedUnivariateSpline(fitX,fitY, k=1)
                popt, pcov = sp.optimize.curve_fit(fitF, fitX, f(fitX))
                topo_leveled[y]=fitY-(popt[0]*fitX+popt[1])
        
        colormap="afmhot"#self.m_colorBarBox.currentText()
        if(self.m_scaleMetric.isChecked()):
            self.ax_topo.axis([0,w,0,h])
            XYmap=self.ax_topo.pcolormesh(np.linspace(0,w,xPx),np.linspace(0,h,yPx),self.topo,cmap=colormap)
            if(lineFit):
                axes2.axis([0,w,0,h])
                XYmap2=axes2.pcolormesh(np.linspace(0,w,xPx),np.linspace(0,h,yPx),topo_leveled,cmap=colormap)
        else:
            self.ax_topo.axis([0,xPx,0,yPx])
            XYmap=self.ax_topo.pcolormesh(self.topo,cmap=colormap)
            if(lineFit):
                axes2.axis([0,xPx,0,yPx])
                XYmap2=axes2.pcolormesh(topo_leveled,cmap=colormap)
        #Colorbar stuff
        cbar = fig.colorbar(XYmap, shrink=.9, pad=0.05, aspect=15)
        cbar.ax.yaxis.set_ticks_position('both')
        cbar.ax.tick_params(axis='y', direction='in')
        if(lineFit):
            cbar2 = fig2.colorbar(XYmap2, shrink=.9, pad=0.05, aspect=15)
            cbar2.ax.yaxis.set_ticks_position('both')
            cbar2.ax.tick_params(axis='y', direction='in')
        #Keep the fig in memory
        if(lineFit):
            self.fig_topo=fig2
        
    def drawLineTopo(self):
        """ Draws the topography read while opening a Line Spectro """
        lineFit=True
        #Get parameters
        w=self.m_params["xL"]
        # h is 0 (line spectro) 
        xPx=self.m_params["xPx"]
        # yPx is 1 (line spectro)
        
        #Set up figure
        if(self.m_vLineBox.isChecked()):
            fig=pylab.figure()
            self.ax_topo=fig.add_subplot(1,1,1)
            self.fig_topo=fig
        
        #Linear fit of the topo
        fitX=np.arange(0,xPx)
        fitF = lambda z,a,b: a*z+b
        #Line spectro so there is only data for y=0
        y=0
        fitY=self.topo[y]
        f = sp.interpolate.InterpolatedUnivariateSpline(fitX,fitY, k=1)
        popt, pcov = sp.optimize.curve_fit(fitF, fitX, f(fitX))
        topo_leveled=fitY-(popt[0]*fitX+popt[1])
        
        if(self.m_scaleMetric.isChecked()):        
            self.ax_topo.plot(np.linspace(0,w,xPx),topo_leveled,label="With line leveling")
            #self.ax_topo.plot(np.linspace(0,w,xPx),self.topo[y]-popt[1],label="Without line leveling")
            self.ax_topo.set_xlim(0,w)
            self.ax_topo.set_xlabel("Distance (nm)")
            self.ax_topo.set_ylabel("Z (nm)")
        else:
            self.ax_topo.set_xlim(0,xPx)
            self.ax_topo.plot(topo_leveled,label="With line leveling")
            self.ax_topo.plot(self.topo[y],label="Without line leveling")
            self.ax_topo.set_ylabel("Z (nm)")
        pylab.legend(loc=0)
        
### Updating methods. Usually called by signals
    def updateAvgVariables(self):
        """ Slot called by the toggling of the average box. Toggling on the box clears everything to start a new averaging. Toggling off averages and plots the data stored by picking spectra """
        if(self.dataLoaded):
            #If toggled on, put everything to zero to be ready to store data
            if(self.m_avgBox.isChecked()):
                self.tot_data=np.zeros(shape=self.m_params["zPt"])
                self.nAvgSpectra=0
            #If toggled off, plot the data stored
            else:
                if(self.nAvgSpectra==0): return
                dataToPlot=self.tot_data/(self.nAvgSpectra)
                self.drawSpectrum(dataToPlot,str(self.nAvgSpectra)+" spectra averaged")
                
    def updateAboveValue(self,value):
        self.m_aboveBar.setValue(value)
        if(self.dataLoaded): self.m_aboveLine.setText(str(self.mapMax-value*(self.mapMax-self.mapMin)/100))
        
    def updateBelowValue(self,value):
        self.m_belowBar.setValue(value)
        if(self.dataLoaded): self.m_belowLine.setText(str(self.mapMin+value*(self.mapMax-self.mapMin)/100))
        
    def updateMap(self):
        self.drawXYMap(self.m_voltageBox.value())
        self.updateAboveValue(self.m_aboveBar.value())
        self.updateBelowValue(self.m_belowBar.value())         
                
    def updateToPointedX(self,event):
        """ Slot called when clicking on the spectrum window. Updates the map according to the position of the cursor when clicked """
        if(event.xdata!=None and event.ydata!=None and self.dataLoaded and self.toolbar_spec._active is None):
            #If the scale is in volts, need to divide by dV to have the index
            if(self.m_scaleVoltage.isChecked()):
                pointedIndex=int((event.xdata-self.m_params["vStart"])/self.m_params["dV"])
            else:
                pointedIndex=int(event.xdata)
            #Update the box
            self.m_voltageBox.setValue(pointedIndex)
            #The map updates itself because the change of value of the voltage box calls the drawXYMap method
                
    def updateVoltageBox(self,Vmin,Vmax,zPt):
        """ Method called by updateWidgets. Sets the values of the voltage box """
        #self.m_voltageBox.setMinimum(Vmin)
        #self.m_voltageBox.setMaximum(Vmax)
        #self.m_voltageBox.setSingleStep(dV)
        self.m_voltageBox.setMinimum(0)
        self.m_voltageBox.setMaximum(zPt-1)
        self.m_voltageBox.setSingleStep(1)
        
    def updateWidgets(self):
        """ Slot called after the reading of the CITS. Sets the values combo box (voltage and channels) and draws the map """
        self.m_channelBox.clear()
        self.m_channelBox.addItems(self.channelList)
        self.drawXYMap(self.m_voltageBox.value())
        self.updateVoltageBox(self.m_params["vStart"],self.m_params["vEnd"],self.m_params["zPt"])
        self.updateAboveValue(self.m_aboveBar.value())
        self.updateBelowValue(self.m_belowBar.value())
        
### Methods related to spectra    

    def averageSpectrum(self,xi,xf,yi,yf):
        """ Averages the spectra contained in the rectangle drawn by the points (xi,yi);(xf,yi);(xf,yf) and (xi,yf) """
        if(self.dataLoaded):
            zPt=self.m_params["zPt"]
            chan=self.m_channelBox.currentIndex()
            avg_data=np.zeros(shape=zPt)
            if(yf<yi):
                t=yf
                yf=yi
                yi=t
            if(xf<xi):
                t=xf
                xf=xi
                xi=t
            N=(yf-yi)*(xf-xi) #Number of spectra averaged
            for y in range(yi,yf):
                for x in range(xi,xf):
                        avg_data+=(self.m_data[chan][y][x]/N)
            self.m_mapWidget.draw()
            self.drawSpectrum(avg_data,"Average")
            
    def averageSpectrumWithValues(self):
        """ Averages the spectra according their values at a certain voltage """
        if(self.dataLoaded):
            voltage=self.m_voltageBox.value()
            chan=self.m_channelBox.currentIndex()
            viewSelected=self.m_viewSelectedBox.isChecked()
            zPt=self.m_params["zPt"]
            avg_data_aboveV=np.zeros(shape=zPt)
            avg_data_belowV=np.zeros(shape=zPt)
            #midpoint=(self.mapMax+self.mapMin)/2
            midpoint2=(self.mapMax-self.mapMin)
            limit_aboveV=self.mapMax-self.m_aboveBar.value()*midpoint2/100
            limit_belowV=self.mapMin+self.m_belowBar.value()*midpoint2/100
            if(limit_aboveV<limit_belowV):
                print("Above and below spectra intersect !")
                return
            N_aboveV=0
            N_belowV=0
            xPx=self.m_params["xPx"]
            yPx=self.m_params["yPx"]
            for y in range(0,yPx):
                for x in range(0,xPx):
                    currentValue=self.m_data[chan][y][x][voltage]
                    if(currentValue>limit_aboveV):
                        avg_data_aboveV+=self.m_data[chan][y][x]
                        N_aboveV+=1
                        if(viewSelected): self.addToPtsClicked(x,y,self.getSpectrumColor(self.nSpectraDrawn))
                    elif(currentValue<limit_belowV):
                        avg_data_belowV+=self.m_data[chan][y][x]
                        N_belowV+=1
                        if(viewSelected): self.addToPtsClicked(x,y,self.getSpectrumColor(self.nSpectraDrawn+1))
            if(N_aboveV!=0): 
                avg_data_aboveV/=N_aboveV
                self.drawSpectrum(avg_data_aboveV,"Average above "+str(limit_aboveV)+" ("+str(N_aboveV)+" spectra averaged)")
            if(N_belowV!=0): 
                avg_data_belowV/=N_belowV
                self.drawSpectrum(avg_data_belowV,"Average below "+str(limit_belowV)+" ("+str(N_belowV)+" spectra averaged)")
            self.drawPtsClicked()
            
    def clearSpectrum(self):
        """ Clears the spectrum window """
        self.ax_spec.clear()
        self.nSpectraDrawn=0
        self.voltageLine=0
        self.clearPtsClicked()
        self.m_specWidget.draw()
            
    def drawSpectrum(self,dataToPlot,label):
        """ Method called each time a spectrum needs to be plotted. Takes care of the derivative and different scales stuff and updates the window """
        finalLabel=label+" - "+str(self.m_channelBox.currentText())
        shift=self.nSpectraDrawn*self.m_shiftBox.value()
        if(self.dataLoaded and dataToPlot.size!=0):
            dV=self.m_params["dV"]
            deriv=sp.signal.savgol_filter(dataToPlot, 5, 2, deriv=1, delta=dV)
            self.lastSpectrum=dataToPlot
            if(self.m_scaleVoltage.isChecked()):
                vStart=self.m_params["vStart"]
                vEnd=self.m_params["vEnd"]
                self.ax_spec.plot(np.arange(vStart,vEnd,dV),shift+dataToPlot,label=finalLabel)
                self.nSpectraDrawn=self.nSpectraDrawn+1
                if(self.m_derivBox.isChecked()):
                    self.ax_spec.plot(np.arange(vStart,vEnd,dV),deriv,label="Derivative of "+finalLabel)
                    self.nSpectraDrawn=self.nSpectraDrawn+1
            else:
                self.ax_spec.plot(shift+dataToPlot,label=finalLabel)
                self.nSpectraDrawn=self.nSpectraDrawn+1
                if(self.m_derivBox.isChecked()):
                    self.ax_spec.plot(deriv,label="Derivative of "+finalLabel)
                    self.nSpectraDrawn=self.nSpectraDrawn+1
        self.ax_spec.legend(loc=0)
        self.m_specWidget.draw()
        
    def fitSpectrum(self):
        if(self.dataLoaded and self.lastSpectrum.size!=0):
            fit_func = lambda v,a,b: a*v+b
            dV=self.m_params["dV"]
            if(self.m_fitCustomCheckbox.isChecked()):
                limitL=float(self.m_fitLowerBox.text())
                limitU=float(self.m_fitUpperBox.text())
                vStart=min(limitL,limitU)
                vEnd=max(limitL,limitU)    
                zPt=self.m_params["zPt"]
                xArray=np.arange(0,zPt)*dV+self.m_params["vStart"]
                #Take the portion to fit
                mask1=xArray>vStart
                mask2=xArray<vEnd
                temp=self.lastSpectrum[mask1]
                dataToFit=temp[mask2]
            else:
                vStart=self.m_params["vStart"]
                vEnd=self.m_params["vEnd"]
                dataToFit=self.lastSpectrum
            X=np.arange(vStart,vEnd,dV)
            popt, pcov = sp.optimize.curve_fit(fit_func, X, dataToFit)
            slope=popt[0]
            coef=popt[1]
            print("Linear fit gives a slope of "+str(slope)+" and a coef of "+str(coef))
            self.ax_spec.plot(X,slope*X+coef,label="Linear fit of last plot")
            self.ax_spec.legend(loc=0)
            self.m_specWidget.draw()
        
    def getSpectrumColor(self,n):
        i=n%len(self.spectrumColor)
        return self.spectrumColor[i]
        
    def launchAvgSpectrum(self):
        """ Slot called to average the spectra of the whole map """
        if(self.dataLoaded):
            xPx=self.m_params["xPx"]
            yPx=self.m_params["yPx"]
            self.averageSpectrum(0,xPx,0,yPx)
            
    def pickSpectrum(self,event):
        """ Method called when a press-and-release event is done at the same location of the map. If the average box is checked, it will keep the data in storage to average it later. Otherwise it plots the spectrum in the corresponding widget """
        if(event.xdata!=None and event.ydata!=None and self.dataLoaded):
            color=self.getSpectrumColor(self.nSpectraDrawn)
            PixelX=int(event.xdata)
            PixelY=int(event.ydata)
            chan=self.m_channelBox.currentIndex()
            if("Slope" in self.channelList[chan]):
                zPt=self.m_params["zPt"]
                zPts=np.arange(0,zPt)
                dataLogCurrent=np.zeros(shape=(zPt))
                dataLine=self.m_data[chan][PixelY][PixelX][0]*self.m_params["dV"]*zPts+self.m_data[chan+1][PixelY][PixelX][0]                
                cutOff=False
                for z in zPts:
                    i=self.m_data[0][PixelY][PixelX][z]
                    if(i<0.01 or cutOff):
                        dataLogCurrent[z]=dataLine[z]
                    else:
                        dataLogCurrent[z]=np.log(i)
                self.drawSpectrum(dataLogCurrent,"Log of current at ["+str(PixelX)+","+str(PixelY)+"]")
                self.drawSpectrum(dataLine,"Linear fit of the log at "+str(PixelX)+","+str(PixelY)+"]")
            else:
                self.m_mapWidget.draw()
                #Add data to the total data to average if in average mod
                if(self.m_avgBox.isChecked()):
                    self.tot_data+=self.m_data[chan][PixelY][PixelX]
                    self.nAvgSpectra+=1
                    #color='white'
                #Plot the data with the desired scale (Volts or index) if in normal mode
                else:
                    dataToPlot=self.m_data[chan][PixelY][PixelX]
                    self.drawSpectrum(dataToPlot,"["+str(PixelX)+","+str(PixelY)+"]")
            self.addToPtsClicked(PixelX,PixelY,color)
            self.drawPtsClicked()
            
### Methods related to the clicks on the map
            
    def onPressOnMap(self,event):
        """ Slot called when a press event is detected. Keeps the pressed location in integer variables """
        if(event.xdata!=None and event.ydata!=None and self.dataLoaded and self.toolbar_map._active is None):
            self.origin_x=int(event.xdata)
            self.origin_y=int(event.ydata)
            if(event.button==1):
                self.motionConnection = self.m_mapWidget.mpl_connect('motion_notify_event', self.drawLine)
                self.lines=self.ax_map.plot([self.origin_x,event.xdata],[self.origin_y,event.ydata],'k--')
            else:
                self.motionConnection = self.m_mapWidget.mpl_connect('motion_notify_event', self.drawRectangle)
                self.rectangle=self.ax_map.add_patch(patches.Rectangle((self.origin_x, self.origin_y),0,0,color=self.getSpectrumColor(self.nSpectraDrawn)))
                
            
    def drawLine(self,event):
        if(event.xdata!=None and event.ydata!=None):
            self.lines[0].set_xdata([self.origin_x+0.5,int(event.xdata)+0.5])
            self.lines[0].set_ydata([self.origin_y+0.5,int(event.ydata)+0.5])
            self.m_mapWidget.draw()
            
    def drawRectangle(self,event):
        if(event.xdata!=None and event.ydata!=None):
            self.rectangle.set_width(int(event.xdata)-self.origin_x)
            self.rectangle.set_height(int(event.ydata)-self.origin_y)
            self.m_mapWidget.draw()
            
            
    def onReleaseOnMap(self,event):
        """ Slot called when a release event is detected. If it happens at the same location of the press, calls the pickSpectrum method. Otherwise, it means that a line was drawn so it makes a cut of the spectra along the line. """
        if(self.dataLoaded and self.toolbar_map._active is None):
            if(event.xdata!=None and event.ydata!=None):
                xf=int(event.xdata)
                yf=int(event.ydata)
            xi=self.origin_x
            yi=self.origin_y
            #If left-click : either a line was drawn or a spectrum picked
            if(event.button==1):
                self.m_mapWidget.mpl_disconnect(self.motionConnection)
                if(event.xdata==None or event.ydata==None):
                     self.lines.pop(0).remove()
                # Cut along the XY line if a line is traced (X or Y different)
                elif(xf!=xi or yf!=yi):
                    self.cutAlongLine(xi,xf,yi,yf)
                #Pick spectrum otherwise
                else:
                    self.lines.pop(0).remove()
                    self.pickSpectrum(event)
            #If right-click : either a rectangle was drawn or the center of the rectangle to average was picked
            else:
                self.m_mapWidget.mpl_disconnect(self.motionConnection)
                if(event.xdata!=None and event.ydata!=None):
                    if(xf!=xi or yf!=yi):
                        self.averageSpectrum(xi,xf,yi,yf)
                    else:
                        n=self.m_rcAvgBox.value()
                        self.rectangle.set_x(max(0,xi-n))
                        self.rectangle.set_y(max(0,yi-n))
                        self.rectangle.set_width(min(self.m_params["xPx"],xf+n)-max(0,xi-n))
                        self.rectangle.set_height(min(self.m_params["yPx"],yf+n)-max(0,yi-n))
                        self.averageSpectrum(max(0,xi-n),min(self.m_params["xPx"],xf+n),max(0,yi-n),min(self.m_params["yPx"],yf+n))
            self.m_mapWidget.draw()


    def addToPtsClicked(self,x,y,color):
        """ Method called when a click was detected on the map. The coordinates are saved in the pts_clicked list with a corresponding color """
        # Shifting the coordinates by (0.5,0.5) so that the dot is plotted at the center of the square
        self.pts_clicked[0].append(x+0.5)
        self.pts_clicked[1].append(y+0.5)
        #color_cycle = self.ax_spec._get_lines.color_cycle
        self.pts_clicked[2].append(color)
        
    def drawPtsClicked(self):
        """ Method that draws all the points saved in pts_clicked. Usually called when the map is drawn or redrawn """
        pts_x=np.array(self.pts_clicked[0])
        pts_y=np.array(self.pts_clicked[1])
        pts_c=np.array(self.pts_clicked[2])
        self.fig_topo.axes[0].scatter(pts_x*self.m_params["xL"]/self.m_params["xPx"],pts_y*self.m_params["yL"]/self.m_params["yPx"],color=pts_c)
        self.fig_topo.canvas.draw()
        self.ax_map.scatter(pts_x,pts_y,color=pts_c)
        self.m_mapWidget.draw()
        
    def clearPtsClicked(self):
        """ Method that clears all saved points and then redraws the map to reflect the change """
        self.pts_clicked=[[],[],[]]
        self.drawXYMap(self.m_voltageBox.value())
        if(self.dataLoaded):
            pass
            #pylab.close("all")
            #self.drawTopo()
        
    def cutAlongLine(self,xi,xf,yi,yf):
        #If the line is vertical, the equation is x=xi with y varying from yi to yf 
        if(xf==xi):
            if(yi>yf): y_plot=np.arange(yf,yi+1)[::-1]
            else: y_plot=np.arange(yi,yf+1)
            x_plot=np.full(shape=y_plot.size,fill_value=xi)
        else:
            #Simple algorithm for cuts
            simple=True
            if(simple):
                #If the line is not vertical, determine its equation y=k*x+c
                k=float(yf-yi)/(xf-xi)
                c=yi-k*xi                
                #Check if there is more y or more x to have to most precise arrangment
                if(abs(xf-xi)>abs(yf-yi)):
                    if(xi>xf): x_plot=np.arange(xf,xi+1)[::-1]
                    else: x_plot=np.arange(xi,xf+1)
                    y_plot=k*x_plot+c
                else:
                    if(yi>yf): y_plot=np.arange(yf,yi+1)[::-1]
                    else: y_plot=np.arange(yi,yf+1)
                    x_plot=(y_plot-c)/k
            # Bresenham algorithm
            else:
                x_plot_p=[]
                y_plot_p=[]
                dx=abs(xi-xf)
                dy=abs(yi-yf)
                x0=xi
                y0=yi
                if(xi<xf): sx=1
                else: sx=-1
                if(yi<yf): sy=1
                else: sy=-1
                err = dx-dy
                while(True):
                    x_plot_p.append(x0)
                    y_plot_p.append(y0)
                    if(x0>=xf and y0>=yf): break
                    e2=err*2
                    if(e2>-dy):
                        err -= dy
                        x0 += sx
                    if(e2<dx):
                        err += dx
                        y0 +=sy
                x_plot=np.array(x_plot_p)
                y_plot=np.array(y_plot_p)
        self.cutPlot(x_plot,y_plot)
            
            
    def cutPlot(self,x_plot,y_plot):
        #Build the data to plot with v as Y and z (number of pixels gone through) as X
        zPt=self.m_params['zPt']
        voltages=np.arange(0,zPt)
        chan=self.m_channelBox.currentIndex()
        z_plot=np.arange(0,x_plot.size)
        dataToPlot=np.ndarray(shape=(zPt,z_plot.size))
        xi=x_plot[0]
        yi=y_plot[0]
        zf=z_plot[-1]
        # Variables needed to compute the metric distances
        dx=self.m_params["xL"]/self.m_params["xPx"]
        dy=self.m_params["yL"]/self.m_params["yPx"]
        metricDistances=np.sqrt((dx*(x_plot-xi))**2+(dy*(y_plot-yi))**2)
        # Bool to keep trace of the selected spectra
        viewSelected=self.m_viewSelectedBox.isChecked()
        # Matlab convention : Y (v) first then X (z)
        # Plot the built map in a new figure
        fig=pylab.figure()
        ax=fig.add_subplot(1,1,1)
        ax.set_title(self.mapName.split(".")[0]+" - Cut "+str(fig.number))
        self.ax_map.text(xi+0.5,yi+0.5,str(fig.number))
        if(self.m_waterfallButton.isChecked()):
            if(self.m_scaleVoltage.isChecked()): 
                voltages=self.m_params["vStart"]+voltages*self.m_params["dV"]
            for z in z_plot:
                xc=int(x_plot[zf-z-1])
                yc=int(y_plot[zf-z-1])
                spectrum=self.m_data[chan][yc][xc]
                offset = (zf-z)*self.m_shiftBox.value()
                ax.plot(voltages,spectrum+offset, 'k', lw=2, zorder=(z+1)*2)
                ax.fill_between(voltages, spectrum+offset, offset, facecolor='w', lw=0, zorder=(z+1)*2-1)
            ax.set_xlim([voltages[0],voltages[-1]])
        else:
            for v in voltages:
                for z in z_plot:
                    xc=int(x_plot[z])
                    yc=int(y_plot[z])
                    dataToPlot[v][z]=self.m_data[chan][yc][xc][v]
                    if(viewSelected): self.addToPtsClicked(xc,yc,color='yellow')
            ax.set_ylabel("Voltage index")
            fig.subplots_adjust(left=0.125,right=0.95,bottom=0.15,top=0.92)
            # Change the scales if needed
            if(self.m_scaleVoltage.isChecked()): 
                voltages=self.m_params["vStart"]+voltages*self.m_params["dV"]
                ax.set_ylabel("Bias (V)")
            if(self.m_scaleMetric.isChecked()):
                mapData=pylab.pcolormesh(metricDistances,voltages,dataToPlot,cmap=self.m_colorBarBox.currentText())
                ax.axis([metricDistances[0],metricDistances[-1],voltages[0],voltages[-1]])
                ax.set_xlabel("Distance (nm)")
            else:
                mapData=pylab.pcolormesh(z_plot,voltages,dataToPlot,cmap=self.m_colorBarBox.currentText())
                ax.axis([z_plot[0],z_plot[-1],voltages[0],voltages[-1]])
                ax.set_xlabel("Pixels")
            #Colorbar set up
            cbar = fig.colorbar(mapData, shrink=.9, pad=0.05, aspect=15)
            cbar.ax.yaxis.set_ticks_position('both')
            cbar.ax.tick_params(axis='y', direction='in')
            if(self.m_cbarCustomCheckbox.isChecked()): mapData.set_clim(float(self.m_cbarLowerBox.text()),float(self.m_cbarUpperBox.text()))
        self.drawPtsClicked()
        
    def launchBigCut(self):
        if(self.dataLoaded):
            self.ax_map.plot([0.5,self.m_params["xPx"]-0.5],[0.5,self.m_params["yPx"]-0.5],'k--')
            self.cutAlongLine(0,self.m_params["xPx"]-1,0,self.m_params["yPx"]-1)
                
            
### Methods related to the map
        
    def clearMap(self):
        """ Unloads the map and clears the map window """
        self.m_mapWidget.figure.clear()
        self.clearPtsClicked()
        self.m_data=np.ndarray([])
        self.dataLoaded=False
        
    def drawXYMap(self,voltage):
        """ Calls the getMapData function and draws the result in the map window with the approriate formatting. Called at each change of voltage """
        #Start everything anew
        fig_map=self.m_mapWidget.figure
        fig_map.clear()
        if(self.dataLoaded):
            xPx=self.m_params["xPx"]
            yPx=self.m_params["yPx"]
            #zPt=self.m_params["zPt"]
            if(xPx==yPx):
                self.ax_map = fig_map.add_subplot(1,1,1,aspect='equal')
            else:
                self.ax_map = fig_map.add_subplot(1,1,1)
            self.ax_map.hold(True)
            fig_map.subplots_adjust(left=0.125,right=0.95,bottom=0.15,top=0.92)
            #Get the data of the map and draw it
            mapData,self.mapMin,self.mapMax=self.getMapData(voltage)
            #Use metric dimensions if the corresponding box is checked (work in progress)
            if(self.m_scaleMetric.isChecked() and False):
                xL=self.m_params["xL"]
                yL=self.m_params["yL"]
                dx=xL/xPx
                dy=yL/yPx
                x_m=np.arange(0,xL,dx)
                y_m=np.arange(0,yL,dy)
                XYmap=self.ax_map.pcolormesh(x_m,y_m,mapData,cmap=self.m_colorBarBox.currentText())
                self.ax_map.axis([0,xL,0,yL])
            #Else, use pixels
            else:
                XYmap=self.ax_map.pcolormesh(mapData,cmap=self.m_colorBarBox.currentText())
                self.ax_map.axis([0,xPx,0,yPx])
            #Set title
            chan=self.m_channelBox.currentText()
            #print(chan)
            self.ax_map.set_title(self.mapName+" - "+chan+"\n"
            +"V="+str(self.m_params["vStart"]+voltage*self.m_params["dV"]))
            #Colorbar stuff
            cbar = fig_map.colorbar(XYmap, shrink=.9, pad=0.05, aspect=15)
            cbar.ax.yaxis.set_ticks_position('both')
            cbar.ax.tick_params(axis='y', direction='in')
            # Image color scale is adjusted to the data:
            if(self.m_cbarCustomCheckbox.isChecked()): XYmap.set_clim(float(self.m_cbarLowerBox.text()),float(self.m_cbarUpperBox.text()))
            else: XYmap.set_clim(self.mapMin,self.mapMax)
            #Plot a dashed line at X=voltage if asked
            if(self.m_vLineBox.isChecked()):
                self.drawVoltageLine(voltage)
            #Draws the points saved in pts_clicked
            self.drawPtsClicked()
            #Not really useful as the draw is called in drawPtsClicked but it is more logical to keep it.
            self.m_mapWidget.draw()
        
        
    def getMapData(self,v):
        """ Returns an array constructed from the data loaded that can be used to display a map at fixed voltage """
        xPx=self.m_params["xPx"]
        yPx=self.m_params["yPx"]
        mapData=np.ndarray(shape=(yPx,xPx))
        chan=self.m_channelBox.currentIndex()
        for y in range(0,yPx):
            for x in range(0,xPx):
                val=self.m_data[chan][y][x][v]
                mapData[y][x]=val
                if(x==0 and y==0):
                    valMin=val
                    valMax=val
                else:
                    if(val<valMin): valMin=val
                    elif(val>valMax): valMax=val
        return mapData,valMin,valMax
        
### Methods related to the voltage guide line in the spectra window
        
    def drawVoltageLine(self,voltage):
        self.clearVoltageLine()
        if(self.dataLoaded):
            #Get the current voltage : real voltage if the scale box is checked, voltage index otherwise
            if(self.m_scaleVoltage.isChecked): currentV=self.m_params["vStart"]+voltage*self.m_params["dV"]
            else: currentV=voltage
            #Plot the dashed line
            self.voltageLine=self.ax_spec.axvline(currentV,color='k',linestyle='--')
            self.m_specWidget.draw()
        
    def clearVoltageLine(self):
        if(self.voltageLine):
            self.ax_spec.lines.pop(self.ax_spec.lines.index(self.voltageLine))
            self.m_specWidget.draw()
        self.voltageLine=0
        
### Post-processing methods
    def addChannel(self,newChannelData,channelName):
        #Gets parameters for reshaping
        nChannels=len(self.channelList)
        xPx=self.m_params["xPx"]
        yPx=self.m_params["yPx"]
        zPt=self.m_params["zPt"]
        try:
            newData=np.zeros(shape=(nChannels+1,yPx,xPx,zPt))
        except MemoryError:
            print("The data is too big ! Or the memory too small...")
            return
        newData[0:nChannels]=self.m_data
        newData[nChannels]=newChannelData
        self.m_data=newData
        del newData
        #Updates the channel list
        self.channelList.append(channelName)
        self.m_channelBox.clear()
        self.m_channelBox.addItems(self.channelList)
        
    def extractOutOfPhase(self,numChanR,numChanPhi):
        #Phase : 9V = pi
        outOfPhase=-self.m_data[numChanR]*np.cos(self.m_data[numChanPhi]*np.pi/9)
        #Add the channel to the data
        self.addChannel(outOfPhase,self.channelList[numChanR]+"cos("+self.channelList[numChanPhi]+")")
        
    def extractDerivative(self,numChanToDeriv):
        dV=self.m_params["dV"]
        yPx=self.m_params["yPx"]
        xPx=self.m_params["xPx"]
        zPt=self.m_params["zPt"]
        derivData=np.zeros(shape=(yPx,xPx,zPt))
        for y in range(0,yPx):
            for x in range(0,xPx):
                derivData[y][x]=sp.signal.savgol_filter(self.m_data[numChanToDeriv][y][x], 9, 2, deriv=1, delta=dV)
        #Add the channel to the data
        self.addChannel(derivData,"Derivative of "+self.channelList[numChanToDeriv])
        
    def extractFFT(self,numChanToFFT,axisOfFFT):
        yPx=self.m_params["yPx"]
        xPx=self.m_params["xPx"]
        zPt=self.m_params["zPt"]
        FFTData=np.zeros(shape=(yPx,xPx,zPt))
        FFTData=np.fft(self.m_data[numChanToFFT],axis=axisOfFFT)
        #Add the channel to the data
        self.addChannel(FFTData,"FFT of "+self.channelList[numChanToFFT])
        
    def computeAngle(Dmoire,k=True):
        return 2*np.arcsin(0.246/(2*Dmoire))*180/np.pi
        
    def extractSlope(self,cutOffValue,numChanToFit):
        """ Do a linear fit of the data in the asked channel and add the slope and the coef found as channels (usually called for zSpectros) """
        yPx=self.m_params["yPx"]
        xPx=self.m_params["xPx"]
        zPt=self.m_params["zPt"]
        dZ=self.m_params["dV"]
        zg=np.zeros(shape=(yPx,xPx,zPt))
        slopeData=np.zeros(shape=(yPx,xPx,zPt))
        coefData=np.zeros(shape=(yPx,xPx,zPt))
        fit_func = lambda v,a,b: a*v+b
        xArray=np.arange(0,zPt)*dZ
        for y in range(0,yPx):
            for x in range(0,xPx):
                    rawData=self.m_data[numChanToFit][y][x]
                    mask=rawData>cutOffValue
                    xArrayFiltered=xArray[mask]
                    dataFiltered=np.log(rawData[mask])
                    popt, pcov = sp.optimize.curve_fit(fit_func, xArrayFiltered, dataFiltered)
                    zg[y][x]=1/20.5*np.log(rawData)+np.arange(0,zPt*dZ,dZ)
                    for z in range(0,zPt):
                        slopeData[y][x][z]=popt[0]
                        coefData[y][x][z]=popt[1]
        #Add the created channel to the data
        self.addChannel(slopeData,"Slope by linear fit of "+self.channelList[numChanToFit])
        #self.addChannel(coefData,"Coef by linear fit of "+self.channelList[numChanToFit])
        self.addChannel(zg,"Zg")
        
    def magicLinearCut(self,i):
        #Launch a linear cut along x with y=0
        xPx=self.m_params["xPx"]
        w=self.m_params["xL"]
        x_plot=np.arange(0,xPx)
        #Build the data to plot with v as Y and z (number of pixels gone through) as X
        zPt=self.m_params['zPt']
        voltages=np.arange(0,zPt)
        chan=self.m_channelBox.currentIndex()
        dataToPlot=np.ndarray(shape=(zPt,x_plot.size))
        # Variables needed to compute the metric distances
        dx=w/xPx
        metricDistances=dx*(x_plot)
        # Matlab convention : Y (v) first then X (z)
        valMax=0
        for v in voltages:
            for x in x_plot:
                val=self.m_data[1][0][x][v]
                dataToPlot[v][x]=val
                if(valMax==0):
                    valMax=val
                if(val>valMax):
                    valMax=val
        # Plot the built map in a new figure
        fig=pylab.figure(figsize=(8,10))
        ax=fig.add_subplot(2,1,1,anchor="W")
        ax.set_title(self.mapName.split(".")[0])
        ax.set_ylabel("Voltage index")
        #fig.subplots_adjust(left=0.125,right=0.95,bottom=0.15,top=0.92)
        # Change the scales if needed
        if(self.m_scaleVoltage.isChecked()): 
            voltages=self.m_params["vStart"]+voltages*self.m_params["dV"]
            ax.set_ylabel("Bias (V)")
            if(self.m_scaleMetric.isChecked()):
                mapData=pylab.pcolormesh(metricDistances,voltages,dataToPlot,cmap=self.m_colorBarBox.currentText())
                ax.axis([metricDistances[0],metricDistances[-1],voltages[0],voltages[-1]])
                ax.set_xlabel("Distance (nm)")
            else:
                mapData=pylab.pcolormesh(x_plot,voltages,dataToPlot,cmap=self.m_colorBarBox.currentText())
                ax.axis([z_plot[0],z_plot[-1],voltages[0],voltages[-1]])
                ax.set_xlabel("Pixels")
            #Colorbar set up
            cbar = fig.colorbar(mapData, shrink=.9, pad=0.05, aspect=15)
            cbar.ax.yaxis.set_ticks_position('both')
            cbar.ax.tick_params(axis='y', direction='in')
            if(self.m_cbarCustomCheckbox.isChecked()): mapData.set_clim(float(self.m_cbarLowerBox.text()),float(self.m_cbarUpperBox.text()))
            else: mapData.set_clim(0,0.8*valMax)
                #cur=self.mapName.replace(",",".")
                #if("nA" in cur):                
                #    mapData.set_clim(0,float(cur[0:-2]))
                #else:
                #    mapData.set_clim(0,float(cur[0:-2])/1000)
            
        #Draw the topo in a subplot
        fitX=np.arange(0,xPx)
        fitF = lambda z,a,b: a*z+b  
        y=0
        fitY=self.topo[y]
        f = sp.interpolate.InterpolatedUnivariateSpline(fitX,fitY, k=1)
        popt, pcov = sp.optimize.curve_fit(fitF, fitX, f(fitX))
        test=fitY-(popt[0]*fitX+popt[1])
        print(y,popt[0],popt[1],pcov[0],pcov[1])
        ax_topo=fig.add_subplot(2,1,2,anchor="W")
        cbar = fig.colorbar(mapData, ax=ax_topo, shrink=.9, pad=0.05, aspect=15)               
               
        ax_topo.plot(np.linspace(0,w,xPx),test,label="With line leveling")
        ax_topo.plot(np.linspace(0,w,xPx),self.topo[y]-popt[1],label="Without line leveling")
        ax_topo.set_xlim(0,w-dx)
        ax_topo.set_xlabel("Distance (nm)")
        ax_topo.set_ylabel("Z (nm)")
        ax_topo.set_ylim(-0.6,0.8)
        pylab.legend(loc=0)        
        
        fig.savefig("E:\\PhD\\Experiments\\STM\\Lc12\\Depouillements\\Depouillement_LS_20B\\Vertical\\"+self.mapName+".png")
                
    def magicFunction(self):
        """ Function dedicated to the execution of bits of code on the fly by pression the 'magic button' """
        i=0
        #pylab.figure()
        path="E:\\PhD\\Experiments\\STM\\Lc12\\LS\\Line_Spectro_20\\"
        for fich in listdir(path):
            print fich
            self.dataLoaded=self.readCitsBin(path+fich)
            pylab.close()
            cur=fich.split('_')[-1].split('.')[0]
            self.mapName=cur
            
            #self.fig_topo.suptitle(cur)
            #self.ax_topo.set_title(self.mapName.split(".")[0].split("-")[-1])
            #self.ax_topo.set_ylim(-.1,0.4)
            #self.fig_topo.savefig("E:\\PhD\\Experiments\\STM\\Lc12\\Depouillements\\Depouillement_LS_4B\\"+self.mapName+"_Z"+".png")
            #pylab.close()
            
            # Coolwarm cmap
            #self.m_colorBarBox.setCurrentIndex(self.m_colorBarBox.findText("coolwarm"))
            
            self.magicLinearCut(i)
            pylab.close()
            i=i+1
            
            #xPx=self.m_params["xPx"]
            #yPx=self.m_params["yPx"]
            #zPt=self.m_params["zPt"]
            #vS=self.m_params["vStart"]
            #vE=self.m_params["vEnd"]
            #dV=self.m_params["dV"]
            #avg_data=np.zeros(shape=zPt)
            #N=xPx*yPx
            #for y in range(0,yPx):
            #    for x in range(0,xPx):
            #            avg_data+=(self.m_data[1][y][x]/N)
            #pylab.plot(np.arange(vS,vE,dV),(avg_data+i*0.3),label=cur)
            #i=i+1
        #pylab.legend()
     
### Complentary functions
     
def plotMapDistrib(data1d,valMax):
    dx=valMax/100
    distrib=np.zeros(shape=101)
    for val in data1d:
        index=int(val//dx)
        distrib[index]=distrib[index]+1
    pylab.figure()
    pylab.plot(distrib)
        
def getDataFromTxtFile(fileName,header=False):
    data=np.genfromtxt(fileName,delimiter="\t",unpack=True)
    return data
    
def diffSingularite(voltage):
    dataS1=getDataFromTxtFile("E:\\PhD\\Experiments\\STM\\Lc12\\Depouillements\\Depouillement_LS_1B\\"+voltage+"\\Singularite1.csv")
    f=sp.interpolate.InterpolatedUnivariateSpline(dataS1[0],dataS1[1], k=1)
    dataS2=getDataFromTxtFile("E:\\PhD\\Experiments\\STM\\Lc12\\Depouillements\\Depouillement_LS_1B\\"+voltage+"\\Singularite2.csv")
    g=sp.interpolate.InterpolatedUnivariateSpline(dataS2[0],dataS2[1], k=1)
    x=np.linspace(0,40,500)
    print("s")
    plot(x,g(x)-f(x),label=voltage)
    return
    
def angleFromSingulariteDiff(dE,adv=False):
    if(adv):
        sol1=2*np.arcsin(8*10**(-3)*(dE/0.108+2+np.sqrt((2+dE/0.108)**2+36)))
        sol2=2*np.arcsin(8*10**(-3)*(dE/0.108+2-np.sqrt((2+dE/0.108)**2+36)))
        return 180*sol1/np.pi,180*sol2/np.pi
    else:
        sol=2*arcsin(16*10**(-3)*(dE/0.108+2))
        return 180*sol/np.pi