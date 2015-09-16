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
import scipy.signal as sp
import matplotlib.animation as animation

class CitsWidget(QDockWidget, Ui_CitsWidget):
### Building methods
    def __init__(self,parent):
        """ Builds the widget with parent widget in arg """
        QDockWidget.__init__(self,parent)
        #Set up the user interface from Designer.
        self.setupUi(self)
        self.setAutoFillBackground(True)
        self.m_data=np.ndarray([])
        self.tot_data=np.ndarray([])
        self.pts_clicked=[[],[],[]]
        self.nAvgSpectra=0
        self.m_params={}
        self.mapName=""
        #Set up figures
        self.m_fwdButton.setChecked(True)
        self.toolbar_map = NavigationToolbar(self.m_mapWidget,self)
        self.toolbar_spec = NavigationToolbar(self.m_specWidget,self)
        self.map_layout.insertWidget(0,self.toolbar_map)
        self.spec_layout.insertWidget(0,self.toolbar_spec)
        self.fig_spec=self.m_specWidget.figure
        self.axes_spec=self.fig_spec.add_subplot(1,1,1)
        self.axes_spec.hold(True)
        self.fig_spec.subplots_adjust(left=0.125,right=0.95,bottom=0.15,top=0.92)
        #Boolean that is True if a map is laoded
        self.dataLoaded=False
        self.connect()
        #self.readTopo("H:\\Experiments\\STM\\Lc0 (Si-260)\\4K\\Spectro ASCII\\TestZ.txt")
    
    def connect(self):
        """ Connects all the signals. Only called in the constructor """
        self.m_openButton.clicked.connect(self.loadCITS)
        self.m_bwdButton.toggled.connect(self.updateWidgets)
        self.m_voltageBox.valueChanged.connect(self.drawXYMap)
        self.m_mapWidget.mpl_connect('button_press_event',self.pickSpectrum)
        self.m_specWidget.mpl_connect('button_press_event',self.updateToPointedX)
        self.m_clearSpec.clicked.connect(self.clearSpectrum)
        self.m_scaleVoltage.toggled.connect(self.clearSpectrum)
        self.m_avgSpec.clicked.connect(self.launchAvgSpectrum)
        self.m_avgBox.toggled.connect(self.updateAvgVariables)
        
### Reading and loading CITS methods
        
    def loadCITS(self):
        """ Slot that launches the reading of the CITS by asking the path of the file """
        filename=QFileDialog.getOpenFileName(self,"Choose a CITS Ascii file","H:\\Experiments\\STM\\Lc0 (Si-260)\\4K\\Spectro ASCII","Ascii file (*.asc);;Text file (*.txt)")
        self.mapName=osp.basename(filename).split()[0]
        print(self.mapName)
        self.dataLoaded=self.readCITS(filename)
        if(self.dataLoaded): self.updateWidgets()
        
    def readCITS(self,filepath):
        """ Reads a CITS file and stores all the parameters"""
        f=open(filepath)
        divider=1
        while(True):
            #Read the parameters of the map until "Start of Data"
            line=f.readline()
            if("x-pixels" in line):
                xPx=int(line.split()[-1])
            if("y-pixels" in line):
                yPx=int(line.split()[-1])
            if("x-length" in line):
                xL=float(line.split()[-1])
            if("y-length" in line):
                yL=float(line.split()[-1])
            if("z-points" in line):
                zPt=int(line.split()[-1])
            if("Device_1_Start" in line):
                vStart=float(line.split()[-2])
            if("Device_1_End" in line):
                vEnd=float(line.split()[-2])
            if("divider" in line):
                divider=int(line.split()[-1])
            if("Start of Data" in line):
                break
        # Matlab convention : columns first then rows hence [y][x]
        self.m_data=np.zeros(shape=(2,yPx,xPx,zPt/2))
        for y in range(0,yPx):
            for x in range(0,xPx):
                #The line read is an array containing the dI/dV (or I(V)) values indexed by the voltage index
                #String to remove the newline at the end and split to transform the string in a list
                #0 means Backwards=False so forward data, 1 means Backwards=True so backward data
                data_list=f.readline().strip().split()
                self.m_data[0][y][x]=data_list[0:zPt/2]
                #No need to reverse the backward data as it was from Vmin to Vmax in the file as the fwd data
                #self.m_data[1][y][x]=(data_list[zPt/2:zPt])[::-1]
                self.m_data[1][y][x]=(data_list[zPt/2:zPt])
        f.close()
        #Store the parameters in a dictonnary to use them later
        self.m_params={"xPx":xPx,"yPx":yPx,"xL":xL,"yL":yL,"zPt":zPt,"vStart":vStart/divider,"vEnd":vEnd/divider,"dV":2*abs(vEnd-vStart)/(divider*zPt),"Div":divider}
        if(divider!=1):
            print "A divider of "+str(divider)+" was found and applied"
        return True
        
### Reading and loading topo images methods (not finished yet)
        
    def readTopo(self,filepath):
        """ Reads a topography file """
        f=open(filepath)
        topo_data=[[]]
        for line in f:
            #Treat the headers differently"
            if(line[0]=='#'):
                if("Width" in line):
                    w=float(line.split()[-2])
                if("Height" in line):
                    h=float(line.split()[-2])
            else:
                topo_data.append(line.strip().split())
        f.close()
        fig_map=self.m_mapWidget.figure
        fig_map.clear()
        ax_map = fig_map.add_subplot(1,1,1,aspect='equal')
        fig_map.subplots_adjust(left=0.125,right=0.95,bottom=0.15,top=0.92)
        XYmap=ax_map.pcolormesh(np.ndarray(topo_data))
        self.m_mapWidget.draw()
        return True
        
### Updating methods. Usually called by signals
   
    def updateAvgVariables(self):
        """ Slot called by the toggling of the average box. Toggling on the box clears everything to start a new averaging. Toggling off averages and plots the data stored by picking spectra """
        if(self.dataLoaded):
            #If toggled on, put everything to zero to be ready to store data
            if(self.m_avgBox.isChecked()):
                self.tot_data=np.zeros(shape=(self.m_params["zPt"]/2))
                self.nAvgSpectra=0
            #If toggled off, plot the data stored
            else:
                if(self.nAvgSpectra==0): return
                dataToPlot=self.tot_data/(self.nAvgSpectra)
                self.drawSpectrum(dataToPlot,str(self.nAvgSpectra)+" spectra averaged")
                
    def updateToPointedX(self,event):
        """ Slot called when clicking on the spectrum window. Updates the map according to the position of the cursor when clicked """
        if(event.xdata!=None and event.ydata!=None and self.dataLoaded):
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
        self.m_voltageBox.setMaximum(zPt/2-1)
        self.m_voltageBox.setSingleStep(1)
        
    def updateWidgets(self):
        """ Slot called after the reading of the CITS. Sets the values of the voltage box and draws the map """
        self.updateVoltageBox(self.m_params["vStart"],self.m_params["vEnd"],self.m_params["zPt"])
        self.drawXYMap(self.m_voltageBox.value())
        
### Methods related to spectra    

    def averageSpectrum(self,Xs,Ys,Xe,Ye):
        """ Averages the spectra contained in the rectangle drawn by the points (Xs,Ys);(Xe,Ys);(Xe;Ye) and (Xs;Ye) """
        if(self.dataLoaded):
            zPt=self.m_params["zPt"]
            bwd=self.m_bwdButton.isChecked()
            avg_data=np.zeros(shape=(zPt/2))
            N=(Ye-Ys)*(Xe-Xs) #Number of spectra averaged
            for y in range(Ys,Ye):
                for x in range(Xs,Xe):
                    for v in range(0,zPt/2):
                        avg_data[v]+=(self.m_data[bwd][y][x][v]/N)
            self.drawSpectrum(avg_data,"Avg")
            
    def clearSpectrum(self):
        """ Clears the spectrum window """
        self.axes_spec.clear()
        self.clearPtsClicked()
        self.m_specWidget.draw()
            
    def drawSpectrum(self,dataToPlot,label):
        """ Method called each time a spectrum needs to be plotted. Takes care of the derivative and different scales stuff and updates the window """
        bwd=self.m_bwdButton.isChecked()
        fb=[" (Fwd)"," (Bwd)"]
        dV=self.m_params["dV"]
        
        if(self.dataLoaded and dataToPlot.size!=0):
            deriv=sp.savgol_filter(dataToPlot, 5, 2, deriv=1, delta=dV)*(-10**9)
            if(self.m_scaleVoltage.isChecked()):
                vStart=self.m_params["vStart"]
                vEnd=self.m_params["vEnd"]
                self.axes_spec.plot(np.arange(vStart,vEnd,dV),dataToPlot,label=label+fb[bwd])
                if(self.m_derivBox.isChecked()):
                    self.axes_spec.plot(np.arange(vStart,vEnd,dV),deriv,label="Derivative of "+label+fb[bwd])
            else:
                self.axes_spec.plot(dataToPlot,label=label+fb[bwd])
                if(self.m_derivBox.isChecked()):
                    self.axes_spec.plot(deriv,label="Derivative of "+label+fb[bwd])
        
        self.axes_spec.legend(loc=0)
        self.m_specWidget.draw()
        
    def launchAvgSpectrum(self):
        """ Slot called to average the spectra of the whole map """
        if(self.dataLoaded):
            xPx=self.m_params["xPx"]
            yPx=self.m_params["yPx"]
            self.averageSpectrum(0,0,xPx,yPx)
            
    def pickSpectrum(self,event):
        """ Slot called when clicking on the map. If the average box is checked, it will keep the data in storage to average it later. Otherwise it plots the spectrum in the corresponding widget """
        if(event.xdata!=None and event.ydata!=None and self.dataLoaded):
            color='black'
            PixelX=int(event.xdata)
            PixelY=int(event.ydata)
            bwd=self.m_bwdButton.isChecked()
            self.m_mapWidget.draw()
            #Add data to the total data to average if in average mod
            if(self.m_avgBox.isChecked()):
                self.tot_data+=self.m_data[bwd][PixelY][PixelX]
                self.nAvgSpectra+=1
                color='white'
            #Plot the data with the desired scale (Volts or index) if in normal mode
            else:
                dataToPlot=self.m_data[bwd][PixelY][PixelX]
                self.drawSpectrum(dataToPlot,"["+str(PixelX)+","+str(PixelY)+"]")
            self.addToPtsClicked(PixelX+0.5,PixelY+0.5,color)
            self.drawPtsClicked()
            
### Methods related to the indicators (points clicked) on the map
    
    def addToPtsClicked(self,x,y,color):
        self.pts_clicked[0].append(x)
        self.pts_clicked[1].append(y)
        #color_cycle = self.axes_spec._get_lines.color_cycle
        self.pts_clicked[2].append(color)
        
    def drawPtsClicked(self):
        pts_x=np.array(self.pts_clicked[0])
        pts_y=np.array(self.pts_clicked[1])
        pts_c=np.array(self.pts_clicked[2])
        self.ax_map.scatter(pts_x,pts_y,color=pts_c)
        self.m_mapWidget.draw()
        
    def clearPtsClicked(self):
        self.pts_clicked=[[],[],[]]
        self.drawXYMap(self.m_voltageBox.value())
            
### Methods related to the map
        
    def clearMap(self):
        """ Unloads the map and clears the map window """
        self.m_mapWidget.figure.clear()
        self.dataLoaded=False
        
    def drawXYMap(self,voltage):
        """ Calls the getMapData function and draws the result in the map window with the approriate formatting. Called at each change of voltage """
        #Start everything anew
        fig_map=self.m_mapWidget.figure
        fig_map.clear()
        self.ax_map = fig_map.add_subplot(1,1,1,aspect='equal')
        self.ax_map.hold(True)
        fig_map.subplots_adjust(left=0.125,right=0.95,bottom=0.15,top=0.92)
        xPx=self.m_params["xPx"]
        yPx=self.m_params["yPx"]
        zPt=self.m_params["zPt"]
        #Get the data of the map and draw it
        mapData,mapMin,mapMax=self.getMapData(voltage)
        XYmap=self.ax_map.pcolormesh(mapData)
        #XYmap=ax_map.imshow(mapData,origin="lower")
        #Set axis limits and title
        self.ax_map.axis([0,xPx,yPx,0])
        self.ax_map.set_title("CITS Map\n"
        +str(xPx)+"x"+str(yPx)+" pixels - "+str(zPt/2)+" points\n"
        +"V="+str(self.m_params["vStart"]+voltage*self.m_params["dV"]))
        #Colorbar stuff
        cbar = fig_map.colorbar(XYmap, shrink=.9, pad=0.05, aspect=15)
        cbar.ax.yaxis.set_ticks_position('both')
        cbar.ax.tick_params(axis='y', direction='in')
        # Image color scale is adjusted to the data:
        XYmap.set_clim(mapMin,mapMax)
        #Plot a dashed line at X=voltage
        #self.axes_spec.plot([voltage,voltage],[mapMin,mapMax],'k--')
        self.drawPtsClicked()
        self.m_mapWidget.draw()
        
    def getMapData(self,v,bwd=False):
        """ Returns an array constructed from the data loaded that can be used to display a map at fixed voltage """
        xPx=self.m_params["xPx"]
        yPx=self.m_params["yPx"]
        mapData=np.ndarray(shape=(yPx,xPx))
        bwd=self.m_bwdButton.isChecked()
        for y in range(0,yPx):
            for x in range(0,xPx):
                val=self.m_data[bwd][y][x][v]
                mapData[y][x]=val
                if(x==0 and y==0):
                    valMin=val
                    valMax=val
                else:
                    if(val<valMin): valMin=val
                    elif(val>valMax): valMax=val
        return mapData,valMin,valMax