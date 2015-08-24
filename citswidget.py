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
import matplotlib.animation as animation

class CitsWidget(QDockWidget, Ui_CitsWidget):
    def __init__(self,parent):
        QDockWidget.__init__(self,parent)
        # Set up the user interface from Designer.
        self.setupUi(self)
        self.setAutoFillBackground(True)
        self.m_data=np.ndarray([])
        self.tot_data=np.ndarray([])
        self.nAvgSpectra=0
        self.m_params={}
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
        
    def launchAni(self):
        line_ani = animation.FuncAnimation(self.m_mapWidget.figure, self.drawXYMap, 25, fargs=(range(0,50)),interval=500, blit=True)
        
    def loadCITS(self):
        filename=QFileDialog.getOpenFileName(self,"Choose a CITS Ascii file","H:\\Experiments\\STM","Ascii file (*.asc);;Text file (*.txt)")
        #filename="H:\\testCITS.asc"
        self.dataLoaded=self.readCITS(filename)
        if(self.dataLoaded): self.updateWidgets()
        
    def updateVoltageBox(self,Vmin,Vmax,zPt):
        #self.m_voltageBox.setMinimum(Vmin)
        #self.m_voltageBox.setMaximum(Vmax)
        #self.m_voltageBox.setSingleStep(dV)
        self.m_voltageBox.setMinimum(0)
        self.m_voltageBox.setMaximum(zPt/2-1)
        self.m_voltageBox.setSingleStep(1)
        
    def updateAvgVariables(self):
        if(self.dataLoaded):
            #If toggled on, put everything to zero to be ready to store data
            if(self.m_avgBox.isChecked()):
                self.tot_data=np.zeros(shape=(self.m_params["zPt"]/2))
                self.nAvgSpectra=0
            #If toggled off, plot the data stored
            else:
                if(self.nAvgSpectra==0): return
                dataToPlot=self.tot_data/(self.nAvgSpectra)
                if(self.m_scaleVoltage.isChecked()):
                    vStart=self.m_params["vStart"]
                    vEnd=self.m_params["vEnd"]
                    dV=self.m_params["dV"]
                    self.axes_spec.plot(np.arange(vStart,vEnd,dV),dataToPlot,label=str(self.nAvgSpectra)+" spectra averaged")
                else:
                    self.axes_spec.plot(dataToPlot,label=str(self.nAvgSpectra)+" spectra averaged")
                self.axes_spec.legend(loc=0)
                self.m_specWidget.draw()
            
        
    def updateWidgets(self):
        self.updateVoltageBox(self.m_params["vStart"],self.m_params["vEnd"],self.m_params["zPt"])
        self.drawXYMap(self.m_voltageBox.value())
        
    def connect(self):
        self.m_openButton.clicked.connect(self.loadCITS)
        self.m_bwdButton.toggled.connect(self.updateWidgets)
        self.m_voltageBox.valueChanged.connect(self.drawXYMap)
        self.m_mapWidget.mpl_connect('button_press_event',self.drawSpectrum)
        self.m_clearSpec.clicked.connect(self.clearSpectrum)
        self.m_scaleVoltage.toggled.connect(self.clearSpectrum)
        self.m_avgSpec.clicked.connect(self.launchAvgSpectrum)
        self.m_avgBox.toggled.connect(self.updateAvgVariables)
        #self.m_testButton.clicked.connect(self.launchAni)
        
    def drawSpectrum(self,event):
        if(event.xdata!=None and event.ydata!=None and self.dataLoaded):
            PixelX=int(event.xdata)
            PixelY=int(event.ydata)
            bwd=self.m_bwdButton.isChecked()
            fb=["_Fwd","_Bwd"]
            #Add data to the total data to average if needed
            if(self.m_avgBox.isChecked()):
                zPt=self.m_params["zPt"]
                #for v in range(0,zPt/2):
                self.tot_data+=self.m_data[bwd][PixelY][PixelX]
                self.nAvgSpectra+=1
                print(self.nAvgSpectra)
                #dataToPlot=self.tot_data/self.nAvgSpectra
            #Plot the data with the desired scale (Volts or index)
            dataToPlot=self.m_data[bwd][PixelY][PixelX]
            if(self.m_scaleVoltage.isChecked()):
                vStart=self.m_params["vStart"]
                vEnd=self.m_params["vEnd"]
                dV=self.m_params["dV"]
                self.axes_spec.plot(np.arange(vStart,vEnd,dV),dataToPlot,label="["+str(PixelX)+","+str(PixelY)+"]"+fb[bwd])
            else:
                self.axes_spec.plot(dataToPlot,label="["+str(PixelX)+","+str(PixelY)+"]"+fb[bwd])
            self.axes_spec.legend(loc=0)
            self.m_specWidget.draw()
            
    def launchAvgSpectrum(self):
        if(self.dataLoaded):
            xPx=self.m_params["xPx"]
            yPx=self.m_params["yPx"]
            self.averageSpectrum(0,0,xPx,yPx)
            
    def averageSpectrum(self,Xs,Ys,Xe,Ye):
        if(self.dataLoaded):
            zPt=self.m_params["zPt"]
            bwd=self.m_bwdButton.isChecked()
            fb=["_Fwd","_Bwd"]
            avg_data=np.zeros(shape=(zPt/2))
            N=(Ye-Ys)*(Xe-Xs) #Number of spectra averaged
            for y in range(Ys,Ye):
                for x in range(Xs,Xe):
                    for v in range(0,zPt/2):
                        avg_data[v]+=(self.m_data[bwd][y][x][v]/N)
            if(self.m_scaleVoltage.isChecked()):
                vStart=self.m_params["vStart"]
                vEnd=self.m_params["vEnd"]
                dV=self.m_params["dV"]
                self.axes_spec.plot(np.arange(vStart,vEnd,dV),avg_data,label="Avg"+fb[bwd])
            else:
                self.axes_spec.plot(avg_data,label="Avg"+fb[bwd])
            self.axes_spec.legend(loc=0)
            self.m_specWidget.draw()
        
            
    def clearSpectrum(self):
        self.axes_spec.clear()
        self.m_specWidget.draw()
        
    def clearMap(self):
        self.m_mapWidget.figure.clear()
        self.dataLoaded=False
        
    def drawXYMap(self,voltage):
        #Some operations
        fig_map=self.m_mapWidget.figure
        fig_map.clear()
        ax_map = fig_map.add_subplot(1,1,1)
        fig_map.subplots_adjust(left=0.125,right=0.95,bottom=0.15,top=0.92)
        xPx=self.m_params["xPx"]
        yPx=self.m_params["yPx"]
        zPt=self.m_params["zPt"]
        #Get the data of the map and draw it
        mapData,mapMin,mapMax=self.getMapData(voltage)
        XYmap=ax_map.pcolormesh(mapData)
        #XYmap=ax_map.imshow(mapData,origin="lower")
        #Set axis limits and title
        ax_map.axis([0,xPx,yPx,0])
        ax_map.set_title("CITS Map\n"
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
        self.m_mapWidget.draw()
        
    def readCITS(self,filepath):
        """ Reads a CITS file and store the parameters"""
        f=open(filepath)
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
            if("Start of Data" in line):
                break
        # Matlab convention : columns first then rows hence [y][x]
        self.m_data=np.zeros(shape=(2,yPx,xPx,zPt/2))
        for y in range(0,yPx):
            for x in range(0,xPx):
                #The line read is an array containing the dI/dV values indexed by the voltage index
                #String to remove the newline at the end and split to transform the string in a list
                #0 means Backwards=False so forward data, 1 means Backwards=True so backward data
                data_list=f.readline().strip().split()
                self.m_data[0][y][x]=data_list[0:zPt/2]
                #No need to reverse the backward data as it was from Vmin to Vmax in the file as the fwd data
                #self.m_data[1][y][x]=(data_list[zPt/2:zPt])[::-1]
                self.m_data[1][y][x]=(data_list[zPt/2:zPt])
        f.close()
        #Store the parameters in a dictonnary to use them later
        self.m_params={"xPx":xPx,"yPx":yPx,"xL":xL,"yL":yL,"zPt":zPt,"vStart":vStart,"vEnd":vEnd,"dV":2*abs(vEnd-vStart)/zPt}
        return True