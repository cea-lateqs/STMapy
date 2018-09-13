# -*- coding: utf-8 -*-
"""
Created on Mon Aug 03 12:04:19 2015

@author: LH242250
"""
import matplotlib
matplotlib.use('qt5agg')
from ui_citswidget import Ui_CitsWidget
from matplotlib.backends.backend_qt5 import NavigationToolbar2QT as NavigationToolbar
import numpy as np
import os.path as osp
import scipy as sp
import scipy.interpolate
import scipy.optimize
import scipy.signal
import matplotlib.pyplot as pyplot
import matplotlib.backend_bases
import struct
import Common.functions as fc
import PyQt5.QtWidgets as QtWidgets
from shape import Shape
from PyQt5.QtCore import QCoreApplication

# noinspection PyPep8Naming
class CitsWidget(QtWidgets.QMainWindow, Ui_CitsWidget):
    ### Building methods
    def __init__(self, parent):
        """ Builds the widget with parent widget in arg """
        QtWidgets.QMainWindow.__init__(self, parent)
        # matplotlib.pyplot.style.use('def')
        # Set up the user interface from Designer.
        self.setupUi(self)
        self.setAutoFillBackground(True)
        self.m_data = np.ndarray([])
        self.tot_data = np.ndarray([])
        self.nAvgSpectra = 0
        self.m_params = {}
        self.channelList = []
        self.mapName = ""
        # Set up figures
        self.toolbar_map = NavigationToolbar(self.m_mapWidget, self)
        self.toolbar_spec = NavigationToolbar(self.m_specWidget, self)
        self.map_layout.insertWidget(0, self.toolbar_map)
        self.spec_layout.insertWidget(0, self.toolbar_spec)
        self.fig_spec = self.m_specWidget.figure
        self.ax_spec = self.fig_spec.add_subplot(111)
        self.fig_spec.subplots_adjust(left=0.125, right=0.95, bottom=0.15, top=0.92)
        # Variables linked to clicks on map
        self.shapes_clicked = []
        self.voltageLine = 0
        # Colormaps
        self.m_colorBarBox.addItems(matplotlib.pyplot.colormaps())
        self.m_colorBarBox.setCurrentIndex(145)
        # Boolean that is True if a map is loaded
        self.dataLoaded = False
        # Other parameters used after map loading
        self.mapType = ""
        self.wdir = "/home/huderl/Bureau/Experiments/STM"
        self.fig_topo = 0
        self.topo = []
        # Connect all widgets
        self.connect()
        self.nSpectraDrawn = 0
        self.spectrumColor = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
        self.lastSpectrum = []
        # Set up layouts
        self.m_avgWidget.hide()
        self.m_cbarWidget.hide()
        self.m_fitUpperLabel.hide()
        self.m_fitLowerLabel.hide()
        self.m_fitUpperBox.hide()
        self.m_fitLowerBox.hide()
        # ­Calls the loading method at launch
        self.askCits()

    def connect(self):
        """ Connects all the signals. Only called in the constructor """
        self.m_openButton.clicked.connect(self.askCits)
        self.m_topoButton.clicked.connect(self.drawTopo)
        self.m_channelBox.currentIndexChanged.connect(self.updateMap)
        self.m_colorBarBox.currentIndexChanged.connect(self.updateMap)
        self.m_voltageBox.valueChanged.connect(self.drawXYMap)
        self.m_derivBox.stateChanged.connect(self.m_derivNBox.setEnabled)
        self.m_mapWidget.mpl_connect('button_press_event', self.onPressOnMap)
        self.m_mapWidget.mpl_connect('button_release_event', self.onReleaseOnMap)
        self.m_specWidget.mpl_connect('button_press_event', self.updateToPointedX)
        self.m_clearSpec.clicked.connect(self.clearSpectrum)
        self.m_scaleVoltage.toggled.connect(self.clearSpectrum)
        # self.m_scaleMetric.toggled.connect(self.updateMap)
        self.m_averageCitsButton.clicked.connect(self.avgSpectrasX)
        self.m_wholeLengthCutButton.clicked.connect(self.launchBigCut)
        self.m_normButton.clicked.connect(self.normalizeCurrentChannel)
        self.m_avgSpec.clicked.connect(self.launchAvgSpectrum)
        self.m_avgBox.toggled.connect(self.updateAvgVariables)
        self.m_vLineBox.toggled.connect(self.clearVoltageLine)
        self.m_avgVButton.clicked.connect(self.averageSpectrumWithValues)
        # Averaging with respect to values
        self.m_avgCheckBox.toggled.connect(self.m_avgWidget.setVisible)
        self.m_aboveBox.valueChanged.connect(self.updateAboveValue)
        self.m_belowBox.valueChanged.connect(self.updateBelowValue)
        # Cbar custom limits
        self.m_cbarCheckBox.toggled.connect(self.m_cbarWidget.setVisible)
        self.m_cbarCustomCheckbox.stateChanged.connect(self.m_cbarUpperBox.setEnabled)
        self.m_cbarCustomCheckbox.stateChanged.connect(self.m_cbarLowerBox.setEnabled)
        # Fit spectrum widget
        self.m_fitSpec.clicked.connect(self.fitSpectrum)
        self.m_fitCustomCheckbox.toggled.connect(self.m_fitUpperBox.setVisible)
        self.m_fitCustomCheckbox.toggled.connect(self.m_fitLowerBox.setVisible)
        self.m_fitCustomCheckbox.toggled.connect(self.m_fitUpperLabel.setVisible)
        self.m_fitCustomCheckbox.toggled.connect(self.m_fitLowerLabel.setVisible)

    ### Reading and loading CITS methods
    def askCits(self):
        """ Slot that only takes care of opening a file dialog for the user to select one or several CITS - Returns the CITS paths """
        cits_names_and_ext = QtWidgets.QFileDialog.getOpenFileNames(self,
                                                                    "Choose a CITS file to read or several to average",
                                                                    self.wdir,
                                                                    "3D binary file (*.3ds);;Ascii file (*.asc);;Text file (*.txt)")
        # getOpenFilesNames retunrs a tuple with the Cits_names as first element and extension as second. We just need the names.
        self.loadCits(cits_names_and_ext[0])

    def loadCits(self, cits_names):
        """ Slot that launches the reading of the CITS given in arguments. Having several CITS will prompt their averaging but they have to be of the same dimensions"""
        n_cits = len(cits_names)
        print(cits_names)
        if n_cits == 0:
            return
        first = True
        for cits in cits_names:
            print(cits)
            extension = cits.split('.')[-1]
            if extension == "asc":
                self.clearMap()
                self.mapType = "Omicron"
                self.dataLoaded = self.readCitsAscii(cits)
            elif extension == "3ds":
                self.clearMap()
                self.mapType = "Nanonis"
                self.dataLoaded = self.readCitsBin(cits)
            elif extension == "txt":
                self.readTopo(cits)
            else:
                print("Extension not recognized")
                self.m_statusBar.showMessage("Extension not recognized")
                return
            # After reading, check if the data was read correctly and update the working directory and the map name
            if self.dataLoaded:
                self.wdir = osp.dirname(cits)
                self.mapName = osp.basename(cits)
                print(self.mapName + " read as a " + self.mapType + " map")
                self.m_statusBar.showMessage(self.mapName + " read as a " + self.mapType + " map")
            else:
                print('Problem while reading ' + cits)
                self.m_statusBar.showMessage('Problem while reading ' + cits)
                return
            # If only one Cits was selected, there is no need to run the averaging code so return after drawing the topo and updating of the widgets
            if n_cits == 1:
                self.drawTopo()
                self.updateWidgets()
                return
            else:  # Else begin the averaging
                if first:  # If this was the first Cits to average, create the mean_data array to store the avergage
                    mean_data = self.m_data / n_cits
                    first = False
                else:  # Else, continue adding to the mean_data
                    mean_data += self.m_data / n_cits
        # If everthing went well and if there was several CITS chosen, clear the map and set the data to mean_data.
        self.mapName = "Average of " + str(n_cits) + " CITS"
        self.clearMap()
        self.dataLoaded = True
        self.m_data = mean_data
        self.updateWidgets()
        return

    def readCitsAscii(self, filepath):
        """ Reads an Ascii CITS file (Omicron) and stores all the parameters"""
        f = open(filepath)
        divider = 1
        unit = 1
        header_end_not_found = True
        for line in f:
            # Read the parameters of the map until "Start of Data"
            # Pixel dimensions in X
            if "x-pixels" in line:
                xPx = int(line.split()[-1])
            # Pixel dimensions in Y
            elif "y-pixels" in line:
                yPx = int(line.split()[-1])
            # Metric dimensions in X
            elif "x-length" in line:
                xL = float(line.split()[-1])
            # Metric dimensions in Y
            elif "y-length" in line:
                yL = float(line.split()[-1])
            # Number of points IN TOTAL therefore the division per 2 to have the number of points per channel
            elif "z-points" in line:
                zPt = int(line.split()[-1])
                # There is zPt/2 points for fwd and zPt/2 points for bwd
                zPt = zPt // 2
            # Starting voltage of the spectro
            elif "Device_1_Start" in line:
                vStart = round(float(line.split()[-2]), 6)
            # Final voltage of the spectro
            elif "Device_1_End" in line:
                vEnd = round(float(line.split()[-2]), 6)
            # Any eventual divider
            elif "divider" in line:
                divider = float(line.split()[-1])
            # Convert nV in V
            elif "value-unit = nV" in line:
                unit = 10 ** (-9)
            elif "Start of Data" in line:
                header_end_not_found = False
                break

        if header_end_not_found:
            print("Problem while reading the file : could not find ':HEADER END:' in file")
            f.close()
            return False
        # Matlab convention : columns first then rows hence [y][x]
        # In Omicron CITS, there is only two channels : fwd and bwd so it is read as such
        self.channelList = ["Data [Fwd]", "Data [Bwd]"]
        self.m_data = np.zeros(shape=(2, yPx, xPx, zPt))
        for y in range(0, yPx):
            for x in range(0, xPx):
                # The line read is an array containing the dI/dV (or I(V)) values indexed by the voltage index
                # Strip to remove the newline at the end and split to transform the string in a list
                data_list = f.readline().strip().split()
                # Forward data
                self.m_data[0][y][x] = data_list[0:zPt]
                self.m_data[0][y][x] = self.m_data[0][y][x] * unit
                # No need to reverse the backward data as it was from Vmin to Vmax in the file as the fwd data
                # Backward data
                self.m_data[1][y][x] = (data_list[zPt:2 * zPt])
                self.m_data[1][y][x] = self.m_data[1][y][x] * unit
        f.close()
        # Store the parameters in a dictonnary to use them later
        self.m_params = {"xPx": xPx, "yPx": yPx, "xL": xL, "yL": yL, "zPt": zPt, "vStart": vStart / divider,
                         "vEnd": vEnd / divider, "dV": abs(vEnd - vStart) / (divider * zPt)}
        if divider != 1:
            print("A divider of " + str(divider) + " was found and applied")

        # Check if a topo file exists and read it if yes
        topopath = osp.join(osp.dirname(filepath), 'Topo.txt')
        if osp.exists(topopath):
            self.readTopo(topopath)
        return True

    def readCitsBin(self, filepath):
        """ Reads a binary CITS file (Nanonis) and stores all the parameters"""
        # The divider is already taken into account by Nanonis during the experiment so no need to process it again*
        f = open(filepath, "rb")
        zSpectro = False
        half = False
        # Read the header of the map until its end ("HEADER_END")
        header_end_not_found = True
        for line in f:
            # Header lines can be treated as regular strings
            line = line.decode("utf-8")
            # Pixel dimensions
            if "Grid dim" in line:
                splitted_line = line.split('"')[1].split()
                xPx = int(splitted_line[0])
                yPx = int(splitted_line[-1])
            # Center coordinates and metric dimensions in nm (Grid settings also contains other data)
            elif "Grid settings" in line:
                xC = float(line.split(";")[0].split("=")[-1]) * (10 ** 9)
                yC = float(line.split(";")[1]) * (10 ** 9)
                xL = float(line.split(";")[-3]) * (10 ** 9)
                yL = float(line.split(";")[-2]) * (10 ** 9)
            elif "Sweep Signal" in line:
                if line.split('"')[1] == 'Z (m)':
                    zSpectro = True
            # Number of points per channel
            elif "Points" in line:
                zPt = int(line.split('=')[-1])
            # Channels recorded
            elif "Channels" in line:
                self.channelList = line.split('"')[1].split(';')
                nChannels = len(self.channelList)
            # Experiment parameters. Not used for now, only the number is recorded to skip the corresponding bytes afterwards
            elif "Experiment parameters" in line:
                nbExpParams = len(line.split(';'))
            # End of the header
            elif ":HEADER_END:" in line:
                header_end_not_found = False
                break

        if header_end_not_found:
            print("Problem while reading the file : could not find ':HEADER END:' in file")
            f.close()
            return False

        # Reading vStart and vEnd (floats of 4 bytes each)
        try:
            reading = struct.unpack('>' + 'f' * 2, f.read(8))
        except struct.error:
            print("Problem while reading the file : number of bytes to read different than what was expected")
            f.close()
            return False
        # If it is a Z-Spectroscopy, put the Z boundaries in nm
        if zSpectro:
            vStart = round(reading[0] * 10 ** 9, 6)
            vEnd = round(reading[1] * 10 ** 9, 6)
        else:
            vStart = round(reading[0], 6)
            vEnd = round(reading[1], 6)
        # Reading experiment parameters (nbExpParams*4 bytes)
        # f.read(nbExpParams*4)
        # Matlab convention : columns first then rows hence [y][x]
        try:
            self.topo = np.zeros(shape=(yPx, xPx))
            self.m_data = np.zeros(shape=(nChannels, yPx, xPx, zPt))
        except MemoryError:
            print("The data is too big ! Or the memory too small...\nI will take half of the channels...\n")
            half = True
        # If the first alloc didn't work, try to halve it
        if half:
            try:
                self.topo = np.zeros(shape=(yPx, xPx))
                self.m_data = np.zeros(shape=(nChannels / 2, yPx, xPx, zPt))
            except MemoryError:
                print("The data is REALLY too big ! Or the memory REALLY too small...\nI give up...\n")
                f.close()
                QtWidgets.QMessageBox.critical(self, 'Oops !',
                                               "The data is REALLY too big ! Or the memory REALLY too small...\nI give up...\n")
                return False
        # Format string for unpacking zPt big-endians floats ('>f')
        fmtString = '>' + 'f' * zPt
        # zPt floats to read of 4 bytes each
        bytesToRead = 4 * zPt
        y = 0
        while y < yPx:
            x = 0
            while x < xPx:
                chan = 0
                b = f.read(nbExpParams * 4)
                try:
                    self.topo[y][x] = struct.unpack('>' + 'f' * nbExpParams, b)[2]
                except struct.error:
                    print(
                        "Problem while reading the topo : number of bytes to read different than what was expected at " + str(
                            x) + " " + str(y))
                while chan < nChannels:
                    # Each channel is written successively by sequences of 4*zPt bytes. I then read these bytes and unpack them as big-endians floats ('>f')
                    b = f.read(bytesToRead)
                    try:
                        if not half or chan < nChannels / 2:
                            self.m_data[chan][y][x] = struct.unpack(fmtString, b)
                    except struct.error:
                        print(
                            "Problem while reading the file : number of bytes to read different than what was expected at" + str(
                                x) + " " + str(y) + " " + str(chan))
                        # Set chan,x,y to exit the loop
                        chan = nChannels
                        x = xPx
                        y = yPx
                    chan = chan + 1
                # After each loop over channels, a new "experiment" begins so I need to skip the vStart, vEnd and experiments parameters floats that are written once again before the channels
                f.read(8)
                # f.read(8+nbExpParams*4)
                x = x + 1
            y = y + 1
        f.close()
        if half:
            self.channelList = self.channelList[0:nChannels / 2]
        # Store the parameters in a dictonnary to use them later
        dV = (vEnd - vStart) / zPt
        self.m_params = {"xPx": xPx, "yPx": yPx, "xC": xC, "yC": yC, "xL": xL, "yL": yL, "zPt": zPt, "vStart": vStart,
                         "vEnd": vEnd, "dV": dV}
        print(self.m_params)
        # self.m_statusBar.showMessage(self.m_params)
        # Convert currents in nA
        for i in range(0, len(self.channelList)):
            chan = self.channelList[i]
            if "(A)" in chan:
                self.m_data[i] = np.abs(self.m_data[i]) * 10 ** 9
                self.channelList[i] = chan.replace("(A)", "(nA)")
        # Convert topo in nm
        self.topo = self.topo * 10 ** 9
        # Level topo
        self.topo = self.levelTopo()
        # Test
        if zSpectro:
            self.extractSlope(0.01, 0)
        # Post-processing
        # self.extractOutOfPhase(1,2)
        # pylab.close()
        return True

    ### Reading and loading topo images methods

    def readTopo(self, filepath):
        """ Reads a topography file (in test) """
        f = open(filepath)
        topo_data = []
        w = 0
        h = 0
        for line in f:
            # Treat the headers differently
            if line[0] == '#':
                if "Width" in line:
                    w = float(line.split()[-2])
                if "Height" in line:
                    h = float(line.split()[-2])
            else:
                topo_data.append(line.strip().split())
        f.close()
        self.topo = (np.asfarray(topo_data))
        return True

    #        #Set up figure
    #        xPx=len(topo_data[0])
    #        yPx=len(topo_data)
    #        self.fig_topo=pylab.figure()
    #        #Connect the close handling
    #        self.fig_topo.canvas.mpl_connect('close_event',self.handleClosingTopo)
    #        self.ax_topo=self.fig_topo.add_subplot(1,1,1,aspect=float(yPx)/xPx)
    #        self.fig_topo.subplots_adjust(left=0.125,right=0.95,bottom=0.15,top=0.92)
    #        topo_rescaled=(np.asfarray(topo_data))*10**(-9)
    #        if(w!=0 and h!=0):
    #            dx=w/xPx
    #            dy=h/yPx
    #            self.ax_topo.axis([0,w,0,h])
    #            self.ax_topo.invert_yaxis()
    #            # Plot data (y first (matlab convention)
    #            # pcolormesh takes *vertices* in arguments so the X (Y) array need to be from 0 to W (H) INCLUDED
    #            XYmap=self.ax_topo.pcolormesh(np.arange(0,w+dx,dx),np.arange(0,h+dy,dy),topo_rescaled,cmap=self.m_colorBarBox.currentText())
    #        else:
    #            self.ax_topo.axis([0,xPx,0,yPx])
    #            #self.ax_topo.invert_yaxis()
    #            XYmap=self.ax_topo.pcolormesh(topo_rescaled,cmap=self.m_colorBarBox.currentText())
    #        #Colorbar stuff
    #        cbar = self.fig_topo.colorbar(XYmap, shrink=.9, pad=0.05, aspect=15)
    #        cbar.ax.yaxis.set_ticks_position('both')
    #        cbar.ax.tick_params(axis='y', direction='in')
    #        #self.addChannel(data,"Z (nm)")

    def levelTopo(self):
        yPx = self.m_params["yPx"]
        xPx = self.m_params["xPx"]
        # Numpy array to save the leveled topo
        topo_leveled = np.zeros(shape=(yPx, xPx))
        fitX = np.arange(0, xPx)

        def fitF(z, a, b):
            return a * z + b

        for y in range(0, yPx):
            fitY = self.topo[y]
            f = sp.interpolate.InterpolatedUnivariateSpline(fitX, fitY, k=1)
            popt, pcov = sp.optimize.curve_fit(fitF, fitX, f(fitX))
            topo_leveled[y] = fitY - (popt[0] * fitX + popt[1])
        # Return the leveled topo
        return topo_leveled

    def drawTopo(self):
        """ Draws the topography read while opening the CITS."""
        if len(self.topo) != 0:
            yPx = self.m_params["yPx"]
            # If yPx==1, it is a Line Spectro so I need to call the specific method to draw the topo
            if yPx == 1:
                self.drawLineTopo()
                return

            lineFit = False
            # Line fitting if necessary
            if lineFit:
                self.topo = self.levelTopo()
            # Get other parameters
            w = self.m_params["xL"]
            h = self.m_params["yL"]
            xPx = len(self.topo[0])
            yPx = len(self.topo)

            # Set up the figure for the plot
            if self.fig_topo == 0:
                self.fig_topo = pyplot.figure()
            else:
                self.fig_topo.clear()
            self.ax_topo = self.fig_topo.add_subplot(1, 1, 1, aspect=float(yPx) / xPx)
            self.fig_topo.subplots_adjust(left=0.125, right=0.95, bottom=0.15, top=0.92)
            # Connect the close handling
            self.fig_topo.canvas.mpl_connect('close_event', self.handleClosingTopo)
            # Put the appropriate title
            if lineFit:
                self.fig_topo.suptitle("Leveled topo (line fit)")
            else:
                self.fig_topo.suptitle("Raw topo data")
            # Choose the colormap
            colormap = "afmhot"  # self.m_colorBarBox.currentText()

            if self.m_scaleMetric.isChecked():
                # If the map is an Omicron one, I have to invert the y-axis
                if self.mapType == "Omicron":
                    self.ax_topo.axis([0, w, h, 0])
                else:
                    self.ax_topo.axis([0, w, 0, h])
                # pcolormesh takes *vertices* in arguments so the X (Y) array need to be from 0 to W (H) INCLUDED 
                XYmap = self.ax_topo.pcolormesh(np.linspace(0, w, xPx + 1), np.linspace(0, h, yPx + 1), self.topo,
                                                cmap=colormap)
            else:
                # If the map is an Omicron one, I have to invert the y-axis
                if self.mapType == "Omicron":
                    self.ax_topo.axis([0, xPx, yPx, 0])
                else:
                    self.ax_topo.axis([0, xPx, 0, yPx])
                XYmap = self.ax_topo.pcolormesh(self.topo, cmap=colormap)
            # Colorbar stuff
            cbar = self.fig_topo.colorbar(XYmap, shrink=.9, pad=0.05, aspect=15)
            cbar.ax.yaxis.set_ticks_position('both')
            cbar.ax.tick_params(axis='y', direction='in')
            self.fig_topo.canvas.draw()

    def drawLineTopo(self):
        """ Draws the topography read while opening a Line Spectro """
        # Get parameters
        w = self.m_params["xL"]
        # h is 0 (line spectro) 
        xPx = self.m_params["xPx"]
        # yPx is 1 (line spectro)

        # Set up figure
        fig = pyplot.figure()
        self.ax_topo = fig.add_subplot(1, 1, 1)
        self.fig_topo = fig
        # Connect the close handling
        self.fig_topo.canvas.mpl_connect('close_event', self.handleClosingTopo)

        if self.m_scaleMetric.isChecked():
            self.ax_topo.plot(np.linspace(0, w, xPx), self.topo[0], label="Without line leveling")
            self.ax_topo.plot(np.linspace(0, w, xPx), self.levelTopo()[0], label="With line leveling")
            self.ax_topo.set_xlim(0, w)
            self.ax_topo.set_xlabel("Distance (nm)")
            self.ax_topo.set_ylabel("Z (nm)")
        else:
            self.ax_topo.set_xlim(0, xPx)
            self.ax_topo.plot(self.topo[0], label="Without line leveling")
            self.ax_topo.plot(self.levelTopo(), label="With line leveling")
            self.ax_topo.set_ylabel("Z (nm)")
        self.ax_topo.legend(loc=0)

    def handleClosingTopo(self, event):
        """ Called when the topo figure is closed - Put back self.fig_topo to 0 to indicate that no topo figure exists """
        print('Topo closing')
        self.fig_topo = 0
        return

    ### Updating methods. Usually called by signals
    def updateAvgVariables(self):
        """ Slot called by the toggling of the average box. Toggling on the box clears everything to start a new averaging. Toggling off averages and plots the data stored by picking spectra """
        if self.dataLoaded:
            # If toggled on, put everything to zero to be ready to store data
            if self.m_avgBox.isChecked():
                self.tot_data = np.zeros(shape=self.m_params["zPt"])
                self.nAvgSpectra = 0
            # If toggled off, plot the data stored
            else:
                if self.nAvgSpectra == 0: return
                dataToPlot = self.tot_data / self.nAvgSpectra
                self.drawSpectrum(dataToPlot, str(self.nAvgSpectra) + " spectra averaged")

    def updateAboveValue(self, value):
        """ Slot called when the above slider is changed. Changes the value of the textbox to reflect the change """
        self.m_aboveBar.setValue(value)
        if self.dataLoaded:
            self.m_aboveLine.setText(str(self.mapMax - value * (self.mapMax - self.mapMin) / 100))

    def updateBelowValue(self, value):
        """ Slot called when the below slider is changed. Changes the value of the textbox to reflect the change """
        self.m_belowBar.setValue(value)
        if self.dataLoaded:
            self.m_belowLine.setText(str(self.mapMin + value * (self.mapMax - self.mapMin) / 100))

    def updateMap(self):
        """ Updates the map by redrawing it. Updates also the above and below sliders """
        self.drawXYMap(self.m_voltageBox.value())
        self.updateAboveValue(self.m_aboveBar.value())
        self.updateBelowValue(self.m_belowBar.value())

    def updateToPointedX(self, event):
        """ Slot called when clicking on the spectrum window. Updates the map according to the position of the cursor when clicked """
        if event.xdata is not None and event.ydata is not None and self.dataLoaded and self.toolbar_spec._active is None:
            # If the scale is in volts, need to divide by dV to have the index
            if self.m_scaleVoltage.isChecked():
                pointedIndex = int((event.xdata - self.m_params["vStart"]) / self.m_params["dV"])
            else:
                pointedIndex = int(event.xdata)
            # Update the box
            self.m_voltageBox.setValue(pointedIndex)
            # The map updates itself because the change of value of the voltage box calls the drawXYMap method

    def updateVoltageBox(self, Vmin, Vmax, zPt):
        """ Method called by updateWidgets. Sets the values of the voltage box """
        # self.m_voltageBox.setMinimum(Vmin)
        # self.m_voltageBox.setMaximum(Vmax)
        # self.m_voltageBox.setSingleStep(dV)
        self.m_voltageBox.setMinimum(0)
        self.m_voltageBox.setMaximum(zPt - 1)
        self.m_voltageBox.setSingleStep(1)

    def updateWidgets(self):
        """ Slot called after the reading of the CITS. Sets the values combo box (voltage and channels) and draws the map """
        self.updateVoltageBox(self.m_params["vStart"], self.m_params["vEnd"], self.m_params["zPt"])
        self.m_channelBox.clear()
        self.m_channelBox.addItems(self.channelList)
        self.drawXYMap(self.m_voltageBox.value())
        self.updateAboveValue(self.m_aboveBar.value())
        self.updateBelowValue(self.m_belowBar.value())
        self.m_CitsAvgBox.setMaximum(self.m_params["xPx"])

    ### Methods related to spectra
    def normalizeSpectrum(self):
        """ Normalizes the spectra of displayed index """
        chan = self.m_channelBox.currentIndex()
        for x in range(self.m_params['xPx']):
            for y in range(self.m_params['yPx']):
                self.m_data[chan][x][y] = self.normalizeDOS(self.m_data[chan][x][y], self.m_params['zPt'])

    def normalizeDOS(self, dos, dos_length):
        mean_l = np.mean(dos[0:dos_length / 4])
        mean_r = np.mean(dos[3 * dos_length / 4:dos_length])
        mean = (mean_l + mean_r) / 2
        return dos / mean

    def averageSpectrum(self, xi, xf, yi, yf):
        """ Averages the spectra contained in the rectangle drawn by the points (xi,yi);(xf,yi);(xf,yf) and (xi,yf) """
        if self.dataLoaded:
            zPt = self.m_params["zPt"]
            chan = self.m_channelBox.currentIndex()
            avg_data = np.zeros(shape=zPt)
            if yf < yi:
                t = yf
                yf = yi
                yi = t
            if xf < xi:
                t = xf
                xf = xi
                xi = t
            N = (yf - yi) * (xf - xi)  # Number of spectra averaged
            for y in range(yi, yf):
                for x in range(xi, xf):
                    avg_data += (self.m_data[chan][y][x] / N)
            xavg = xi + (xf - xi) / 2
            yavg = yi + (yf - yi) / 2
            self.drawSpectrum(avg_data, "Average around " + "(" + str(xavg) + "," + str(yavg) + ")", xavg, yavg)
            return avg_data

    def averageSpectrumWithValues(self):
        """ Averages the spectra according their values at a certain voltage """
        if self.dataLoaded:
            voltage = self.m_voltageBox.value()
            chan = self.m_channelBox.currentIndex()
            viewSelected = self.m_viewSelectedBox.isChecked()
            zPt = self.m_params["zPt"]
            avg_data_aboveV = np.zeros(shape=zPt)
            avg_data_belowV = np.zeros(shape=zPt)
            # midpoint=(self.mapMax+self.mapMin)/2
            midpoint2 = (self.mapMax - self.mapMin)
            limit_aboveV = self.mapMax - self.m_aboveBar.value() * midpoint2 / 100
            limit_belowV = self.mapMin + self.m_belowBar.value() * midpoint2 / 100
            if limit_aboveV < limit_belowV:
                print("Above and below spectra intersect !")
                self.m_statusBar.showMessage("Above and below spectra intersect !")
                return
            N_aboveV = 0
            N_belowV = 0
            xPx = self.m_params["xPx"]
            yPx = self.m_params["yPx"]
            xPts = []
            yPts = []
            cPts = []
            for y in range(0, yPx):
                for x in range(0, xPx):
                    currentValue = self.m_data[chan][y][x][voltage]
                    if currentValue > limit_aboveV:
                        avg_data_aboveV += self.m_data[chan][y][x]
                        N_aboveV += 1
                        if viewSelected:
                            xPts.append(x)
                            yPts.append(y)
                            cPts.append(self.getSpectrumColor(self.nSpectraDrawn))
                    elif currentValue < limit_belowV:
                        avg_data_belowV += self.m_data[chan][y][x]
                        N_belowV += 1
                        if viewSelected:
                            xPts.append(x)
                            yPts.append(y)
                            cPts.append(self.getSpectrumColor(self.nSpectraDrawn + 1))
            if N_aboveV != 0:
                avg_data_aboveV /= N_aboveV
                self.drawSpectrum(avg_data_aboveV,
                                  "Average above " + str(limit_aboveV) + " (" + str(N_aboveV) + " spectra averaged)")
            if N_belowV != 0:
                avg_data_belowV /= N_belowV
                self.drawSpectrum(avg_data_belowV,
                                  "Average below " + str(limit_belowV) + " (" + str(N_belowV) + " spectra averaged)")
            if viewSelected: self.ax_map.plot(xPts, yPts, c=cPts, marker='o', ls='None')

    def clearSpectrum(self):
        """ Clears the spectrum window """
        self.ax_spec.clear()
        self.nSpectraDrawn = 0
        self.voltageLine = 0
        self.clearShapesClicked()
        # self.drawTopo()
        self.m_specWidget.draw()

    def drawSpectrum(self, dataToPlot, label, x=-1, y=-1):
        """ Method called each time a spectrum needs to be plotted. Takes care of the derivative and different scales stuff and updates the window """
        finalLabel = label + " - " + str(self.m_channelBox.currentText())
        shiftX = str(self.m_shiftXBox.text()).lower()
        if shiftX == "topo":
            shiftX = self.topo[y][x]
        else:
            shiftX = float(shiftX)
        shiftY = self.nSpectraDrawn * float(self.m_shiftYBox.text())
        if self.dataLoaded and dataToPlot.size != 0:
            dV = self.m_params["dV"]
            deriv = np.fabs(sp.signal.savgol_filter(dataToPlot, self.m_derivNBox.value(), 3, deriv=1, delta=dV))

            if self.m_logBox.isChecked():
                dataToPlot = np.log(dataToPlot)

            self.lastSpectrum = [dataToPlot, finalLabel]
            if self.m_scaleVoltage.isChecked():
                vStart = self.m_params["vStart"]
                vEnd = self.m_params["vEnd"]
                zPt = self.m_params["zPt"]
                V = np.arange(vStart, vEnd, dV) + shiftX
                # Check consistency of V array with the number of points.
                if len(V) != zPt:
                    print('Round-off error while computing the voltage array: dV ({}) might be too small. Computing from number of points instead.'.format(dV))
                    # If this fails, compute from the number of points to be sure to have the same dim as dataToPlot.
                    V = np.linspace(vStart, vEnd, zPt) + shiftX
                self.ax_spec.plot(V, shiftY + dataToPlot, label=finalLabel, c=self.getSpectrumColor(self.nSpectraDrawn))
                if self.m_derivBox.isChecked():
                    self.ax_spec.plot(V, shiftY + deriv, label="Derivative of " + finalLabel, marker='o', markersize=3.,
                                      c=self.getSpectrumColor(self.nSpectraDrawn))
                self.nSpectraDrawn = self.nSpectraDrawn + 1
            else:
                self.ax_spec.plot(shiftY + dataToPlot, label=finalLabel, c=self.getSpectrumColor(self.nSpectraDrawn))
                if self.m_derivBox.isChecked():
                    self.ax_spec.plot(shiftY + deriv, label="Derivative of " + finalLabel, marker='o', markersize=3.,
                                      c=self.getSpectrumColor(self.nSpectraDrawn))
                self.nSpectraDrawn = self.nSpectraDrawn + 1
        if not self.m_legendBox.isChecked():
            self.ax_spec.legend(loc=0)
        self.m_specWidget.draw()

    def fitSpectrum(self):
        """ Linear fitting of the last spectrum plotted that is stored in lastSpectrum """
        if self.dataLoaded and self.lastSpectrum[0].size != 0:
            def fit_func(v, a, b):
                return a * v + b
            dV = self.m_params["dV"]
            if self.m_fitCustomCheckbox.isChecked():
                limitL = float(self.m_fitLowerBox.text())
                limitU = float(self.m_fitUpperBox.text())
                vStart = min(limitL, limitU)
                vEnd = max(limitL, limitU)
                zPt = self.m_params["zPt"]
                xArray = np.arange(0, zPt) * dV + self.m_params["vStart"]
                # Take the portion to fit
                mask1 = xArray > vStart
                mask2 = xArray < vEnd
                temp = self.lastSpectrum[0][mask1]
                dataToFit = temp[mask2]
            else:
                vStart = self.m_params["vStart"]
                vEnd = self.m_params["vEnd"]
                dataToFit = self.lastSpectrum[0]
            X = np.arange(vStart, vEnd, dV)
            popt, pcov = sp.optimize.curve_fit(fit_func, X, dataToFit)
            slope = popt[0]
            coef = popt[1]
            print("Linear fit gives a slope of " + str(slope) + " and a coef of " + str(coef))
            self.ax_spec.plot(X, slope * X + coef, label="Linear fit of " + self.lastSpectrum[1],
                              color=self.getSpectrumColor(self.nSpectraDrawn - 1), linestyle="--")
            self.ax_spec.legend(loc=0)
            self.m_specWidget.draw()

    def getSpectrumColor(self, n):
        """ Returns the color corresponding to a spectrum of given index according to spectrumColor """
        i = n % len(self.spectrumColor)
        return self.spectrumColor[i]

    def launchAvgSpectrum(self):
        """ Slot called to average the spectra of the whole map """
        if self.dataLoaded:
            xPx = self.m_params["xPx"]
            yPx = self.m_params["yPx"]
            self.averageSpectrum(0, xPx, 0, yPx)

    def pickSpectrum(self, event):
        """ Method called when a press-and-release event is done at the same location of the map. If the average box is checked, it will keep the data in storage to average it later. Otherwise it plots the spectrum in the corresponding widget """
        if event.xdata is not None and event.ydata is not None and self.dataLoaded:
            PixelX = int(event.xdata)
            PixelY = int(event.ydata)
            chan = self.m_channelBox.currentIndex()
            if "Slope" in self.channelList[chan]:
                zPt = self.m_params["zPt"]
                zPts = np.arange(0, zPt)
                dataLogCurrent = np.zeros(shape=zPt)
                dataLine = self.m_data[chan][PixelY][PixelX][0] * self.m_params["dV"] * zPts + \
                           self.m_data[chan + 1][PixelY][PixelX][0]
                cutOff = False
                for z in zPts:
                    i = self.m_data[0][PixelY][PixelX][z]
                    if i < 0.01 or cutOff:
                        dataLogCurrent[z] = dataLine[z]
                    else:
                        dataLogCurrent[z] = np.log(i)
                self.drawSpectrum(dataLogCurrent, "Log of current at [" + str(PixelX) + "," + str(PixelY) + "]", PixelX, PixelY)
                # self.drawSpectrum(dataLine,"Linear fit of the log at "+str(PixelX)+","+str(PixelY)+"]")
            else:
                self.m_mapWidget.draw()
                # Add data to the total data to average if in average mod
                if self.m_avgBox.isChecked():
                    self.tot_data += self.m_data[chan][PixelY][PixelX]
                    self.nAvgSpectra += 1
                    # color='white'
                # Plot the data with the desired scale (Volts or index) if in normal mode
                else:
                    dataToPlot = self.m_data[chan][PixelY][PixelX]
                    self.drawSpectrum(dataToPlot, "[" + str(PixelX) + "," + str(PixelY) + "]", PixelX, PixelY)

    ### Methods related to the clicks on the map
    def onPressOnMap(self, event):
        """ Slot called when a press event is detected. Creates a Shape object that will be dynamically updated when the mouse moves """
        if (event.xdata is not None and event.ydata is not None) and (self.dataLoaded and self.toolbar_map._active is None):
            if self.m_scaleMetric.isChecked():
                ratioX = self.m_params['xL'] / self.m_params['xPx']
                ratioY = self.m_params['yL'] / self.m_params['yPx']
            else:
                ratioX = 1
                ratioY = 1
            self.currentShape = Shape(event, self.m_mapWidget.figure, self.fig_topo,
                                      self.getSpectrumColor(self.nSpectraDrawn), ratioX, ratioY)
            self.motionConnection = self.m_mapWidget.mpl_connect('motion_notify_event', self.currentShape.update)
            print('PRESS')

    def onReleaseOnMap(self, event):
        """ Slot called when a release event is detected. Disconnects the updating of the currentShape and launch the appropriate method depending on which button was pressed and where it was released """
        if self.dataLoaded and self.toolbar_map._active is None:
            if event.xdata is not None and event.ydata is not None:
                # Disconnects the updating of the current Shape and gets its coordinates
                self.m_mapWidget.mpl_disconnect(self.motionConnection)
                xi = self.currentShape.xi
                yi = self.currentShape.yi
                xf = self.currentShape.xf
                yf = self.currentShape.yf
                # If left-click : either a line was drawn or a spectrum picked
                if event.button == 1:
                    if event.xdata is None or event.ydata is None:
                        self.currentShape.remove()
                    # Cut along the XY line if a line is traced (X or Y different)
                    elif xf != xi or yf != yi:
                        self.cutAlongLine(xi, xf, yi, yf)
                    # Pick spectrum otherwise and change the line shape to a point
                    else:
                        self.currentShape.changeToPt()
                        self.pickSpectrum(event)
                # If right-click : either a rectangle was drawn or the center of the rectangle to average was picked
                else:
                    if event.xdata is not None and event.ydata is not None:
                        if xf != xi or yf != yi:
                            self.averageSpectrum(xi, xf, yi, yf)
                        # If X=Y, we need to force the updating of the Shape so it is drawn around the X,Y point and not starting at X,Y
                        else:
                            n = self.m_rcAvgBox.value()
                            self.currentShape.forceUpdate(max(0, xi - n), max(0, yi - n),
                                                          min(self.m_params["xPx"], xf + n),
                                                          min(self.m_params["yPx"], yf + n))
                            self.averageSpectrum(max(0, xi - n), min(self.m_params["xPx"], xf + n), max(0, yi - n),
                                                 min(self.m_params["yPx"], yf + n))
                # Add the current Shape to the list of clicked Shapes
                self.addToShapesClicked(self.currentShape)
            print('RELEASE')

    def addToShapesClicked(self, shape):
        """ Method called when a release was detected on the map. The Shape is saved in the shapes_clicked list """
        self.shapes_clicked.append(shape)

    def drawShapesClicked(self, fig_map, fig_topo, recreate):
        """ Method that draws all the Shapes saved in shapes_clicked or recreates them if the XYmap was cleared"""
        for shape in self.shapes_clicked:
            if recreate:
                shape.recreate(fig_map, fig_topo)
            else:
                shape.draw()

    def clearShapesClicked(self):
        """ Method that clears all saved Shapes and then redraws the map to reflect the change. No need to call Shape.remove as the XYmap will be cleared """
        self.shapes_clicked = []
        self.drawXYMap(self.m_voltageBox.value())

    def cutAlongLine(self, xi, xf, yi, yf):
        """ Method that finds the positions of the pixels forming the line from (xi,yi) to (xf,yf). The plotting of the spectra is done in cutPlot called at the end """
        # If the line is vertical, the equation is x=xi with y varying from yi to yf
        self.m_statusBar.showMessage(
            'Cut from ' + '(' + str(xi) + ',' + str(yi) + ')' + ' to ' + '(' + str(xf) + ',' + str(yf) + ') in ' + str(
                self.m_channelBox.currentText()))
        if xf == xi:
            if yi > yf:
                y_plot = np.flipud(np.arange(yf, yi + 1))
            else:
                y_plot = np.arange(yi, yf + 1)
            x_plot = np.full(shape=y_plot.size, fill_value=xi)
        else:
            # Simple algorithm for cuts
            simple = True
            if simple:
                # If the line is not vertical, determine its equation y=k*x+c
                k = float(yf - yi) / (xf - xi)
                c = yi - k * xi
                # Check if there is more y or more x to have to most precise arrangment
                if abs(xf - xi) > abs(yf - yi):
                    if xi > xf:
                        x_plot = np.flipud(np.arange(xf, xi + 1))
                    else:
                        x_plot = np.arange(xi, xf + 1)
                    y_plot = k * x_plot + c
                else:
                    if yi > yf:
                        y_plot = np.flipud(np.arange(yf, yi + 1))
                    else:
                        y_plot = np.arange(yi, yf + 1)
                    x_plot = (y_plot - c) / k
            # Bresenham algorithm
            else:
                x_plot_p = []
                y_plot_p = []
                dx = abs(xi - xf)
                dy = abs(yi - yf)
                x0 = xi
                y0 = yi
                if xi < xf:
                    sx = 1
                else:
                    sx = -1
                if yi < yf:
                    sy = 1
                else:
                    sy = -1
                err = dx - dy
                while True:
                    x_plot_p.append(x0)
                    y_plot_p.append(y0)
                    if x0 >= xf and y0 >= yf: break
                    e2 = err * 2
                    if e2 > -dy:
                        err -= dy
                        x0 += sx
                    if e2 < dx:
                        err += dx
                        y0 += sy
                x_plot = np.array(x_plot_p)
                y_plot = np.array(y_plot_p)
        return self.cutPlot(x_plot, y_plot)

    def cutPlot(self, x_plot, y_plot):
        """ Method called by cutAlongLine. Plots the spectra for the pixels of positions (x_plot[i],y_plot[i]) """
        # Build the data to plot with v as Y and z (number of pixels gone through) as X
        zPt = self.m_params['zPt']
        voltages = np.arange(0, zPt)
        chan = self.m_channelBox.currentIndex()
        z_plot = np.arange(0, x_plot.size)
        dataToPlot = np.ndarray(shape=(zPt, z_plot.size))
        xi = x_plot[0]
        yi = y_plot[0]
        zi = z_plot[0]
        zf = z_plot[-1]
        # Variables needed to compute the metric distances
        dx = self.m_params["xL"] / self.m_params["xPx"]
        dy = self.m_params["yL"] / self.m_params["yPx"]
        metricDistances = np.sqrt((dx * (x_plot - xi)) ** 2 + (dy * (y_plot - yi)) ** 2)
        # Bool to keep trace of the selected spectra
        viewSelected = self.m_viewSelectedBox.isChecked()
        # Matlab convention : Y (v) first then X (z)
        # Plot the built map in a new figure
        fig = pyplot.figure()
        # fig.canvas.mplconnect()
        ax = fig.add_subplot(1, 1, 1)
        ax.set_title(self.mapName.split(".")[0] + " - Cut " + str(fig.number))
        self.ax_map.text(xi + 0.5, yi + 0.5, str(fig.number))
        if self.m_waterfallButton.isChecked():
            if self.m_scaleVoltage.isChecked():
                voltages = self.m_params["vStart"] + voltages * self.m_params["dV"]
            for z in z_plot:
                xc = int(x_plot[z])
                yc = int(y_plot[z])
                spectrum = self.m_data[chan][yc][xc]
                offset = (zf - z) * float(self.m_shiftYBox.text())
                ax.plot(voltages, spectrum + offset, 'k', zorder=(z + 1) * 2)
                # Uncomment this to enable white filling under the curves
                ax.fill_between(voltages, spectrum + offset, offset, facecolor='w', lw=0, zorder=(z + 1) * 2 - 1)
                # if(z==zi): last_sp=np.amin(spectrum)
                # ax.fill_between(voltages, last_sp, spectrum+offset, facecolor='w', zorder=2*z)
                # last_sp=spectrum+offset
            ax.set_xlim([voltages[0], voltages[-1]])
        else:
            for v in voltages:
                for z in z_plot:
                    xc = int(x_plot[z])
                    yc = int(y_plot[z])
                    dataToPlot[v][z] = self.m_data[chan][yc][xc][v]  # /self.m_data[chan][yc][xc][0]
                    if viewSelected: self.addToPtsClicked(xc, yc, color='yellow')
            ax.set_ylabel("Voltage index")
            fig.subplots_adjust(left=0.125, right=0.95, bottom=0.15, top=0.92)
            # Pcolormesh takes vertices as arguments so need to add the last vertex to have the last quad plotted
            #            voltages=np.append(voltages,voltages[-1]+1)
            #            metricDistances=np.append(metricDistances,metricDistances[-1]+(metricDistances[1]-metricDistances[0]))
            #            z_plot=np.append(z_plot,z_plot[-1]+1)
            # Change the scales if needed
            if self.m_scaleVoltage.isChecked():
                voltages = self.m_params["vStart"] + voltages * self.m_params["dV"]
                ax.set_ylabel("Bias (V)")
            if self.m_scaleMetric.isChecked():
                mapData = pyplot.pcolormesh(metricDistances, voltages, dataToPlot, cmap=self.m_colorBarBox.currentText())
                ax.axis([metricDistances[0], metricDistances[-1], voltages[0], voltages[-1]])
                ax.set_xlabel("Distance (nm)")
            else:
                mapData = pyplot.pcolormesh(z_plot, voltages, dataToPlot, cmap=self.m_colorBarBox.currentText())
                ax.axis([z_plot[0], z_plot[-1], voltages[0], voltages[-1]])
                ax.set_xlabel("Pixels")
            # Colorbar set up
            cbar = fig.colorbar(mapData, shrink=.9, pad=0.05, aspect=15)
            cbar.ax.yaxis.set_ticks_position('both')
            cbar.ax.tick_params(axis='y', direction='in')
            if self.m_cbarCustomCheckbox.isChecked():
                mapData.set_clim(float(self.m_cbarLowerBox.text()), float(self.m_cbarUpperBox.text()))
            else:
                mapData.set_clim(0)
            return metricDistances, voltages, dataToPlot

    def launchBigCut(self):
        """ Launches a cut in the big diagonal """
        if self.dataLoaded:
            if self.m_scaleMetric.isChecked():
                ratioX = self.m_params['xL'] / self.m_params['xPx']
                ratioY = self.m_params['yL'] / self.m_params['yPx']
            else:
                ratioX = 1
                ratioY = 1
            # Simulates a left-click press at (0,0) and a release at xPx-1,yPx-1
            simEvent = matplotlib.backend_bases.MouseEvent('button_press_event', self.m_mapWidget.figure.canvas, 0, 0,
                                                           button=1)
            simEvent.xdata = 0
            simEvent.ydata = 0
            self.currentShape = Shape(simEvent, self.m_mapWidget.figure, self.fig_topo,
                                      self.getSpectrumColor(self.nSpectraDrawn), ratioX, ratioY)
            simEvent.xdata = self.m_params['xPx'] - 1
            simEvent.ydata = self.m_params['yPx'] - 1
            self.currentShape.update(simEvent)
            X, Y, Z = self.cutAlongLine(0, self.m_params["xPx"] - 1, 0, self.m_params["yPx"] - 1)
            self.m_mapWidget.draw()
            return X, Y, Z

    ### Methods related to the map
    def clearMap(self):
        """ Unloads the map and clears the map window """
        self.m_mapWidget.figure.clear()
        self.clearShapesClicked()
        self.m_data = np.ndarray([])
        # self.m_params.clear()
        self.mapType = ""
        self.dataLoaded = False

    def drawXYMap(self, voltage):
        """ Calls the getMapData function and draws the result in the map window with the approriate formatting. Called at each change of voltage """
        # Start everything anew
        fig_map = self.m_mapWidget.figure
        fig_map.clear()
        if self.dataLoaded:
            xPx = self.m_params["xPx"]
            yPx = self.m_params["yPx"]
            # zPt=self.m_params["zPt"]
            if xPx == yPx:
                self.ax_map = fig_map.add_subplot(1, 1, 1, aspect='equal')
            else:
                self.ax_map = fig_map.add_subplot(1, 1, 1)
            # self.ax_map.hold(True)
            fig_map.subplots_adjust(left=0.125, right=0.95, bottom=0.15, top=0.92)
            # Get the data of the map and draw it
            mapData, self.mapMin, self.mapMax = self.getMapData(voltage)
            # Use metric dimensions if the corresponding box is checked (work in progress)
            if self.m_scaleMetric.isChecked() and False:
                xL = self.m_params["xL"]
                yL = self.m_params["yL"]
                x_m = np.linspace(0, xL, xPx)
                y_m = np.linspace(0, yL, yPx)
                XYmap = self.ax_map.pcolormesh(x_m, y_m, mapData, cmap=self.m_colorBarBox.currentText())
                self.ax_map.axis([0, xL, 0, yL])
            # Else, use pixels
            else:
                XYmap = self.ax_map.pcolormesh(mapData, cmap=self.m_colorBarBox.currentText())
                # If the map is an Omicron one, I have to invert the y-axis
                if self.mapType == "Omicron":
                    self.ax_map.axis([0, xPx, yPx, 0])
                else:
                    self.ax_map.axis([0, xPx, 0, yPx])
            # Set title
            chan = self.m_channelBox.currentText()
            # print(chan)
            self.ax_map.set_title(self.mapName + " - " + chan + "\n"
                                  + "V=" + str(self.m_params["vStart"] + voltage * self.m_params["dV"]))
            # Colorbar stuff
            cbar = fig_map.colorbar(XYmap, shrink=.9, pad=0.05, aspect=15)
            cbar.ax.yaxis.set_ticks_position('both')
            cbar.ax.tick_params(axis='y', direction='in')
            # Image color scale is adjusted to the data:
            if self.m_cbarCustomCheckbox.isChecked():
                XYmap.set_clim(float(self.m_cbarLowerBox.text()), float(self.m_cbarUpperBox.text()))
            else:
                XYmap.set_clim(self.mapMin, self.mapMax)
            # Plot a dashed line at X=voltage if asked
            if self.m_vLineBox.isChecked():
                self.drawVoltageLine(voltage)
            # Recreate the shapes points saved in shapes_clicked
            self.drawShapesClicked(fig_map, self.fig_topo, True)
            self.m_mapWidget.draw()

    def getMapData(self, v):
        """ Returns an array built from the data loaded that can be used to display a map at fixed voltage """
        xPx = self.m_params["xPx"]
        yPx = self.m_params["yPx"]
        mapData = np.ndarray(shape=(yPx, xPx))
        chan = self.m_channelBox.currentIndex()
        valMin = np.inf
        valMax = -np.inf
        for y in range(0, yPx):
            for x in range(0, xPx):
                val = self.m_data[chan][y][x][v]
                mapData[y][x] = val
                if val < valMin:
                    valMin = val
                elif val > valMax:
                    valMax = val
        return mapData, valMin, valMax

    ### Methods related to the voltage guide line in the spectra window

    def drawVoltageLine(self, voltage):
        """ Draws the vertical line at the given voltage """
        self.clearVoltageLine()
        if self.dataLoaded:
            # Get the current voltage : real voltage if the scale box is checked, voltage index otherwise
            if self.m_scaleVoltage.isChecked:
                currentV = self.m_params["vStart"] + voltage * self.m_params["dV"]
            else:
                currentV = voltage
            # Plot the dashed line
            self.voltageLine = self.ax_spec.axvline(currentV, color='k', linestyle='--')
            self.m_specWidget.draw()

    def clearVoltageLine(self):
        """ Removes the vertical voltage line """
        if self.voltageLine:
            self.ax_spec.lines.pop(self.ax_spec.lines.index(self.voltageLine))
            self.m_specWidget.draw()
        self.voltageLine = 0

    ### Post-processing methods
    def addChannel(self, newChannelData, channelName):
        """ Adds a channel to the whole data array """
        # Gets parameters for reshaping
        nChannels = len(self.channelList)
        xPx = self.m_params["xPx"]
        yPx = self.m_params["yPx"]
        zPt = self.m_params["zPt"]
        try:
            newData = np.zeros(shape=(nChannels + 1, yPx, xPx, zPt))
        except MemoryError:
            print("The data is too big ! Or the memory too small...")
            return
        newData[0:nChannels] = self.m_data
        newData[nChannels] = newChannelData
        self.m_data = newData
        del newData
        # Updates the channel list
        self.channelList.append(channelName)
        self.m_channelBox.clear()
        self.m_channelBox.addItems(self.channelList)

    def avgSpectrasX(self):
        """ Slot that averages the spectras in the X direction of the loaded CITS and replaces the loaded CITS by the result """
        if self.dataLoaded:  # Check if a CITS was loaded
            # Get the needed params
            Navg = self.m_CitsAvgBox.value()
            xPx = self.m_params['xPx']
            yPx = self.m_params['yPx']
            zPt = self.m_params['zPt']
            # Try the allocation
            try:
                new_data = np.zeros(shape=(len(self.channelList), yPx, xPx / Navg, zPt))
            except MemoryError:
                print('Not enough memory to do this operation')
                return
            # Average in X for each channel and Y
            for chan in range(0, len(self.channelList)):
                for y in range(0, yPx):
                    # To average, for each x, add every spectrum between x and x+Navg (averaging window) divided by Navg.
                    for x in range(0, xPx, Navg):
                        if (
                                x + Navg > xPx):  # If the averaging window is not contained in the CITS data, stop averaging. The last spectras of this window will then will dismissed.
                            break
                        else:  # Else, average by adding the spectras of the avergaging window and dividing by Navg
                            spectra = np.zeros(zPt)
                            for i in range(0, Navg):
                                spectra = spectra + self.m_data[chan][y][x + i] / Navg
                                # Store the result in new_data
                                new_data[chan][y][x / Navg] = spectra
            # If eveything went well, clear the current map and replace it by the created data and change xPx
            self.mapName = "Average_" + str(Navg)
            self.clearMap()
            self.m_params['xPx'] = xPx / Navg
            self.dataLoaded = True
            self.m_data = new_data
            self.updateWidgets()

    def extractOutOfPhase(self, numChanR, numChanPhi):
        """ Extracts the out of phase component and adds it to the data """
        # Phase : 9V = pi
        outOfPhase = -self.m_data[numChanR] * np.cos(self.m_data[numChanPhi] * np.pi / 9)
        # Add the channel to the data
        self.addChannel(outOfPhase, self.channelList[numChanR] + "cos(" + self.channelList[numChanPhi] + ")")

    def extractDerivative(self, numChanToDeriv):
        """ Extracts the derivative of a channel and adds it to the data """
        dV = self.m_params["dV"]
        yPx = self.m_params["yPx"]
        xPx = self.m_params["xPx"]
        zPt = self.m_params["zPt"]
        derivData = np.zeros(shape=(yPx, xPx, zPt))
        for y in range(0, yPx):
            for x in range(0, xPx):
                derivData[y][x] = sp.signal.savgol_filter(self.m_data[numChanToDeriv][y][x], 9, 2, deriv=1, delta=dV)
        # Add the channel to the data
        self.addChannel(derivData, "Derivative of " + self.channelList[numChanToDeriv])

    def extractFFT(self, numChanToFFT, axisOfFFT):
        """ Extracts the FFT of a channel and adds it to the data """
        yPx = self.m_params["yPx"]
        xPx = self.m_params["xPx"]
        zPt = self.m_params["zPt"]
        FFTData = np.fft.fft(self.m_data[numChanToFFT], axis=axisOfFFT)
        # Add the channel to the data
        self.addChannel(FFTData, "FFT of " + self.channelList[numChanToFFT])

    def computeAngle(self, Dmoire, k=True):
        """ Computes the twist angle of a given graphene moiré of period Dmoire """
        return 2 * np.arcsin(0.246 / (2 * Dmoire)) * 180 / np.pi

    def extractSlope(self, cutOffValue, numChanToFit):
        """ Do a linear fit of the data in the asked channel and add the slope and the coef found as channels (usually called for zSpectros) """
        yPx = self.m_params["yPx"]
        xPx = self.m_params["xPx"]
        zPt = self.m_params["zPt"]
        dZ = self.m_params["dV"]
        zg = np.zeros(shape=(yPx, xPx, zPt))
        slopeData = np.zeros(shape=(yPx, xPx, zPt))
        coefData = np.zeros(shape=(yPx, xPx, zPt))
        fit_func = lambda v, a, b: a * v + b
        xArray = np.arange(0, zPt) * dZ
        for y in range(0, yPx):
            for x in range(0, xPx):
                rawData = self.m_data[numChanToFit][y][x]
                mask = rawData > cutOffValue
                xArrayFiltered = xArray[mask]
                dataFiltered = np.log(rawData[mask])
                popt, pcov = sp.optimize.curve_fit(fit_func, xArrayFiltered, dataFiltered)
                zg[y][x] = 1 / 20.5 * np.log(rawData) + np.arange(0, zPt * dZ, dZ) + self.topo[y][x]
                for z in range(0, zPt):
                    slopeData[y][x][z] = popt[0]
                    coefData[y][x][z] = popt[1]
        # Add the created channel to the data
        self.addChannel(slopeData, "Slope by linear fit of " + self.channelList[numChanToFit])
        self.addChannel(coefData, "Coef by linear fit of " + self.channelList[numChanToFit])
        self.addChannel(zg, "Zg")

    def normalizeCurrentChannel(self):
        if self.dataLoaded:
            yPx = self.m_params["yPx"]
            xPx = self.m_params["xPx"]
            zPt = self.m_params["zPt"]
            chan = self.m_channelBox.currentIndex()
            normData = np.zeros(shape=(yPx, xPx, zPt))
            for y in range(0, yPx):
                for x in range(0, xPx):
                    normData[y][x] = fc.normalizeLDOS(self.m_data[chan][y][x], 10)
            self.addChannel(normData, "Normalized " + self.channelList[chan])
