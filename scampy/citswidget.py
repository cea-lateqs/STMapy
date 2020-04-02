# -*- coding: utf-8 -*-

import logging
import os.path
import numpy as np
import scipy as sp
import scipy.optimize
import scipy.signal
import scipy.io
import matplotlib
from matplotlib.backends.backend_qt5 import NavigationToolbar2QT as NavigationToolbar
from matplotlib import pyplot
from PyQt5 import QtWidgets, uic
from .shape import generateShape, changeToDot
from .reads import readCitsAscii, readCits3dsBin, readCitsSm4Bin
from .processing import levelTopo, normalizeDOS, linear_fit_function, findPixelsOnLine

# Set explictly the backend to Qt for consistency with pyqt.
matplotlib.use("qt5agg")


# noinspection PyPep8Naming
class CitsWidget(QtWidgets.QMainWindow):
    # %% Building methods
    def __init__(self, parent):
        """ Builds the widget with parent widget in arg """
        QtWidgets.QMainWindow.__init__(self, parent)
        uic.loadUi(os.path.join(os.path.dirname(__file__), "ui_citswidget.ui"), self)
        # Set up the user interface from Designer.
        self.setAutoFillBackground(True)
        # Initiate arrays
        self.m_data = np.ndarray([])
        self.tot_data = np.ndarray([])
        self.nAvgSpectra = 0
        self.m_params = {}
        self.channelList = []
        self.mapName = ""
        # Set up figures
        self.toolbar_map = NavigationToolbar(self.m_mapWidget, self)
        self.toolbar_spec = NavigationToolbar(self.m_specWidget, self)
        self.m_saveCsvButton = QtWidgets.QPushButton(
            "Save as CSV"
        )  # Add csv save button
        self.toolbar_spec.addWidget(self.m_saveCsvButton)
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
        # Boolean that is True if a map is loaded
        self.dataLoaded = False
        # Read config to set wdir and matplotlib stylesheet
        self.wdir = ""
        self.readConfig()
        # Set default colormap
        index = self.m_colorBarBox.findText(self.default_cmap)
        if index >= 0:
            self.m_colorBarBox.setCurrentIndex(index)
        else:
            logging.warning(
                "{0} is not a valid colormap name. Cannot be set as default.".format(
                    self.default_cmap
                )
            )
        # Other parameters used after map loading
        self.mapType = ""
        self.fig_topo = 0
        self.topo = []
        self.topo_colormap = "afmhot"
        # Connect all widgets
        self.connect()
        self.nSpectraDrawn = 0
        self.spectrumColor = ["b", "g", "r", "c", "m", "y", "k"]
        self.lastSpectrum = []
        # Set up layouts
        self.m_avgWidget.hide()
        self.m_cbarWidget.hide()
        self.m_fitUpperLabel.hide()
        self.m_fitLowerLabel.hide()
        self.m_fitUpperBox.hide()
        self.m_fitLowerBox.hide()
        # Check 'Colorbar settings' by default
        self.m_cbarCheckBox.setChecked(True)
        # Calls the loading method at launch
        if self.autoload:
            self.askCits()

    def connect(self):
        """ Connects all the signals. Only called in the constructor """
        self.m_openButton.clicked.connect(self.askCits)
        self.m_topoButton.clicked.connect(self.drawTopo)
        self.m_channelBox.currentIndexChanged.connect(self.updateMap)
        self.m_colorBarBox.currentIndexChanged.connect(self.updateMap)
        self.m_voltageBox.valueChanged.connect(self.drawXYMap)
        self.m_derivBox.stateChanged.connect(self.m_derivNBox.setEnabled)
        self.m_mapWidget.mpl_connect("button_press_event", self.onPressOnMap)
        self.m_mapWidget.mpl_connect("button_release_event", self.onReleaseOnMap)
        self.m_specWidget.mpl_connect("button_press_event", self.updateToPointedX)
        self.m_saveCsvButton.clicked.connect(self.saveSpectra)
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

    # %% Reading and loading CITS methods
    def askCits(self):
        """ Slot that only takes care of opening a file dialog
        for the user to select one or several CITS - Returns the CITS paths """
        cits_names_and_ext = QtWidgets.QFileDialog.getOpenFileNames(
            self,
            "Choose a CITS file to read or several to average",
            self.wdir,
            " RHK file (*.sm4);;Matlab file (*.mat);;3D binary file (*.3ds);;Ascii file (*.asc);;Text file (*.txt)",
        )
        # getOpenFilesNames retunrs a tuple with Cits_names as first element
        # and extension as second. We just need the names.
        self.loadCits(cits_names_and_ext[0])

    def loadCits(self, cits_names):
        """ Slot that launches the reading of the CITS given in arguments.
        Having several CITS will prompt their averaging
        but they have to be of the same dimensions"""
        n_cits = len(cits_names)
        logging.debug("CITS to load: {0}".format(cits_names))
        if n_cits == 0:
            return
        first = True
        for cits in cits_names:
            logging.info("Loading {0}...".format(cits))
            extension = cits.split(".")[-1]
            if extension == "asc":
                self.clearMap()
                self.mapType = "Omicron"
                self.topo, self.m_data, self.channelList, self.m_params = readCitsAscii(
                    cits
                )
                self.dataLoaded = True
            elif extension == "3ds":
                self.clearMap()
                self.mapType = "Nanonis"
                # The zSpectro bool has to be checked in readCits3dsBin
                zSpectro = False
                if zSpectro:
                    (
                        self.topo,
                        self.m_data,
                        self.channelList,
                        self.m_params,
                        slopeData,
                        slopeDataName,
                        coefData,
                        coefDataName,
                        zg,
                    ) = readCits3dsBin(cits, zSpectro)
                    self.addChannel(slopeData, slopeDataName)
                    self.addChannel(coefData, coefDataName)
                    self.addChannel(zg, "Zg")
                else:
                    (
                        self.topo,
                        self.m_data,
                        self.channelList,
                        self.m_params,
                    ) = readCits3dsBin(cits, zSpectro)
                self.dataLoaded = True
            elif extension == "txt":
                self.readTopo(cits)
            elif extension == "" or extension == "mat":
                self.clearMap()
                self.mapType = "Sm4 to .mat"
                self.dataLoaded = self.loadCitsSm4(cits)
            elif extension == "sm4":
                self.clearMap()
                self.mapType = "Sm4"
                (
                    self.topo,
                    self.m_data,
                    self.channelList,
                    self.m_params,
                    average,
                ) = readCitsSm4Bin(cits)
                addAverage = True
                if addAverage:
                    self.addChannel(average, "average")
                self.dataLoaded = True
            else:
                logging.error(
                    "Extension {0} not supported ! Only asc, 3ds, txt and mat are valid.".format(
                        extension
                    )
                )
                return
            # After reading, check if the data was read correctly and update the working directory and the map name
            if self.dataLoaded:
                self.wdir = os.path.dirname(cits)
                self.mapName = os.path.basename(cits)
                logging.info(self.mapName + " read as a " + self.mapType + " map")
            else:
                logging.error("Problem while reading " + cits)
                return
            # If only one Cits was selected, there is no need to run the averaging code so return after drawing the topo
            # and updating of the widgets
            if n_cits == 1:
                self.updateWidgets()
                self.drawTopo()
                return
            else:  # Else begin the averaging
                if (
                    first
                ):  # If this was the first Cits to average, create the mean_data array to store the avergage
                    mean_data = self.m_data / n_cits
                    first = False
                else:  # Else, continue adding to the mean_data
                    mean_data += self.m_data / n_cits
        # If everthing went well and if there was several CITS chosen,
        # clear the map and set the data to mean_data.
        self.mapName = "Average of {0} CITS".format(n_cits)
        self.clearMap()
        self.dataLoaded = True
        self.m_data = mean_data
        self.updateWidgets()
        return

    def readConfig(self):
        config = {}
        with open("config.txt") as f:
            for line in f:
                (key, val) = line.split()
                config[key] = val
        if "working_directory" in config.keys():
            if os.path.exists(config["working_directory"]):
                self.wdir = config["working_directory"]
            else:
                logging.warning(
                    "{} is not a valid path ! Check your config file.".format(
                        config["working_directory"]
                    )
                )

        if "matplotlib_stylesheet" in config.keys():
            try:
                matplotlib.pyplot.style.use(config["matplotlib_stylesheet"])
            except IOError:
                logging.warning(
                    "{} was not found in the .matplotlib folder. Using default parameters for matplotlib...".format(
                        config["matplotlib_stylesheet"]
                    )
                )

        if "autoload" in config.keys():
            autoload = config["autoload"].lower()
            if "no" in autoload or "false" in autoload:
                self.autoload = False
            else:
                self.autoload = True
        else:
            self.autoload = False

        if "default_cmap" in config.keys():
            self.default_cmap = config["default_cmap"]
        else:
            self.default_cmap = "magma_r"

    # %% Reading and loading topo images methods

    def drawTopo(self, lineFit=True):
        """ Draws the topography read while opening the CITS."""
        if self.topo.ndim == 2:
            # Get parameters
            xPx, yPx = self.topo.shape
            # If m_params["yPx"] == 1, it is a Line Spectro so I need to call the specific
            # method to draw the topo
            if self.m_params["yPx"] == 1:
                self.drawLineTopo()
                return

            # Line fitting if necessary
            if lineFit:
                self.topo = levelTopo(self.topo)
                figtitle = "Leveled topo (line fit)"
            else:
                figtitle = "Raw topo data"
            # Set up the figure for the plot
            if self.fig_topo == 0:
                self.fig_topo = pyplot.figure()
            else:
                self.fig_topo.clear()
            self.ax_topo = self.fig_topo.add_subplot(1, 1, 1, aspect=yPx / xPx)
            self.fig_topo.subplots_adjust(left=0.125, right=0.95, bottom=0.15, top=0.92)
            # Connect the close handling
            self.fig_topo.canvas.mpl_connect("close_event", self.handleClosingTopo)
            # Put the appropriate title
            self.fig_topo.suptitle(figtitle)
            # Choose the colormap
            colormap = self.topo_colormap

            if self.m_scaleMetric.isChecked():
                logging.debug("Scale metric box is checked")
                max_x = self.m_params["xL"]
                max_y = self.m_params["yL"]
            else:
                max_x = xPx
                max_y = yPx
                logging.debug("Scale metric box is checked")
            # If the map is an Omicron one, I have to invert the y-axis
            if self.mapType == "Omicron":
                self.ax_topo.axis([0, max_x, max_y, 0])
                # pcolormesh takes *vertices* in arguments
                # so the X (Y) array need to be from 0 to W (H) INCLUDED
                XYmap = self.ax_topo.pcolormesh(
                    np.linspace(0, max_x, xPx + 1),
                    np.linspace(0, max_y, yPx + 1),
                    self.topo,
                    cmap=colormap,
                )
            elif self.mapType == "Sm4 to .mat":
                XYmap = self.ax_topo.imshow(self.topo, extent=[0, max_x, 0, max_y])
                # We plot our spec location points.
                locations = True
                if locations:
                    logging.debug("Spectrum Locations will be printed")
                    patch = self.m_params["Patch"]
                    for m in range(len(patch)):  # iterate over the number of locations
                        self.ax_topo.add_patch(patch[m])
            else:
                self.ax_topo.axis([0, max_x, 0, max_y])
                # pcolormesh takes *vertices* in arguments so the X (Y) array need to be from 0 to W (H) INCLUDED
                XYmap = self.ax_topo.pcolormesh(
                    np.linspace(0, max_x, xPx + 1),
                    np.linspace(0, max_y, yPx + 1),
                    self.topo,
                    cmap=colormap,
                )

            # Colorbar stuff
            cbar = self.fig_topo.colorbar(XYmap, shrink=0.9, pad=0.05, aspect=15)
            cbar.ax.yaxis.set_ticks_position("both")
            cbar.ax.tick_params(axis="y", direction="in")
            self.fig_topo.canvas.draw()
            self.fig_topo.show()

    def drawLineTopo(self):
        """ Draws the topography read while opening a Line Spectro """
        xPx = self.m_params["xPx"]
        # yPx is 1 (line spectro)

        # Set up figure
        fig = pyplot.figure()
        self.ax_topo = fig.add_subplot(1, 1, 1)
        self.fig_topo = fig
        # Connect the close handling
        self.fig_topo.canvas.mpl_connect("close_event", self.handleClosingTopo)

        if self.m_scaleMetric.isChecked():
            max_x = self.m_params["xL"]
            x_label = "Distance (nm)"
        else:
            max_x = xPx
            x_label = ""
        self.ax_topo.plot(
            np.linspace(0, max_x, xPx), self.topo[0], label="Without line leveling"
        )
        self.ax_topo.plot(
            np.linspace(0, max_x, xPx),
            levelTopo(self.topo)[0],
            label="With line leveling",
        )
        self.ax_topo.set_xlim(0, max_x)
        self.ax_topo.set_xlabel(x_label)
        self.ax_topo.set_ylabel("Z (nm)")
        self.ax_topo.legend(loc=0)
        self.fig_topo.show()

    def handleClosingTopo(self, event):
        """ Called when the topo figure is closed.
        Put back self.fig_topo to 0 to indicate that no topo figure exists.
        """
        logging.debug("Topo closing")
        self.fig_topo = 0
        return

    # %% Updating methods. Usually called by signals
    def updateAvgVariables(self):
        """ Slot called by the toggling of the average box.
        Toggling on the box clears everything to start a new averaging.
        Toggling off averages and plots the data stored by picking spectra """
        if self.dataLoaded:
            # If toggled on, put everything to zero to be ready to store data
            if self.m_avgBox.isChecked():
                self.tot_data = np.zeros(shape=self.m_params["zPt"])
                self.nAvgSpectra = 0
            # If toggled off, plot the data stored
            else:
                if self.nAvgSpectra == 0:
                    return
                dataToPlot = self.tot_data / self.nAvgSpectra
                self.drawSpectrum(
                    dataToPlot, str(self.nAvgSpectra) + " spectra averaged"
                )

    def updateAboveValue(self, value):
        """ Slot called when the above slider is changed.
        Changes the value of the textbox to reflect the change """
        self.m_aboveBar.setValue(value)
        if self.dataLoaded:
            self.m_aboveLine.setText(
                str(self.mapMax - value * (self.mapMax - self.mapMin) / 100)
            )

    def updateBelowValue(self, value):
        """ Slot called when the below slider is changed.
        Changes the value of the textbox to reflect the change """
        self.m_belowBar.setValue(value)
        if self.dataLoaded:
            self.m_belowLine.setText(
                str(self.mapMin + value * (self.mapMax - self.mapMin) / 100)
            )

    def updateMap(self):
        """ Updates the map by redrawing it.
        Updates also the above and below sliders """
        self.drawXYMap(self.m_voltageBox.value())
        self.updateAboveValue(self.m_aboveBar.value())
        self.updateBelowValue(self.m_belowBar.value())

    def updateToPointedX(self, event):
        """ Slot called when clicking on the spectrum window.
        Updates map according to the position of the cursor when clicked """
        if (
            event.xdata is not None
            and event.ydata is not None
            and self.dataLoaded
            and self.toolbar_spec._active is None
        ):
            # If the scale is in volts, need to divide by dV to have the index
            if self.m_scaleVoltage.isChecked():
                pointedIndex = int(
                    (event.xdata - self.m_params["vStart"]) / self.m_params["dV"]
                )
            else:
                pointedIndex = int(event.xdata)
            # Update the box
            self.m_voltageBox.setValue(pointedIndex)
            # The map updates itself because the change of value of the voltage box calls the drawXYMap method

    def updateVoltageBox(self, Vmin, Vmax, zPt):
        """ Method called by updateWidgets.
        Sets the values of the voltage box """
        # self.m_voltageBox.setMinimum(Vmin)
        # self.m_voltageBox.setMaximum(Vmax)
        # self.m_voltageBox.setSingleStep(dV)
        self.m_voltageBox.setMinimum(0)
        self.m_voltageBox.setMaximum(zPt - 1)
        self.m_voltageBox.setSingleStep(1)

    def updateWidgets(self):
        """ Slot called after the reading of the CITS.
        Sets the values combo box (voltage and channels) and draws the map """
        self.updateVoltageBox(
            self.m_params["vStart"], self.m_params["vEnd"], self.m_params["zPt"]
        )
        self.m_channelBox.clear()
        self.m_channelBox.addItems(self.channelList)
        self.drawXYMap(self.m_voltageBox.value())
        self.updateAboveValue(self.m_aboveBar.value())
        self.updateBelowValue(self.m_belowBar.value())
        self.m_CitsAvgBox.setMaximum(self.m_params["xPx"])
        logging.debug("Widgets updated !")

    # %% Methods related to spectra
    def normalizeSpectrum(self):
        """ Normalizes the spectra of displayed index """
        chan = self.m_channelBox.currentIndex()
        # TODO: can be vectorised
        for y in range(self.m_params["xPx"]):
            for x in range(self.m_params["yPx"]):
                self.m_data[chan, y, x] = normalizeDOS(
                    self.m_data[chan, y, x], self.m_params["zPt"]
                )

    def averageSpectrum(self, xi, xf, yi, yf):
        """ Averages the spectra contained in the rectangle drawn
        by the points (xi,yi);(xf,yi);(xf,yf) and (xi,yf) """
        if self.dataLoaded:
            chan = self.m_channelBox.currentIndex()
            if yf < yi:
                yi, yf = yf, yi
            if xf < xi:
                xi, xf = xf, xi
            avg_data = np.mean(self.m_data[chan, yi:yf, xi:xf], axis=(0, 1))
            xavg = xi + (xf - xi) // 2
            yavg = yi + (yf - yi) // 2
            self.drawSpectrum(
                avg_data, "Average around ({0},{1})".format(xavg, yavg), xavg, yavg,
            )
            return avg_data

    def averageSpectrumWithValues(self):
        """ Averages spectra according their values at a certain voltage """
        if self.dataLoaded:
            voltage = self.m_voltageBox.value()
            chan = self.m_channelBox.currentIndex()
            viewSelected = self.m_viewSelectedBox.isChecked()
            # midpoint=(self.mapMax+self.mapMin)/2
            midpoint2 = self.mapMax - self.mapMin
            limit_aboveV = self.mapMax - self.m_aboveBar.value() * midpoint2 / 100
            limit_belowV = self.mapMin + self.m_belowBar.value() * midpoint2 / 100
            if limit_aboveV < limit_belowV:
                logging.error("Above and below spectra intersect !")
                return
            isAbove = self.m_data[chan, :, :, voltage] > limit_aboveV
            isBelow = self.m_data[chan, :, :, voltage] < limit_belowV
            aboveValues = self.m_data[chan, :, :, :][isAbove]
            belowValues = self.m_data[chan, :, :, :][isBelow]

            N_aboveV = len(aboveValues)
            if N_aboveV != 0:
                self.drawSpectrum(
                    np.mean(aboveValues, axis=0),
                    "Average above {0} ({1} spectra averaged)".format(
                        limit_aboveV, N_aboveV
                    ),
                )
            N_belowV = len(belowValues)
            if N_belowV != 0:
                self.drawSpectrum(
                    np.mean(belowValues, axis=0),
                    "Average below {0} ({1} spectra averaged)".format(
                        limit_belowV, N_belowV
                    ),
                )
            if viewSelected:
                # TODO: to reimplement with np.indices or np.where
                # self.ax_map.plot(xPts, yPts, c=cPts, marker="o", ls="None")
                pass

    def clearSpectrum(self):
        """ Clears the spectrum window """
        self.ax_spec.clear()
        self.nSpectraDrawn = 0
        self.voltageLine = 0
        self.clearShapesClicked()
        # self.drawTopo()
        self.m_specWidget.draw()

    def drawSpectrum(self, dataToPlot, label, x=-1, y=-1):
        """ Method called each time a spectrum needs to be plotted.
        Takes care of the derivative and different scales stuff and updates the window.
        """
        shiftX = str(self.m_shiftXBox.text()).lower()
        if shiftX == "topo":
            shiftX = self.topo[y, x]
        else:
            shiftX = float(shiftX)
        shiftY = self.nSpectraDrawn * float(self.m_shiftYBox.text())
        if self.dataLoaded and dataToPlot.size != 0:
            if self.m_logBox.isChecked():
                dataToPlot = np.log(dataToPlot)

            finalLabel = "{0} - {1}".format(label, self.m_channelBox.currentText())
            # TODO: remove lastSpectrum with a ref to lines plotted in the spectrum widget
            self.lastSpectrum = [dataToPlot, finalLabel]
            if self.m_scaleVoltage.isChecked():
                dV = self.m_params["dV"]
                vStart = self.m_params["vStart"]
                vEnd = self.m_params["vEnd"]
                zPt = self.m_params["zPt"]
                voltages = np.arange(vStart, vEnd, dV) + shiftX
                # Check consistency of voltages array with the number of points.
                if len(voltages) != zPt:
                    logging.warning(
                        "Round-off error while computing the voltage array: dV ({}) might be too small."
                        "Computing from number of points instead.".format(dV)
                    )
                    # If this fails, compute from the number of points to be sure to have the same dim as dataToPlot.
                    voltages = np.linspace(vStart, vEnd, zPt) + shiftX
                self.ax_spec.plot(
                    voltages,
                    shiftY + dataToPlot,
                    label=finalLabel,
                    c=self.getSpectrumColor(self.nSpectraDrawn),
                )
                if self.m_derivBox.isChecked():
                    deriv = sp.signal.savgol_filter(
                        dataToPlot, self.m_derivNBox.value(), 3, deriv=1, delta=dV
                    )
                    self.ax_spec.plot(
                        voltages,
                        shiftY + deriv,
                        label="Derivative of " + finalLabel,
                        marker="o",
                        markersize=3.0,
                        c=self.getSpectrumColor(self.nSpectraDrawn),
                    )
                self.nSpectraDrawn = self.nSpectraDrawn + 1
            else:
                self.ax_spec.plot(
                    shiftY + dataToPlot,
                    label=finalLabel,
                    c=self.getSpectrumColor(self.nSpectraDrawn),
                )
                if self.m_derivBox.isChecked():
                    deriv = sp.signal.savgol_filter(
                        dataToPlot, self.m_derivNBox.value(), 3, deriv=1, delta=dV
                    )
                    self.ax_spec.plot(
                        shiftY + deriv,
                        label="Derivative of " + finalLabel,
                        marker="o",
                        markersize=3.0,
                        c=self.getSpectrumColor(self.nSpectraDrawn),
                    )
                self.nSpectraDrawn = self.nSpectraDrawn + 1
        if not self.m_legendBox.isChecked():
            self.ax_spec.legend(loc=0)
        self.m_specWidget.draw()

    def fitSpectrum(self):
        """ Linear fitting of the last spectrum plotted
        that is stored in lastSpectrum """
        if self.dataLoaded and self.lastSpectrum[0].size != 0:

            dV = self.m_params["dV"]
            if self.m_fitCustomCheckbox.isChecked():
                limitL = float(self.m_fitLowerBox.text())
                limitU = float(self.m_fitUpperBox.text())
                vStart = min(limitL, limitU)
                vEnd = max(limitL, limitU)
                zPt = self.m_params["zPt"]
                voltages = np.arange(zPt) * dV + self.m_params["vStart"]
                # Take the portion to fit
                mask1 = voltages > vStart
                mask2 = voltages < vEnd
                dataToFit = self.lastSpectrum[0][np.logical_and(mask1, mask2)]
            else:
                vStart = self.m_params["vStart"]
                vEnd = self.m_params["vEnd"]
                dataToFit = self.lastSpectrum[0]
            X = np.arange(vStart, vEnd, dV)
            popt, pcov = sp.optimize.curve_fit(linear_fit_function, X, dataToFit)
            slope = popt[0]
            coef = popt[1]
            logging.info(
                "Linear fit gives a slope of {} and a coef of {}".format(slope, coef)
            )
            self.ax_spec.plot(
                X,
                slope * X + coef,
                label="Linear fit of " + self.lastSpectrum[1],
                color=self.getSpectrumColor(self.nSpectraDrawn - 1),
                linestyle="--",
            )
            self.ax_spec.legend(loc=0)
            self.m_specWidget.draw()

    def getSpectrumColor(self, n):
        """ Returns the color corresponding to a spectrum of given index
        according to spectrumColor """
        # TODO: this should be removed in favour of matplotlib color_cycle or an iterator
        i = n % len(self.spectrumColor)
        return self.spectrumColor[i]

    def launchAvgSpectrum(self):
        """ Slot called to average the spectra of the whole map """
        if self.dataLoaded:
            xPx = self.m_params["xPx"]
            yPx = self.m_params["yPx"]
            self.averageSpectrum(0, xPx, 0, yPx)

    def pickSpectrum(self, event):
        """ Method called when a press-and-release event is done
        at the same location of the map. If the average box is checked,
        it will keep the data in storage to average it later.
        Otherwise it plots the spectrum in the corresponding widget """
        if event.xdata is not None and event.ydata is not None and self.dataLoaded:
            PixelX = int(event.xdata)
            PixelY = int(event.ydata)
            chan = self.m_channelBox.currentIndex()
            # TODO: the slope part should be put in a function
            if "Slope" in self.channelList[chan]:
                dataLogCurrent = np.where(
                    self.m_data[0, PixelY, PixelX] < 0.01,
                    np.log(self.m_data[0][PixelY][PixelX]),
                    self.m_data[0, PixelY, PixelX],
                )
                self.drawSpectrum(
                    dataLogCurrent,
                    "Log of current at [{}, {}]".format(PixelX, PixelY),
                    PixelX,
                    PixelY,
                )
                # self.drawSpectrum(dataLine,"Linear fit of the log at "+str(PixelX)+","+str(PixelY)+"]")
            else:
                self.m_mapWidget.draw()
                # Add data to the total data to average if in average mod
                if self.m_avgBox.isChecked():
                    self.tot_data += self.m_data[chan, PixelY, PixelX]
                    self.nAvgSpectra += 1
                # Plot the data with the desired scale (Volts or index) if in normal mode
                else:
                    self.drawSpectrum(
                        self.m_data[chan, PixelY, PixelX],
                        "[{}, {}]".format(PixelX, PixelY),
                        PixelX,
                        PixelY,
                    )

    def saveSpectra(self):
        """
        Saves the plotted spectra *as-is* as CSV data
        """
        output_path, fmt = QtWidgets.QFileDialog.getSaveFileName(
            self, "Choose a path to save the CSV file", self.wdir, "CSV (*.csv)"
        )
        header = ""
        output = []
        for line in self.ax_spec.lines:
            # If first line, add the X-axis data
            if len(output) == 0:
                header += "X-Axis,"
                output.append(line._x)
            output.append(line._y)
            # Replace ',' in the labels to avoid problems when exporting as CSV
            header += line._label.replace(",", "/") + ","
        # Crop last comma after looping on lines
        header = header[:-1]
        # Save if there is an output
        if len(output) != 0:
            # Tranposing the output array to get columns format corresponding to the header
            np.savetxt(
                output_path,
                np.transpose(np.array(output)),
                delimiter=",",
                header=header,
            )
            logging.info("Finished exporting csv at {}".format(output_path))

    # %% Methods related to the clicks on the map
    def onPressOnMap(self, event):
        """ Slot called when a press event is detected. Creates a Shape object
        that will be dynamically updated when the mouse moves """
        if (event.xdata is not None and event.ydata is not None) and (
            self.dataLoaded and self.toolbar_map._active is None
        ):
            # TODO: the ratio could be put in attributes
            if self.m_scaleMetric.isChecked():
                ratioX = self.m_params["xL"] / self.m_params["xPx"]
                ratioY = self.m_params["yL"] / self.m_params["yPx"]
            else:
                ratioX = 1
                ratioY = 1
            self.currentShape = generateShape(
                event,
                self.m_mapWidget.figure,
                self.fig_topo,
                self.getSpectrumColor(self.nSpectraDrawn),
                ratioX,
                ratioY,
            )
            self.motionConnection = self.m_mapWidget.mpl_connect(
                "motion_notify_event", self.currentShape.update
            )
            logging.debug("PRESS")

    def onReleaseOnMap(self, event):
        """ Slot called when a release event is detected.
        Disconnects the updating of the currentShape and
        launches the appropriate method depending on which button
        was pressed and where it was released """
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
                    # Cut along the XY line if a line is traced (X or Y different)
                    if xf != xi or yf != yi:
                        self.cutAlongLine(xi, xf, yi, yf)
                    # Pick spectrum otherwise and change the line shape to a point
                    else:
                        self.currentShape = changeToDot(self.currentShape)
                        self.pickSpectrum(event)
                # If right-click : either a rectangle was drawn or the center of the rectangle to average was picked
                else:
                    if xf != xi or yf != yi:
                        self.averageSpectrum(xi, xf, yi, yf)
                    # If X=Y, we need to force the updating of the Shape so it is drawn around the X,Y point
                    # and not starting at X,Y (TODO: This can be refactored)
                    else:
                        n = self.m_rcAvgBox.value()
                        self.currentShape.forceUpdate(
                            max(0, xi - n),
                            max(0, yi - n),
                            min(self.m_params["xPx"], xf + n),
                            min(self.m_params["yPx"], yf + n),
                        )
                        self.averageSpectrum(
                            max(0, xi - n),
                            min(self.m_params["xPx"], xf + n),
                            max(0, yi - n),
                            min(self.m_params["yPx"], yf + n),
                        )
                # Add the current Shape to the list of clicked Shapes
                self.addToShapesClicked(self.currentShape)
            logging.debug("RELEASE")

    def addToShapesClicked(self, shape):
        """ Method called when a release was detected on the map.
        The Shape is saved in the shapes_clicked list """
        self.shapes_clicked.append(shape)

    def drawShapesClicked(self, fig_map, fig_topo, recreate):
        """ Method that draws all the Shapes saved in shapes_clicked
        or recreates them if the XYmap was cleared"""
        for shape in self.shapes_clicked:
            if recreate:
                shape.recreate(fig_map, fig_topo)
            else:
                shape.draw()

    def clearShapesClicked(self):
        """ Method that clears all saved Shapes and
        then redraws the mapto reflect the change.
        No need to call Shape.remove as the XYmap will be cleared """
        self.shapes_clicked = []
        self.drawXYMap(self.m_voltageBox.value())

    def cutAlongLine(self, xi, xf, yi, yf):
        """ Method that finds the positions of the pixels forming
        the line from (xi,yi) to (xf,yf).
        The plotting of the spectra is done in cutPlot called at the end """
        # If the line is vertical, the equation is x=xi with y varying from yi to yf
        logging.info(
            "Cut from ({}, {}) to ({}, {}) in {}".format(
                xi, yi, xf, yf, self.m_channelBox.currentText()
            )
        )
        x_plot, y_plot = findPixelsOnLine(xi, xf, yi, yf)
        return self.cutPlot(x_plot, y_plot)

    def cutPlot(self, x_plot, y_plot):
        """ Method called by cutAlongLine.
        Plots the spectra for the pixels of positions (x_plot[i],y_plot[i]) """
        # Build the data to plot with v as Y and z (number of pixels gone through) as X
        zPt = self.m_params["zPt"]
        voltages = np.arange(zPt)
        chan = self.m_channelBox.currentIndex()
        z_plot = np.arange(x_plot.size)
        dataToPlot = np.ndarray(shape=(zPt, z_plot.size))
        xi = x_plot[0]
        yi = y_plot[0]
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
        ax = fig.add_subplot(1, 1, 1)
        ax.set_title(self.mapName.split(".")[0] + " - Cut " + str(fig.number))
        self.ax_map.text(xi + 0.5, yi + 0.5, str(fig.number))
        # TODO: The waterfall and the 2D should be put in seperate functions.
        if self.m_waterfallButton.isChecked():
            if self.m_scaleVoltage.isChecked():
                voltages = self.m_params["vStart"] + voltages * self.m_params["dV"]
            for z in z_plot:
                xc = int(x_plot[z])
                yc = int(y_plot[z])
                spectrum = self.m_data[chan, yc, xc]
                offset = (zf - z) * float(self.m_shiftYBox.text())
                ax.plot(voltages, spectrum + offset, "k", zorder=(z + 1) * 2)
                # White filling under the curves
                ax.fill_between(
                    voltages,
                    spectrum + offset,
                    offset,
                    facecolor="w",
                    lw=0,
                    zorder=(z + 1) * 2 - 1,
                )
            ax.set_xlim([voltages[0], voltages[-1]])
        else:
            for z in z_plot:
                xc = int(x_plot[z])
                yc = int(y_plot[z])
                dataToPlot[:, z] = self.m_data[chan, yc, xc, :]
                if viewSelected:
                    self.addToPtsClicked(xc, yc, color="yellow")
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
                mapData = pyplot.pcolormesh(
                    metricDistances,
                    voltages,
                    dataToPlot,
                    cmap=self.m_colorBarBox.currentText(),
                )
                ax.axis(
                    [metricDistances[0], metricDistances[-1], voltages[0], voltages[-1]]
                )
                ax.set_xlabel("Distance (nm)")
            else:
                mapData = pyplot.pcolormesh(
                    z_plot, voltages, dataToPlot, cmap=self.m_colorBarBox.currentText()
                )
                ax.axis([z_plot[0], z_plot[-1], voltages[0], voltages[-1]])
                ax.set_xlabel("Pixels")
            # Colorbar set up
            cbar = fig.colorbar(mapData, shrink=0.9, pad=0.05, aspect=15)
            cbar.ax.yaxis.set_ticks_position("both")
            cbar.ax.tick_params(axis="y", direction="in")
            if self.m_cbarCustomCheckbox.isChecked():
                mapData.set_clim(
                    float(self.m_cbarLowerBox.text()), float(self.m_cbarUpperBox.text())
                )
            else:
                mapData.set_clim(0)
        fig.show()
        return metricDistances, voltages, dataToPlot

    def launchBigCut(self):
        """ Launches a cut in the big diagonal """
        if self.dataLoaded:
            if self.m_scaleMetric.isChecked():
                ratioX = self.m_params["xL"] / self.m_params["xPx"]
                ratioY = self.m_params["yL"] / self.m_params["yPx"]
            else:
                ratioX = 1
                ratioY = 1
            # Simulates a left-click press at (0,0) and a release at xPx-1,yPx-1
            simEvent = matplotlib.backend_bases.MouseEvent(
                "button_press_event", self.m_mapWidget.figure.canvas, 0, 0, button=1
            )
            simEvent.xdata = 0
            simEvent.ydata = 0
            self.currentShape = generateShape(
                simEvent,
                self.m_mapWidget.figure,
                self.fig_topo,
                self.getSpectrumColor(self.nSpectraDrawn),
                ratioX,
                ratioY,
            )
            simEvent.xdata = self.m_params["xPx"] - 1
            simEvent.ydata = self.m_params["yPx"] - 1
            self.currentShape.update(simEvent)
            X, Y, Z = self.cutAlongLine(
                0, self.m_params["xPx"] - 1, 0, self.m_params["yPx"] - 1
            )
            self.m_mapWidget.draw()
            return X, Y, Z

    # %% Methods related to the map
    def clearMap(self):
        """ Unloads the map and clears the map window """
        self.m_mapWidget.figure.clear()
        self.clearShapesClicked()
        self.m_data = np.ndarray([])
        # self.m_params.clear()
        self.mapType = ""
        self.dataLoaded = False

    def drawXYMap(self, voltage):
        """ Calls the getMapData function and draws the result
        in the map window with the approriate formatting.
        Called at each change of voltage """
        # Start everything anew
        fig_map = self.m_mapWidget.figure
        fig_map.clear()
        if self.dataLoaded:
            xPx = self.m_params["xPx"]
            yPx = self.m_params["yPx"]
            if xPx == yPx:
                self.ax_map = fig_map.add_subplot(1, 1, 1, aspect="equal")
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
                XYmap = self.ax_map.pcolormesh(
                    x_m, y_m, mapData, cmap=self.m_colorBarBox.currentText()
                )
                self.ax_map.axis([0, xL, 0, yL])
            # Else, use pixels
            else:
                XYmap = self.ax_map.pcolormesh(
                    mapData, cmap=self.m_colorBarBox.currentText()
                )
                # If the map is an Omicron one, I have to invert the y-axis
                if self.mapType == "Omicron":
                    self.ax_map.axis([0, xPx, yPx, 0])
                else:
                    self.ax_map.axis([0, xPx, 0, yPx])
            # Set title
            chan = self.m_channelBox.currentText()
            self.ax_map.set_title(
                self.mapName
                + " - "
                + chan
                + "\n"
                + "V="
                + str(self.m_params["vStart"] + voltage * self.m_params["dV"])
            )
            # Colorbar stuff
            cbar = fig_map.colorbar(XYmap, shrink=0.9, pad=0.05, aspect=15)
            cbar.ax.yaxis.set_ticks_position("both")
            cbar.ax.tick_params(axis="y", direction="in")
            # Image color scale is adjusted to the data:
            if self.m_cbarCustomCheckbox.isChecked():
                XYmap.set_clim(
                    float(self.m_cbarLowerBox.text()), float(self.m_cbarUpperBox.text())
                )
            else:
                XYmap.set_clim(self.mapMin, self.mapMax)
            # Plot a dashed line at X=voltage if asked
            if self.m_vLineBox.isChecked():
                self.drawVoltageLine(voltage)
            # Recreate the shapes points saved in shapes_clicked
            self.drawShapesClicked(fig_map, self.fig_topo, True)
            self.m_mapWidget.draw()

    def getMapData(self, v):
        """ Returns an array built from the data loaded that can be used to
        display a map at fixed voltage """
        mapData = self.m_data[self.m_channelBox.currentIndex(), :, :, v]
        return mapData, np.min(mapData), np.max(mapData)

    # %% Methods related to the voltage guide line in the spectra window

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
            self.voltageLine = self.ax_spec.axvline(currentV, color="k", linestyle="--")
            self.m_specWidget.draw()

    def clearVoltageLine(self):
        """ Removes the vertical voltage line """
        if self.voltageLine:
            self.ax_spec.lines.pop(self.ax_spec.lines.index(self.voltageLine))
            self.m_specWidget.draw()
        self.voltageLine = 0

    # Post-processing methods
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
            logging.error(
                "The data is too big ! Or the memory too small... Aborting addChannel !"
            )
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
        """ Slot that averages the spectras in the X direction of the loaded
        CITS and replaces the loaded CITS by the result """
        if self.dataLoaded:  # Check if a CITS was loaded
            # Get the needed params
            Navg = int(self.m_CitsAvgBox.value())
            xPx = self.m_params["xPx"]
            yPx = self.m_params["yPx"]
            zPt = self.m_params["zPt"]
            # Try the allocation
            try:
                new_data = np.zeros(
                    shape=(len(self.channelList), yPx, xPx // Navg, zPt)
                )
            except MemoryError:
                logging.error("Not enough memory to average spectra !")
                return
            # Average in X for each channel and Y
            for i_chan in range(len(self.channelList)):
                for y in range(yPx):
                    # To average, for each x, add every spectrum between x and x+Navg divided by Navg.
                    for x in range(0, xPx, Navg):
                        # If the averaging window is not contained in the CITS data, stop averaging.
                        # The last spectra of this window will then will dismissed.
                        if x + Navg > xPx:
                            break
                        # Else, average by adding the spectra of the avergaging window and dividing by Navg
                        else:
                            spectra = np.zeros(zPt)
                            for i in range(int(Navg)):
                                spectra = spectra + self.m_data[i_chan, y, x + i] / Navg
                                # Store the result in new_data
                                new_data[i_chan, y, x // Navg] = spectra
            # If eveything went well, clear the current map and replace it by
            # the created data and change xPx
            self.clearMap()
            self.mapName = "Average_" + str(Navg)
            self.m_params["xPx"] = xPx / Navg
            self.dataLoaded = True
            self.m_data = new_data
            self.updateWidgets()

    def extractDerivative(self, numChanToDeriv):
        """ Extracts the derivative of a channel and adds it to the data """
        if self.dataLoaded:
            derivData = sp.signal.savgol_filter(
                self.m_data[numChanToDeriv],
                9,
                2,
                deriv=1,
                delta=self.m_params["dV"],
                axis=-1,
            )
            # Add the channel to the data
            self.addChannel(
                derivData, "Derivative of " + self.channelList[numChanToDeriv]
            )

    def extractFFT(self, numChanToFFT, axisOfFFT):
        """ Extracts the FFT of a channel and adds it to the data """
        if self.dataLoaded:
            FFTData = np.fft.fft(self.m_data[numChanToFFT], axis=axisOfFFT)
            # Add the channel to the data
            self.addChannel(FFTData, "FFT of " + self.channelList[numChanToFFT])

    def normalizeCurrentChannel(self):
        if self.dataLoaded:
            yPx = self.m_params["yPx"]
            xPx = self.m_params["xPx"]
            zPt = self.m_params["zPt"]
            chan = self.m_channelBox.currentIndex()
            normData = np.zeros(shape=(yPx, xPx, zPt))
            for y in range(yPx):
                for x in range(xPx):
                    normData[y][x] = self.normalizeDOS(self.m_data[chan][y][x], 10)
            self.addChannel(normData, "Normalized " + self.channelList[chan])
