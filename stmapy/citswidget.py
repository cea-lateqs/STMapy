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
from .shape import generateShape, changeToDot, Dot, Line
from .reads import readCitsAscii, readCits3dsBin, readCitsSm4Bin, setUpConfig
from .processing import (
    normalizeDOS,
    linearFitFunction,
    findPixelsOnLine,
    directionAverageCITS,
    FeenstraNormalization,
    FeenstraNormalization_log,
    derivate_IV,
)
from .plotting import (
    visualize_statistics, 
    visualize_statistics_topo, 
    waterfallPlot, 
    lineTopoPlot, 
    imageTopoPlot
)
from PIL import Image, ImageFile
from scipy.signal import find_peaks

#Set style
matplotlib.pyplot.style.use(os.path.join(os.path.dirname(__file__),'stmapyStyle.mplstyle')) 

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
        # TODO: CITS should be a dedicated object
        self.cits_data = np.ndarray([])
        self.cits_add_data = np.ndarray([])
        self.cits_params = {}
        self.channelList = []
        self.cits_name = ""
        self.summed_spectra = np.ndarray([])
        self.nb_summed_spectra = 0
        # Set up map figure
        self.xy_map = None
        self.toolbar_map = NavigationToolbar(self.ui_mapWidget, self)
        # Set up spectrum figure
        self.ax_spec = self.ui_specWidget.figure.add_subplot(111)
        self.ui_specWidget.figure.subplots_adjust()
        self.voltageLine = None
        self.toolbar_spec = NavigationToolbar(self.ui_specWidget, self)
        # Add csv save button
        self.ui_saveCsvButton = QtWidgets.QPushButton("Save as CSV")
        self.toolbar_spec.addWidget(self.ui_saveCsvButton)
        self.map_layout.insertWidget(0, self.toolbar_map)
        self.spec_layout.insertWidget(0, self.toolbar_spec)
        # Variables linked to clicks on map
        self.shapes_clicked = []
        # Colormaps
        self.ui_colorBarBox.addItems(matplotlib.pyplot.colormaps())
        # Boolean that is True if a map is loaded
        self.dataLoaded = False
        # Boolean that is True if a map is loaded
        self.dataAdded = False
        # Parse config
        self.config = setUpConfig(
            os.path.join(os.path.dirname(__file__), "config.json")
        )
        self.wdir = self.config["working_directory"]
        self.ui_colorBarBox.setCurrentIndex(
            self.ui_colorBarBox.findText(self.config["default_cmap"])
        )
        # Other parameters used after map loading
        self.mapType = ""
        self.fig_topo = None
        self.topo = []
        # Connect all widgets
        self.connect()
        self.nSpectraDrawn = 0
        self.lastSpectrum = []
        # Set up layouts
        self.ui_avgWidget.hide()
        self.ui_cbarWidget.hide()
        self.ui_fitUpperLabel.hide()
        self.ui_fitLowerLabel.hide()
        self.ui_fitUpperBox.hide()
        self.ui_fitLowerBox.hide()
        # Calls the loading method at launch
        if self.config["autoload"]:
            self.askCits()

    def connect(self):
        """ Connects all the signals. Only called in the constructor """
        self.ui_openButton.clicked.connect(self.askCits)
        self.ui_topoButton.clicked.connect(self.drawTopo)
        self.ui_addCITSButton.clicked.connect(self.addCits)
        self.ui_FFTButton.clicked.connect(self.calculateFFT)
        self.ui_derivateButton.clicked.connect(self.derivateCurrentChannel)
        self.ui_FeenstraBox.clicked.connect(self.FeenstraCits)
        self.ui_channelBox.currentIndexChanged.connect(self.updateMap)
        self.ui_colorBarBox.currentTextChanged.connect(self.setColorMap)
        self.ui_indexBox.valueChanged.connect(self.drawXYMap)
        self.ui_derivBox.stateChanged.connect(self.ui_derivNBox.setEnabled)
        self.ui_mapWidget.mpl_connect("button_press_event", self.onPressOnMap)
        self.ui_mapWidget.mpl_connect("button_release_event", self.onReleaseOnMap)
        self.ui_specWidget.mpl_connect("button_press_event", self.updateToPointedX)
        self.ui_saveCsvButton.clicked.connect(self.saveSpectra)
        self.ui_clearSpec.clicked.connect(self.clearSpectrum)
        self.ui_scaleVoltage.toggled.connect(self.clearSpectrum)
        self.ui_averageCitsButton.clicked.connect(self.averageCITSOnAxis)
        self.ui_wholeLengthCutButton.clicked.connect(self.launchBigCut)
        self.ui_normButton.clicked.connect(self.normalizeCurrentChannel)
        self.ui_avgSpec.clicked.connect(self.launchAvgSpectrum)
        self.ui_avgBox.toggled.connect(self.updateAvgVariables)
        self.ui_vLineBox.toggled.connect(self.toggleVoltageLine)
        self.ui_avgVButton.clicked.connect(self.averageSpectrumWithValues)
        # Averaging with respect to values
        self.ui_avgCheckBox.toggled.connect(self.ui_avgWidget.setVisible)
        self.ui_aboveBox.valueChanged.connect(self.updateAboveValue)
        self.ui_belowBox.valueChanged.connect(self.updateBelowValue)
        # Cbar custom limits
        self.ui_cbarCheckBox.toggled.connect(self.ui_cbarWidget.setVisible)
        self.ui_cbarCustomCheckbox.toggled.connect(self.setColorMapLims)
        self.ui_cbarCustomCheckbox.toggled.connect(self.ui_cbarUpperBox.setEnabled)
        self.ui_cbarCustomCheckbox.toggled.connect(self.ui_cbarLowerBox.setEnabled)
        self.ui_cbarUpperBox.textChanged.connect(self.setColorMapMax)
        self.ui_cbarLowerBox.textChanged.connect(self.setColorMapMin)
        # Fit spectrum widget
        self.ui_fitSpec.clicked.connect(self.fitSpectrum)
        self.ui_fitCustomCheckbox.toggled.connect(self.ui_fitUpperBox.setVisible)
        self.ui_fitCustomCheckbox.toggled.connect(self.ui_fitLowerBox.setVisible)
        self.ui_fitCustomCheckbox.toggled.connect(self.ui_fitUpperLabel.setVisible)
        self.ui_fitCustomCheckbox.toggled.connect(self.ui_fitLowerLabel.setVisible)
        self.ui_MakeGifButton.clicked.connect(self.make_gif)


    @property
    def metric_ratios(self):
        ratioX = 1
        ratioY = 1
        if self.ui_scaleMetric.isChecked():
            try:
                ratioX = self.cits_params["xL"] / self.cits_params["xPx"]
                ratioY = self.cits_params["yL"] / self.cits_params["yPx"]
            except KeyError:
                logging.warning(
                    "Metric ratios could not be computed. Is the CITS loaded with all parameters ?"
                )
        return {"x": ratioX, "y": ratioY}    
        #TODO : the plotting on topo is incorrect when npx(CITS) neq npx(Topo)

    # %% Reading and loading CITS methods
    def askCits(self):
        """ Slot that only takes care of opening a file dialog
        for the user to select one or several CITS - Returns the CITS paths """
        cits_names_and_ext = QtWidgets.QFileDialog.getOpenFileNames(
            self,
            "Choose a CITS file to read or several to average",
            self.wdir,
            "Ascii file (*.asc);;RHK file (*.sm4);;3D binary file (*.3ds);;Text file (*.txt);;All files (*)",
        )
        # getOpenFilesNames returns a tuple with Cits_names as first element
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
                (
                    self.topo,
                    self.cits_data,
                    self.channelList,
                    self.cits_params,
                ) = readCitsAscii(cits)
                self.dataLoaded = True
            elif extension == "3ds":
                self.clearMap()
                self.mapType = "Nanonis"
                (
                    self.topo,
                    self.cits_data,
                    self.channelList,
                    self.cits_params,
                    zSpectroData,
                ) = readCits3dsBin(cits)
                if zSpectroData is not None:
                    self.addChannel(
                        zSpectroData["slopeData"],
                        "Slope by linear fit of {}".format(zSpectroData["zChannel"]),
                    )
                    self.addChannel(
                        zSpectroData["coefData"],
                        "Coef by linear fit of {}".format(zSpectroData["zChannel"]),
                    )
                    self.addChannel(zSpectroData["Zg"], "Zg")
                self.dataLoaded = True
            elif extension == "txt":
                self.readTopo(cits)
            elif extension == "sm4":
                self.clearMap()
                self.mapType = "Sm4"
                (
                    self.topo,
                    self.cits_data,
                    self.channelList,
                    self.cits_params,
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
                self.cits_name = os.path.basename(cits)
                logging.info(self.cits_name + " read as a " + self.mapType + " map")
            else:
                logging.error("Problem while reading " + cits)
                return
            # If only one Cits was selected, there is no need to run the averaging code
            # so return after drawing the topo and updating of the widgets
            if n_cits == 1:
                self.updateWidgets()
                self.drawTopo()
                return
            else:  # Else begin the averaging
                # If this was the first Cits to average, create the mean_data array to store the average
                if first:
                    mean_data = self.cits_data / n_cits
                    first = False
                else:  # Else, continue adding to the mean_data
                    mean_data += self.cits_data / n_cits
        # If everthing went well and if there was several CITS chosen,
        # clear the map and set the data to mean_data.
        self.cits_name = "Average of {0} CITS".format(n_cits)
        self.clearMap()
        self.dataLoaded = True
        self.cits_data = mean_data
        self.updateWidgets()
        return
    
    # %% Reading and loading topo images methods

    def drawTopo(self):
        """ Draws the topography read while opening the CITS."""
        assert self.topo.ndim == 2
        yPx, xPx = self.topo.shape
        level_algo = self.config["level_topo"]

        # Set up the figure for the plot
        if self.fig_topo is None:
            self.fig_topo = pyplot.figure()
            self.fig_topo.canvas.mpl_connect("close_event", self.handleClosingTopo)
        else:
            self.fig_topo.clear()

        if self.ui_scaleMetric.isChecked():
            logging.debug("Scale metric box is checked")
            max_x = self.cits_params["xL"]
            max_y = self.cits_params["yL"]
            ax_label = "Distance (nm)"
        else:
            max_y, max_x = yPx, xPx
            logging.debug("Scale metric box is checked")
            ax_label = "Pixels"

        if yPx == 1:
            self.ax_topo = lineTopoPlot(self.fig_topo, max_x, self.topo)
            self.ax_topo.set_ylabel("Z (nm)")
        else:
            self.ax_topo = imageTopoPlot(
                self.fig_topo,
                max_x,
                max_y,
                self.topo,
                self.config["topo_cmap"],
                level_algo,
            )
            if self.mapType == "Omicron":
                self.ax_topo.axis([0, max_x, max_y, 0])
            else:
                self.ax_topo.axis([0, max_x, 0, max_y])
            self.ax_topo.set_ylabel(ax_label)
        self.ax_topo.set_xlabel(ax_label)

        self.fig_topo.subplots_adjust()
        self.fig_topo.canvas.draw()
        self.fig_topo.show()

    def handleClosingTopo(self, event):
        """ Called when the topo figure is closed.
        Put back self.fig_topo to 0 to indicate that no topo figure exists.
        """
        logging.debug("Topo closing")
        self.fig_topo = None

#%% Adding channels
    def addCits(self):
        """ Slot that launches the reading of the I(V) CITS 
        that corresponds to the preniously opened single dIdV CITS.
        TODO : Only implemented for ascii files"""  
        if self.dataLoaded :
            IVpath = QtWidgets.QFileDialog.getOpenFileNames(
                self,
                "Choose the I(V) CITS file correponding to the alreday loaded dIdV CITS file",
                self.wdir,
                "Ascii file (*.asc)",
            )
            # getOpenFilesNames returns a tuple with Cits_names as first element
            # and extension as second. We just need the names.
            (
                topo,
                cits_add_data,
                IV_channelList,
                cits_params,
            ) = readCitsAscii(IVpath[0][0])
            self.cits_data = np.concatenate((self.cits_data, cits_add_data), axis = 0)
            self.dataAdded = True
            IV_channelList = ['IV_'+x for x in IV_channelList]
            self.channelList = self.channelList + IV_channelList
            self.ui_channelBox.addItems(IV_channelList)
        else:
            logging.error(
                "Please first load single dIdV CITS"
                )
        return

    def normalizeCurrentChannel(self):
        """ Adds the normalized current data as a new channel """
        if self.dataLoaded:
            chan = self.ui_channelBox.currentIndex()
            self.addChannel(
                normalizeDOS(self.cits_data[chan]),
                "Normalized " + self.channelList[chan],
            )
            self.ui_channelBox.setCurrentIndex(len(self.channelList) - 1)

    def derivateCurrentChannel(self):
        """ Adds the derivative of current data as a new channel """
        if self.dataLoaded:
            chan = self.ui_channelBox.currentIndex()
            self.addChannel(
                derivate_IV(self.cits_data[chan], self.cits_params["dV"]),
                "Derivated " + self.channelList[chan],
            )
            self.ui_channelBox.setCurrentIndex(len(self.channelList) - 1)

    def calculateFFT(self):
        """ Adds the 2D FFT data as a new channel """
        if self.dataLoaded:
            chan = self.ui_channelBox.currentIndex()
            self.addChannel(
                np.sqrt((np.abs(np.fft.fftshift(
                    np.fft.fft2(self.cits_data[chan], axes=(- 3, - 2)), axes=(- 3, - 2)
                        )
                 )
                )),
                "FFT of " + self.channelList[chan],
            )
            self.ui_channelBox.setCurrentIndex(len(self.channelList) - 1)

    def FeenstraCits(self):
          """ Slot that launches the calculation of Feenstra CITS 
          corresponding to the previously opened single dIdV and I(V) CITS.
          TODO : Only implemented for ascii files"""  
          if self.dataLoaded and self.dataAdded :
              self.cits_params["LIsensi"] = self.ui_LIBox.value()
              self.cits_params["AmpMod"] = self.ui_AmpBox.value()
              self.cits_params["Gain"] = self.ui_GainBox.value()
              cits_feenstra = FeenstraNormalization(self.cits_data, self.cits_params, self.channelList)
              # cits_feenstra = FeenstraNormalization_log(self.cits_data, self.cits_params, self.channelList)
              Feenstra_channelList = ['Feenstra_Data [Fwd]', 'Feenstra_Data [Bwd]']
              self.channelList = self.channelList + Feenstra_channelList
              self.ui_channelBox.addItems(Feenstra_channelList)
              self.cits_data = np.concatenate((self.cits_data, cits_feenstra), axis = 0)
          else:
              logging.error(
                  "Please first load single dIdV CITS and its corresponding I(V) CITS"
                  )
          return

    # %% Updating methods. Usually called by signals
    def updateAvgVariables(self):
        """ Slot called by the toggling of the average box.
        Toggling on the box clears everything to start a new averaging.
        Toggling off averages and plots the data stored by picking spectra """
        if self.dataLoaded:
            # If toggled on, put everything to zero to be ready to store data
            if self.ui_avgBox.isChecked():
                self.summed_spectra = np.zeros(shape=self.cits_params["zPt"])
                self.nb_summed_spectra = 0
            # If toggled off, plot the data stored
            else:
                if self.nb_summed_spectra == 0:
                    return
                self.drawSpectrum(
                    self.summed_spectra / self.nb_summed_spectra,
                    str(self.nb_summed_spectra) + " spectra averaged",
                )

    def updateAboveValue(self, value):
        """ Slot called when the above slider is changed.
        Changes the value of the textbox to reflect the change """
        self.ui_aboveBar.setValue(value)
        if self.dataLoaded:
            self.ui_aboveLine.setText(
                str(self.mapMax - value * (self.mapMax - self.mapMin) / 100)
            )

    def updateBelowValue(self, value):
        """ Slot called when the below slider is changed.
        Changes the value of the textbox to reflect the change """
        self.ui_belowBar.setValue(value)
        if self.dataLoaded:
            self.ui_belowLine.setText(
                str(self.mapMin + value * (self.mapMax - self.mapMin) / 100)
            )

    def updateMap(self):
        """ Updates the map by redrawing it.
        Updates also the above and below sliders """
        self.drawXYMap(self.ui_indexBox.value())
        self.updateAboveValue(self.ui_aboveBar.value())
        self.updateBelowValue(self.ui_belowBar.value())

    def updateToPointedX(self, event):
        """ Slot called when clicking on the spectrum window.
        Updates map according to the position of the cursor when clicked """
        if (
            event.xdata is not None
            and event.ydata is not None
            and self.dataLoaded
            and not self.toolbar_spec.mode.value
        ):
            # If the scale is in volts, need to divide by dV to have the index
            if self.ui_scaleVoltage.isChecked():
                pointedIndex = int(
                    (event.xdata - self.cits_params["vStart"]) / self.cits_params["dV"]
                )
            else:
                pointedIndex = int(event.xdata)
            # Update the box
            self.ui_indexBox.setValue(pointedIndex)
            # The map updates itself because the change of value of the voltage box calls the drawXYMap method

    def updateIndexBox(self):
        """ Method called by updateWidgets.
        Sets the values of the index box """
        self.ui_indexBox.setMinimum(0)
        self.ui_indexBox.setMaximum(self.cits_params["zPt"] - 1)
        self.ui_indexBox.setSingleStep(1)

    def updateWidgets(self):
        """ Slot called after the reading of the CITS.
        Sets the values combo box (voltage and channels) and draws the map """
        self.updateIndexBox()
        self.ui_channelBox.clear()
        self.ui_channelBox.addItems(self.channelList)
        self.drawXYMap(self.ui_indexBox.value())
        self.updateAboveValue(self.ui_aboveBar.value())
        self.updateBelowValue(self.ui_belowBar.value())
        self.ui_CitsAvgBox.setMaximum(self.cits_params["xPx"])
        logging.debug("Widgets updated !")

    # %% Methods related to spectra
    def averageSpectrum(self, xi, xf, yi, yf):
        """ Averages the spectra contained in the rectangle drawn
        by the points (xi,yi);(xf,yi);(xf,yf) and (xi,yf) """
        if self.dataLoaded:
            chan = self.ui_channelBox.currentIndex()
            if yf < yi:
                yi, yf = yf, yi
            if xf < xi:
                xi, xf = xf, xi
            avg_data = np.mean(self.cits_data[chan, yi:yf, xi:xf], axis=(0, 1))
            xavg = xi + (xf - xi) // 2
            yavg = yi + (yf - yi) // 2
            self.drawSpectrum(
                avg_data, "Average around ({0},{1})".format(xavg, yavg), xavg, yavg,
            )
            return avg_data

    def averageSpectrumWithValues(self):
        """ Averages spectra according their values at a certain voltage """
        if self.dataLoaded:
            voltage_index = self.ui_indexBox.value()
            chan = self.ui_channelBox.currentIndex()
            viewSelected = self.ui_viewSelectedBox.isChecked()
            # midpoint=(self.mapMax+self.mapMin)/2
            midpoint2 = self.mapMax - self.mapMin
            limit_aboveV = self.mapMax - self.ui_aboveBar.value() * midpoint2 / 100
            limit_belowV = self.mapMin + self.ui_belowBar.value() * midpoint2 / 100
            if limit_aboveV < limit_belowV:
                logging.error("Above and below spectra intersect !")
                return
            isAbove = self.cits_data[chan, :, :, voltage_index] > limit_aboveV
            isBelow = self.cits_data[chan, :, :, voltage_index] < limit_belowV
            aboveValues = self.cits_data[chan, :, :, :][isAbove]
            belowValues = self.cits_data[chan, :, :, :][isBelow]

            N_aboveV = len(aboveValues)
            if N_aboveV != 0:
                self.drawSpectrum(
                    np.mean(aboveValues, axis=0),
                    "Average above {0} ({1} spectra averaged)".format(
                        limit_aboveV, N_aboveV
                    ),
                )
                if viewSelected:
                    pixel_pos = np.argwhere(isAbove).T + 0.5
                    color = self.getSpectrumColor(self.nSpectraDrawn - 1)
                    self.ax_map.plot(
                        pixel_pos[0], pixel_pos[1], c=color, marker="o", ls="None",
                    )
                    self.ui_mapWidget.draw()
            N_belowV = len(belowValues)
            if N_belowV != 0:
                self.drawSpectrum(
                    np.mean(belowValues, axis=0),
                    "Average below {0} ({1} spectra averaged)".format(
                        limit_belowV, N_belowV
                    ),
                )
                if viewSelected:
                    pixel_pos = np.argwhere(isBelow).T + 0.5
                    color = self.getSpectrumColor(self.nSpectraDrawn - 1)
                    self.ax_map.plot(
                        pixel_pos[0], pixel_pos[1], c=color, marker="o", ls="None",
                    )
                    self.ui_mapWidget.draw()

    def clearSpectrum(self):
        """ Clears the spectrum window """
        self.ax_spec.clear()
        self.nSpectraDrawn = 0
        self.voltageLine = None
        self.clearShapesClicked()
        self.ui_specWidget.draw()

    def drawSpectrum(self, spectrum_data, label, x=-1, y=-1):
        """ Method called each time a spectrum needs to be plotted.
        Takes care of the derivative and different scales stuff and updates the window.
        """
        shiftX = str(self.ui_shiftXBox.text()).lower()
        if shiftX == "topo":
            shiftX = self.topo[y, x]
        else:
            shiftX = float(shiftX)
        shiftY = self.nSpectraDrawn * float(self.ui_shiftYBox.text())
        if self.dataLoaded and spectrum_data.size != 0:
            if self.ui_logBox.isChecked():
                data_to_plot = np.log(spectrum_data)
            else:
                data_to_plot = spectrum_data

            finalLabel = "{0} - {1}".format(label, self.ui_channelBox.currentText())
            # TODO: remove lastSpectrum with a ref to lines plotted in the spectrum widget
            self.lastSpectrum = [data_to_plot, finalLabel]
            if self.ui_scaleVoltage.isChecked():
                dV = self.cits_params["dV"]
                vStart = self.cits_params["vStart"]
                vEnd = self.cits_params["vEnd"]
                zPt = self.cits_params["zPt"]
                voltages = np.arange(vStart, vEnd, dV) + shiftX
                self.ax_spec.set_xlabel('$V_{bias}$ (eV)')
                self.ax_spec.set_ylabel('LDOS (arb. units)')
                # Check consistency of voltages array with the number of points.
                if len(voltages) != zPt:
                    logging.warning(
                        "Round-off error while computing the voltage array: dV ({}) might be too small."
                        "Computing from number of points instead.".format(dV)
                    )
                    # If this fails, compute from the number of points to be sure to have the same dim as data_to_plot.
                    voltages = np.linspace(vStart, vEnd, zPt) + shiftX
                self.ax_spec.plot(
                    voltages,
                    shiftY + data_to_plot,
                    label=finalLabel,
                    c=self.getSpectrumColor(self.nSpectraDrawn),
                )
                if self.ui_derivBox.isChecked():
                    deriv = sp.signal.savgol_filter(
                        data_to_plot, self.ui_derivNBox.value(), 3, deriv=1, delta=dV
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
                    shiftY + data_to_plot,
                    label=finalLabel,
                    c=self.getSpectrumColor(self.nSpectraDrawn),
                )
                if self.ui_derivBox.isChecked():
                    deriv = sp.signal.savgol_filter(
                        data_to_plot, self.ui_derivNBox.value(), 3, deriv=1, delta=dV
                    )
                    self.ax_spec.plot(
                        shiftY + deriv,
                        label="Derivative of " + finalLabel,
                        marker="o",
                        markersize=3.0,
                        c=self.getSpectrumColor(self.nSpectraDrawn),
                    )
                self.nSpectraDrawn = self.nSpectraDrawn + 1
        if not self.ui_legendBox.isChecked():
            self.ax_spec.legend(loc=0)
        self.ui_specWidget.draw()

    def fitSpectrum(self):
        """ Linear fitting of the last spectrum plotted
        that is stored in lastSpectrum """
        if self.dataLoaded and self.lastSpectrum[0].size != 0:

            dV = self.cits_params["dV"]
            if self.ui_fitCustomCheckbox.isChecked():
                limitL = float(self.ui_fitLowerBox.text())
                limitU = float(self.ui_fitUpperBox.text())
                vStart = min(limitL, limitU)
                vEnd = max(limitL, limitU)
                zPt = self.cits_params["zPt"]
                voltages = np.arange(zPt) * dV + self.cits_params["vStart"]
                # Take the portion to fit
                range_mask = np.logical_and(voltages > vStart, voltages < vEnd)
                dataToFit = self.lastSpectrum[0][range_mask]
            else:
                vStart = self.cits_params["vStart"]
                vEnd = self.cits_params["vEnd"]
                dataToFit = self.lastSpectrum[0]
            X = np.arange(vStart, vEnd, dV)
            popt, pcov = sp.optimize.curve_fit(linearFitFunction, X, dataToFit)
            slope, coef = popt
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
            self.ui_specWidget.draw()

    def getSpectrumColor(self, n):
        """ Returns the color corresponding to a spectrum of given index
        according to spectrumColor """
        return matplotlib.rcParams[
            'axes.prop_cycle'
            ].by_key()['color'][
                n % len(matplotlib.rcParams['axes.prop_cycle'].by_key()['color'])
                ]

    def launchAvgSpectrum(self):
        """ Slot called to average the spectra of the whole map """
        if self.dataLoaded:
            self.averageSpectrum(0, self.cits_params["xPx"], 0, self.cits_params["yPx"])

    def pickSpectrum(self, event):
        """ Method called when a press-and-release event is done
        at the same location of the map. If the average box is checked,
        it will keep the data in storage to average it later.
        Otherwise it plots the spectrum in the corresponding widget """
        if self.dataLoaded and event.xdata is not None and event.ydata is not None:
            PixelX = int(event.xdata/self.metric_ratios['x'])
            PixelY = int(event.ydata/self.metric_ratios['y'])
            chan = self.ui_channelBox.currentIndex()
            # TODO: the slope part should be put in a function
            if "Slope" in self.channelList[chan]:
                dataLogCurrent = np.where(
                    self.cits_data[0, PixelY, PixelX] < 0.01,
                    np.log(self.cits_data[0][PixelY][PixelX]),
                    self.cits_data[0, PixelY, PixelX],
                )
                self.drawSpectrum(
                    dataLogCurrent,
                    "Log of current at [{}, {}]".format(PixelX, PixelY),
                    PixelX,
                    PixelY,
                )
                # self.drawSpectrum(dataLine,"Linear fit of the log at "+str(PixelX)+","+str(PixelY)+"]")
            else:
                self.ui_mapWidget.draw()
                # Add data to the total data to average if in average mod
                if self.ui_avgBox.isChecked():
                    self.summed_spectra += self.cits_data[chan, PixelY, PixelX]
                    self.nb_summed_spectra += 1
                # Plot the data with the desired scale (Volts or index) if in normal mode
                else:
                    self.drawSpectrum(
                        self.cits_data[chan, PixelY, PixelX],
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
            self.dataLoaded and not self.toolbar_map.mode.value
        ):
            print(event, type(event))
            self.currentShape = generateShape(
                event,               # np.array([PixelX,PixelY]),
                self.ui_mapWidget.figure,
                self.fig_topo,
                self.getSpectrumColor(self.nSpectraDrawn),
                self.metric_ratios["x"],
                self.metric_ratios["y"],
            )
            self.motionConnection = self.ui_mapWidget.mpl_connect(
                "motion_notify_event", self.currentShape.update
            )
            logging.debug("PRESS")

    def onReleaseOnMap(self, event):
        """ Slot called when a release event is detected.
        Disconnects the updating of the currentShape and
        launches the appropriate method depending on which button
        was pressed and where it was released """
        if self.dataLoaded and not self.toolbar_map.mode.value:

            if event.xdata is not None and event.ydata is not None:
                square=False
                # Disconnects the updating of the current Shape and gets its coordinates
                self.ui_mapWidget.mpl_disconnect(self.motionConnection)
                xi = self.currentShape.xi
                yi = self.currentShape.yi
                xf = self.currentShape.xf
                yf = self.currentShape.yf
                PixelY = int(event.ydata/self.metric_ratios['y'])
                # If left-click : either a line was drawn or a spectrum picked
                if event.button == 1:
                    # Cut along the XY line if a line is traced (X or Y different)
                    if xf != xi or yf != yi:
                        self.cutAlongLine(int(xi/self.metric_ratios['x']), 
                                          int(xf/self.metric_ratios['x']), 
                                          int(yi/self.metric_ratios['y']), 
                                          int(yf/self.metric_ratios['y']))
                    # Pick spectrum otherwise and change the line shape to a point
                    else:
                        self.currentShape = changeToDot(self.currentShape)
                        self.pickSpectrum(event)
                # If right-click : either a rectangle was drawn or the center of the rectangle to average was picked
                else:
                    square=True
                    if xf != xi or yf != yi:
                        self.averageSpectrum(int(xi/self.metric_ratios['x']), 
                                             int(xf/self.metric_ratios['x']), 
                                             int(yi/self.metric_ratios['y']), 
                                             int(yf/self.metric_ratios['y']))
                    # If X=Y, we need to force the updating of the Shape so it is drawn around the X,Y point
                    # and not starting at X,Y (TODO: This can be refactored)
                    else:
                        n = self.ui_rcAvgBox.value()
                        self.currentShape.forceUpdate(
                            max(0, xi - n),
                            max(0, yi - n),
                            min(self.cits_params["xPx"], xf + n),
                            min(self.cits_params["yPx"], yf + n),
                        )
                        self.averageSpectrum(
                            max(0, xi - n),
                            min(self.cits_params["xPx"], xf + n),
                            max(0, yi - n),
                            min(self.cits_params["yPx"], yf + n),
                        )
                
                #if asked, find peaks that define gap and correlate them
                if self.ui_StatisticsGapBox.isChecked():
                    self.find_peaks_statistics([xi,xf], [yi,yf], square)
                    
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
        self.drawXYMap(self.ui_indexBox.value())

    def cutAlongLine(self, xi, xf, yi, yf):
        """ Method that finds the positions of the pixels forming
        the line from (xi,yi) to (xf,yf).
        The plotting of the spectra is done in cutPlot called at the end """
        # If the line is vertical, the equation is x=xi with y varying from yi to yf
        logging.info(
            "Cut from ({}, {}) to ({}, {}) in {}".format(
                xi, yi, xf, yf, self.ui_channelBox.currentText()
            )
        )
        x_plot, y_plot = findPixelsOnLine(xi, xf, yi, yf)
        if self.ui_viewSelectedBox.isChecked():
            for x, y in zip(x_plot, y_plot):
                self.addToShapesClicked(
                    Dot(
                        x,
                        y,
                        self.ui_mapWidget.figure,
                        self.fig_topo,
                        "yellow",
                        self.metric_ratios["x"],
                        self.metric_ratios["y"],
                    )
                )
        return self.cutPlot(x_plot, y_plot)

    def cutPlot(self, x_plot, y_plot):
        """ 
        Method called by cutAlongLine.
        Plots the spectra for the pixels of positions (x_plot[i],y_plot[i]) 
        if asked, plots FFT of cut
        """
        # Build the data to plot with v as Y and z (number of pixels gone through) as X
        zPt = self.cits_params["zPt"]
        voltage_indices = np.arange(zPt)
        chan = self.ui_channelBox.currentIndex()
        xi = x_plot[0]
        yi = y_plot[0]
        # Matlab convention : Y (v) first then X (z)
        # Plot the built map in a new figure
        fig = pyplot.figure()
        ax = fig.add_subplot(1, 1, 1)
        self.ax_map.text(xi + 0.5, yi + 0.5, str(fig.number))
        self.ui_mapWidget.draw()

        voltage_array = (
            self.cits_params["vStart"] + voltage_indices * self.cits_params["dV"]
            if self.ui_scaleVoltage.isChecked()
            else voltage_indices
        )

        if self.ui_waterfallButton.isChecked():
            waterfallPlot(
                ax,
                voltage_array,
                self.cits_data[chan, y_plot, x_plot],
                offset=float(self.ui_shiftYBox.text()),
            )
        else:
            # Change the scales if needed
            if self.ui_scaleVoltage.isChecked():
                ax.set_ylabel("Bias (V)")
            else:
                ax.set_ylabel("Voltage index")
            ax.set_ylim([voltage_array[0], voltage_array[-1]])

            if self.ui_scaleMetric.isChecked():
                dx = self.metric_ratios['x']
                dy = self.metric_ratios['y']
                if self.ui_plotFFTBox.isChecked():
                    d=np.sqrt((dx * (x_plot[-1] - xi)) ** 2 + (dy * (y_plot[-1] - yi)) ** 2)
                    x_array = 2*np.pi/d*np.arange(
                        -len(x_plot)/(2), len(x_plot)/(2),1
                    )
                    print(len(x_plot), len(x_array))
                    ax.set_xlabel("K Distance ($nm^{-1}$)")
                else:
                    x_array = np.sqrt((dx * (x_plot - xi)) ** 2 + (dy * (y_plot - yi)) ** 2)
                    ax.set_xlabel("Distance (nm)")
            else:
                x_array = np.arange(len(x_plot))
                ax.set_xlabel("Pixels")
                print(len(x_plot), len(x_array))
            ax.set_xlim([x_array[0], x_array[-1]])

            #display FFT if asked 
            if self.ui_plotFFTBox.isChecked():
                ax.set_title("{} - FFT of Cut {}".format(self.ui_channelBox.currentText(), fig.number))
                datafft = np.array([np.abs(np.fft.fftshift(
                    np.fft.fft(self.cits_data[chan, y_plot, x_plot, i])
                )) for i in range(np.shape(self.cits_data)[-1])])
                #work in progress to optimize plotting : Log Norm too long ?
                if self.ui_cbarCustomCheckbox.isChecked():
                    mapData = pyplot.pcolormesh(
                    x_array,
                    voltage_array,
                    np.sqrt(datafft),
                    cmap=self.ui_colorBarBox.currentText()
                    # norm=matplotlib.colors.LogNorm(vmin=float(self.ui_cbarLowerBox.text()), 
                    #                       vmax=float(self.ui_cbarUpperBox.text()))
                      )
                else:
                    mapData = pyplot.pcolormesh(
                    x_array,
                    voltage_array,
                    datafft,
                    cmap=self.ui_colorBarBox.currentText(),
                    norm=matplotlib.colors.LogNorm(vmin=5),
                    )
            else:
                ax.set_title("{} - Cut {}".format(self.ui_channelBox.currentText(), fig.number))
                mapData = pyplot.pcolormesh(
                    x_array,
                    voltage_array,
                    self.cits_data[chan, y_plot, x_plot, :].T,
                    cmap=self.ui_colorBarBox.currentText()
                )
            
            # Colorbar set up
            fig.subplots_adjust(left=0.125, right=0.95, bottom=0.15, top=0.92)
            cbar = fig.colorbar(mapData, shrink=0.9, pad=0.05, aspect=15)
            cbar.ax.yaxis.set_ticks_position("both")
            cbar.ax.tick_params(axis="y", direction="in")
            cbar.set_label('LDOS (arb. units)', rotation=90)
            if self.ui_cbarCustomCheckbox.isChecked():
                mapData.set_clim(
                    float(self.ui_cbarLowerBox.text()),
                    float(self.ui_cbarUpperBox.text()),
                )

            else:
                mapData.set_clim(0)
        fig.show()

        

    def launchBigCut(self):
        """ Launches a cut in the big diagonal """
        if self.dataLoaded:
            self.cutAlongLine(
                0, self.cits_params["xPx"] - 1, 0, self.cits_params["yPx"] - 1
            )
            # Draw the line
            self.currentShape = Line(
                0,
                0,
                self.ui_mapWidget.figure,
                self.fig_topo,
                self.getSpectrumColor(self.nSpectraDrawn),
                self.metric_ratios["x"],
                self.metric_ratios["y"],
                xf=self.cits_params["xPx"] - 1,
                yf=self.cits_params["yPx"] - 1,
            )
            self.currentShape.draw()
            self.ui_mapWidget.draw()

    # %% Methods related to the map
    def clearMap(self):
        """ Unloads the map and clears the map window """
        self.ui_mapWidget.figure.clear()
        self.clearShapesClicked()
        self.cits_data = np.ndarray([])
        # self.cits_params.clear()
        self.mapType = ""
        self.dataLoaded = False

    def createXYMap(self, map_data):
        """ Calls the getMapData function and draws the result
        in the map window with the approriate formatting. """
        fig_map = self.ui_mapWidget.figure
        xPx = self.cits_params["xPx"]
        yPx = self.cits_params["yPx"]
        self.ax_map = self.ui_mapWidget.figure.add_subplot(111)
        self.ui_mapWidget.figure.subplots_adjust()
        if self.ui_ForceAspect_Box.isChecked() :
            self.ax_map.set_aspect("equal")
        else:
            self.ax_map.set_aspect("equal" if xPx == yPx else "auto")
        # Use metric dimensions if the corresponding box is checked
        if self.ui_scaleMetric.isChecked() :
            xL = self.cits_params["xL"]
            yL = self.cits_params["yL"]
            x_m = np.linspace(0, xL, xPx+1)
            y_m = np.linspace(0, yL, yPx+1)
            XYmap = self.ax_map.pcolormesh(
                x_m, y_m, map_data, cmap=self.ui_colorBarBox.currentText(),
             )
            self.ax_map.invert_yaxis()
            self.ax_map.set_ylabel('Y (nm)')
            self.ax_map.set_xlabel('X (nm)')
        # Else, use pixels
        else:
            XYmap = self.ax_map.pcolormesh(
                map_data, cmap=self.ui_colorBarBox.currentText()
            )
            # If the map is an Omicron one, I have to invert the y-axis
            if self.mapType == "Omicron":
                self.ax_map.axis([0, xPx, yPx, 0])
            else:
                self.ax_map.axis([0, xPx, 0, yPx])
        # Colorbar stuff
        cbar = fig_map.colorbar(XYmap, shrink=0.9, pad=0.05, aspect=15, ax=self.ax_map)
        cbar.set_label('LDOS (arb. units)', rotation=90)
        # Recreate the shapes points saved in shapes_clicked
        self.drawShapesClicked(fig_map, self.fig_topo, True)
        return XYmap

    def drawXYMap(self, voltage_index):
        """ Redraws the XYMap
        Called at each change of voltage """
        if self.dataLoaded:
            # Get the data of the map and draw it
            mapData, self.mapMin, self.mapMax = self.getMapData(voltage_index)
            # if self.xy_map is None:
            self.ui_mapWidget.figure.clear()
            self.xy_map = self.createXYMap(mapData)
            # Set title
            self.ax_map.set_title(
                self.cits_name
                + " - "
                + self.ui_channelBox.currentText()
                + "\n"
                + "V="
                + str(
                    self.cits_params["vStart"] + voltage_index * self.cits_params["dV"]
                )
                + " (V)"
            )
            self.setColorMapLims()
            self.ui_mapWidget.draw()
            # Plot a dashed line at X=voltage if asked
            if self.ui_vLineBox.isChecked():
                self.drawVoltageLine(voltage_index)

    def getMapData(self, v):
        """ Returns an array built from the data loaded that can be used to
        display a map at fixed voltage """
        mapData = self.cits_data[self.ui_channelBox.currentIndex(), :, :, v]
        return mapData, np.min(mapData), np.max(mapData)

    def setColorMap(self, color_map):
        if self.xy_map is not None:
            self.xy_map.set_cmap(color_map)
            self.ui_mapWidget.draw()

    def setColorMapMin(self, min_value_text):
        if self.xy_map is not None:
            min_value = self.mapMin
            if self.ui_cbarCustomCheckbox.isChecked():
                try:
                    min_value = float(min_value_text)
                except ValueError:
                    pass
            self.xy_map.set_clim(vmin=min_value)
            self.ui_mapWidget.draw()

    def setColorMapMax(self, max_value_text):
        if self.xy_map is not None:
            max_value = self.mapMax
            if self.ui_cbarCustomCheckbox.isChecked():
                try:
                    max_value = float(max_value_text)
                except ValueError:
                    pass
            self.xy_map.set_clim(vmax=max_value)
            self.ui_mapWidget.draw()

    def setColorMapLims(self):
        self.setColorMapMin(self.ui_cbarLowerBox.text())
        self.setColorMapMax(self.ui_cbarUpperBox.text())

    # %% Methods related to the voltage guide line in the spectra window

    def toggleVoltageLine(self, checked):
        """ Slot attached to the toggle of the voltageLine checkbox. """
        if checked:
            self.drawVoltageLine(self.ui_indexBox.value())
        else:
            self.clearVoltageLine()

    def createVoltageLine(self, voltage):
        """ Creates the voltage line. Called in drawVoltageLine. """
        # Plot the dashed line
        self.voltageLine = self.ax_spec.axvline(voltage, color="k", linestyle="--")

    def drawVoltageLine(self, voltage_index):
        """ Draws the vertical line at the given voltage index. Takes care of the creation if needed. """
        # self.clearVoltageLine()
        if self.dataLoaded:
            # Get the current voltage : real voltage if the scale box is checked, voltage index otherwise
            if self.ui_scaleVoltage.isChecked:
                voltage = (
                    self.cits_params["vStart"] + voltage_index * self.cits_params["dV"]
                )
            else:
                voltage = voltage_index
            if self.voltageLine is None:
                self.createVoltageLine(voltage)
            else:
                self.voltageLine.set_xdata((voltage, voltage))
            self.ui_specWidget.draw()

    def clearVoltageLine(self):
        """ Removes the vertical voltage line. """
        if self.voltageLine is not None:
            self.ax_spec.lines.pop(self.ax_spec.lines.index(self.voltageLine))
            self.voltageLine = None
            self.ui_specWidget.draw()

    # Post-processing methods
    def addChannel(self, new_channel_data, channel_name):
        """ Adds a channel to the whole data array """
        try:
            self.cits_data = np.concatenate(
                (self.cits_data, new_channel_data[np.newaxis, ...]), axis=0
            )
        except MemoryError:
            logging.error(
                "The data is too big ! Or the memory too small... Aborting addChannel !"
            )
            return
        # Updates the channel list
        self.channelList.append(channel_name)
        self.ui_channelBox.clear()
        self.ui_channelBox.addItems(self.channelList)

    def averageCITSOnAxis(self, direction="x"):
        """ Slot that averages the spectras in one direction of the loaded
        CITS and replaces the loaded CITS by the result """
        if self.dataLoaded:  # Check if a CITS was loaded
            # Get the needed params
            Navg = int(self.ui_CitsAvgBox.value())
            try:
                new_data = directionAverageCITS(self.cits_data, Navg, direction)
            except MemoryError:
                logging.error("Not enough memory to average spectra !")
                return

            # If eveything went well, clear the current map and replace it by
            # the created data.
            self.clearMap()
            self.cits_name = "Average_on_{}_every_{}_spectra".format(direction, Navg)
            self.dataLoaded = True
            self.cits_data = new_data
            self.cits_params["yPx"] = new_data.shape[1]
            self.cits_params["xPx"] = new_data.shape[2]
            self.updateWidgets()

    def extractDerivative(self, chan_index_to_deriv):
        """ Extracts the derivative of a channel and adds it to the data """
        if self.dataLoaded:
            derivData = sp.signal.savgol_filter(
                self.cits_data[chan_index_to_deriv],
                9,
                2,
                deriv=1,
                delta=self.cits_params["dV"],
                axis=-1,
            )
            # Add the channel to the data
            self.addChannel(
                derivData, "Derivative of " + self.channelList[chan_index_to_deriv]
            )
            self.ui_channelBox.setCurrentIndex(len(self.channelList) - 1)

    def extractFFT(self, chan_index_to_fft, axis_to_fft):
        """ Extracts the FFT of a channel and adds it to the data """
        if self.dataLoaded:
            FFTData = np.fft.fft(self.cits_data[chan_index_to_fft], axis=axis_to_fft)
            # Add the channel to the data
            self.addChannel(FFTData, "FFT of " + self.channelList[chan_index_to_fft])
            self.ui_channelBox.setCurrentIndex(len(self.channelList) - 1)

#%% Extract gif

    def make_gif(self):
        n = self.ui_gifstartBox.value()
        m = self.ui_gifstopBox.value()
        step = self.ui_gifstepBox.value()
        xPx = self.cits_params["xPx"]
        yPx = self.cits_params["yPx"]
        assert m>n
        frames = []
        chan = self.ui_channelBox.currentIndex()

        for i in np.arange(n,m, step):
            figif = pyplot.figure()
            ax_gif = figif.add_subplot(111)
            figif.subplots_adjust()
            # Set title
            pyplot.title(
                "$V_{BIAS}$="
                + str(
                    round(
                        self.cits_params["vStart"] + i * self.cits_params["dV"]
                        ,3)
                )
                + " V"
            )
            if self.ui_ForceAspect_Box.isChecked():
                ax_gif.set_aspect("equal")
            # Use metric dimensions if the corresponding box is checked (work in progress)
            if self.ui_scaleMetric.isChecked() :
                xL = self.cits_params["xL"]
                yL = self.cits_params["yL"]
                x_m = np.linspace(0, xL, xPx+1)
                y_m = np.linspace(0, yL, yPx+1)
                XYmap = ax_gif.pcolormesh(
                    x_m, y_m, self.cits_data[chan][:,:,i],
                    cmap=self.ui_colorBarBox.currentText()
                )
                ax_gif.invert_yaxis()
                ax_gif.set_ylabel('Y (nm)')
                ax_gif.set_xlabel('X (nm)')
            # Else, use pixels
            else:
                XYmap = ax_gif.pcolormesh(
                    self.cits_data[chan][:,:,i], cmap=self.ui_colorBarBox.currentText()
                )
                # If the map is an Omicron one, invert the y-axis
                if self.mapType == "Omicron":
                    pyplot.axis([0, xPx, yPx, 0])
                else:
                    pyplot.axis([0, xPx, 0, yPx])
            # Colorbar
            cbar = pyplot.colorbar(XYmap, shrink=0.9, pad=0.05, aspect=15)
            cbar.set_label('LDOS (arb. units)', rotation=90)
            if self.ui_cbarCustomCheckbox.isChecked():
                XYmap.set_clim(
                    float(self.ui_cbarLowerBox.text()),
                    float(self.ui_cbarUpperBox.text()),
                )
            else:
                XYmap.set_clim(0)

            pyplot.savefig(self.wdir+'/giftemp{}.png'.format(i))

            frame = Image.open(self.wdir+'/giftemp{}.png'.format(i))
            frames.append(frame)

            pyplot.close()

            os.remove(self.wdir+'/giftemp{}.png'.format(i)) 

        frame_one = frames[0]

        frame_one.save(self.wdir+'/'+
                       self.cits_name+'channel '+self.channelList[chan]+'.gif',
                       format="GIF", append_images=frames,
                        save_all=True, duration=500, loop=0, dpi=50)

#%% methods related to correlation between peak position and gap value
    def find_peaks_statistics(self, x, y,square=False):
        DEgap = 0.026 #V
        centergap = -0.1 #V
        chan = self.ui_channelBox.currentIndex()

        voltage_array = np.arange(
            int((np.abs(self.cits_params["vStart"]-centergap)-2*DEgap)/self.cits_params["dV"]),
            int((np.abs(self.cits_params["vStart"]-centergap)+2*DEgap)/self.cits_params["dV"])
            )

        print(len(x))
        
        peaks_index = []
        topo_array = []
        if square:
            average=False #work in progress
            if average:
                n_avg = int(self.ui_CitsAvgBox.value())
                self.averageSpectrum()
                avg_data = np.mean(self.cits_data[chan, y[0]:y[1], x[0]:x[1]], axis=(0, 1))

            else:
                ldos_array = self.cits_data[chan, y[0]:y[1], x[0]:x[1],voltage_array[0]-1:voltage_array[-1]]
                print(np.shape(ldos_array), np.shape(self.cits_data[chan]), np.shape(self.topo))
    
                # find peaks
                nTx = np.shape(self.topo)[1]
                nTy = np.shape(self.topo)[0]
                nSx = np.shape(self.cits_data[chan])[1]
                nSy = np.shape(self.cits_data[chan])[0]
                rx = nTx/nSx
                ry = nTy/nSy
                print(int(round(y[0]*ry)),int(round(x[0]*rx)), int(round(y[1]*ry)), int(x[1]*rx))
                print(np.shape(self.topo[int(round(y[0]*ry)):int(round(y[1]*ry)), 
                                         int(round(x[0]*rx)):int(round(x[1]*rx))]))
                for i, ldos in enumerate([
                        ldos_array[n,m,:] 
                        for n in range(np.shape(ldos_array)[0]) 
                        for m in range(np.shape(ldos_array)[1])
                        ]):
                    print(i)
                    peaks, properties = find_peaks(ldos, prominence=0.5, width=0.01)
                    if len(peaks)==2 :
                        peaks_index.append(peaks)
                        if self.fig_topo:
                            n = int(round(y[0]*ry))+int(round((i//(np.shape(ldos_array)[0])-1)*ry))
                            m = int(round(x[0]*rx))+int(round((i%(np.shape(ldos_array)[0]))*rx))
                            topo_array.append(self.topo[n,m])

        else:
            x_plot, y_plot = findPixelsOnLine(x[0], x[1], y[0], y[1])
            ldos_array = self.cits_data[chan, y_plot, x_plot, voltage_array[0]-1:voltage_array[-1]]
            print(np.shape(ldos_array))

            figTest = 1
            if figTest and not(square):
                # Plot the in a new figure
                fig = pyplot.figure()
                ax = fig.add_subplot(1, 1, 1)  
                waterfallPlot(
                    ax,
                    voltage_array*self.cits_params["dV"],
                    ldos_array, #self.cits_data[chan, y_plot, x_plot],
                    offset=float(self.ui_shiftYBox.text()),
                )
                
                if self.ui_scaleVoltage.isChecked():
                    ax.set_xlabel("Bias (V)")
                else:
                    ax.set_xlabel("Voltage index")
    
            # find peaks
            for i, ldos in enumerate(ldos_array):
                peaks, properties = find_peaks(ldos, prominence=0.4, width=0.01)
                if len(peaks)==2 :
                    peaks_index.append(peaks)
                    if figTest:
                        ax.plot((voltage_array[0]+peaks)*self.cits_params['dV'], 
                                ldos[peaks]+i*float(self.ui_shiftYBox.text()), 'x', color='r')
                # print(properties["prominences"], properties["widths"])
        
        print (np.shape(peaks_index))
        fig_visualize = pyplot.figure()
        if self.fig_topo is None:
            visualize_statistics(
                fig_visualize, 
                self.cits_params["vStart"]+(voltage_array[0]+peaks_index)*self.cits_params['dV']
                )
        else:
            visualize_statistics_topo(
                fig_visualize, 
                self.cits_params["vStart"]+(voltage_array[0]+peaks_index)*self.cits_params['dV'],
                topo_array
                )        
        
        return (voltage_array[0]+peaks_index)*self.cits_params['dV']
        
