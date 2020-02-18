# -*- coding: utf-8 -*-

# Set explictly the backend to Qt for consistency with pyqt.
import matplotlib
matplotlib.use('qt5agg')
from scampy.ui_citswidget import Ui_CitsWidget
from matplotlib.backends.backend_qt5 import NavigationToolbar2QT as NavigationToolbar
import numpy as np
import os.path
import scipy as sp
import scipy.interpolate
import scipy.optimize
import scipy.signal
import scipy.io
from matplotlib import pyplot as pyplot
from matplotlib.patches import Circle
import matplotlib.backend_bases
import struct
import PyQt5.QtWidgets as QtWidgets
from scampy.shape import generateShape, changeToDot

# noinspection PyPep8Naming
class CitsWidget(QtWidgets.QMainWindow, Ui_CitsWidget):
    #%% Building methods
    def __init__(self, parent):
        """ Builds the widget with parent widget in arg """
        QtWidgets.QMainWindow.__init__(self, parent)
        # Set up the user interface from Designer.
        self.setupUi(self)
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
        self.m_saveCsvButton = QtWidgets.QPushButton('Save as CSV') # Add csv save button
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
        # Other parameters used after map loading
        self.mapType = ""
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
        self.m_mapWidget.mpl_connect('button_press_event', self.onPressOnMap)
        self.m_mapWidget.mpl_connect('button_release_event', self.onReleaseOnMap)
        self.m_specWidget.mpl_connect('button_press_event', self.updateToPointedX)
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


    #%% Reading and loading CITS methods
    def askCits(self):
        """ Slot that only takes care of opening a file dialog for the user to select one or several CITS - Returns the CITS paths """
        cits_names_and_ext = QtWidgets.QFileDialog.getOpenFileNames(self,
                                                                    "Choose a CITS file to read or several to average",
                                                                    self.wdir,
                                                                    " RHK file (*.sm4);;Matlab file (*.mat);;3D binary file (*.3ds);;Ascii file (*.asc);;Text file (*.txt)")
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
            elif extension == "" or extension == "mat":
                self.clearMap()
                self.mapType = "Sm4 to .mat"
                self.dataLoaded = self.loadCitsSm4(cits)
            elif extension == "sm4":
                self.clearMap()
                print('sm4')
                self.mapType = "Sm4"
                self.dataLoaded = self.readCitsSm4Bin(cits)
            else:
                print("Extension not recognized")
                self.m_statusBar.showMessage("Extension not recognized")
                return
            # After reading, check if the data was read correctly and update the working directory and the map name
            if self.dataLoaded:
                self.wdir = os.path.dirname(cits)
                self.mapName = os.path.basename(cits)
                print(self.mapName + " read as a " + self.mapType + " map")
                self.m_statusBar.showMessage(self.mapName + " read as a " + self.mapType + " map")
            else:
                print('Problem while reading ' + cits)
                self.m_statusBar.showMessage('Problem while reading ' + cits)
                return
            # If only one Cits was selected, there is no need to run the averaging code so return after drawing the topo and updating of the widgets
            print(n_cits)
            if n_cits == 1:
                self.updateWidgets()
                print('updated')
                self.drawTopo()
                print('drawn')
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

    def readConfig(self):
        config = {}
        with open("config.txt") as f:
            for line in f:
                (key, val) = line.split()
                config[key] = val
        if 'working_directory' in config.keys():
            if os.path.exists(config['working_directory']):
                self.wdir = config['working_directory']
            else:
                print("{} is not a valid path ! Check your config file.".format(config['working_directory']))

        if 'matplotlib_stylesheet' in config.keys():
            try:
                matplotlib.pyplot.style.use(config['matplotlib_stylesheet'])
            except IOError:
                print("{} was not found in the .matplotlib folder. Using default parameters for matplotlib...".format(config['matplotlib_stylesheet']))

        if 'autoload' in config.keys():
            autoload = config['autoload'].lower()
            if 'no' in autoload or 'false' in autoload:
                self.autoload = False
            else:
                self.autoload = True
        else:
            self.autoload = False

        if 'default_cmap' in config.keys():
            self.default_cmap = config['default_cmap']
        else:
            self.default_cmap = 'magma_r'

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
        topopath = os.path.join(os.path.dirname(filepath), 'Topo.txt')
        if os.path.exists(topopath):
            self.readTopo(topopath)
        return True
    
    def loadCitsSm4(self, filepath):
        zSpectro = False
        #file loading : see sm4_m_reading
        loadedData = scipy.io.loadmat(filepath)
        Spatial = loadedData['Spatial']
        Spectral = loadedData['Spectral']
        Image = Spatial['TopoData'][0, 0]
        FImage = Image[0, 0]# Image forward
        BImage = Image[0, 1]#Image bacward
        spectraType = Spectral['type'][0, 0]#see scipy.org : scipy.io.loadmat
        if spectraType == 'Point':
            print('You didnt load a CITS. Use sm4_reader to read this data')
 
        #find out how many measurement locations were taken along a line :
        #each measurement is taken at a coordinate, but several measurements (repetitions) are taken on the same spot. Eg if repetitions = 4, 2 measurements, one forward, one backward (saved in right direction)
        numOfMeasurements = np.shape(Spectral['xCoord'][0, 0])[0]
        repetitions = 0
        pointdiff = 0
        repindex = 0
        while pointdiff == 0:  #step through data points until there is a difference between the current (x,y) point and the next (x1,y1)point
            pointdiff = (Spectral['xCoord'][0,0][repindex+1]- Spectral['xCoord'][0,0][repindex]) + (Spectral['yCoord'][0,0][repindex+1]- Spectral['yCoord'][0,0][repindex]) #adds the x and y difference values together.  if this is anything other than 0, we have found the limit of the repetitions
            repetitions = repetitions+1
            repindex = repindex+1
        numberOfPlots = numOfMeasurements/repetitions #this is the number of distinct plotting locations
        print(numOfMeasurements)
        #Load spectral data
        try:
            SpectralData_y = Spectral['dIdV_Line_Data'][0, 0]
        except ValueError:
            print('The file seems to contain no dIdV_Line_Data. Use sm4_reader to read this data or check the .mat file loading/Saving')
        SpectralData_x = Spectral['xdata'][0, 0] * 1000#mV
        # Center coordinates and metric dimensions in nm
        xL = Spatial['width'][0, 0] * (10 ** 9)
        yL = Spatial['height'][0, 0] * (10 ** 9)
        xC = Spatial['xoffset'][0, 0] * (10 ** 9)
        yC = Spatial['yoffset'][0, 0] * (10 ** 9)
        #repetition_index = 1 #0(forw),1(back),2(forw)... according to the number of repetitions.
        #size spatial data
        xPx = int((Spatial['lines'][0, 0]))
        yPx = int((Spatial['points'][0, 0]))
        #size spectral data ( not necessarily the same in RHK )
        xSpec = int(np.sqrt(numberOfPlots))
        ySpec = int(np.sqrt(numberOfPlots))
        zPt = int(len(SpectralData_y[:, 0]))#nbre of points in spec data
#        x_m = np.linspace(0, xL*1e+9, xPx)#?
#        y_m = np.linspace(0, yL*1e+9, yPx)
        try:
            self.topo = np.zeros(shape=(xPx, yPx))
            self.m_data = np.zeros(shape=(repetitions, ySpec, xSpec, zPt))
            print(np.shape(self.m_data))
        except MemoryError:
            print("The data is too big ! Or the memory too small...")
            return False
            
        #in Spectraldata_y, dIdV info corresponds to the spec data saved from left to right and increasing y(downwards in RHK), with same nbre of repetions at each spot
        for r in range(repetitions):#even : forward, odd : backwards
            for y in range(ySpec):#len(SpectralData_y[:,0])):
                for x in range(xSpec):
                    self.m_data[r][y][x] = SpectralData_y[:, (xSpec*y+x)*repetitions+r]
        
        patch = []
        for m in range(0, numOfMeasurements, repetitions): #iterate over the number of locations
                            patch.append(Circle((-xC + xL/2 + Spectral['xCoord'][0,0][m]*1e+9,- yC + yL/2 + Spectral['yCoord'][0,0][m]*1e+9), yL/5000,facecolor='r',edgecolor='None'))

        self.channelList = ['Data {}'.format(i) for i in range(repetitions)]

        self.m_params = {"xPx": xSpec, "yPx": ySpec, "xL": xL, "yL": yL, "zPt": zPt,
                         "vStart": SpectralData_x[0],"vEnd": SpectralData_x[-1], "dV": abs(SpectralData_x[-1] - SpectralData_x[0])/zPt,
                         "Patch": patch}

        self.topo = FImage
        print(np.shape(self.topo))
        
        #!!! create average Data :
        average = True
        if average:
            average = np.zeros(shape=(ySpec, xSpec, zPt))
            for y in range(ySpec):#len(SpectralData_y[:,0])):
                for x in range(xSpec):
                    for r in range(repetitions):
                        average[y][x] += SpectralData_y[:, (xSpec*y+x)*repetitions+r]/repetitions
            self.addChannel(average, "average")
        return True

    def string(self, array):
        array = list(map(lambda x: chr(x), array))
        array = ("".join(array))
        return array

    def readCitsSm4Bin(self, filepath):
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
        Signature = self.string(np.fromfile(f, dtype=np.uint16, count=18))
        print(Signature)
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

        # Read and record each Pages Headers and Data
        PageHeader = []
        TextStrings = []
        Spatial = []
        Spectral = []
        SpectralInfo = []
        SpatialInfo = []
        topocount = 0
        currentcount = 0
        dIdV_Line_Speccount = 0
        for j in range(PageIndexHeader_PageCount):
            #f.seek(nbytes,0 = a partir du début) va au byte n + 1
            f.seek(PageIndex[j]['ObjectList'][0]['Offset'], 0)

        # Read Page Header    
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
            # Skip flags and reserved data
            f.seek(1+3+60, 1);
            # Read the Object List
            for i in range(PageHeader[j]['ObjectListCount']):
                dict1 = {'ObjectID': int(np.fromfile(f, dtype=np.uint32, count=1)[0]),
                         'Offset': int(np.fromfile(f, dtype=np.uint32, count=1)[0]),
                         'Size': int(np.fromfile(f, dtype=np.uint32, count=1)[0])}
                PageHeader[j]['ObjectList'].append(dict1)
    
            # PageheaderObjectList
            c = int(np.fromfile(f, dtype=np.uint16, count=1)[0])
            d = {'strLabel': self.string(np.fromfile(f, dtype=np.uint16, count=c))} # eg "current image"
            count = int(np.fromfile(f, dtype=np.uint16, count=1)[0])
            d['strSystemText'] = self.string(np.fromfile(f, dtype=np.uint16, count=count))
    
            count = int(np.fromfile(f, dtype=np.uint16, count=1)[0])
            d['strSessionText'] = self.string(np.fromfile(f, dtype=np.uint16, count=count))
    
            count = int(np.fromfile(f, dtype=np.uint16, count=1)[0])
            d['strUserText'] = self.string(np.fromfile(f, dtype=np.uint16, count=count))
    
            count = int(np.fromfile(f, dtype=np.uint16, count=1)[0])
            d['strPath'] = self.string(np.fromfile(f, dtype=np.uint16, count=count))
    
            count = int(np.fromfile(f, dtype=np.uint16, count=1)[0])
            d['strDate'] = self.string(np.fromfile(f, dtype=np.uint16, count=count))
    
            count = int(np.fromfile(f, dtype=np.uint16, count=1)[0])  # DAQ time
            d['strTime'] = self.string(np.fromfile(f, dtype=np.uint16, count=count))
    
            count = int(np.fromfile(f, dtype=np.uint16, count=1)[0])  # physical units of x axis
            d['strXUnits'] = self.string(np.fromfile(f, dtype=np.uint16, count=count))
    
            count = int(np.fromfile(f, dtype=np.uint16, count=1)[0])
            d['strYUnits'] = self.string(np.fromfile(f, dtype=np.uint16, count=count))
    
            count = int(np.fromfile(f, dtype=np.uint16, count=1)[0])
            d['strZUnits'] = self.string(np.fromfile(f, dtype=np.uint16, count=count))
    
            count = int(np.fromfile(f, dtype=np.uint16, count=1)[0])
            d['strYLabel'] = self.string(np.fromfile(f, dtype=np.uint16, count=count))
    
            count = int(np.fromfile(f, dtype=np.uint16, count=1)[0])
            d['strStatusChannelText'] = self.string(np.fromfile(f, dtype=np.uint16, count=count))
    
            count = int(np.fromfile(f, dtype=np.uint16, count=1)[0])
            # contains last saved line count for an image data page
            d['strCompletedLineCount'] = self.string(np.fromfile(f, dtype=np.uint16, count=count))
    
            count = int(np.fromfile(f, dtype=np.uint16, count=1)[0])
            # Oversampling count for image data pages
            d['strOverSamplingCount'] = self.string(np.fromfile(f, dtype=np.uint16, count=count))
    
            count = int(np.fromfile(f, dtype=np.uint16, count=1)[0])
            # voltage at which the sliced image is created from the spectra page.  empty if not a sliced image
            d['strSlicedVoltage'] = self.string(np.fromfile(f, dtype=np.uint16, count=count))
    
            count = int(np.fromfile(f, dtype=np.uint16, count=1)[0])
            # PLLPro status text: blank, master or user
            d['strPLLProStatus'] = self.string(np.fromfile(f, dtype=np.uint16, count=count))
    
            count = int(np.fromfile(f, dtype=np.uint16, count=1)[0])
            # ZPI controller item's set-point unit
            d['strSetpointUnit'] = self.string(np.fromfile(f, dtype=np.uint16, count=count))
    
            count = int(np.fromfile(f, dtype=np.uint16, count=1)[0])
            # stores value of CH1 and CH2 if they are in hardware space
            d['strCHDriveValues'] = self.string(np.fromfile(f, dtype=np.uint16, count=count))

            TextStrings.append(d)

            # Page Header – Sequential Data Page - Get Datas
            Data = []
            f.seek(PageIndex[j]['ObjectList'][1]['Offset'],0)
            Data.append(np.fromfile(f, dtype=np.uint32, count=int(round(PageIndex[j]['ObjectList'][1]['Size']/4))))
#            print(np.shape(Data))
            #/4 because total data size is divided by the number of bytes that use each 'long' data
            ScaleData = PageHeader[j]['ZOffset']+Data[0]*PageHeader[j]['ZScale']
            ScaleData = np.reshape(ScaleData,(PageHeader[j]['Width'],PageHeader[j]['Height']), order="F")
            #order Fortran = "F" to match readCITSsm4File function
            print(np.shape(ScaleData))

            ###################### Spatial Data
            # is it label topo ? # cf p.28 of SM4 DATA FILE FORMAT V5.pdf
            print(TextStrings[j]['strLabel'])
            if TextStrings[j]['strLabel'] == 'Topography' and PageHeader[j]['PageType_DataSource'] == 1:
                topocount += 1
                Spatialpagenumber = j  # records last spatial page; to read the rest of the spatial data later
    # !!! rotate ? transpose ?
                Spatial.append({'TopoData': np.rot90(ScaleData,3)})
                SpatialInfo.append({'TopoUnit': TextStrings[j]['strZUnits']})
                if topocount == 1:
                    FImage = Spatial[-1]['TopoData']  # Image forward
                    print(np.shape(FImage))
                elif topocount ==2:
                    BImage = Spatial[-1]['TopoData']  # Image bacward
                else:
                    print('there is more topo data than expected')
            # is it spatial Current data ?
            elif TextStrings[j]['strLabel'] == 'Current' and PageHeader[j]['PageType_DataSource'] == 2:
                currentcount += 1
                Spatial.append({'CurrentData': np.rot90(ScaleData,3)})  
                Spatial[-1]['CurrentUnit'] = TextStrings[j]['strZUnits']
                Spatialpagenumber = j
                if currentcount == 1:
                    FImage_I = Spatial[-1]['CurrentData']  # Image forward
                elif currentcount ==2:
                    BImage_I = Spatial[-1]['CurrentData']  # Image bacward
                else:
                    print('there is more current data than expected')
    
            ###################### Spectral Data - can be Point or Line
            # Is it Spectral Point(38) ? Not taken in charge currently
            elif PageHeader[j]['PageType_DataSource'] == 38:
                print('You didnt load a CITS. Use sm4_reader to read this data')
                
            # CITS is recarded as a line of LIA
            # Is it Spectral Line(16) LIA ?
            elif TextStrings[j]['strLabel'] == 'LIA' and PageHeader[j]['PageType_DataSource'] == 16:
                dIdV_Line_Speccount += 1;
                Spectral.append({'dIdV_Line_Data': ScaleData})
                Linespectrapagenumber = j
                if dIdV_Line_Speccount == 1 :
                    SpectralData_y = Spectral[-1]['dIdV_Line_Data']
                    #Get spectra coordinates offset
                    for objectNbre in range(PageHeader[j]['ObjectListCount']):
                        # 7 == Tip track Info Header
                        if PageHeader[j]['ObjectList'][objectNbre]['ObjectID'] == 7:
                            TipTrackInfo_offset = PageHeader[j]['ObjectList'][objectNbre]['Offset']
                            TipTrackInfo_Size = PageHeader[j]['ObjectList'][objectNbre]['Size']
                        # 8 == Tip track Data
                        elif PageHeader[j]['ObjectList'][objectNbre]['ObjectID'] == 8:
                            TipTrackData_offset = PageHeader[j]['ObjectList'][objectNbre]['Offset']
                            TipTrackData_Size = PageHeader[j]['ObjectList'][objectNbre]['Size']
                else:
                    print('there is more spectral data than expected')
    
            # We could add 'Spaectral line Current', or 'PLL Amplitude' ;
            # 'PLL Phase' ; 'dF' ; 'PLL Drive'spectra, AFM, ... in the same way
            else:
                print('Data Type '+str(TextStrings[j]['strLabel'])
                      + ' or PageType_DataSource '
                      + str(PageHeader[j]['PageType_DataSource'])
                      + ' not found. Check spelling or add it.')
    
        ###################### Get spectra coordinates
        # Nbre of measurements taken along a line :
        nbreScans = np.shape(SpectralData_y)[1]
        print(nbreScans)
        xCoord = np.zeros(nbreScans)
        yCoord = np.zeros(nbreScans)
        f.seek(TipTrackData_offset,0)  # Go to beggining of header
        for i in range(nbreScans):
            f.seek(4, 1)  # skip start time
            a = np.float(np.fromfile(f, dtype=np.float32, count=1)[0])
    #        print(a)
            xCoord[i] = a
            yCoord[i] = np.float(np.fromfile(f, dtype=np.float32, count=1)[0])
            f.seek(16, 1)  # skip dx dy xcumul ycumul fields
#        print(xCoord[0])
#        print((PageHeader[Linespectrapagenumber]['Width']))

        SpectralData_x = PageHeader[Linespectrapagenumber]['XOffset'] + PageHeader[Linespectrapagenumber]['XScale'] * np.array(list(range(0, PageHeader[Linespectrapagenumber]['Width']))) * 1000.0#mV
#        print(np.shape(SpectralData_x))

        # each measurement is taken at a coordinate,
        # but several measurements (repetitions) are taken on the same spot.
        # Eg if repetitions = 4, 2, one forward, one backward (saved in right direction)
        repetitions = 0
        pointdiff = 0
        repindex = 0
        # step through data points until there is a difference between the current (x,y) point and the next (x1,y1) point
        while pointdiff == 0: 
            pointdiff = (xCoord[repindex+1]- xCoord[repindex]) + (xCoord[repindex+1]- xCoord[repindex]) #adds the x and y difference values together.  if this is anything other than 0, we have found the limit of the repetitions
            repetitions = repetitions+1
            repindex = repindex+1
        numberOfPlots = nbreScans/repetitions  # number of distinct plotting locations
        # Center coordinates and metric dimensions in nm
        xL = abs(PageHeader[Spatialpagenumber]['XScale']*PageHeader[Spatialpagenumber]['Width']) * (10 ** 9)
        yL = abs(PageHeader[Spatialpagenumber]['YScale']*PageHeader[Spatialpagenumber]['Height']) * (10 ** 9)
        xC = PageHeader[Spatialpagenumber]['XOffset'] * (10 ** 9)
        yC = PageHeader[Spatialpagenumber]['YOffset'] * (10 ** 9) * (10 ** 9)
        # repetition_index = 1 #0(forw),1(back),2(forw)...
        # size spatial data
        xPx = int(PageHeader[Spatialpagenumber]['Height'])
        yPx = int(PageHeader[Spatialpagenumber]['Width'])
        # size spectral data ( not necessarily the same in RHK )
        xSpec = int(np.sqrt(numberOfPlots))
        ySpec = int(np.sqrt(numberOfPlots))
        zPt = np.shape(SpectralData_y)[0]  #  nbre of points in spec data
    #        x_m = np.linspace(0, xL*1e+9, xPx)#?
    #        y_m = np.linspace(0, yL*1e+9, yPx)
        try:
            self.topo = np.zeros(shape=(xPx, yPx))
            self.m_data = np.zeros(shape=(repetitions, ySpec, xSpec, zPt))
            print(np.shape(self.m_data))
        except MemoryError:
            print("The data is too big ! Or the memory too small...")
            return False

        # in Spectraldata_y, dIdV info corresponds to the spec data saved
        # from left to right and increasing y(downwards in RHK),
        # with same nbre of repetions at each spot
        for r in range(repetitions):  # even : forward, odd : backwards
            for y in range(ySpec):
                for x in range(xSpec):
                    self.m_data[r][y][x] = SpectralData_y[:, (xSpec*y+x)*repetitions+r]

        patch = []
        for m in range(0, int(numberOfPlots), repetitions):  # iterate over the number of locations
            patch.append(Circle((-xC + xL/2 + xCoord[m]*1e+9, - yC + yL/2 + yCoord[m]*1e+9), yL/5000, facecolor='r', edgecolor='None'))

        self.channelList = ['Data {}'.format(i) for i in range(repetitions)]

        self.m_params = {"xPx": xSpec, "yPx": ySpec, "xL": xL, "yL": yL,
                         "zPt": zPt, "vStart": SpectralData_x[0],
                         "vEnd": SpectralData_x[-1],
                         "dV": abs(SpectralData_x[-1] - SpectralData_x[0])/zPt,
                         "Patch": patch}

        self.topo = FImage
        print(np.shape(self.topo))

        # create average Data :
        average = 1
        if average:
            average = np.zeros(shape=(ySpec, xSpec, zPt))
            for y in range(ySpec):#len(SpectralData_y[:,0])):
                for x in range(xSpec):
                    for r in range(repetitions):
                        average[y][x] += SpectralData_y[:, (xSpec*y+x)*repetitions+r]/repetitions
            self.addChannel(average, "average")
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
        for y in range(yPx):
            for x in range(xPx):
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
        f.close()
        if half:
            self.channelList = self.channelList[0:nChannels // 2]
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
        return True


#%% Reading and loading topo images methods    
    def readTopo(self, filepath):#used for txt files
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
        yPx, xPx = self.topo.shape
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
            # Get parameters
            w = self.m_params["xL"]
            h = self.m_params["yL"]
            xPx = len(self.topo[0])
            yPx = len(self.topo[:,0])
            yspec = self.m_params["yPx"]
            # If yspec==1, it is a Line Spectro so I need to call the specific method to draw the topo
            if yspec == 1:
                print('ySpec = 1')
                self.drawLineTopo()
                return
        
#!!! decide linefit
            lineFit = True
            # Line fitting if necessary
            if lineFit:
                self.topo = self.levelTopo()
            # Set up the figure for the plot
            if self.fig_topo == 0:
                print('self fig topo = 0')
                self.fig_topo = pyplot.figure()
            else:
                print('self fig topo =! 0')
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
                print( 'scale metric is checked' )
                # If the map is an Omicron one, I have to invert the y-axis
                if self.mapType == "Omicron":
                    self.ax_topo.axis([0, w, h, 0])
                    # pcolormesh takes *vertices* in arguments so the X (Y) array need to be from 0 to W (H) INCLUDED 
                    XYmap = self.ax_topo.pcolormesh(np.linspace(0, w, xPx + 1), np.linspace(0, h, yPx + 1), self.topo,
                                                    cmap=colormap)
                if self.mapType == "Sm4 to .mat" : 
                    XYmap = self.ax_topo.imshow(self.topo,extent=[0,w,0,h]);
                                #We plot our spec location points.
                    locations = True
                    if locations == True:
                        print('Spectrum Locations will be printed')
                        patch = self.m_params['Patch']
                        for m in range (0,len(patch)): #iterate over the number of locations
                            self.ax_topo.add_patch(patch[m]);

                else:
                    self.ax_topo.axis([0, w, 0, h])
                        # pcolormesh takes *vertices* in arguments so the X (Y) array need to be from 0 to W (H) INCLUDED 
                    XYmap = self.ax_topo.pcolormesh(np.linspace(0, w, xPx + 1), np.linspace(0, h, yPx + 1), self.topo,
                                                    cmap=colormap)
                
            else:
                # If the map is an Omicron one, I have to invert the y-axis
                if self.mapType == "Omicron":
                    self.ax_topo.axis([0, xPx, yPx, 0])
                    #if self.maptype == "Sm4" : ax.imshow(FImage,extent=[0,fxscale*1e+9,0,fyscale*1e+9], vmin = (imagerms - 5*imagestd) , vmax = (imagerms + 5*imagestd));#
                else:
                    print(self.maptype)
                    self.ax_topo.axis([0, xPx, 0, yPx])
                XYmap = self.ax_topo.pcolormesh(self.topo, cmap=colormap)
            # Colorbar stuff
            cbar = self.fig_topo.colorbar(XYmap, shrink=.9, pad=0.05, aspect=15)
            cbar.ax.yaxis.set_ticks_position('both')
            cbar.ax.tick_params(axis='y', direction='in')
            self.fig_topo.canvas.draw()
            self.fig_topo.show()

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
        self.fig_topo.show()

    def handleClosingTopo(self, event):
        """ Called when the topo figure is closed - Put back self.fig_topo to 0 to indicate that no topo figure exists """
        print('Topo closing')
        self.fig_topo = 0
        return

#%% Updating methods. Usually called by signals
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
        self.m_voltageBox.setMaximum(zPt-1)
        print(zPt-1)
        self.m_voltageBox.setSingleStep(1)

    def updateWidgets(self):
        """ Slot called after the reading of the CITS. Sets the values combo box (voltage and channels) and draws the map """
        self.updateVoltageBox(self.m_params["vStart"], self.m_params["vEnd"], self.m_params["zPt"])
        print('done0')
        self.m_channelBox.clear()
        self.m_channelBox.addItems(self.channelList)
        print('done1!')
        self.drawXYMap(self.m_voltageBox.value())
        print('done2')
        self.updateAboveValue(self.m_aboveBar.value())
        self.updateBelowValue(self.m_belowBar.value())
        self.m_CitsAvgBox.setMaximum(self.m_params["xPx"])
        print('done!')

#%% Methods related to spectra
    def normalizeSpectrum(self):
        """ Normalizes the spectra of displayed index """
        chan = self.m_channelBox.currentIndex()
        for x in range(self.m_params['xPx']):
            for y in range(self.m_params['yPx']):
                self.m_data[chan][x][y] = self.normalizeDOS(self.m_data[chan][x][y], self.m_params['zPt'])

    def normalizeDOS(self, dos, dos_length):
        mean_l = np.mean(dos[0:dos_length // 4])
        mean_r = np.mean(dos[3 * dos_length // 4:dos_length])
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

    def saveSpectra(self):
        """
        Saves the plotted spectra *as-is* as CSV data
        """
        output_path, fmt = QtWidgets.QFileDialog.getSaveFileName(self, "Choose a path to save the CSV file", self.wdir, "CSV (*.csv)")
        header = ''
        output = []
        for line in self.ax_spec.lines:
            # If first line, add the X-axis data
            if len(output) == 0:
                header += 'X-Axis,'
                output.append(line._x)
            output.append(line._y)
            # Replace ',' in the labels to avoid problems when exporting as CSV
            header += line._label.replace(',', '/') + ','
        # Crop last comma after looping on lines
        header = header[:-1]
        # Save if there is an output
        if len(output) != 0:
            # Tranposing the output array to get columns format corresponding to the header
            np.savetxt(output_path, np.transpose(np.array(output)), delimiter=',', header=header)
            print('Finished exporting csv at {}'.format(output_path))


#%% Methods related to the clicks on the map
    def onPressOnMap(self, event):
        """ Slot called when a press event is detected. Creates a Shape object that will be dynamically updated when the mouse moves """
        if (event.xdata is not None and event.ydata is not None) and (self.dataLoaded and self.toolbar_map._active is None):
            if self.m_scaleMetric.isChecked():
                ratioX = self.m_params['xL'] / self.m_params['xPx']
                ratioY = self.m_params['yL'] / self.m_params['yPx']
            else:
                ratioX = 1
                ratioY = 1
            self.currentShape = generateShape(event, self.m_mapWidget.figure, self.fig_topo,
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
                        self.currentShape = changeToDot(self.currentShape)
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
            fig.show()
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
            self.currentShape = generateShape(simEvent, self.m_mapWidget.figure, self.fig_topo,
                                      self.getSpectrumColor(self.nSpectraDrawn), ratioX, ratioY)
            simEvent.xdata = self.m_params['xPx'] - 1
            simEvent.ydata = self.m_params['yPx'] - 1
            self.currentShape.update(simEvent)
            X, Y, Z = self.cutAlongLine(0, self.m_params["xPx"] - 1, 0, self.m_params["yPx"] - 1)
            self.m_mapWidget.draw()
            return X, Y, Z

#%% Methods related to the map
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

#%% Methods related to the voltage guide line in the spectra window

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
            newData = np.zeros(shape=(nChannels+1, yPx, xPx, zPt))
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
        if self.dataLoaded == True :  # Check if a CITS was loaded
            # Get the needed params
            Navg = self.m_CitsAvgBox.value()
            xPx = self.m_params['xPx']
            yPx = self.m_params['yPx']
            zPt = self.m_params['zPt']
            # Try the allocation
            try:
                new_data = np.zeros(shape=(len(self.channelList), yPx, int(xPx / Navg), zPt))
                print (np.shape(new_data))
            except MemoryError:
                print('Not enough memory to do this operation')
                return
            # Average in X for each channel and Y
            for chan in range(0, len(self.channelList)):
                for y in range(0, yPx):
                    # To average, for each x, add every spectrum between x and x+Navg (averaging window) divided by Navg.
                    for x in range(0, xPx, int(Navg)):
                        if (
                                x + Navg > xPx):  # If the averaging window is not contained in the CITS data, stop averaging. The last spectras of this window will then will dismissed.
                            break
                        else:  # Else, average by adding the spectras of the avergaging window and dividing by Navg
                            spectra = np.zeros(zPt)
                            for i in range(0, int(Navg)):
                                spectra = spectra + self.m_data[chan][y][x + i] / Navg
                                # Store the result in new_data
                                new_data[chan][y][int(x / Navg)] = spectra
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
                    normData[y][x] = self.normalizeDOS(self.m_data[chan][y][x], 10)
            self.addChannel(normData, "Normalized " + self.channelList[chan])
