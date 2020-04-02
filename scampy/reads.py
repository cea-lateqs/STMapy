#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" Reading functions. """
import logging
import numpy as np
import struct
import os.path
import PyQt5.QtWidgets as QtWidgets
from scampy.processing import extractSlope, levelTopo, stringify


def readCitsAscii(filepath):
    """ Reads an Ascii CITS file (Omicron) and stores all the parameters"""
    with open(filepath, "r") as f:
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
            logging.error(
                "Problem while reading the file : could not find ':HEADER END:' in file"
            )
            return False
        # Matlab convention : columns first then rows hence [y][x]
        # In Omicron CITS, there is only two channels : fwd and bwd so it is read as such
        channelList = ["Data [Fwd]", "Data [Bwd]"]
        m_data = np.zeros(shape=(2, yPx, xPx, zPt))
        for y in range(yPx):
            for x in range(xPx):
                # The line read is an array containing the dI/dV (or I(V)) values indexed by the voltage index
                # Strip to remove the newline at the end and split to transform the string in a list
                data_list = f.readline().strip().split()
                # Forward data
                m_data[0][y][x] = np.float64(data_list[0:zPt]) * unit
                # No need to reverse the backward data as it was from Vmin to Vmax in the file as the fwd data
                # Backward data
                m_data[1][y][x] = np.float64(data_list[zPt : 2 * zPt]) * unit
    # Store the parameters in a dictonnary to use them later
    m_params = {
        "xPx": xPx,
        "yPx": yPx,
        "xL": xL,
        "yL": yL,
        "zPt": zPt,
        "vStart": vStart / divider,
        "vEnd": vEnd / divider,
        "dV": abs(vEnd - vStart) / (divider * zPt),
    }
    if divider != 1:
        logging.info("A divider of " + str(divider) + " was found and applied")

    # Check if a topo file exists and read it if yes
    topopath = os.path.join(os.path.dirname(filepath), "Topo.txt")
    topo = readTopo(topopath) if os.path.exists(topopath) else None
    return topo, m_data, channelList, m_params


def readTopo(filepath):
    """ Reads a topography file (in test). Used for txt files. """
    with open(filepath, "r") as f:
        topo_data = []
        width, height = None, None
        for line in f:
            # Treat the headers differently
            if line[0] == "#":
                if "Width" in line:
                    width = float(line.split()[-2])
                if "Height" in line:
                    height = float(line.split()[-2])
            else:
                topo_data.append(line.strip().split())
    if width is None or height is None:
        raise IOError("No width/height read in topo file {} !".format(filepath))
    return np.asfarray(topo_data)


def readCits3dsBin(filepath, zSpectro):
    # The divider is already taken into account by Nanonis during the experiment so no need to process it again*
    half = False
    # Read the header of the map until its end ("HEADER_END")
    header_end_not_found = True
    with open(filepath, "rb") as f:
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
                if line.split('"')[1] == "Z (m)":
                    zSpectro = True
            # Number of points per channel
            elif "Points" in line:
                zPt = int(line.split("=")[-1])
            # Channels recorded
            elif "Channels" in line:
                channelList = line.split('"')[1].split(";")
                nChannels = len(channelList)
            # Experiment parameters. Not used for now, only the number is recorded to skip the corresponding bytes afterwards
            elif "Experiment parameters" in line:
                nbExpParams = len(line.split(";"))
            # End of the header
            elif ":HEADER_END:" in line:
                header_end_not_found = False
                break

        if header_end_not_found:
            logging.error(
                "Problem while reading the file : could not find ':HEADER END:' in file"
            )
            f.close()
            return False
        # Reading vStart and vEnd (floats of 4 bytes each)
        try:
            reading = struct.unpack(">" + "f" * 2, f.read(8))
        except struct.error:
            logging.error(
                "Problem while reading the file : number of bytes to read different than what was expected"
            )
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
            topo = np.zeros(shape=(yPx, xPx))
            m_data = np.zeros(shape=(nChannels, yPx, xPx, zPt))
        except MemoryError:
            logging.warning(
                "The data is too big ! Or the memory too small...\nI will take half of the channels...\n"
            )
            half = True
        # If the first alloc didn't work, try to halve it
        if half:
            try:
                topo = np.zeros(shape=(yPx, xPx))
                m_data = np.zeros(shape=(nChannels / 2, yPx, xPx, zPt))
            except MemoryError:
                logging.error(
                    "The data is REALLY too big ! Or the memory REALLY too small...\nI give up...\n"
                )
                f.close()
                QtWidgets.QMessageBox.critical(
                    "Oops !",
                    "The data is REALLY too big ! Or the memory REALLY too small...\nI give up...\n",
                )
                return False
        # Format string for unpacking zPt big-endians floats ('>f')
        fmtString = ">" + "f" * zPt
        # zPt floats to read of 4 bytes each
        bytesToRead = 4 * zPt
        for y in range(yPx):
            for x in range(xPx):
                chan = 0
                b = f.read(nbExpParams * 4)
                try:
                    topo[y][x] = struct.unpack(">" + "f" * nbExpParams, b)[2]
                except struct.error:
                    logging.error(
                        "Problem while reading the topo : number of bytes to read different than what was expected at "
                        + str(x)
                        + " "
                        + str(y)
                    )
                while chan < nChannels:
                    # Each channel is written successively by sequences of 4*zPt bytes. I then read these bytes and unpack them as big-endians floats ('>f')
                    b = f.read(bytesToRead)
                    try:
                        if not half or chan < nChannels / 2:
                            m_data[chan][y][x] = struct.unpack(fmtString, b)
                    except struct.error:
                        logging.error(
                            "Problem while reading the file : number of bytes to read different than what was expected at"
                            + str(x)
                            + " "
                            + str(y)
                            + " "
                            + str(chan)
                        )
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
        channelList = channelList[0 : nChannels // 2]
    # Store the parameters in a dictonnary to use them later
    dV = (vEnd - vStart) / zPt
    m_params = {
        "xPx": xPx,
        "yPx": yPx,
        "xC": xC,
        "yC": yC,
        "xL": xL,
        "yL": yL,
        "zPt": zPt,
        "vStart": vStart,
        "vEnd": vEnd,
        "dV": dV,
    }

    # self.m_statusBar.showMessage(self.m_params)
    # Convert currents in nA
    for i in range(len(channelList)):
        chan = channelList[i]
        if "(A)" in chan:
            m_data[i] = np.abs(m_data[i]) * 10 ** 9
            channelList[i] = chan.replace("(A)", "(nA)")
    # Convert topo in nm
    topo = topo * 10 ** 9
    # Level topo
    topo = levelTopo(topo)
    # Test
    if zSpectro:
        slopeData, slopeDataName, coefData, coefDataName, zg = extractSlope(
            topo, m_data, m_params, channelList, 0.01, 0
        )
        # return data and data to be added by addchannel
        return (
            topo,
            m_data,
            channelList,
            m_params,
            slopeData,
            slopeDataName,
            coefData,
            coefDataName,
            zg,
        )

    else:
        return topo, m_data, channelList, m_params


def sm4readFileHeader(f, ObjectIDCode):

    Header_size = int(np.fromfile(f, dtype=np.uint16, count=1)[0])
    Signature = stringify(np.fromfile(f, dtype=np.uint16, count=18))
    Total_Pagecount = int(np.fromfile(f, dtype=np.uint32, count=1)[0])
    ObjectListCount = int(np.fromfile(f, dtype=np.uint32, count=1)[0])
    ObjectFieldSize = int(np.fromfile(f, dtype=np.uint32, count=1)[0])
    Reserved = int(np.fromfile(f, dtype=np.uint32, count=2)[0])

    # iterate over the known objects from file header
    #        ObjectlistName = []
    #        ObjectlistOffset = []
    #        ObjectlistSize = []
    #        for i in range(ObjectListCount):
    #            ObjectlistName.append(ObjectIDCode[np.fromfile(f, dtype=np.uint32, count=1)[0]])
    #            ObjectlistOffset.append(np.fromfile(f, dtype=np.uint32, count=1))
    #            ObjectlistSize.append(np.fromfile(f, dtype=np.uint32, count=1))
    f.seek(ObjectListCount * 3 * 4, 1)

    # Read and record the Page Index Header
    PageIndexHeader_PageCount = int(np.fromfile(f, dtype=np.uint32, count=1)[0])
    # the number of pages in the Page Index Array
    PageIndexHeader_ObjectListCount = int(np.fromfile(f, dtype=np.uint32, count=1)[0])
    # count of objects stored after Page Index Header
    # currently there is just one: Page Index Array
    PageIndexHeader_Reserved = int(np.fromfile(f, dtype=np.uint32, count=2)[0])

    # Read and record the Page Index Array
    PageIndexHeader_ObjectID = int(np.fromfile(f, dtype=np.uint32, count=1)[0])
    PageIndexHeader_Offset = int(np.fromfile(f, dtype=np.uint32, count=1)[0])
    PageIndexHeader_Size = int(np.fromfile(f, dtype=np.uint32, count=1)[0])

    # Get info on each page
    PageIndex = []
    for pageindex in range(PageIndexHeader_PageCount):
        PageIndex.append(
            {
                "PageID": np.fromfile(f, dtype=np.uint16, count=8)[0],
                "PageDataType": int(np.fromfile(f, dtype=np.uint32, count=1)[0]),
                "PageSourceType": int(np.fromfile(f, dtype=np.uint32, count=1)[0]),
                "ObjectListCount": int(np.fromfile(f, dtype=np.uint32, count=1)[0]),
                "MinorVersion": int(np.fromfile(f, dtype=np.uint32, count=1)[0]),
                "ObjectList": [],
            }
        )
        for i in range(PageIndex[pageindex]["ObjectListCount"]):
            dict1 = {
                "ObjectID": ObjectIDCode[
                    int(np.fromfile(f, dtype=np.uint32, count=1)[0])
                ],
                "Offset": int(np.fromfile(f, dtype=np.uint32, count=1)[0]),
                "Size": int(np.fromfile(f, dtype=np.uint32, count=1)[0]),
            }
            PageIndex[pageindex]["ObjectList"].append(dict1)

    return PageIndexHeader_PageCount, PageIndex


def sm4readpageHeader(
    f, PageHeader, TextStrings, PageIndex, PageIndexHeader_PageCount, j
):
    # f.seek(nbytes,0 = a partir du début) va au byte n + 1
    f.seek(PageIndex[j]["ObjectList"][0]["Offset"], 0)

    # Read Page Header
    PageHeader.append(
        {
            "FieldSize": int(np.fromfile(f, dtype=np.uint16, count=1)[0]),
            "StringCount": int(np.fromfile(f, dtype=np.uint16, count=1)[0]),
            "PageType_DataSource": int(np.fromfile(f, dtype=np.uint32, count=1)[0]),
            "DataSubSource": int(np.fromfile(f, dtype=np.uint32, count=1)[0]),
            "LineType": int(np.fromfile(f, dtype=np.uint32, count=1)[0]),
            "XCorner": int(np.fromfile(f, dtype=np.uint32, count=1)[0]),
            "YCorner": int(np.fromfile(f, dtype=np.uint32, count=1)[0]),
            "Width": int(np.fromfile(f, dtype=np.uint32, count=1)[0]),
            "Height": int(np.fromfile(f, dtype=np.uint32, count=1)[0]),
            "ImageType": int(np.fromfile(f, dtype=np.uint32, count=1)[0]),
            "ScanDirection": int(np.fromfile(f, dtype=np.uint32, count=1)[0]),
            "GroupID": int(np.fromfile(f, dtype=np.uint32, count=1)[0]),
            "PageDataSize": int(np.fromfile(f, dtype=np.uint32, count=1)[0]),
            "MinZValue": int(np.fromfile(f, dtype=np.uint32, count=1)[0]),
            "MaxZValue": int(np.fromfile(f, dtype=np.uint32, count=1)[0]),
            "XScale": np.float(
                np.fromfile(f, dtype=np.float32, count=1)[0]
            ),  # single = 4 bytes Floating-point numbers
            "YScale": np.float(np.fromfile(f, dtype=np.float32, count=1)[0]),
            "ZScale": np.float(np.fromfile(f, dtype=np.float32, count=1)[0]),
            "XYScale": np.float(np.fromfile(f, dtype=np.float32, count=1)[0]),
            "XOffset": np.float(np.fromfile(f, dtype=np.float32, count=1)[0]),
            "YOffset": np.float(np.fromfile(f, dtype=np.float32, count=1)[0]),
            "ZOffset": np.float(np.fromfile(f, dtype=np.float32, count=1)[0]),
            "Period": np.float(np.fromfile(f, dtype=np.float32, count=1)[0]),
            "Bias": np.float(np.fromfile(f, dtype=np.float32, count=1)[0]),
            "Current": np.float(np.fromfile(f, dtype=np.float32, count=1)[0]),
            "Angle": np.float(np.fromfile(f, dtype=np.float32, count=1)[0]),
            "ColorInfoListCount": int(np.fromfile(f, dtype=np.uint32, count=1)[0]),
            "GridXSize": int(np.fromfile(f, dtype=np.uint32, count=1)[0]),
            "GridYSize": int(np.fromfile(f, dtype=np.uint32, count=1)[0]),
            "ObjectListCount": int(np.fromfile(f, dtype=np.uint32, count=1)[0]),
            "ObjectList": [],
        }
    )
    # Skip flags and reserved data
    f.seek(1 + 3 + 60, 1)
    # Read the Object List
    for iloop in range(PageHeader[j]["ObjectListCount"]):
        dict1 = {
            "ObjectID": int(np.fromfile(f, dtype=np.uint32, count=1)[0]),
            "Offset": int(np.fromfile(f, dtype=np.uint32, count=1)[0]),
            "Size": int(np.fromfile(f, dtype=np.uint32, count=1)[0]),
        }
        PageHeader[j]["ObjectList"].append(dict1)

    # read text describing file
    c = int(np.fromfile(f, dtype=np.uint16, count=1)[0])
    d = {"strLabel": stringify(np.fromfile(f, dtype=np.uint16, count=c))}
    count = int(np.fromfile(f, dtype=np.uint16, count=1)[0])
    d["strSystemText"] = stringify(np.fromfile(f, dtype=np.uint16, count=count))

    count = int(np.fromfile(f, dtype=np.uint16, count=1)[0])
    d["strSessionText"] = stringify(np.fromfile(f, dtype=np.uint16, count=count))

    count = int(np.fromfile(f, dtype=np.uint16, count=1)[0])
    d["strUserText"] = stringify(np.fromfile(f, dtype=np.uint16, count=count))

    count = int(np.fromfile(f, dtype=np.uint16, count=1)[0])
    d["strPath"] = stringify(np.fromfile(f, dtype=np.uint16, count=count))

    count = int(np.fromfile(f, dtype=np.uint16, count=1)[0])
    d["strDate"] = stringify(np.fromfile(f, dtype=np.uint16, count=count))

    count = int(np.fromfile(f, dtype=np.uint16, count=1)[0])  # DAQ time
    d["strTime"] = stringify(np.fromfile(f, dtype=np.uint16, count=count))

    count = int(np.fromfile(f, dtype=np.uint16, count=1)[0])  # physical units of x axis
    d["strXUnits"] = stringify(np.fromfile(f, dtype=np.uint16, count=count))

    count = int(np.fromfile(f, dtype=np.uint16, count=1)[0])
    d["strYUnits"] = stringify(np.fromfile(f, dtype=np.uint16, count=count))

    count = int(np.fromfile(f, dtype=np.uint16, count=1)[0])
    d["strZUnits"] = stringify(np.fromfile(f, dtype=np.uint16, count=count))

    count = int(np.fromfile(f, dtype=np.uint16, count=1)[0])
    d["strYLabel"] = stringify(np.fromfile(f, dtype=np.uint16, count=count))

    count = int(np.fromfile(f, dtype=np.uint16, count=1)[0])
    d["strStatusChannelText"] = stringify(np.fromfile(f, dtype=np.uint16, count=count))

    count = int(np.fromfile(f, dtype=np.uint16, count=1)[0])
    # contains last saved line count for an image data page
    d["strCompletedLineCount"] = stringify(np.fromfile(f, dtype=np.uint16, count=count))

    count = int(np.fromfile(f, dtype=np.uint16, count=1)[0])
    # Oversampling count for image data pages
    d["strOverSamplingCount"] = stringify(np.fromfile(f, dtype=np.uint16, count=count))

    count = int(np.fromfile(f, dtype=np.uint16, count=1)[0])
    # voltage at which the sliced image is created from the spectra page.
    # empty if not a sliced image
    d["strSlicedVoltage"] = stringify(np.fromfile(f, dtype=np.uint16, count=count))

    count = int(np.fromfile(f, dtype=np.uint16, count=1)[0])
    # PLLPro status text: blank, master or user
    d["strPLLProStatus"] = stringify(np.fromfile(f, dtype=np.uint16, count=count))

    count = int(np.fromfile(f, dtype=np.uint16, count=1)[0])
    # ZPI controller item's set-point unit
    d["strSetpointUnit"] = stringify(np.fromfile(f, dtype=np.uint16, count=count))

    count = int(np.fromfile(f, dtype=np.uint16, count=1)[0])
    # stores value of CH1 and CH2 if they are in hardware space
    d["strCHDriveValues"] = stringify(np.fromfile(f, dtype=np.uint16, count=count))

    TextStrings.append(d)

    return PageHeader, TextStrings


def readCitsSm4Bin(filepath):
    ObjectIDCode = [
        "Undefined",
        "Page Index Header",
        "Page Index Array",
        "Page Header",
        "Page Data",
        "Image Drift Header",
        "Image Drift",
        "Spec Drift Header",
        "Spec Drift Data (with X,Y coordinates)",
        "Color Info",
        "String data",
        "Tip Track Header",
        "Tip Track Data",
        "PRM",
        "Thumbnail",
        "PRM Header",
        "Thumbnail Header",
        "API Info",
        "History Info",
        "Piezo Sensitivity",
        "Frequency Sweep Data",
        "Scan Processor Info",
        "PLL Info",
        "CH1 Drive Info",
        "CH2 Drive Info",
        "Lockin0 Info",
        "Lockin1 Info",
        "ZPI Info",
        "KPI Info",
        "Aux PI Info",
        "Low-pass Filter0 Info",
        "Low-pass Filter1 Info",
    ]
    with open(filepath, "rb") as f:
        # f = open(filepath, "rb")

        # Get main info from File Header
        PageIndexHeader_PageCount, PageIndex = sm4readFileHeader(f, ObjectIDCode)

        # Read and record each Pages Headers and Data
        PageHeader = []
        TextStrings = []
        Spatial = []
        Spectral = []
        SpatialInfo = []
        topocount = 0
        currentcount = 0
        dIdV_Line_Speccount = 0

        for j in range(PageIndexHeader_PageCount):

            PageHeader, TextStrings = sm4readpageHeader(
                f, PageHeader, TextStrings, PageIndex, PageIndexHeader_PageCount, j
            )

            # Page Header – Sequential Data Page - Get Datas
            Data = []
            f.seek(PageIndex[j]["ObjectList"][1]["Offset"], 0)
            Data.append(
                np.fromfile(
                    f,
                    dtype=np.int32,
                    count=int(round(PageIndex[j]["ObjectList"][1]["Size"] / 4)),
                )
            )
            # /4 because total data size is divided by the number of bytes that use each 'long' data
            #!!! This takes too much time ?
            ScaleData = [
                x * PageHeader[j]["ZScale"] + PageHeader[j]["ZOffset"] for x in Data[-1]
            ]
            ScaleData = np.reshape(
                ScaleData, (PageHeader[j]["Width"], PageHeader[j]["Height"]), order="F"
            )
            # order Fortran = "F" to match readCITSsm4File function

            ###################### Spatial Data
            # is it label topo ? # cf p.28 of SM4 DATA FILE FORMAT V5.pdf
            logging.debug("Page {0} contains {1}".format(j, TextStrings[j]["strLabel"]))
            if (
                TextStrings[j]["strLabel"] == "Topography"
                and PageHeader[j]["PageType_DataSource"] == 1
            ):
                topocount += 1
                Spatialpagenumber = j
                # records last spatial page; to read the rest of the spatial data later
                Spatial.append({"TopoData": np.rot90(ScaleData, 3)})
                SpatialInfo.append({"TopoUnit": TextStrings[j]["strZUnits"]})
                if topocount == 1:
                    FImage = Spatial[-1]["TopoData"]  # Image forward
                elif topocount == 2:
                    BImage = Spatial[-1]["TopoData"]  # Image bacward
                else:
                    logging.warning("There is more topo data than expected")
            # is it spatial Current data ?
            elif (
                TextStrings[j]["strLabel"] == "Current"
                and PageHeader[j]["PageType_DataSource"] == 2
            ):
                currentcount += 1
                Spatial.append({"CurrentData": np.rot90(ScaleData, 3)})
                Spatial[-1]["CurrentUnit"] = TextStrings[j]["strZUnits"]
                Spatialpagenumber = j
                if currentcount == 1:
                    FImage_I = Spatial[-1]["CurrentData"]  # Image forward
                elif currentcount == 2:
                    BImage_I = Spatial[-1]["CurrentData"]  # Image bacward
                else:
                    logging.warning("There is more current data than expected")

            ###################### Spectral Data - can be Point or Line
            # Is it Spectral Point(38) ? Not taken in charge currently
            elif PageHeader[j]["PageType_DataSource"] == 38:
                logging.error("You didnt load a CITS. Use sm4_reader to read this data")

            # CITS is recarded as a line of LIA
            # Is it Spectral Line(16) LIA ?
            elif (
                TextStrings[j]["strLabel"] == "LIA"
                and PageHeader[j]["PageType_DataSource"] == 16
            ):
                dIdV_Line_Speccount += 1
                Spectral.append({"dIdV_Line_Data": ScaleData})
                Linespectrapagenumber = j
                if dIdV_Line_Speccount == 1:
                    SpectralData_y = Spectral[-1]["dIdV_Line_Data"]
                    # Get spectra coordinates offset
                    for objectNbre in range(PageHeader[j]["ObjectListCount"]):
                        # 8 == Tip track Data
                        if PageHeader[j]["ObjectList"][objectNbre]["ObjectID"] == 8:
                            TipTrackData_offset = PageHeader[j]["ObjectList"][
                                objectNbre
                            ]["Offset"]
                            TipTrackData_Size = PageHeader[j]["ObjectList"][objectNbre][
                                "Size"
                            ]
                else:
                    logging.warning("There is more spectral data than expected")

            # We could add 'Spaectral line Current', or 'PLL Amplitude' ;
            # 'PLL Phase' ; 'dF' ; 'PLL Drive'spectra, AFM, ... in the same way
            else:
                logging.warning(
                    "Data Type "
                    + str(TextStrings[j]["strLabel"])
                    + " or PageType_DataSource "
                    + str(PageHeader[j]["PageType_DataSource"])
                    + " not found. Check spelling or add it."
                )

        ###################### Get spectra coordinates
        # Nbre of measurements taken along a line :
        nbreScans = np.shape(SpectralData_y)[1]
        xCoord = np.zeros(nbreScans)
        yCoord = np.zeros(nbreScans)
        f.seek(TipTrackData_offset, 0)  # Go to beggining of header
        for i in range(nbreScans):
            f.seek(4, 1)  # skip start time
            xCoord[i] = np.float(np.fromfile(f, dtype=np.float32, count=1)[0])
            yCoord[i] = np.float(np.fromfile(f, dtype=np.float32, count=1)[0])
            f.seek(16, 1)  # skip dx dy xcumul ycumul fields

    SpectralData_x = (
        PageHeader[Linespectrapagenumber]["XOffset"]
        + PageHeader[Linespectrapagenumber]["XScale"]
        * np.array(list(range(PageHeader[Linespectrapagenumber]["Width"])))
    ) * 1000.0  # mV

    # each measurement is taken at a coordinate,
    # but several measurements (repetitions) are taken on the same spot.
    # Eg if repetitions = 4, 2, one forward, one backward (saved in right direction)
    repetitions = 0
    pointdiff = 0
    repindex = 0
    # step through data points until there is a difference between the
    # current (x,y) point and the next (x1,y1) point
    while pointdiff == 0:
        pointdiff = (xCoord[repindex + 1] - xCoord[repindex]) + (
            xCoord[repindex + 1] - xCoord[repindex]
        )
        # adds the x and y difference values together.  if this is anything other than 0, we have found the limit of the repetitions
        repetitions = repetitions + 1
        repindex = repindex + 1
    numberOfPlots = nbreScans / repetitions  # number of distinct plotting locations
    # Center coordinates and metric dimensions in nm
    xL = abs(
        PageHeader[Spatialpagenumber]["XScale"] * PageHeader[Spatialpagenumber]["Width"]
    ) * (10 ** 9)
    yL = abs(
        PageHeader[Spatialpagenumber]["YScale"]
        * PageHeader[Spatialpagenumber]["Height"]
    ) * (10 ** 9)

    # size spatial data
    xPx = int(PageHeader[Spatialpagenumber]["Height"])
    yPx = int(PageHeader[Spatialpagenumber]["Width"])
    # size spectral data ( not necessarily the same in RHK )
    xSpec = int(np.sqrt(numberOfPlots))
    ySpec = int(np.sqrt(numberOfPlots))
    zPt = np.shape(SpectralData_y)[0]  # nbre of points in spec data

    try:
        topo = np.zeros(shape=(xPx, yPx))
        m_data = np.zeros(shape=(repetitions, ySpec, xSpec, zPt))
    except MemoryError:
        logging.error("The data is too big ! Or the memory too small...")
        return False

    # in Spectraldata_y, dIdV info corresponds to the spec data saved
    # from spatial left to right and increasing y(downwards in RHK),
    # with same nbre of repetions at each spot
    for r in range(int(repetitions)):  # even : forward, odd : backwards
        for y in range(ySpec):
            for x in range(xSpec):
                m_data[r][y][x] = SpectralData_y[:, (xSpec * y + x) * repetitions + r]

    channelList = ["Data {}".format(i) for i in range(int(repetitions))]

    m_params = {
        "xPx": xSpec,
        "yPx": ySpec,
        "xL": xL,
        "yL": yL,
        "zPt": zPt,
        "vStart": SpectralData_x[0],
        "vEnd": SpectralData_x[-1],
        "dV": abs(SpectralData_x[-1] - SpectralData_x[0]) / zPt,
    }

    topo = FImage

    # create average Data :
    average = np.zeros(shape=(ySpec, xSpec, zPt))
    for y in range(ySpec):  # len(SpectralData_y[:,0])):
        for x in range(xSpec):
            for r in range(repetitions):
                average[y][x] += (
                    SpectralData_y[:, (xSpec * y + x) * repetitions + r] / repetitions
                )

    return topo, m_data, channelList, m_params, average
