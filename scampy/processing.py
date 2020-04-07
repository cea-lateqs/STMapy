# -*- coding: utf-8 -*-
"""
Data processing functions
"""
import numpy as np
import scipy as sp
import scipy.interpolate
import scipy.optimize


def directionAverageCITS(cits_data, Navg, direction):
    """ Averages the CITS in one direction. """
    nChannels, yPx, xPx, zPt = cits_data.shape

    if direction == "y":
        new_data = np.zeros(shape=(nChannels, yPx // Navg, xPx, zPt))
        for y in range(0, yPx, Navg):
            new_data[:, y // Navg, :, :] = cits_data[:, y : y + Navg, :, :].mean(axis=1)
    else:
        new_data = np.zeros(shape=(nChannels, yPx, xPx // Navg, zPt))
        for x in range(0, xPx, Navg):
            new_data[:, :, x // Navg, :] = cits_data[:, :, x : x + Navg, :].mean(axis=2)

    return new_data


def computeAngle(self, Dmoire):
    """ Computes twist angle of a given graphene moiré of period Dmoire. """
    return 2 * np.arcsin(0.246 / (2 * Dmoire)) * 180 / np.pi


def normalizeDOS(dos, dos_length):
    mean_l = np.mean(dos[: dos_length // 4])
    mean_r = np.mean(dos[3 * dos_length // 4 :])
    mean = (mean_l + mean_r) / 2
    return dos / mean


def linearFitFunction(x, a, b):
    """ Basic function for linear fitting. """
    return a * x + b


def levelTopo(topo):
    yPx, xPx = topo.shape
    # Numpy array to save the leveled topo
    topo_leveled = np.zeros(shape=(yPx, xPx))
    fitX = np.arange(xPx)

    for y in range(yPx):
        fitY = topo[y]
        f = sp.interpolate.InterpolatedUnivariateSpline(fitX, fitY, k=1)
        popt, pcov = sp.optimize.curve_fit(linearFitFunction, fitX, f(fitX))
        topo_leveled[y] = fitY - (popt[0] * fitX + popt[1])
    # Return the leveled topo
    return topo_leveled


def extractSlope(topo, dataToFit, dZ, cutOffValue):
    """ Do a linear fit of the data and returns the slope and the coef found.
    Called for zSpectros.
    """
    yPx, xPx, zPt = dataToFit.shape
    zg = np.zeros(shape=(yPx, xPx, zPt))
    slopeData = np.zeros(shape=(yPx, xPx, zPt))
    coefData = np.zeros(shape=(yPx, xPx, zPt))
    xArray = np.arange(zPt) * dZ
    for y in range(yPx):
        for x in range(xPx):
            rawData = dataToFit[y, x]
            mask = rawData > cutOffValue
            xArrayFiltered = xArray[mask]
            dataFiltered = np.log(rawData[mask])
            popt, pcov = sp.optimize.curve_fit(
                linearFitFunction, xArrayFiltered, dataFiltered
            )
            zg[y, x] = 1 / 20.5 * np.log(rawData) + np.arange(zPt * dZ, dZ) + topo[y, x]
            slopeData[y, x, :] = popt[0]
            coefData[y, x, :] = popt[1]

    # return data to be added by "addchannel" protocol in citswidget
    return slopeData, coefData, zg


def findPixelsOnLine(xi, xf, yi, yf, use_bresenham=False):
    """ Finds the pixels on the line from (xi, yi) to (xf, yf) included.
    Returns integer arrays.
    """
    # First treat the vertical line case
    if xf == xi:
        y_plot = np.arange(min(yi, yf), max(yi, yf) + 1)
        x_plot = np.full(shape=y_plot.size, fill_value=xi)
    elif use_bresenham:  # Bresenham algo for non-vertical lines
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
            if x0 >= xf and y0 >= yf:
                break
            e2 = err * 2
            if e2 > -dy:
                err -= dy
                x0 += sx
            if e2 < dx:
                err += dx
                y0 += sy
        x_plot = np.array(x_plot_p)
        y_plot = np.array(y_plot_p)
    else:  # Simple algo for non-vertical lines
        # Determine its equation y=k*x+c
        k = float(yf - yi) / (xf - xi)
        c = yi - k * xi
        # Check if there is more y or more x to have to most precise arrangment
        if abs(xf - xi) > abs(yf - yi):
            x_plot = np.arange(min(xi, xf), max(xi, xf) + 1)
            y_plot = k * x_plot + c
        else:
            y_plot = np.arange(min(yi, yf), max(yi, yf) + 1)
            x_plot = (y_plot - c) / k
    return x_plot.astype(int), y_plot.astype(int)


def stringify(array):
    return "".join([chr(x) for x in array])
