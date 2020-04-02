#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Data processing functions
"""
import numpy as np
import scipy as sp
import scipy.interpolate
import scipy.optimize


def linear_fit_function(x, a, b):
    return a * x + b


def levelTopo(topo):
    yPx, xPx = topo.shape
    # Numpy array to save the leveled topo
    topo_leveled = np.zeros(shape=(yPx, xPx))
    fitX = np.arange(0, xPx)

    def fitF(z, a, b):
        return a * z + b

    for y in range(0, yPx):
        fitY = topo[y]
        f = sp.interpolate.InterpolatedUnivariateSpline(fitX, fitY, k=1)
        popt, pcov = sp.optimize.curve_fit(linear_fit_function, fitX, f(fitX))
        topo_leveled[y] = fitY - (popt[0] * fitX + popt[1])
    # Return the leveled topo
    return topo_leveled


def extractSlope(topo, m_data, m_params, channelList, cutOffValue, numChanToFit):
    """ Do a linear fit of the data in the asked channel and add the slope and the coef found as channels (usually called for zSpectros) """
    yPx = m_params["yPx"]
    xPx = m_params["xPx"]
    zPt = m_params["zPt"]
    dZ = m_params["dV"]
    zg = np.zeros(shape=(yPx, xPx, zPt))
    slopeData = np.zeros(shape=(yPx, xPx, zPt))
    coefData = np.zeros(shape=(yPx, xPx, zPt))
    xArray = np.arange(zPt) * dZ
    for y in range(yPx):
        for x in range(xPx):
            rawData = m_data[numChanToFit][y][x]
            mask = rawData > cutOffValue
            xArrayFiltered = xArray[mask]
            dataFiltered = np.log(rawData[mask])
            popt, pcov = sp.optimize.curve_fit(
                linear_fit_function, xArrayFiltered, dataFiltered
            )
            zg[y][x] = 1 / 20.5 * np.log(rawData) + np.arange(zPt * dZ, dZ) + topo[y][x]
            for z in range(zPt):
                slopeData[y][x][z] = popt[0]
                coefData[y][x][z] = popt[1]
    # Add the created channel to the data
    slopeDataName = "Slope by linear fit of " + channelList[numChanToFit]
    coefDataName = "Coef by linear fit of " + channelList[numChanToFit]

    # return data to be added by "addchannel" protocol in citswidget
    return slopeData, slopeDataName, coefData, coefDataName, zg


def stringify(array):
    return "".join([chr(x) for x in array])
