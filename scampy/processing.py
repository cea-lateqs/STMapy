# -*- coding: utf-8 -*-
"""
Data processing functions
"""
import numpy as np
import scipy as sp
import scipy.interpolate
import scipy.optimize


def directionAverageCITS(cits_data, nb_to_avg, direction):
    """ Averages the CITS in one direction. """
    nb_channels, yPx, xPx, zPt = cits_data.shape

    if direction == "y":
        new_data = np.zeros(shape=(nb_channels, yPx // nb_to_avg, xPx, zPt))
        for y in range(0, yPx, nb_to_avg):
            new_data[:, y // nb_to_avg, :, :] = cits_data[
                :, y : y + nb_to_avg, :, :
            ].mean(axis=1)
    else:
        new_data = np.zeros(shape=(nb_channels, yPx, xPx // nb_to_avg, zPt))
        for x in range(0, xPx, nb_to_avg):
            new_data[:, :, x // nb_to_avg, :] = cits_data[
                :, :, x : x + nb_to_avg, :
            ].mean(axis=2)

    return new_data


def normalizeDOS(dos, dos_length=None):
    if dos_length is None:
        dos_length = dos.shape[-1]
    mean_l = np.mean(dos[..., : dos_length // 4], axis=-1)
    mean_r = np.mean(dos[..., 3 * dos_length // 4 :], axis=-1)
    mean = (mean_l + mean_r) / 2
    return dos / mean[..., np.newaxis]


def linearFitFunction(x, a, b):
    """ Basic function for linear fitting. """
    return a * x + b


def levelTopo(topo):
    yPx, xPx = topo.shape
    # Numpy array to save the leveled topo
    topo_leveled = np.zeros_like(topo)
    fit_x = np.arange(xPx)

    # TODO: Can be vectorized ?
    for y in range(yPx):
        fit_y = topo[y]
        f = sp.interpolate.InterpolatedUnivariateSpline(fit_x, fit_y, k=1)
        popt, pcov = sp.optimize.curve_fit(linearFitFunction, fit_x, f(fit_y))
        topo_leveled[y] = fit_y - (popt[0] * fit_x + popt[1])
    # Return the leveled topo
    return topo_leveled


def extractSlope(topo, data_to_fit, delta_z, cut_off_value):
    """ Do a linear fit of the data and returns the slope and the coef found.
    Called for zSpectros.
    """
    yPx, xPx, zPt = data_to_fit.shape

    slope_data = np.zeros(shape=(yPx, xPx, zPt))
    coef_data = np.zeros(shape=(yPx, xPx, zPt))
    altitudes = np.arange(zPt) * delta_z

    cut_off_filter = data_to_fit > cut_off_value
    # TODO: Can be vectorized ?
    for y in range(yPx):
        for x in range(xPx):
            filtered_altitudes = altitudes[cut_off_filter[y][x]]
            filtered_data = data_to_fit[y, x, cut_off_filter[y][x]]
            popt, pcov = sp.optimize.curve_fit(
                linearFitFunction, filtered_altitudes, filtered_data
            )
            slope_data[y, x, :] = popt[0]
            coef_data[y, x, :] = popt[1]
    zg = 1 / 20.5 * np.log(data_to_fit) + altitudes + topo[y, x, np.newaxis]

    # return data to be added by "addchannel" protocol in citswidget
    return slope_data, coef_data, zg


def findPixelsOnLine(xi, xf, yi, yf, use_bresenham=False):
    """ Finds the pixels on the line from (xi, yi) to (xf, yf) included.
    Returns integer arrays.
    """
    # First treat the vertical line case
    if xf == xi:
        if yi < yf:
            y_plot = np.arange(yi, yf + 1)
        else:
            y_plot = np.arange(yf, yi + 1)[::-1]
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
            if xi < xf:
                x_plot = np.arange(xi, xf + 1)
            else:
                x_plot = np.arange(xf, xi + 1)[::-1]
            y_plot = k * x_plot + c
        else:
            if yi < yf:
                y_plot = np.arange(yi, yf + 1)
            else:
                y_plot = np.arange(yf, yi + 1)[::-1]
            x_plot = (y_plot - c) / k
    return x_plot.astype(int), y_plot.astype(int)


def stringify(array):
    return "".join([chr(x) for x in array])
