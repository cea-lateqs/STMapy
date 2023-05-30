# -*- coding: utf-8 -*-
"""
Data processing functions
"""
import numpy as np
import scipy as sp
import scipy.interpolate
import scipy.optimize


def directionAverageCITS(cits_data, nb_to_avg, direction):
    """Averages the CITS in one direction."""
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


def normalizeDOS(dos, norm_length=None):
    """
    Normalize a DOS by dividing the dos by a mean of left-side DOS
    and right-side DOS of length norm_length along the last axis.
    Default: norm_length = dos.shape[-1] // 4.
    """
    if norm_length is None:
        norm_length = dos.shape[-1] // 4
    mean_l = np.mean(dos[..., :norm_length], axis=-1)
    mean_r = np.mean(dos[..., -norm_length:], axis=-1)
    mean = (mean_l + mean_r) / 2
    return dos / mean[..., np.newaxis]


def linearFitFunction(x, a, b):
    """Basic function for linear fitting."""
    return a * x + b


def planeFitFunction(yx, a, b, c):
    """Basic function for plane fitting."""
    y, x = yx
    return a * x + b * y + c


def levelTopoLine(topo):
    """Fit and substract a line background for each line of the topo."""
    yPx, xPx = topo.shape
    topo_leveled = np.zeros_like(topo)
    fit_x = np.arange(xPx)

    for y in range(yPx):
        fit_y = topo[y]
        f = sp.interpolate.InterpolatedUnivariateSpline(fit_x, fit_y, k=1)
        popt, pcov = sp.optimize.curve_fit(linearFitFunction, fit_x, f(fit_x))
        topo_leveled[y] = fit_y - linearFitFunction(fit_x, popt[0], popt[1])

    return topo_leveled


def levelTopoPlane(topo):
    """Fit and substract by a plane the topo."""
    Y, X = np.indices(topo.shape)
    fit_yx = np.vstack((Y.ravel(), X.ravel()))
    popt, pcov = sp.optimize.curve_fit(planeFitFunction, fit_yx, topo.ravel())
    return topo - planeFitFunction(fit_yx, *popt).reshape(topo.shape)


def extractSlope(topo, data_to_fit, delta_z, cut_off_value):
    """Do a linear fit of the data and returns the slope and the coef found.
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
    """Finds the pixels on the line from (xi, yi) to (xf, yf) included.
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


def FeenstraNormalization(data, params, channelList):
    """
    Performs conductance normalization as described in
    R. M. Feenstra et al., Surface Science 181, (1986)
    https://www.sciencedirect.com/science/article/pii/0039602887901701
    """
    dIdVindex = channelList.index("Data [Fwd]")
    dIdVFwd = data[dIdVindex, :, :, :]
    dIdVBwd = data[dIdVindex + 1, :, :, :]
    IVindex = channelList.index("IV_Data [Fwd]")
    IVFwd = abs(data[IVindex, :, :, :])  # +50*10**(-3) # unit : nA from matrix files
    IVBwd = abs(
        data[IVindex + 1, :, :, :]
    )  # +50*10**(-3) # unit : nA from matrix files
    # IVFwd = np.where(IVFwd<0.00001,1,IVFwd)
    # IVBwd = np.where(IVBwd<0.00001,1,IVBwd)
    # Vbias = np.array([np.log(abs(params["vStart"] + i * params["dV"]))
    #                     if not(abs(params["vStart"] + i * params["dV"])==0)
    #                     else 1 for i in range(np.shape(IVFwd)[-1])])
    Vbias = np.array(
        [abs(params["vStart"] + i * params["dV"]) for i in range(np.shape(IVFwd)[-1])]
    )
    Vbias = np.concatenate([Vbias] * np.shape(dIdVFwd)[-2], axis=0)
    Vbias = np.concatenate([Vbias] * np.shape(dIdVFwd)[-3], axis=0)
    Vbias = np.reshape(Vbias, np.shape(dIdVFwd))
    Norm_conductanceFwd = (
        dIdVFwd
        * params["LIsensi"]
        * Vbias
        / (10 * params["AmpMod"] * params["Gain"] * 10 ** (-6) * IVFwd)
    )
    Norm_conductanceBwd = (
        dIdVBwd
        * params["LIsensi"]
        * Vbias
        / (10 * params["AmpMod"] * params["Gain"] * 10 ** (-6) * IVBwd)
    )
    Norm_conductance = np.concatenate(([Norm_conductanceFwd], [Norm_conductanceBwd]))
    return Norm_conductance


def savgolderivative(
    data,
    window,
    polyorder,
    shape,
    dV,
    deriv=0,
):
    d = [
        sp.signal.savgol_filter(
            data[i, j, :],
            window_length=window,
            polyorder=polyorder,
            deriv=deriv,
            delta=dV,
        )
        for i in range(shape[0])
        for j in range(shape[1])
    ]
    return np.reshape(d, np.shape(data))


def diffderivative(data, shape):
    d = [
        np.hstack((data[i, j, 0], data[i, j, :], data[i, j, -1]))[2:]
        - np.hstack((data[i, j, 0], data[i, j, :], data[i, j, -1]))[:-2]
        for i in range(shape[0])
        for j in range(shape[1])
    ]
    # d = np.reshape(d, (np.shape(data)[0], np.shape(data)[1], np.shape(data)[2]-2))
    d = np.reshape(d, np.shape(data))
    return d


def FeenstraNormalization_log(data, params, channelList):
    """
    Performs conductance normalization as described in
    R. M. Feenstra et al., Surface Science 181, (1986)
    https://www.sciencedirect.com/science/article/pii/0039602887901701
    """
    savgolfilter = 1
    IVindex = channelList.index("IV_Data [Fwd]")
    IVFwd = abs(data[IVindex, :, :, :])
    IVBwd = abs(data[IVindex + 1, :, :, :])
    lnIVFwd = np.array(np.log(IVFwd + 5 * 10 ** (-3)))  # unit : nA from matrix files
    lnIVBwd = np.array(np.log(IVBwd + 5 * 10 ** (-3)))
    lnVbias = np.array(
        [
            np.log(abs(params["vStart"] + i * params["dV"]))
            if not (abs(params["vStart"] + i * params["dV"]) == 0)
            else 1
            for i in range(np.shape(lnIVFwd)[-1])
        ]
    )
    # Vbias = np.array([params["vStart"] + i * params["dV"] if not (("vStart"] + i * params["dV"]) == 0) else for i in range(len(IVFwd))])
    lnVbias = np.concatenate([lnVbias] * np.shape(lnIVFwd)[-2], axis=0)
    lnVbias = np.concatenate([lnVbias] * np.shape(lnIVFwd)[-3], axis=0)
    lnVbias = np.reshape(lnVbias, np.shape(lnIVFwd))
    if savgolfilter:
        dlnIVFwd = savgolderivative(lnIVFwd, 31, 2, np.shape(lnIVFwd[:, :, 0]), deriv=1)
        dlnIVBwd = savgolderivative(lnIVBwd, 31, 2, np.shape(lnIVFwd[:, :, 0]), deriv=1)
        dlnVbias = savgolderivative(lnVbias, 31, 2, np.shape(lnIVFwd[:, :, 0]), deriv=1)
    else:
        dlnIVFwd = diffderivative(lnIVFwd, np.shape(lnIVFwd[:, :, 0]))
        dlnIVBwd = diffderivative(lnIVBwd, np.shape(lnIVFwd[:, :, 0]))
        dlnVbias = diffderivative(lnVbias, np.shape(lnIVFwd[:, :, 0]))
    Norm_conductanceFwd = dlnIVFwd / dlnVbias
    Norm_conductanceBwd = dlnIVBwd / dlnVbias
    Norm_conductance = np.concatenate(([Norm_conductanceFwd], [Norm_conductanceBwd]))
    return Norm_conductance


def derivate_IV(IVdata, dV):
    """
    Calculates derivative of IVdata
    Using either savgol derivative or savgol filtering and by hand derivative
    """
    savgolfilter = 1
    if savgolfilter:
        dIVdata = -savgolderivative(
            IVdata, 11, 2, np.shape(IVdata[:, :, 0]), dV, deriv=1
        )
        # dIVdata = savgolderivative(IVdata, 81, 2, np.shape(IVdata[:,:,0]))
    else:
        dIVdata = (
            -diffderivative(
                savgolderivative(IVdata, 31, 2, np.shape(IVdata[:, :, 0])),
                np.shape(IVdata[:, :, 0]),
            )
            / dV
        )
    return dIVdata
