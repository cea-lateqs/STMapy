#-*- coding: utf-8 -*-
from matplotlib import pyplot
import scipy.interpolate
import numpy as np


# Complementary functions
def plotMapDistrib(data1d, valMax):
    """ ??? """
    dx = valMax / 100
    distrib = np.zeros(shape=101)
    for val in data1d:
        index = int(val // dx)
        distrib[index] = distrib[index] + 1
    pyplot.figure()
    pyplot.plot(distrib)


def getDataFromTxtFile(fileName, header=False):
    data = np.genfromtxt(fileName, delimiter="\t", unpack=True)
    return data


def diffSingularite(voltage):
    dataS1 = getDataFromTxtFile(
        "E:\\PhD\\Experiments\\STM\\Lc12\\Depouillements\\Depouillement_LS_1B\\" + voltage + "\\Singularite1.csv")
    f = scipy.interpolate.InterpolatedUnivariateSpline(dataS1[0], dataS1[1], k=1)
    dataS2 = getDataFromTxtFile(
        "E:\\PhD\\Experiments\\STM\\Lc12\\Depouillements\\Depouillement_LS_1B\\" + voltage + "\\Singularite2.csv")
    g = scipy.interpolate.InterpolatedUnivariateSpline(dataS2[0], dataS2[1], k=1)
    x = np.linspace(0, 40, 500)
    print("s")
    pyplot.plot(x, g(x) - f(x), label=voltage)
    return


def angleFromSingulariteDiff(dE, adv=False):
    """ Returns the angle from the difference of energy of van Hove singularities """
    if adv:
        sol1 = 2 * np.arcsin(8 * 10 ** (-3) * (dE / 0.108 + 2 + np.sqrt((2 + dE / 0.108) ** 2 + 36)))
        sol2 = 2 * np.arcsin(8 * 10 ** (-3) * (dE / 0.108 + 2 - np.sqrt((2 + dE / 0.108) ** 2 + 36)))
        return 180 * sol1 / np.pi, 180 * sol2 / np.pi
    else:
        sol = 2 * np.arcsin(16 * 10 ** (-3) * (dE / 0.108 + 2))
        return 180 * sol / np.pi