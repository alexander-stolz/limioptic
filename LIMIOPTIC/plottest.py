"""
Aus der Erinnerung:

X: Position in mm vom BPM Draht
Y: Strom in muA

Am FN habe ich mit dem BPM Strahldicken bestimmt und mit den Simion- und Limioptic-Simulationen verglichen.
Es war mir damit möglich die FNVB-Fokussierung zu bestimmen (s. Bestimmung der FNVB) und die FN-Quelle in Limioptic einzubauen.
"""


#!/usr/bin/env python

from pylab import show, grid, figure
import sys

print "plottest >> ", sys.argv

messung1X = [42.6, 44, 46, 48, 50, 52, 50, 48, 46, 44, 42, 47, 49, 54, 40]
messung1Y = [4.425, 3.345, 2.22, 1.95, 3.49, 6.285, 3.325, 1.99, 2.22, 3.55, 5.4, 1.895, 2.17, 8.25, 7.52]

messung2X = [60, 58, 56, 52, 59]
messung2Y = [2.63, 2.945, 3.925, 7.05, 2.735]

messung3X = [58, 56, 54, 52, 50, 48]
messung3Y = [3.65, 2.5, 2.25, 3.21, 4.7, 6.05]

messung4X = [0, 2, 4, 6, 8, 10, 12]
messung4Y = [1.92, 2.07, 2.05, 2.23, 2.94, 4.27, 6.645]

messung5X = [0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 9, 11]
messung5Y = [5.3, 5.0, 4.5, 4.05, 3.6, 2.75, 3.65, 4.85, 6.85, 7.6, 7.45, 3.45, 3.1]

lithium1X = [38.0, 38.0, 40.0, 40.0, 42.0, 42.0, 44.0, 44.0, 46.0, 46.0, 48.0, 48.0, 50.0, 50.0, 53.0, 53.0, 56.0, 56.0]
lithium1Y = [2.29, 2.42, 1.93, 1.80, 1.41, 1.44, 1.49, 1.50, 3.32, 3.20, 1.39, 1.46, 1.34, 1.27, 2.55, 2.17, 3.82, 3.78]

lithium2X = [44.4, 44.4, 42.5, 34.5, 34.5, 34.5, 34.5, 38.5, 38.5, 40.7, 40.7, 37.5, 37.5, 39.5, 39.5, 39.0, 39.0, 38.0, 38.0]
lithium2Y = [7.43, 8.06, 3.86, 2.32, 2.51, 2.14, 2.14, 1.97, 1.95, 2.62, 2.62, 2.04, 1.99, 2.22, 2.14, 2.03, 1.99, 1.96, 1.91]

messungC1X = [45.4, 46.5, 48, 50, 52, 54, 45, 43]
messungC1Y = [3.75, 3.4, 2.85, 3.85, 6.3, 9.75, 4.75, 5.9]
messungC1Y = [_ / 1.9 * 2.7 for _ in messungC1Y]

messungC2X = [65, 63, 61, 59, 57, 55, 53, 51, 49]
messungC2Y = [8.7, 9.75, 9, 6.05, 3.2, 1.7, 2.1, 2.6, 2.6]
messungC2Y = [_ / 1.9 * 2.7 for _ in messungC2Y]

messungC3X = [45, 46, 47, 43, 41, 40, 39, 38, 37, 36, 35]
messungC3Y = [10.6, 9.45, 8.65, 7.1, 3.45, 2.4, 2.25, 2.65, 3.8, 4.95, 6.25]
messungC3Y = [_ / 1.9 * 2.7 for _ in messungC3Y]

messungC4X = [45, 47, 43, 41, 42, 46, 44]
messungC4Y = [2.6, 5.6, 2.3, 4.5, 3.25, 3.75, 2.]
messungC4Y = [_ / 1.9 * 2.7 for _ in messungC4Y]

messungC5X = [30, 32, 34, 35, 36, 37, 38, 39, 40]
messungC5Y = [6.7, 7.35, 5.1, 3.5, 2.2, 2.1, 3.3, 5.45, 7.55]
messungC5Y = [_ / 1.9 * 2.7 for _ in messungC5Y]

messungBEL1X = [0, 2, 4, 6, 8, 10, 12]
messungBEL1Y = [1.92, 2.07, 2.05, 2.23, 2.94, 4.27, 6.645]

messungBEL2X = [0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 9, 11]
messungBEL2Y = [5.3, 5, 4.5, 4.05, 3.6, 2.75, 3.65, 4.85, 6.85, 7.6, 7.45, 3.45, 3.1]

messungCEL1X = [14, 12, 10, 8, 6, 9, 11, 13]
messungCEL1Y = [8.4, 4.4, 2.6, 3.35, 4.1, 3.05, 2.65, 5.85]
messungCEL1Y = [_ / 1.9 * 2.7 for _ in messungCEL1Y]

messungCEL2X = [7.99, 7.5, 7, 6.5, 6, 5.5, 5, 4, 3, 2, 1]
messungCEL2Y = [9.05, 6.3, 5.15, 4.85, 4.1, 3.75, 3.45, 2.75, 2.4, 2.25, 1.9]
messungCEL2Y = [_ / 1.9 * 2.7 for _ in messungCEL2Y]

messungCEL3X = [8, 10, 11, 12, 13, 14.5, 9, 7]
messungCEL3Y = [0.32, 0.5, 0.65, 0.79, 0.9, 0.8, 0.45, 0.29]
messungCEL3Y = [_ / 1.9 * 2.7 for _ in messungCEL3Y]

try:
    fig = figure()
    for parameterNr in range(1, len(sys.argv)):
        file = sys.argv[parameterNr]
        dataX = [float(line.split()[0]) for line in open(file)]
        dataY = [float(line.split()[1]) for line in open(file)]

        if file.endswith("messungB1"):
            messungX = messung1X
            messungY = messung1Y
        elif file.endswith("messungB2"):
            messungX = messung2X
            messungY = messung2Y
        elif file.endswith("messungB3"):
            messungX = messung3X
            messungY = messung3Y
        elif file.endswith("messungB4"):
            messungX = messung4X
            messungY = messung4Y
        elif file.endswith("messungB5"):
            messungX = messung5X
            messungY = messung5Y
        elif file.endswith("lithium1"):
            messungX = lithium1X
            messungY = lithium1Y
        elif file.endswith("lithium2"):
            messungX = lithium2X
            messungY = lithium2Y
        elif file.endswith("messungC1"):
            messungX = messungC1X
            messungY = messungC1Y
        elif file.endswith("messungC2"):
            messungX = messungC2X
            messungY = messungC2Y
        elif file.endswith("messungC3"):
            messungX = messungC3X
            messungY = messungC3Y
        elif file.endswith("messungC4"):
            messungX = messungC4X
            messungY = messungC4Y
        elif file.endswith("messungC5"):
            messungX = messungC5X
            messungY = messungC5Y
        elif file.endswith("messungBEL1"):
            messungX = messungBEL1X
            messungY = messungBEL1Y
        elif file.endswith("messungBEL2"):
            messungX = messungBEL2X
            messungY = messungBEL2Y
        elif file.endswith("messungCEL1"):
            messungX = messungCEL1X
            messungY = messungCEL1Y
        elif file.endswith("messungCEL2"):
            messungX = messungCEL2X
            messungY = messungCEL2Y
        elif file.endswith("messungCEL3"):
            messungX = messungCEL3X
            messungY = messungCEL3Y

        ax = fig.add_subplot((len(sys.argv)-1) // 2 + (len(sys.argv) - 1) % 2, 2, parameterNr)
        #ax = fig.add_subplot(1, 1, parameterNr)
        ax.plot(dataX, dataY, "-", messungX, messungY, "o")
        grid()
    show()
except Exception, e:
    print e
