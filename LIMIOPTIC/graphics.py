#!/usr/bin/env python

import threading
import time
import pyqtgraph
from pyqtgraph import QtCore, QtGui
import sys
from math import *


class Plot(threading.Thread):
    def __init__(self, title="", size=(800, 350), grid=(True, True), noline=False):
        threading.Thread.__init__(self)
        # self.useAntialiasing()
        self.win = pyqtgraph.GraphicsWindow(title=title)
        self.win.resize(size[0], size[1])
        self._plots = []
        self._lines = []
        self._x = []
        self._y = []
        self._errors = []
        self._symbol_colors = []
        self._pen_colors = []
        self._plots.append(self.addPlot(noline=noline))
        self._plots[0].showGrid(x=grid[0], y=grid[1])
        self.winos = None
        self.live = False
        self.getOS()
        self._starttime = self.getTime()
        self.start()

    def setLive(self, live=True):
        self.live = live

    def addPlot(self, grid=(True, True), row=None, col=None, noline=False):
        """ returns: plot object """
        self._plots.append(self.win.addPlot(row=row, col=col))
        self._plots[-1].showGrid(x=grid[0], y=grid[1])
        if (not noline):
            self.addLine(plotindex=-1)
        return self._plots[-1]

    def getPlots(self, index=None):
        if index:
            return self._plots[index]
        else:
            return self._plots

    def addLine(self, plot=None, plotindex=0, color=(250, 250, 250), errors=[]):
        """ returns: line object """
        if plot:
            self._lines.append(plot.plot())
        else:
            self._lines.append(self._plots[plotindex].plot())
        self._x.append([])
        self._y.append([])
        self._errors.append(errors)
        self._symbol_colors.append(color)
        return self._lines[-1]

    def getLines(self, index=None):
        if index:
            return self._lines[index]
        else:
            return self._lines

    def addPoint(self, index=0, update=True, *point):
        if len(point) == 1:
            y = point[0]
            if self.live:
                x = self.getTime() - self._starttime
            else:
                x = len(self._x[index]) - 1
        else:
            x = point[0]
            y = point[1]
        self._x[index].append(x)
        self._y[index].append(y)
        if update:
            self.update(index)

    def getOS(self):
        try:
            winver = sys.winver
            self.winos = True
            return "Windows"
        except:
            self.winos = False
            return "Unix-like"

    def getTime(self):
        if self.winos is None:
            self.getOS()
        if self.winos:
            return time.perf_counter()
        else:
            return time.time()

    def getData(self, index=None):
        if not index:
            return self._x, self._y
        else:
            return self._x[index], self._y[index]

    def setData(self, lineindex, data, errors=[]):
        self._x[lineindex] = data[0]
        self._y[lineindex] = data[1]
        self._errors[lineindex] = errors

    def update(self, lineindex=None):
        if lineindex:
            self._lines[lineindex].setData(x=self._x[lineindex], y=self._y[lineindex], pen=None, symbolPen=self._symbol_colors[lineindex])
        else:
            for lineindex in range(len(self._lines)):
                self._lines[lineindex].setData(x=self._x[lineindex], y=self._y[lineindex], pen=None, symbol="o", symbolPen=self._symbol_colors[lineindex], brushSize=1, symbolSize=5)
                # self._lines[lineindex].setData(x=self._x[lineindex], y=self._y[lineindex], pen=self._colors[lineindex])

    def startTimer(self, timeout=1000):
        self.timer = QtCore.QTimer()
        self.timer.timeout.connect(self.update)
        self.timer.start(timeout)

    def findFit(self, lineindex, func=None, params=None):
        from scipy import optimize, exp, sqrt
        from numpy import array
        if params:
            p = params
        else:
            p = [0, 0, 0]
        xdata = array(self.getData(lineindex)[0][0])
        ydata = array(self.getData(lineindex)[1][0])

        ufunc = lambda p, x: p[0] + p[1] * x + p[2] * x**2
        if not func:
            func = ufunc
        erfun = lambda p, x, y: func(p, x) - y
        result = optimize.leastsq(erfun, p, (xdata, ydata))
        return result

    def findPolynomialFit(self, lineindex, num_of_params, shift_x=False):
        p = [0] * num_of_params
        if shift_x:
            p += [0]

        funcString = ""
        for i in range(num_of_params):
            funcString += "+p[{0}]*x**{0}".format(i)
        if shift_x:
            funcString.replace("x", "(x+p[{}])".format(num_of_params))
        funcString = "lambda p, x: " + funcString

        myfunc = eval(funcString)

        return self.findFit(lineindex, func=myfunc, params=p)

    def plotFunction(self, lineindex, func, xmin, xmax, num=500):
        try:
            func.__name__
        except:
            try:
                func = eval("lambda x:" + func)
            except:
                raise Exception("func is no function")
        x = [xmin + (xmax - xmin) * _ / (num - 1.) for _ in range(int(num))]
        y = list(map(func, x))
        self._x[lineindex] += x
        self._y[lineindex] += y
        self.update(lineindex)

    def plotFile(self, plotindex, f):
        try:
            data  = [line.split() for line in open(f)]
            dataX = []
            dataY = []
            for row in data:
                if len(row) > 1:
                    try:
                        dataX.append(float(row[0]))
                        dataY.append(float(row[1]))
                    except Exception as e:
                        # print(e)
                        pass
                else:
                    try:
                        dataY.append(float(row[1]))
                    except Exception as e:
                        # print(e)
                        pass
            # for line in range(0, len(dataY)):
            #     # TODO
            #     # dataY = [float(row[line]) for row in data[headerline:]]
            self.addLine(plotindex=plotindex)
            self._x[-1] = dataX if dataX else list(range(len(dataY)))
            self._y[-1] = dataY

            self.update()
        except Exception as e:
            print(e)

    def useAntialiasing(self, aa=True):
        pyqtgraph.setConfigOptions(antialias=aa)

    def streamFile(self, index=0, f=""):
            try:
                self.addPoint(index, True, self.getTime(), float(open(f, "r").read()))
            except Exception as e:
                print(e)


class Thread(threading.Thread):
    def __init__(self):
        threading.Thread.__init__(self)


def startApplication():
    QtGui.QApplication.instance().exec_()




if __name__ == "__main__":
    if len(sys.argv) > 1:
        myplot = Plot(",".join(sys.argv[1:]))
        myplot.plotFile(0, sys.argv[1])
        for f in sys.argv[2:]:
            myplot.addPlot()
            myplot.plotFile(plotindex=-1, f=f)
        startApplication()
    else:
        myplot = Plot("... ..  ..  . .  .  .    .   fastplot")
        print("possible commands:\n----------------------")
        for key in dir(myplot):
            if key.startswith("_") or key in dir(threading.Thread):
                continue
            print((key + ",",))
        print("quit\n----------------------")
        print("drag and drop filename to plot datafile\nfor some purposes you may access the plotwindow-object explicitly with 'myplot'")
        print("----------------------\n")
        while True:
            lastinput = eval(input("> "))
            if lastinput == "quit":
                break
            try:
                strippedinput = ""
                for char in lastinput:
                    if char == "(":
                        break
                    strippedinput += char
                if strippedinput in dir(myplot):
                    print(eval("myplot." + lastinput))
                else:
                    print(eval(lastinput))
            except Exception as e:
                try:
                    myplot.plotFile(len(myplot.getPlots()) - 1, lastinput)
                except:
                    print(e)