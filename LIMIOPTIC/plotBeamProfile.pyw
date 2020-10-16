#!/usr/bin/env python

import graphics
import threading
import time
import sys
import glob
from numpy import genfromtxt

data = None

class DataGrabber(threading.Thread):
    def __init__(self, datafiles=["particles.dat"]):
        threading.Thread.__init__(self)
        self.datafiles = datafiles
        self.plot = graphics.Plot(noline=True, title="Beam Profile. x -> y")
        self.plot.addLine(color=(255, 255, 255))
        # self.plot.addLine(color=(0, 255, 0))
        for i in range(len(datafiles) - 1):
            _plot = self.plot.addPlot()
            self.plot.addLine(plot=_plot, color=(255, 255, 255))
            # self.plot.addLine(plot=_plot, color=(0, 255, 0))
        self.plot.startTimer(200)

    def run(self):
        self.running = True
        while self.running:
            iplot = 0
            for f in self.datafiles:
                try:
                    # data = open(f, "r").readlines()
                    data = genfromtxt(f)
                    x    = data[:, 0]
                    # dx   = data[:, 1]
                    y    = data[:, 2]
                    # dy   = data[:, 3]
                    self.plot.setData(lineindex=iplot, data=[x, y])
                    iplot += 1
                except:
                    time.sleep(.2)
            # self.plot.update()
            time.sleep(.5)

    def close(self):
        self.running = False


def start(num):
    global data
    data = DataGrabber(["particles.dat"])
    data.start()
    #graphics.startApplication()
    #data.close()

def stop():
    global data
    try:
        data.close()
    except:
        print("failed")


if __name__ == "__main__":
    if len(sys.argv) > 1:
        data = DataGrabber(datafiles=sys.argv[1:])
    else:
        df = glob.glob("particles_*.dat")
        data = DataGrabber(datafiles=df)

    data.start()
    graphics.startApplication()
    data.close()
