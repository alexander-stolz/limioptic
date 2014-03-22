#!/usr/bin/env python

import graphics
import threading
import time
import sys

data = None

class DataGrabber(threading.Thread):
	def __init__(self, datafiles=["particles0.dat"]):
		threading.Thread.__init__(self)
		self.datafiles = datafiles
		self.plot = graphics.Plot(noline=True)
		self.plot.addLine(color=(255, 0, 0))
		self.plot.addLine(color=(0, 255, 0))
		for i in xrange(len(datafiles) - 1):
			_plot, index = self.plot.addPlot()
			self.plot.addLine(plot=_plot, color=(255, 0, 0))
			self.plot.addLine(plot=_plot, color=(0, 255, 0))
		self.plot.startTimer(200)

	def run(self):
		self.running = True
		while self.running:
			iplot = 0
			for f in self.datafiles:
				data = [line.split() for line in open(f)]
				x    = [float(row[0]) for row in data]
				dx   = [float(row[1]) for row in data]
				y    = [float(row[2]) for row in data]
				dy   = [float(row[3]) for row in data]
				self.plot.setData(lineindex=iplot,     data=[x, dx])
				self.plot.setData(lineindex=iplot + 1, data=[y, dy])
				iplot += 2
			#self.plot.update()
			time.sleep(.3)

	def close(self):
		self.running = False


def start(num):
	global data
	data = DataGrabber(["particles{}.dat".format(i) for i in xrange(num)])
	data.start()
	#graphics.startApplication()
	#data.close()

def stop():
	global data
	try:
		data.close()
	except:
		print "failed"


if __name__ == "__main__":
	data = DataGrabber(sys.argv[1:])
	data.start()
	graphics.startApplication()
	data.close()
