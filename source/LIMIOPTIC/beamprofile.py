#!/usr/bin/env python

from pylab import hist, arange, show #, figure, ion
#import time
#import threading


#ion()
#fig = figure()

while True:
    x = [float(line.split()[0]) for line in open("particles.dat")]
    hist(x, bins=arange(-9., 9., 18. / 40.))
    show()