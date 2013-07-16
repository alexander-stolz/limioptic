#!/usr/bin/env python

from pylab import *
from scipy import optimize
import random
from PyQt4 import QtGui as gui
from PyQt4 import QtCore as core
import sys
#from fitting import *
#import pickle

#### INFO Bereich ####
"""
self.Source = [[x, dx, y, dy, dk, dm], ...]
Entspricht Limioptic Format

SRIM:
Zeilen. 0-9 Header, dann

 Ion  Atom   Energy        Depth       Lateral-Position        Atom Direction
 Numb Numb    (eV)          X(A)        Y(A)       Z(A)      Cos(X)  Cos(Y) Cos(Z)
T    1 17 ,2071535E+08   1000177E-02 -,2065E+02 -,4704E+02   ,9999303 -,0104408 -,0055138
T    2 17 ,2055665E+08   1000197E-02 -,4745E+02  ,1586E+02   ,9999768 -,0063336  ,0025058
...
"""


class ImportSource():
    def __init__(self):
        self.mittel = [0.]*6
        self.sigma  = [0.]*6
        self.Source = []
        self.Selection = []
        self.UserInteraction = UserInteraction(self)
        self.foilparameters = {}

    def LoadSource(self, filename, filetype="limioptic"):
        self.SourceFile = filename
        self.mittel = [0.]*6
        self.sigma  = [0.]*6
        self.Source = []
        self.Selection = []
        processedParticleNr = 0
        currentLineNr = 1

        if filetype == "limioptic":
            for line in open(self.SourceFile):
                currentLineNr += 1
                if currentLineNr % 100 == 0: print "\r{} particles processed.".format(currentLineNr),
                try:
                    temp = [float(elem) for elem in line.split()]
                    if len(temp) == 6:
                        self.Source.append(temp)
                    elif len(temp) == 4:
                        temp += [0., 0.]
                        self.Source.append(temp)
                    else:
                        print "length of line is not 4 or 6 in line:", currentLineNr
                        continue
                    processedParticleNr += 1
                    for j in xrange(6):
                        self.mittel[j] += temp[j]
                except Exception, e:
                    print "something went really wrong in line", currentLineNr
                    print "=========================\n", e, "\n=========================\n"

        elif filetype == "SRIM":
            content = open(self.SourceFile).readlines()[12:]
            currentLineNr = 12
            _laenge = len(content) + 12
            for line in content:
                currentLineNr += 1
                if currentLineNr % 100 == 0: print "\r{} of {}".format(currentLineNr, _laenge),
                try:
                    (_energy, _s, _x, _y, _cosS, _cosX, _cosY) = [float(x.replace(",", ".")) for x in line[1:].split()[2:]]
                    self.Source.append([_x*1e-7, _cosX*1e3, _y*1e-7, _cosY*1e3, _energy*1e-6, 0.])
                    processedParticleNr += 1
                    for j in xrange(6):
                        self.mittel[j] += [_x*1e-7, _cosX*1e3, _y*1e-7, _cosY*1e3, _energy*1e-6, 0.][j]
                except:
                    print "\rsomething weired happened in line", currentLineNr, "of", len(content)
            currentLineNr -= 12
            print

        ## Mittelwert berechnen
        self.mittel = [self.mittel[i] / processedParticleNr for i in xrange(6)]

        ## Fehler berechnen
        summeDeltaX = [0.]*6
        for i in xrange(processedParticleNr):
            summeDeltaX = [summeDeltaX[j] + (self.Source[i][j] - self.mittel[j])**2 for j in xrange(6)]
            if i % 100 == 0: print "\r{} of {}".format(i, processedParticleNr),
        self.sigma = [(summeDeltaX[j] / (processedParticleNr - 1))**.5 for j in xrange(6)]

        print "\r{} of {} particles read.".format(len(self.Source), currentLineNr)

        if filetype == "limioptic":
            print "average energy loss is\t{} {}.\t(+/- {})".format(self.mittel[4], "permille", self.sigma[4])
            self.foilparameters["dk"] = self.sigma[4]
        elif filetype == "SRIM":
            print "average remaining energy is\t{} {}.\t(+/- {})".format(self.mittel[4], "MeV", self.sigma[4])
            self.foilparameters["dk"] = self.sigma[4] / self.mittel[4] * 1000.

        #print "average mass loss is\t{} permille.\t(+/- {})".format(self.mittel[5], self.sigma[5])

        ## Energie Normalisieren
        if filetype == "SRIM":
            if self.mittel[4] != 0.:
                for i in xrange(len(self.Source)):
                    self.Source[i][4] = ((self.Source[i][4] / self.mittel[4]) - 1) * 1000.
            else:
                print "there was an error while calculating the relative energy deviations. Setting dK = 0."
                for i in xrange(len(self.Source)):
                    self.Source[i][4] = 0.
        else:
            for i in xrange(len(self.Source)):
                self.Source[i][4] -= self.mittel[4]

    def NormalizeEnergy(self):
        for i in xrange(len(self.Source)):
            self.Source[i][4] -= self.mittel[4]

    def NormalizeMass(self):
        for i in xrange(len(self.Source)):
            self.Source[i][5] -= self.mittel[5]

    def ShowFits(self, firstfit=True):
        art = "x", "x'", "y", "y'"
        ax  = [0]*4
        fig, ((ax[0], ax[1]), (ax[2], ax[3])) = subplots(nrows=2, ncols=2)

        # show the histogram
        for z in xrange(4):
            n, bins, patches = ax[z].hist(
                [row[z] for row in self.Source],
                bins=arange(
                    self.mittel[z] - 3. * self.sigma[z],
                    self.mittel[z] + 3. * self.sigma[z],
                    self.sigma[z] / 10.),
                normed=True)

            binsX = [(bins[i+1] + bins[i]) / 2. for i in xrange(len(bins) - 1)]
            binsY = n

            """
            # gauss-fit.
            # normpdf = (norm)alverteilte (p)ropability (d)ensity (f)unction
            y = normpdf(
                arange(
                    self.mittel[z] - 3 * self.sigma[z],
                    self.mittel[z] + 3 * self.sigma[z],
                    self.sigma[z] / 100.),
                self.mittel[z],
                self.sigma[z])
            ax[z].plot(
                arange(
                    self.mittel[z] - 3 * self.sigma[z],
                    self.mittel[z] + 3 * self.sigma[z],
                    self.sigma[z] / 100.),
                y)
            """

            # gauss-fit
            gauss  = lambda p, x: (2. * pi * p[0] * p[0])**(-.5) * exp(-(x - p[1])**2 / (2 * p[0] * p[0]))
            errfkt = lambda p, x, y: gauss(p, x) - y

            p0 = [self.sigma[z], self.mittel[z]]
            p1, success = optimize.leastsq(errfkt, p0[:], args=(binsX, binsY))

            ax[z].plot(binsX, gauss(p1, binsX), "r-")

            if firstfit:
                self.foilparameters[art[z]] = p1[0]

            """
            # cauchy-fit
            cauchy = lambda p, x: p[0] / (pi * (x - p[1])**2 + (p[0])**2)
            errfkt = lambda p, x, y: cauchy(p, x) - y

            p0 = [self.sigma[z], self.mittel[z]]
            p1, success = optimize.leastsq(errfkt, p0[:], args=(binsX, binsY))

            ax[z].plot(binsX, cauchy(p1, binsX), "g-")
            """

            """
            # voigt-fit
            voigt  = lambda p, x: p[0] * 1. / (1. + ((x - p[2]) / p[1])**2) + (1. - p[0]) * exp(-log(2) * ((x - p[2]) / p[1])**2)
            errfkt = lambda p, x, y: voigt(p, x) - y

            p0 = [.5, self.sigma[z], self.mittel[z]]
            p1, success = optimize.leastsq(errfkt, p0[:], args=(binsX, binsY))

            largeX = arange(self.mittel[z] - 3 * self.sigma[z], self.mittel[z] + 3 * self.sigma[z], self.sigma[z] / 100.)

            ax[z].plot(largeX, voigt(p1, largeX), "y-")
            """

            ax[z].grid(True)
            ax[z].set_title(art[z])

            print "mu {}\t=".format(art[z]), p1[1], "\nsigma {}\t=".format(art[z]), p1[0]

        fig.canvas.set_window_title("Source Beam - {}".format(self.SourceFile.replace("\\", "/")))
        tight_layout()
        show()

    def SelectRunaways(self, f):
        self.Selection = []
        for i in xrange(len(self.Source)):
            x  = self.Source[i][0]
            xx = self.Source[i][1]
            y  = self.Source[i][2]
            yy = self.Source[i][3]
            phi1 = arctan(xx / x)
            phi2 = arctan(yy / y)
            if (((f*self.sigma[0]*cos(phi1))**2+(f*self.sigma[1]*sin(phi1))**2 < x**2+xx**2) or ((f*self.sigma[2]*cos(phi2))**2+(f*self.sigma[3]*sin(phi2))**2 < y**2+yy**2)):
                self.Selection.append(self.Source[i])
                #print phi1, self.sigma[0],(self.sigma[0]*cos(phi1))**2+(self.sigma[1]*sin(phi1))**2 , x**2+xx**2
        self.SaveSource(data=self.Selection)
        print "\nFilter applied. Before:", len(self.Source), "Particles, after:", len(self.Selection), "Particles.\n"

    def SelectRandom(self, num):
        self.Selection = []
        random.seed()
        for i in xrange(num):
            self.Selection.append(self.Source[random.randint(0, len(self.Source) - 1)])
        self.SaveSource(data=self.Selection)
        print "\nFilter applied. Before:", len(self.Source), "Particles, after:", len(self.Selection), "Particles.\n"

    def SelectBorder(self, f, step):
        self.Selection = []
        for i in xrange(0, 360, int(step)):
            radi = radians(i)
            for j in xrange(0, 360, int(step)):
                radj = radians(j)
                self.Selection.append(
                    [
                        f * self.sigma[0] * cos(radi) * cos(radj),
                        f * self.sigma[1] * cos(radi) * sin(radj),
                        f * self.sigma[2] * sin(radi) * cos(radj),
                        f * self.sigma[3] * sin(radi) * sin(radj),
                        self.Source[random.choice(range(len(self.Source)))][4],
                        0.
                    ])
        self.SaveSource(data=self.Selection)
        print "\nFilter applied. Before:", len(self.Source), "Particles, after:", len(self.Selection), "Particles.\n"

    def SaveSource(self, data=None, filename="source.dat"):
        if data is None: data = self.Source
        f = open(filename, "w")
        for line in data:
            print >> f, line[0], line[1], line[2], line[3], line[4], line[5]
        f.close()


class UserInteraction(gui.QDialog):
    def __init__(self, parent):
        self.parent = parent
        gui.QDialog.__init__(self)
        self.setWindowTitle("Filter Particles")

    def ChooseFilter(self):
        vbox = gui.QVBoxLayout()
        self.button = []
        self.button.append(gui.QPushButton("Select Runaways"))
        self.button.append(gui.QPushButton("Select Random Particles"))
        self.button.append(gui.QPushButton("Select Border"))
        self.button.append(gui.QPushButton("Select All Particles"))

        for i in xrange(len(self.button)):
            vbox.addWidget(self.button[i])

        self.connect(self.button[0], core.SIGNAL("clicked()"), self.filter1)
        self.connect(self.button[1], core.SIGNAL("clicked()"), self.filter2)
        self.connect(self.button[2], core.SIGNAL("clicked()"), self.filter3)
        self.connect(self.button[3], core.SIGNAL("clicked()"), self.filter4)

        self.setLayout(vbox)
        self.adjustSize()
        self.show()

    def filter1(self):
        self.parent.SelectRunaways(3.)
        self.parent.Source = self.parent.Selection
        self.parent.ShowFits(firstfit=False)
        self.close()

    def filter2(self):
        self.parent.SelectRandom(1000)
        self.parent.Source = self.parent.Selection
        self.parent.ShowFits(firstfit=False)
        self.close()

    def filter3(self):
        self.parent.SelectBorder(3., 10)
        self.parent.Source = self.parent.Selection
        self.parent.ShowFits(firstfit=False)
        self.close()

    def filter4(self):
        self.parent.Selection = self.parent.Source
        self.close()


if __name__ == "__main__":
    #source = "TRANSMIT.txt"
    #source = "C:\\Users\\astolz\\Dropbox\\uni\\aaaa\\pp_to_alex.out"
    app = gui.QApplication(sys.argv)
    myapp = ImportSource()

    import Tkinter
    import tkFileDialog

    root = Tkinter.Tk()
    root.withdraw()

    options = {}
    options["filetypes"] = [("limioptic", ".out"), ("SRIM", "TRANSMIT.txt")]
    source = tkFileDialog.askopenfile(mode="r", **options).name

    root.quit()
    root.destroy()

    myapp.LoadSource(source, filetype=("SRIM" if source.endswith("TRANSMIT.txt") else "limioptic"))
    #myapp.NormalizeEnergy()
    myapp.ShowFits()
    myapp.UserInteraction.ChooseFilter()
    app.exec_()
    #myapp.Source = myapp.Selection
    #myapp.ShowFits()
    sys.exit()
