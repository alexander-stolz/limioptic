#!/usr/bin/env python

"""
Hier wird hauptsaechlich die GUI definiert.
Der Text wird mit limioptic.ExecText() an limioptic.py uebergeben.
Dort sind die verschiedenen Funktionen definiert.
Die eigentliche Berechnung findet in climioptic.cpp statt.


Licence
=============
The software LIMIOPTIC maintained by Alexander Stolz is
freely available and distributable. However, if you use
it for some work whose results are made public, then you
have to reference it properly.
"""

try:
    QString = unicode
except NameError:
    # Python 3
    QString = str

# for buidling executable
if __name__ == "__main__":
    PY2EXE = False
else:
    PY2EXE = True

print("loading:", end=' ')
print("os", end=' ')
import os

print("sys", end=' ')
import sys

print("limioptic", end=' ')
import limioptic                                                                # ionenoptische berechnungen

print("threading", end=' ')
import threading                                                                # multithreading

print("time", end=' ')
import time                                                                     # debuggen + sleep

# print("re", end=' ')
# import re

# import serial                                                                 # auslesen des potis
if not PY2EXE:
    try:
        import vtk                                                              # grafische ausgabe ueber VTK
        print("vtk", end=' ')
    except:
        print("<vtk NOT FOUND>", end=' ')
        vtk = False
else:
    vtk = False
print("urllib", end=' ')
import urllib.request, urllib.parse, urllib.error                                                                   # zum ueberpruefen auf updates

# print("ams_spicker", end=' ')
# import ams_spicker

print("importsrc", end=' ')
from importsrc import ImportSource

print("syntax_highlighting", end=' ')
import syntax

try:
    import pyqtgraph
    # from pyqtgraph import QtCore, QtGui
    print("PyQtGraph", end=' ')
except:
    print("<PyQtGraph NOT FOUND>", end=' ')
    pyqtgraph = False
    print("Qt", end=' ')
from PyQt5 import QtCore, QtGui, QtWidgets
print("optimize", end=' ')
from scipy import optimize

# print "plotBeamprofile",
# import plotEmittance
print("helper", end=' ')
import helper

print("pickle")
import pickle

from collections import Counter

# print "calculator"
# import function_calculator

#################################################
#################################################


class inputcontrol(QtWidgets.QDialog):
    """ Hier wird das Fenster mit den Schiebereglern zur Variablenmanipulation definiert """
    def __init__(self, mode):
        QtWidgets.QDialog.__init__(self)

        global NumberOfInputs

        self.setWindowIcon(QtGui.QIcon('logo.png'))
        self.mode = mode
        self.changing = False

        if (mode == "qt"):
            self.mode = "qt"
            self.plotwindow = plot_qt(self)    # PyQtGraph
        if (mode == "2d"):
            self.plotwindow = plot_vtk(self)     # VTK
        if (mode == "3d"):
            self.plotwindow = plot_vtk3d(self)     # VTK

        # Ab hier Definition des Layouts
        # self.setWindowFlags(QtCore.Qt.WindowStaysOnTopHint)
        self.setGeometry(100, screen.height() - 300, 650, 1)

        self.setWindowTitle("input control")
        self.vbox = QtWidgets.QVBoxLayout()
        self.layout = QtWidgets.QGridLayout()

        CheckInputNumber = True
        myapp.textedit.moveCursor(QtGui.QTextCursor.Start)
        _NumberOfInputs = -1
        while (CheckInputNumber):
            _NumberOfInputs += 1
            texttofind = QString("INPUT[{}]".format(_NumberOfInputs))
            if not (myapp.textedit.find(texttofind)):
                CheckInputNumber = False
        _NumberOfInputs += 1
        if (_NumberOfInputs > NumberOfInputs):
            NumberOfInputs = _NumberOfInputs

        self.min     = []
        self.input   = []
        self.slider  = []
        self.info    = []
        self.infobox = []

        for i in range(NumberOfInputs):
            self.info.append(QtWidgets.QCheckBox())
            # self.info.append(QtGui.QLabel("#%1.0f" % (i)))
            self.min.append(QtWidgets.QDoubleSpinBox())
            self.slider.append(QtWidgets.QSlider(QtCore.Qt.Horizontal, self))
            self.input.append(QtWidgets.QDoubleSpinBox())
            self.infobox.append(QtWidgets.QLineEdit())
            self.min[i].setRange(-100., 100.)
            self.infobox[i].setPlaceholderText("INPUT[{}]".format(i))
            self.infobox[i].setText(BEZEICHNUNGEN[i])
            self.infobox[i].setToolTip("append > or # to the text and see what happens.")
            self.input[i].setDecimals(4)
            self.slider[i].setOrientation(QtCore.Qt.Horizontal)
            self.slider[i].setRange(0, 500000)
            self.slider[i].setSingleStep(500)
            self.input[i].setSingleStep(0.0001)
            self.input[i].setRange(-100., 100.)
            self.input[i].setPrefix("= ")
            self.min[i].setSingleStep(.01)
            self.input[i].setValue(INPUT[i])

            if (INPUT[i] > 5.):
                self.min[i].setValue(INPUT[i] - 5.)
            self.slider[i].setValue(int(INPUT[i] * 100000.))

            self.layout.addWidget(self.info[i], i, 0)
            self.layout.addWidget(self.min[i], i, 1)
            self.layout.addWidget(self.slider[i], i, 2)
            self.layout.addWidget(self.input[i], i, 3)
            self.layout.addWidget(self.infobox[i], i, 4)

        self.layout.setColumnStretch(2, 100)
        self.vbox.addLayout(self.layout)

        moreinputsbox        = QtWidgets.QGridLayout()
        # self.plusbutton      = QtGui.QPushButton("+ INPUT")
        # self.plusbutton.setToolTip("add a slider. if you don't use a slider, e.g. value = 1 and no text, it will not be saved.")
        # self.minusbutton     = QtGui.QPushButton("- INPUT")
        self.optimize_button = QtWidgets.QPushButton("optimize selected parameters")
        # moreinputsbox.addWidget(self.plusbutton, 0, 2)
        # moreinputsbox.addWidget(self.minusbutton, 0, 0)
        self.NumberOfInputsSpinBox = QtWidgets.QSpinBox()
        self.NumberOfInputsSpinBox.setMinimum(1)
        self.NumberOfInputsSpinBox.setValue(8)
        self.NumberOfInputsSpinBox.setToolTip("how many sliders do you want to see?")
        self.NumberOfInputsSpinBox.setSuffix(" INPUTs")
        moreinputsbox.addWidget(self.optimize_button, 0, 1)
        moreinputsbox.addWidget(self.NumberOfInputsSpinBox, 0, 2)
        moreinputsbox.setColumnStretch(1, 100)
        self.vbox.addLayout(moreinputsbox)

        # self.vbox.addWidget(self.optimize_button)

        if (self.mode == "3d"):
            opacitybox   = QtWidgets.QHBoxLayout()
            self.olabel  = QtWidgets.QLabel("Opacity")
            self.oslider = QtWidgets.QSlider(QtCore.Qt.Horizontal, self)
            self.oslider.setRange(0, 100)
            self.oslider.setValue(OPACITY)
            opacitybox.addWidget(self.olabel)
            opacitybox.addWidget(self.oslider)
            self.vbox.addLayout(opacitybox)

            scalebox     = QtWidgets.QHBoxLayout()
            self.slabel  = QtWidgets.QLabel("Scale(z)")
            self.sslider = QtWidgets.QSlider(QtCore.Qt.Horizontal, self)
            self.sslider.setRange(1, 200)
            self.sslider.setValue(SCALE3D)
            scalebox.addWidget(self.slabel)
            scalebox.addWidget(self.sslider)
            self.vbox.addLayout(scalebox)

        self.setLayout(self.vbox)

        # Die Elementgroessen werden neu angepasst (Buttons, Slider, ...)
        # self.adjustSize()

        # Der Thread des Output-Fensters wird gestartet
        self.plotwindow.start()
        if (self.mode == "qt"):
            self.plotwindow.update()

        # Wichtig, damit Aenderungen von Variablen nicht waerend des Berechnens/Plottens stattfinden
        self.threadlock = threading.Lock()

        for i in range(NumberOfInputs):
            self.slider[i].valueChanged.connect(self.slidertoinput)
            self.slider[i].sliderReleased.connect(self.centerslider)
            self.input[i].valueChanged.connect(self.inputtoslider)
            self.min[i].valueChanged.connect(self.slidertoinput)
            self.infobox[i].textChanged.connect(self.infochange)
        if (self.mode == "3d"):
            self.oslider.valueChanged.connect(self.setopacity)
            self.sslider.connect(self.setscale)
        self.optimize_button.clicked.connect(self.optimizeBeam)
        # self.connect(self.minusbutton, QtCore.SIGNAL("clicked()"), self.minusinput)
        # self.connect(self.plusbutton, QtCore.SIGNAL("clicked()"), self.minusinput)
        self.NumberOfInputsSpinBox.valueChanged.connect(self.changeNumberOfInputs)

        self.show()

        # Solange True findet Kommunikation mit Interface statt.
        self.update_on = True

        # Serielle Kommunikation findet in eigenem Thread statt.
        if (PORT != "NONE"):
            t_readserial = threading.Thread(target=self.readserial, args=())
            t_readserial.start()

        # print threading.enumerate()

    def changeNumberOfInputs(self):
        global NumberOfInputs

        for i in range(NumberOfInputs, self.NumberOfInputsSpinBox.value()):
            self.plusinput(i)
        for i in range(NumberOfInputs, self.NumberOfInputsSpinBox.value(), -1):
            self.minusinput(i - 1)

        NumberOfInputs = self.NumberOfInputsSpinBox.value()

    def minusinput(self, row):
        self.layout.removeWidget(self.info[row])
        self.layout.removeWidget(self.min[row])
        self.layout.removeWidget(self.slider[row])
        self.layout.removeWidget(self.input[row])
        self.layout.removeWidget(self.infobox[row])
        self.info[row].deleteLater()
        self.min[row].deleteLater()
        self.slider[row].deleteLater()
        self.input[row].deleteLater()
        self.infobox[row].deleteLater()
        del self.info[row]
        del self.min[row]
        del self.slider[row]
        del self.input[row]
        del self.infobox[row]
        self.adjustSize()
        self.resize(500, 1)

    def plusinput(self, row):
        # global NumberOfInputs
        # NumberOfInputs += 1

        self.info.append(QtWidgets.QCheckBox())
        self.min.append(QtWidgets.QDoubleSpinBox())
        self.slider.append(QtWidgets.QSlider(QtCore.Qt.Horizontal, self))
        self.input.append(QtWidgets.QDoubleSpinBox())
        self.infobox.append(QtWidgets.QLineEdit())
        self.min[-1].setRange(-100., 100.)
        self.infobox[-1].setPlaceholderText("INPUT[{}]".format(row))
        self.infobox[-1].setText(BEZEICHNUNGEN[row])
        self.input[-1].setDecimals(4)
        self.slider[-1].setOrientation(QtCore.Qt.Horizontal)
        self.slider[-1].setRange(0, 500000)
        self.slider[-1].setSingleStep(500)
        self.input[-1].setSingleStep(0.0001)
        self.input[-1].setRange(-100., 100.)
        self.input[-1].setPrefix("= ")
        self.min[-1].setSingleStep(.01)
        self.input[-1].setValue(INPUT[row])

        if (INPUT[row] > 5.):
            self.min[row].setValue(INPUT[row] - 5.)
        self.slider[row].setValue(int(INPUT[row] * 100000.))

        self.layout.addWidget(self.info[-1], row, 0)
        self.layout.addWidget(self.min[-1], row, 1)
        self.layout.addWidget(self.slider[-1], row, 2)
        self.layout.addWidget(self.input[-1], row, 3)
        self.layout.addWidget(self.infobox[-1], row, 4)

        self.slider[-1].valueChanged.connect(self.slidertoinput)
        self.slider[-1].sliderReleased.connect(self.centerslider)
        self.input[-1].valueChanged.connect(self.inputtoslider)
        self.min[-1].valueChanged.connect(self.slidertoinput)
        self.infobox[-1].textChanged.connect(self.infochange)

    def optimizeBeam(self):
        optIndex = []
        for i in range(NumberOfInputs):
            if self.info[i].isChecked():
                optIndex.append(i)

        beamline = str(myapp.textedit.toPlainText())

        i = 0
        _param_ = []
        for opt in optIndex:
            beamline = beamline.replace("INPUT[{}]".format(opt), "_param_[{}]".format(i))
            beamline = beamline.replace("BeamProfile()", "# BeamProfile()")
            _param_.append(INPUT[opt])
            i += 1

        # for i in xrange(NumberOfInputs):
        #    beamline = beamline.replace("INPUT[{}]".format(i), str(INPUT[i]))

        print("computing...")
        self.optimize_button.setText("optimizing.. please wait.")
        app.processEvents()

        # _beamline = ""
        # for row in beamline.split("\n"):
        #     if row.startswith("Slit("):
        #         a, b, c, d = re.search(r"Slit\((\d+\.*\d*)\s*,\s*(\d+\.*\d*)\s*,\s*(\d+\.*\d*)\s*,(\d+\.*\d*)\)", row).groups()
        #         _beamline += f"Slit({a}, {b}, {c}, {d}, output=False)\n"
        #     else:
        #         _beamline += row + "\n"

        result = optimize.minimize(limioptic.ErrFkt, _param_[:], (beamline.replace("Slit(", "Slit_muted("), INPUT, SourceObj.Source))

        self.optimize_button.setText("optimize selected parameters")

        print(result)

        print(limioptic.optic.GetSpotSize())

        j = 0
        self.changing = True
        for i in optIndex:
            INPUT[i] = result.x[j]
            self.input[i].setValue(result.x[j])
            j += 1
        self.changing = False
        self.inputtoslider()

    def setopacity(self):
        """ Nur fuer die 3D-Ausgabe """
        OPACITY = self.oslider.value()
        self.plotwindow.actor.GetProperty().SetOpacity(OPACITY / 1000.)
        self.plotwindow.render = True

    def setscale(self):
        """ Nur fuer die 3D-Ausgabe """
        global SCALE3D
        SCALE3D = self.sslider.value()
        self.plotwindow.neu()

    def readserial(self):
        """ Serielle Kommunikation mit dem Interface (obsolet) """
        ser = serial.Serial(PORT, 9600)
        time.sleep(2)
        print("serial start")
        while (self.update_on):
            ser.write(1)
            # Interface braucht kurz Zeit zum Antworten
            time.sleep(0.1)
            a = ser.readline().split()
            self.slider[0].setValue(int(a[0]) * 488)
            self.slider[1].setValue(int(a[1]) * 488)
            self.slider[2].setValue(int(a[2]) * 488)
        print("serial end")

    def inputtoslider(self):
        """ Variablenaenderungen werden auf Slider uebertragen. """
        if not self.changing:
            global INPUT
            self.changing = True
            for i in range(NumberOfInputs):
                    INPUT[i] = self.input[i].value()
            for i in range(NumberOfInputs):
                    self.min[i].setValue(INPUT[i] - 2.5)
                    self.slider[i].setValue(250000)

            if (RUNNINGQT):
                self.plotwindow.update(self.calculate())
            if (RUNNING2D):
                self.plotwindow.update = True
            if (RUNNING3D):
                self.plotwindow.neu()

            self.changing = False

    def centerslider(self):
        self.changing = True
        for i in range(NumberOfInputs):
            self.min[i].setValue(INPUT[i] - 2.5)
            self.slider[i].setValue(250000)
        self.changing = False

    def slidertoinput(self):
        """ Aenderung des Sliders wird auf self.INPUT uebertragen. """
        if not self.changing:
            global INPUT
            self.changing = True
            for i in range(NumberOfInputs):
                self.input[i].setValue(self.slider[i].value() / 100000. + self.min[i].value())
            for i in range(NumberOfInputs):
                INPUT[i] = self.input[i].value()

            if (RUNNINGQT):
                self.plotwindow.update(self.calculate())
            if (RUNNING2D):
                self.plotwindow.update = True
            if (RUNNING3D):
                self.plotwindow.neu()

            self.changing = False

    def closeit(self):
        """ Wird aufgerufen, wenn das Outputfenster geschlossen wird. """
        global BEZEICHNUNGEN, RUNNING2D, RUNNING3D, NumberOfInputs

        for i in range(NumberOfInputs):
            BEZEICHNUNGEN[i] = self.infobox[i].text()

        self.update_on = False

        if (self.mode == "qt"):
            self.plotwindow.closeme()
        else:
            self.plotwindow.iren.Disable()
            self.plotwindow.iren.EnableRenderOff()
            self.plotwindow.iren.TerminateApp()
            self.plotwindow.iren.DestroyTimer(self.plotwindow.timer)
            if (self.mode == "2d"):
                RUNNING2D = False
                del self.plotwindow.view
                del self.plotwindow.iren
                del self.plotwindow.chart
            if (self.mode == "3d"):
                RUNNING3D = False
                del self.plotwindow.ren
                del self.plotwindow.actor
                del self.plotwindow.mapper
                del self.plotwindow.iren
                del self.plotwindow.renwin
                del self.plotwindow.mydata
                del self.plotwindow.polylines
                del self.plotwindow.mycells
                del self.plotwindow.writer
                del self.plotwindow.w2iFilter
                for segment in self.plotwindow.segments:
                    del segment
                del self.plotwindow.segments
        # del self.plotwindow
        # plotEmittance.stop()
        # print threading.enumerate()

        print("saving autosave..", end=' ')
        myfile = open(backup_file + ".lim", "w")
        myfile.write(str(myapp.textedit.toPlainText()))
        myfile.close()

        myfile = open(backup_file + ".var", "w")
        for i in range(8):
            print("{} = {}".format(BEZEICHNUNGEN[i], INPUT[i]), file=myfile)
        for i in range(NumberOfInputs - 1, 7, -1):
            if (INPUT[i] == 1.) and (BEZEICHNUNGEN[i] == ""):
                NumberOfInputs -= 1
                continue
            break
        for i in range(8, NumberOfInputs):
                print("{} = {}".format(BEZEICHNUNGEN[i], INPUT[i]), file=myfile)
        myfile.close()

        del self.plotwindow

        time.sleep(.5)
        self.close()
        print("\rsaved to {}".format(backup_file))

    def infochange(self):
        """ Spezialbefehle im "Beschreibung" Feld """
        for i in range(NumberOfInputs):
            if self.infobox[i].text()[-1:] == ">":
                myapp.textedit.setText("{0}\t=\tINPUT[{1}]\t# {2}\n{3}".format(
                    self.infobox[i].text()[:-1], i, INPUT[i], myapp.textedit.toPlainText()))
                self.infobox[i].setText(self.infobox[i].text()[:-1])
            if self.infobox[i].text()[-1:] == "#":
                myapp.textedit.setText("{0}\t=\t{2}\t# INPUT[{1}]\n{3}".format(
                    self.infobox[i].text()[:-1], i, INPUT[i], myapp.textedit.toPlainText()))
                self.infobox[i].setText(self.infobox[i].text()[:-1])

    def calculate(self):
        """ Neuberechnung """
        limioptic.optic.Clear()
        limioptic.geo_s.Reset()
        limioptic.geo_y.Reset()
        limioptic.s = 0.
        try:
            limioptic.ExecText(str(myapp.textedit.toPlainText()), INPUT, SourceObj.Source)
            limioptic.optic.CalculateTrajectories()
        except Exception as e:
            print("\n\nFehler in der Eingabe! ({})".format(limioptic.lastFunction), "\n===============================\n", e, "\n===============================\n\n")
            return

        parts = int(limioptic.optic.GetParticleNum())                    # Anz. Partikel
        segs  = int(limioptic.optic.GetTrajectoriesSize() / parts / 8)   # Anz. Segmente

        # Die Trajektorien liegen als array.array("d", ..) vor
        _xi = [limioptic.GetTrajectory(i, 0) for i in range(parts)]
        _yi = [limioptic.GetTrajectory(i, 2) for i in range(parts)]
        _z  = limioptic.GetTrajectory(0, 6)

        # Umwandlung in float arrays
        xi = [None] * parts
        yi = [None] * parts
        for part in range(parts):
            xi[part] = [None] * segs
            yi[part] = [None] * segs
            for seg in range(segs):
                xi[part][seg] = float(_xi[part][seg])
                yi[part][seg] = float(_yi[part][seg])
        zi = [float(_z[seg]) for seg in range(segs)]

        return xi, yi, zi, segs, parts


class MyQtWindow(pyqtgraph.GraphicsWindow):
    def closeEvent(self, event):
        event.accept()
        myapp.inputwindowqt.closeit()


class plot_qt(threading.Thread):
    """ 2D Plot mit PyQtGraph """
    def __init__(self, parent):
        threading.Thread.__init__(self)
        self.name = "qt plot window thread"
        self.running = True
        self.parent = parent
        if not myapp.menu_plot_bg.isChecked():
            pyqtgraph.setConfigOptions(background="w")
        else:
            pyqtgraph.setConfigOptions(background="k")
        if myapp.menu_output_smoothing.isChecked():
            pyqtgraph.setConfigOptions(antialias=True)
        else:
            pyqtgraph.setConfigOptions(antialias=False)
        self.win = MyQtWindow(title="LIMIOPTIC - Output (2D)")
        self.win.setWindowIcon(QtGui.QIcon('logo.png'))
        self.win.resize(1000, 450)
        self.plot1 = self.win.addPlot()
        if myapp.menu_plot_splitview.isChecked():
            self.plot2 = self.win.addPlot(row=1, col=0)
            self.plot2.showGrid(x=True, y=True)
        self.plot1.showGrid(x=True, y=True)
        self.lineX = self.plot1.plot()
        if myapp.menu_plot_splitview.isChecked():
            self.lineY = self.plot2.plot()
            self.linegeo2 = self.plot2.plot()
        else:
            self.lineY = self.plot1.plot()
        self.linegeo = self.plot1.plot()
        self.labels = []

        if myapp.menu_output_file.isChecked():
            (xi, yi, zi, segs, parts) = self.parent.calculate()
            with open("beam.dat", "w") as f:
                for part in range(parts):
                    print(f"# particle {part}", file=f)
                    print(f"# --------------------------------------", file=f)
                    for seg in range(segs):
                        print(zi[seg], xi[part][seg], yi[part][seg], file=f)

    def run(self):
        while self.running:
            time.sleep(.1)

    def update(self, xxx_todo_changeme=(None, None, None, None, None)):
        (xi, yi, zi, segs, parts) = xxx_todo_changeme
        if xi is None:
            (xi, yi, zi, segs, parts) = self.parent.calculate()

        # viel schneller, als mehrere plots
        x_all = []
        y_all = []
        z_all = []
        # self.plot1.clear()
        for part in range(parts):
            x_all += xi[part] + [0., 0.]
            y_all += yi[part] + [0., 0.]
            z_all += zi + [zi[-1], zi[0]]
            # self.plot1.plot(x=zi, y=xi[part], pen=(255, 0, 0))
            # self.plot1.plot(x=zi, y=yi[part], pen=(0, 255, 0))

        if myapp.menu_plot_x.isChecked():
            self.lineX.setData(x=z_all, y=x_all, pen=(255, 0, 0))
        else:
            self.lineX.setData(x=[], y=[], pen=(255, 0, 0))
        if myapp.menu_plot_y.isChecked():
            self.lineY.setData(x=z_all, y=y_all, pen=(0, 255, 0))
        else:
            self.lineY.setData(x=[], y=[], pen=(255, 0, 0))

        # Labels
        for label in self.labels:
            self.plot1.removeItem(label)
            del label
        self.labels = []
        for label in limioptic.textArray:
            self.labels.append(pyqtgraph.TextItem(label[1], angle=-90, anchor=(0, .5), color=(130, 130, 130)))
            self.labels[-1].setPos(label[0], 56)
            self.plot1.addItem(self.labels[-1])

        iele  = limioptic.geo_s.GetNumberOfTuples()
        geo_s = [limioptic.geo_s.GetValue(i) for i in range(iele)]
        geo_y = [limioptic.geo_y.GetValue(i) for i in range(iele)]
        if myapp.menu_plot_geo.isChecked():
            if not myapp.menu_plot_bg.isChecked():
                self.linegeo.setData(x=geo_s, y=geo_y, pen=(0, 0, 250))
                if myapp.menu_plot_splitview.isChecked():
                    self.linegeo2.setData(x=geo_s, y=geo_y, pen=(0, 0, 250))
            else:
                self.linegeo.setData(x=geo_s, y=geo_y, pen=(170, 170, 170))
                if myapp.menu_plot_splitview.isChecked():
                    self.linegeo2.setData(x=geo_s, y=geo_y, pen=(170, 170, 170))

    def closeme(self):
        print("close qt2 window")
        try:
            del self.win
            del self.lineX
            del self.lineY
            del self.linegeo
            del self.plot1
        except:
            print("error in qt closeme")
        self.running = False


class Segment:
    """ Fuer die 3D-Ausgabe der Geometrie """
    def __init__(self, s1, s2, radius):
        self.source     = vtk.vtkLineSource()
        self.tubefilter = vtk.vtkTubeFilter()
        self.mapper     = vtk.vtkPolyDataMapper()
        self.actor      = vtk.vtkActor()

        self.source.SetPoint1(s1, 0, 0)
        self.source.SetPoint2(s2, 0, 0)

        self.tubefilter.SetInputConnection(self.source.GetOutputPort())
        self.tubefilter.SetRadius(radius)
        self.tubefilter.SetNumberOfSides(50)
        self.tubefilter.Update()

        self.mapper.SetInputConnection(self.tubefilter.GetOutputPort())
        self.actor.SetMapper(self.mapper)

        self.actor.GetProperty().SetOpacity(.25)
        self.actor.GetProperty().SetColor(255, 255, 0)

        # print "added: ", s1, s2, radius

    def __del__(self):
        del self.actor
        del self.mapper
        del self.tubefilter
        del self.source


if vtk:
    class Text3D(vtk.vtkActor):
        """ Labels in der 3D-Ausgabe """
        def __init__(self, pos, text, posY=60.):
            posY = posY / SUPERSCALE3D
            self.vectorlabel = vtk.vtkVectorText()
            self.vectorlabel.SetText(text)
            self.extrusionfilter = vtk.vtkLinearExtrusionFilter()
            self.extrusionfilter.SetInputConnection(self.vectorlabel.GetOutputPort())
            self.extrusionfilter.SetExtrusionTypeToNormalExtrusion()
            self.extrusionfilter.SetVector(0, 0, 1)
            self.extrusionfilter.SetScaleFactor(.5)

            self.mapper = vtk.vtkPolyDataMapper()
            self.mapper.SetInputConnection(self.extrusionfilter.GetOutputPort())

            self.SetMapper(self.mapper)
            self.SetPosition(pos, posY, 0)
            self.SetScale(.2, .2, .2)
            self.RotateZ(45)


class plot_vtk3d(threading.Thread):
    """ 3D-Ausgabe """
    def __init__(self, parent):
        threading.Thread.__init__(self)
        self.parent = parent
        self.render = False
        self.threadlock = threading.Lock()

    def run(self):
        (xi, yi, zi, segs, parts) = self.parent.calculate()

        # Hier landen alle Punkte
        self.mypoints = vtk.vtkPoints()

        # Wir erzeugen nur eine Cell, dh. eine lange Polyline (Geschwindigkeit)
        self.mycells = vtk.vtkCellArray()

        # Eine neue Polygonlinie definieren.
        self.polylines = [vtk.vtkPolyLine() for line in range(parts + 1 + len(zi) + 1)]

        # Alle Punkte in mypoints laden und mit polylines verknuepfen. Jede Polyline kommt in eine Cell.
        for part in range(parts):
            self.polylines[part].GetPointIds().SetNumberOfIds(segs)
            for seg in range(segs):
                self.mypoints.InsertNextPoint(
                    zi[seg] * SCALE3D / SUPERSCALE3D,
                    xi[part][seg] / SUPERSCALE3D,
                    yi[part][seg] / SUPERSCALE3D)
                self.polylines[part].GetPointIds().SetId(seg, part * segs + seg)
            self.mycells.InsertNextCell(self.polylines[part])

        # Sollbahn
        self.polylines[parts].GetPointIds().SetNumberOfIds(2)
        self.mypoints.InsertNextPoint(0., 0., 0.)
        self.mypoints.InsertNextPoint(zi[-1] * SCALE3D / SUPERSCALE3D, 0., 0.)
        self.polylines[parts].GetPointIds().SetId(0, parts * segs)
        self.polylines[parts].GetPointIds().SetId(1, parts * segs + 1)
        self.mycells.InsertNextCell(self.polylines[parts])

        # Marker
        marker = 1
        for seg in zi:
            self.mypoints.InsertNextPoint(seg * SCALE3D / SUPERSCALE3D, -10. / SUPERSCALE3D, 0)
            self.mypoints.InsertNextPoint(seg * SCALE3D / SUPERSCALE3D, +10. / SUPERSCALE3D, 0)
            self.polylines[parts + marker].GetPointIds().SetNumberOfIds(2)
            self.polylines[parts + marker].GetPointIds().SetId(0, parts * segs + 2 * marker)
            self.polylines[parts + marker].GetPointIds().SetId(1, parts * segs + 2 * marker + 1)
            self.mycells.InsertNextCell(self.polylines[parts + marker])
            marker += 1

        # Geometrie
        if (myapp.menu_plot_geo.isChecked()):
            iele  = limioptic.geo_s.GetNumberOfTuples()
            # self.polylines[-1].GetPointIds().SetNumberOfIds(iele)
            self.segments = []
            for point in range(iele / 2):
                s1 = limioptic.geo_s.GetValue(point * 2) * SCALE3D
                s2 = limioptic.geo_s.GetValue((point * 2) + 1) * SCALE3D
                radius = limioptic.geo_y.GetValue(point * 2)
                if ((s2 - s1) / SCALE3D > .1) and (radius != 55.):
                    self.segments.append(Segment(
                        s1 / SUPERSCALE3D,
                        s2 / SUPERSCALE3D,
                        radius / SUPERSCALE3D))
                # self.mypoints.InsertNextPoint(limioptic.geo_s.GetValue(point) * SCALE3D, limioptic.geo_y.GetValue(point), 0.)
                # self.polylines[-1].GetPointIds().SetId(point, parts * segs + 2 * marker + point)

            # self.mycells.InsertNextCell(self.polylines[-1])

        # Datenobjekt erzeugen
        self.mydata = vtk.vtkPolyData()

        # Punkte hinzufuegen
        self.mydata.SetPoints(self.mypoints)

        # Polyline hinzufuegen. Ab nun verknuepft
        self.mydata.SetLines(self.mycells)

        # Der Mapper macht aus Daten mehr
        self.mapper = vtk.vtkPolyDataMapper()

        # Hier bekommt er die Daten
        self.mapper.SetInput(self.mydata)

        # Der Actor platziert den Mapper im Raum
        self.actor = vtk.vtkActor()

        # Hier bekommt er den Mapper
        self.actor.SetMapper(self.mapper)
        self.actor.GetProperty().SetOpacity(OPACITY / 1000.)
        # self.actor.GetProperty().SetColor(0., 0., 0.)

        self.ren = vtk.vtkRenderer()
        # self.ren.SetBackground(.1, .2, .4)
        self.ren.GradientBackgroundOn()
        self.ren.SetBackground(.01, .01, .01)
        self.ren.SetBackground2(.2, .2, .2)
        self.ren.AddActor(self.actor)
        for segment in self.segments:
            self.ren.AddActor(segment.actor)

        # Axen
        # self.axisactor = vtk.vtkCubeAxesActor()
        # self.axisactor.SetXAxisRange(0., zi[-1])
        # self.axisactor.SetYAxisRange(-30., 30)
        # self.axisactor.SetBounds(self.actor.GetBounds())
        # self.axisactor.SetXTitle("s in m")
        # self.axisactor.YAxisVisibilityOff()
        # self.axisactor.ZAxisVisibilityOff()
        # self.axisactor.SetCamera(self.ren.GetActiveCamera())
        # self.axisactor.GetProperty().SetColor(1, 1, 1)
        # self.ren.AddActor(self.axisactor)
        self.ren.ResetCamera()

        # Labels
        labels = []
        for label in limioptic.textArray:
            labels.append(Text3D(label[0] * SCALE3D / SUPERSCALE3D, label[1]))
            self.ren.AddActor(labels[-1])
        sourcelabel = Text3D(-1 * SCALE3D / SUPERSCALE3D, "*", -8.)
        self.ren.AddActor(sourcelabel)

        # Renderfenster
        self.renwin = vtk.vtkRenderWindow()
        self.renwin.AddRenderer(self.ren)
        self.renwin.SetSize(900, 450)

        self.iren = vtk.vtkRenderWindowInteractor()

        self.iren.SetRenderWindow(self.renwin)
        self.iren.SetInteractorStyle(vtk.vtkInteractorStyleTrackballCamera())

        # Text
        # txtFilter = vtk.vtkLinearExtrusionFilter()
        # txtFilter.SetScaleFactor(20.)
        # txtMapper = vtk.vtkPolyDataMapper()
        # txtMapper.SetInputConnection(txtFilter.GetOutputPort())

        # for name in limioptic.textArray
        #    text = vtk.vtkVectorText()
        #    text.SetText(name[1])
        #    txtFilter.SetInputConnection(text.GetOutputPort())
        #    txtActor = vtk.vtkActor()
        #    txtActor.SetMapper(txtMapper)
        #    self.ren.AddActor(txtActor)

        # Exporter
        self.writer = vtk.vtk.vtkPNGWriter()
        self.w2iFilter = vtk.vtkWindowToImageFilter()
        self.w2iFilter.SetInput(self.renwin)
        self.w2iFilter.SetMagnification(10)
        self.writer.SetInput(self.w2iFilter.GetOutput())

        # self.exporter = vtk.vtkGL2PSExporter()
        # self.exporter.SetRenderWindow(self.renwin)
        # self.exporter.SetFileFormatToPS()
        # self.exporter.CompressOff()
        # self.exporter.SetSortToBsp()
        # self.exporter.DrawBackgroundOn()
        # self.exporter.Write3DPropsAsRasterImageOn()

        # Damit der Speicher wieder freigegeben wird, wenn das Fenster geschlossen wird.
        self.iren.AddObserver("ExitEvent", lambda o, e, a=myapp.inputwindow3d: a.closeit())

        # Linien glatter machen
        if (myapp.menu_output_smoothing.isChecked()):
            self.renwin.LineSmoothingOn()

        # Rot-Gruen-Brille (aktivieren mit '3')
        self.renwin.SetStereoTypeToAnaglyph()

        self.iren.Initialize()
        self.renwin.Render()

        self.renwin.SetWindowName("LIMIOPTIC - Output (3D)  |  r: reset view - space: take screenshot")

        # self.animate wird aufgerufen wenn der Timer zuschlaegt (alle 50 ms)
        self.iren.AddObserver("TimerEvent", self.animate)
        self.iren.AddObserver("KeyPressEvent", self.keypress)
        self.timer = self.iren.CreateRepeatingTimer(50)

        # Ab hier laeuft die Schleife bis das Fenster geschlossen wird
        self.iren.Start()

    def neu(self, xxx_todo_changeme1=(None, None, None, None, None)):
        """ Wird aufgerufen wenn eine Variable geaendert wird, oder Strg + H gedrueckt wird """
        (xi, yi, zi, segs, parts) = xxx_todo_changeme1
        if xi is None:
            (xi, yi, zi, segs, parts) = self.parent.calculate()

        for part in range(parts):
                for seg in range(segs):
                    self.mypoints.InsertPoint(
                        part * segs + seg,
                        zi[seg] * SCALE3D / SUPERSCALE3D,
                        xi[part][seg] / SUPERSCALE3D,
                        yi[part][seg] / SUPERSCALE3D)

        self.mypoints.Modified()
        self.render = True

    def animate(self, obj=None, event=None):
        """ Nur animieren, wenn etwas geaendert wurde. Die Funktion wird alle 50 ms aufgerufen """
        # if self.render:
        # self.render = False
        self.threadlock.acquire()
        self.renwin.Render()
        self.threadlock.release()

    def keypress(self, obj=None, event=None):
        key = obj.GetKeySym()
        if key == "space":
            global SCREENSHOTNUMBER
            now = time.localtime()
            self.w2iFilter.Modified()
            self.writer.SetFileName("screenshot_{} ({}.{}.{} at {}.{}).png".format(
                SCREENSHOTNUMBER,
                now[0],
                now[1],
                now[2],
                now[3],
                now[4]))
            self.writer.Write()
            SCREENSHOTNUMBER += 1
            self.renwin.Render()


class plot_vtk(threading.Thread):
    """ Die Ausgabe in 2D """
    def __init__(self, parent):
        self.parent = parent
        threading.Thread.__init__(self)
        self.activeInput = 0
        self.render = False
        self.update = False
        self.threadlock = threading.Lock()

    def run(self):
        # 2d Szene und xy Chart erzeugen.
        self.view = vtk.vtkContextView()
        self.view.GetRenderer().SetBackground(1., 1., 1.)
        if myapp.menu_plot_bg.isChecked():
            self.view.GetRenderer().SetBackground(0., 0., 0.)
        # self.view.GetRenderWindow().SetSize(screen.width() - 375, screen.height() - 40)
        self.view.GetRenderWindow().SetSize(800, 350)

        self.chart = vtk.vtkChartXY()
        self.chart.GetAxis(vtk.vtkAxis.BOTTOM).SetTitle("beamline (m)")
        self.chart.GetAxis(vtk.vtkAxis.LEFT).SetTitle("deviation (mm)")

        self.Ticks  = vtk.vtkDoubleArray()
        self.Labels = vtk.vtkStringArray()

        if len(limioptic.textArray) > 0:
            self.chart.GetAxis(vtk.vtkAxis.BOTTOM).SetBehavior(vtk.vtkAxis.CUSTOM)
            self.Ticks.Reset()
            self.Labels.Reset()
            for name in limioptic.textArray:
                self.Ticks.InsertNextValue(name[0])
                self.Labels.InsertNextValue(name[1])
            self.chart.GetAxis(vtk.vtkAxis.BOTTOM).SetTickPositions(self.Ticks)
            self.chart.GetAxis(vtk.vtkAxis.BOTTOM).SetTickLabels(self.Labels)

        self.view.GetScene().AddItem(self.chart)

        if (myapp.menu_plot_geo.isChecked()):
            limioptic.geo_s.Reset()
            limioptic.geo_y.Reset()
            limioptic.s = 0.
        else:
            limioptic.s = -1.

        # Exporter
        # self.writer = vtk.vtk.vtkPNGWriter()
        # self.w2iFilter = vtk.vtkWindowToImageFilter()
        # self.w2iFilter.SetInput(self.view.GetRenderWindow())
        # self.w2iFilter.SetMagnification(10)
        # self.writer.SetInput(self.w2iFilter.GetOutput())

        # Berechnungen durchfuehren
        (xi, yi, zi, segs, parts) = self.parent.calculate()
        iele = limioptic.GetTrajectory(0, 7)

        if (myapp.menu_output_file.isChecked()):
            ausgabe = open("output_markers.dat", "w")

        # Marker setzen
        if myapp.menu_plot_marker.isChecked():
            self.markertable = vtk.vtkTable()
            self.markersY = vtk.vtkFloatArray()
            self.markersY.SetName("Marker")
            self.markersY.InsertNextValue(-10)
            self.markersY.InsertNextValue(10)
            self.markertable.AddColumn(self.markersY)

            self.markersX = []

            self.linelist = []
            if (len(iele) > 0):
                for i in range(1, len(iele)):
                    if (iele[i] != iele[i - 1]):
                        self.linelist.append(zi[i - 1])
                    self.linelist.append(zi[len(iele) - 1])

            if (myapp.menu_output_file.isChecked()):
                print(0, 0, file=ausgabe)
            for i in range(1, len(self.linelist)):
                self.markersX.append(vtk.vtkFloatArray())
                self.markersX[i - 1].SetName("Marker %1.0f" % (i))
                self.markersX[i - 1].InsertNextValue(self.linelist[i - 1])
                self.markersX[i - 1].InsertNextValue(self.linelist[i - 1])
                if (myapp.menu_output_file.isChecked()):
                    print(self.linelist[i - 1], 0, file=ausgabe)
                    print(self.linelist[i - 1], 10, file=ausgabe)
                    print(self.linelist[i - 1], -10, file=ausgabe)
                    print(self.linelist[i - 1], 0, file=ausgabe)
                self.markertable.AddColumn(self.markersX[i - 1])
                self.line2 = self.chart.AddPlot(0)
                self.line2.SetInput(self.markertable, i, 0)
                self.line2.SetColor(.7, .7, .7)
                self.line2.SetWidth(1.)
            if (myapp.menu_output_file.isChecked()):
                ausgabe.close()

        # Erzeuge Wertetabelle
        self.table = vtk.vtkTable()
        self.arrZ = vtk.vtkFloatArray()
        self.arrZ.SetName("Z Achse")
        self.arrX = []
        self.arrY = []

        self.arrZ.SetNumberOfValues(segs)
        self.arrZ.Reset()
        for seg in range(segs):
            self.arrZ.InsertNextValue(zi[seg])

        self.table.AddColumn(self.arrZ)

        j = 0
        if (plotx):
            if (myapp.menu_output_file.isChecked()):
                ausgabe = open("output_xbeam.dat", "w")
            for j in range(parts):
                self.arrX.append(vtk.vtkFloatArray())
                self.arrX[j].SetName("X-Strahl %1.0f Start: %1.2f" % (j, xi[j][0]))
                for i in range(segs):
                    self.arrX[j].InsertNextValue(xi[j][i])
                    if (myapp.menu_output_file.isChecked()) and not (xi[j][i] == yi[j][i] == 0.):
                        print(zi[i], xi[j][i], file=ausgabe)
                self.table.AddColumn(self.arrX[j])
                if (myapp.menu_output_file.isChecked()):
                    ausgabe.write("\n")

                self.line = self.chart.AddPlot(0)
                self.line.SetInput(self.table, 0, j + 1)
                # self.line.SetColor(255, 0, 0, 255)
                self.line.SetColor(XCOLOR[0], XCOLOR[1], XCOLOR[2], XCOLOR[3])
                self.line.SetWidth(.7)
            if (myapp.menu_output_file.isChecked()):
                ausgabe.close()
        if (ploty):
            if (myapp.menu_output_file.isChecked()):
                ausgabe = open("output_ybeam.dat", "w")
            for l in range(parts):
                self.arrY.append(vtk.vtkFloatArray())
                self.arrY[l].SetName("Y-Strahl %1.0f Start: %1.2f" % (l, yi[l][0]))
                for i in range(segs):
                    self.arrY[l].InsertNextValue(yi[l][i])
                    if (myapp.menu_output_file.isChecked()) and not (xi[l][i] == yi[l][i] == 0.):
                        print(zi[i], yi[l][i], file=ausgabe)
                self.table.AddColumn(self.arrY[l])
                if (myapp.menu_output_file.isChecked()):
                    ausgabe.write("\n")

                self.line = self.chart.AddPlot(0)
                self.line.SetInput(self.table, 0, l + j + xy)
                # self.line.SetColor(0, 255, 0, 255)
                self.line.SetColor(YCOLOR[0], YCOLOR[1], YCOLOR[2], YCOLOR[3])
                self.line.SetWidth(.7)
            if (myapp.menu_output_file.isChecked()):
                ausgabe.close()

        if (myapp.menu_plot_geo.isChecked()):
            limioptic.geolines.AddColumn(limioptic.geo_s)
            limioptic.geolines.AddColumn(limioptic.geo_y)
            self.line3 = self.chart.AddPlot(0)
            self.line3.SetInput(limioptic.geolines, 0, 1)
            if myapp.menu_plot_bg.isChecked():
                self.line3.SetColor(1, 1, 0)
            else:
                self.line3.SetColor(0., 0., 1.)
            self.line3.SetWidth(.7)

        for name in limioptic.textArray:
            self.Ticks.InsertNextValue(name[0])
            self.Labels.InsertNextValue(name[1])

        if (myapp.menu_output_smoothing.isChecked()):
            self.view.GetRenderWindow().LineSmoothingOn()
        self.view.Render()

        self.iren = self.view.GetInteractor()
        self.iren.AddObserver("TimerEvent", self.animate)
        self.iren.AddObserver("KeyPressEvent", self.keypress)
        self.iren.AddObserver("ExitEvent", lambda o, e, a=myapp.inputwindow2d: a.closeit())
        self.timer = self.iren.CreateRepeatingTimer(50)

        self.view.GetRenderWindow().SetWindowName("LIMIOPTIC - Output (2D)  ||  TAB: toggle axis labels | [0..7]: select input | up/down/left/right: change input")

        self.iren.Start()

    def keypress(self, obj=None, event=None):
        key = obj.GetKeySym()
        if key == "Tab":
            if self.chart.GetAxis(vtk.vtkAxis.BOTTOM).GetBehavior() == vtk.vtkAxis.CUSTOM:
                self.chart.GetAxis(vtk.vtkAxis.BOTTOM).SetBehavior(vtk.vtkAxis.AUTO)
            else:
                self.chart.GetAxis(vtk.vtkAxis.BOTTOM).SetBehavior(vtk.vtkAxis.CUSTOM)
                self.Ticks.Reset()
                self.Labels.Reset()
                for name in limioptic.textArray:
                    self.Ticks.InsertNextValue(name[0])
                    self.Labels.InsertNextValue(name[1])
                self.chart.GetAxis(vtk.vtkAxis.BOTTOM).SetTickPositions(self.Ticks)
                self.chart.GetAxis(vtk.vtkAxis.BOTTOM).SetTickLabels(self.Labels)
            self.render = True
        elif key in (str(i) for i in range(8)):
            self.activeInput = int(key)
            self.parent.setWindowTitle("input control (INPUT[{}])".format(self.activeInput))
        elif key in ("Up", "Down", "Left", "Right"):
            self.threadlock.acquire()
            if key == "Up":
                self.parent.min[self.activeInput].setValue(self.parent.min[self.activeInput].value() + .01)
            if key == "Down":
                self.parent.min[self.activeInput].setValue(self.parent.min[self.activeInput].value() - .01)
            if key == "Right":
                self.parent.input[self.activeInput].setValue(self.parent.input[self.activeInput].value() + .0001)
            if key == "Left":
                self.parent.input[self.activeInput].setValue(self.parent.input[self.activeInput].value() - .0001)
            self.threadlock.release()
        # elif key == "space":
        #    global SCREENSHOTNUMBER
        #    now = time.localtime()
        #    self.w2iFilter.Modified()
        #    self.writer.SetFileName("screenshot_{} ({}.{}.{} at {}.{}).png".format(
        #        SCREENSHOTNUMBER,
        #        now[0],
        #        now[1],
        #        now[2],
        #        now[3],
        #        now[4]))
        #    self.writer.Write()
        #    SCREENSHOTNUMBER += 1

    def animate(self, obj=None, event=None):
        # if self.render:
        # self.render = False
        self.threadlock.acquire()
        self.view.Render()
        self.threadlock.release()
        if self.update:
            self.update = False
            self.neu()

    def neu(self, xxx_todo_changeme2=(None, None, None, None, None)):
        (xi, yi, zi, segs, parts) = xxx_todo_changeme2
        self.threadlock.acquire()
        if (myapp.menu_plot_geo.isChecked()):
            limioptic.geo_s.Reset()
            limioptic.geo_y.Reset()
            limioptic.s = 0.
        else:
            limioptic.s = -1.

        try:
            if xi is None:
                (xi, yi, zi, segs, parts) = self.parent.calculate()
        except Exception as e:
            print("\n\ntraceback:", "\n===============================\n", e, "\n===============================\n\n")
            self.threadlock.release()
            return
        iele = limioptic.GetTrajectory(0, 7)

        # erzeuge wertetabelle
        if (plotx):
            for part in range(parts):
                self.arrX[part].Reset()
                for seg in range(segs):
                    self.arrX[part].InsertNextValue(xi[part][seg])
        if (ploty):
            for part in range(parts):
                self.arrY[part].Reset()
                for seg in range(segs):
                    self.arrY[part].InsertNextValue(yi[part][seg])

        self.arrZ.Reset()
        for i in range(segs):
            self.arrZ.InsertNextValue(zi[i])

        # Marker setzen
        if myapp.menu_plot_marker.isChecked():
            self.linelist = []
            if (len(iele) > 0):
                for i in range(1, len(iele)):
                    if (iele[i] != iele[i - 1]):
                        self.linelist.append(zi[i - 1])
                self.linelist.append(zi[len(iele) - 1])

            for i in range(1, len(self.linelist)):
                try:                                                        # error falls neues element hinzugefuegt wurde
                    self.markersX[i - 1].Reset()
                    self.markersX[i - 1].InsertNextValue(self.linelist[i - 1])
                    self.markersX[i - 1].InsertNextValue(self.linelist[i - 1])
                except:                                                     # neu initialisieren
                    # self.iren.GetRenderWindow().Finalize()
                    # self.iren.TerminateApp()
                    # self.parent.closeit()
                    print("you must close the output-window and rerender (Ctrl+G)")
                    msg = QtWidgets.QMessageBox()
                    msg.setText("you must close the output-window and rerender\n(Ctrl+G)")
                    msg.exec_()
                    self.threadlock.release()
                    return
                    self.markersX.append(vtk.vtkFloatArray())
                    self.markersX[i - 1].SetName("Marker {}".format(i))
                    self.markersX[i - 1].InsertNextValue(self.linelist[i - 1])
                    self.markersX[i - 1].InsertNextValue(self.linelist[i - 1])
                    self.markertable.AddColumn(self.markersX[i - 1])
            self.markertable.Modified()

        self.table.Modified()
        if (myapp.menu_plot_geo.isChecked()):
            limioptic.geolines.Modified()
        if myapp.menu_plot_bg.isChecked():
            self.view.GetRenderer().SetBackground(0, 0, 0)
        else:
            self.view.GetRenderer().SetBackground(1, 1, 1)

        if self.chart.GetAxis(vtk.vtkAxis.BOTTOM).GetBehavior() == vtk.vtkAxis.CUSTOM:
            self.Ticks.Reset()
            self.Labels.Reset()
            for name in limioptic.textArray:
                self.Ticks.InsertNextValue(name[0])
                self.Labels.InsertNextValue(name[1])
            self.chart.GetAxis(vtk.vtkAxis.BOTTOM).SetTickPositions(self.Ticks)
            self.chart.GetAxis(vtk.vtkAxis.BOTTOM).SetTickLabels(self.Labels)

        self.threadlock.release()

        self.render = True


#################################################
#################################################
class CQtLimioptic(QtWidgets.QMainWindow):
    """ Hier wird das Hauptfenster definiert in dem die Beamline eingegeben werden kann """
    def __init__(self):
        QtWidgets.QMainWindow.__init__(self)
        self.setGeometry(
            screen.width() - 350,
            30,
            350 - 10,
            screen.height() - 40)
        self.setWindowTitle('LIMIOPTIC')

        # enable the staus bar
        self.statusBar().showMessage('')

        # erzeuge Menu->File->Load Action
        menu_file_load = QtWidgets.QAction('Open', self)
        menu_file_load.setShortcut('Ctrl+O')
        menu_file_load.setStatusTip('Load Data from File')
        menu_file_load.triggered.connect(self.LoadFile)

        menu_file_loadauto = QtWidgets.QAction('Open _save.lim', self)
        menu_file_loadauto.setShortcut('Ctrl+-')
        menu_file_loadauto.setStatusTip('Load last autosave')
        menu_file_loadauto.triggered.connect(self.LoadAutosave)

        # erzeuge Menu->File->SaveAs Action
        menu_file_saveas = QtWidgets.QAction('Save as..', self)
        menu_file_saveas.setStatusTip('Save File as')
        menu_file_saveas.setShortcut('Ctrl+Alt+S')
        menu_file_saveas.triggered.connect(self.SaveFileAs)

        menu_file_save = QtWidgets.QAction('Save', self)
        menu_file_save.setStatusTip('Save File')
        menu_file_save.setShortcut('Ctrl+S')
        menu_file_save.triggered.connect(self.SaveFile)

        # erzeuge Menu->File->Exit Action
        menu_file_exit = QtWidgets.QAction(QtGui.QIcon('icon/exit.png'), 'Exit', self)
        menu_file_exit.setShortcut('Ctrl+Q')
        menu_file_exit.setStatusTip('Exit Program')
        menu_file_exit.triggered.connect(self.close)

        # translate 2d
        menu_translate_2d = QtWidgets.QAction('2D (VTK)', self)
        menu_translate_2d.setShortcut('Ctrl+G')
        menu_translate_2d.setStatusTip('Translate Text to GUI 2D')
        # if not vtk:
        menu_translate_2d.setEnabled(False)
        menu_translate_2d.triggered.connect(self.plot2d)
        # translate 3d
        menu_translate_3d = QtWidgets.QAction('3D (VTK)', self)
        menu_translate_3d.setShortcut('Ctrl+H')
        menu_translate_3d.setStatusTip('Translate Text to GUI 3D')
        # if not vtk:
        menu_translate_3d.setEnabled(False)
        menu_translate_3d.triggered.connect(self.plot3d)
        # translate qt 2d
        menu_translate_qt = QtWidgets.QAction('Beam trajectories', self)
        menu_translate_qt.setShortcut('Ctrl+F')
        menu_translate_qt.setStatusTip('Translate Text to GUI 2D')
        if not pyqtgraph:
            menu_translate_qt.setEnabled(False)
        menu_translate_qt.triggered.connect(self.plotqt)
        
        menu_plot_beamprofiles = QtWidgets.QAction('Beam profiles', self)
        # menu_plot_emittance.setShortcut('Ctrl+F')
        # menu_plot_emittance.setStatusTip('Translate Text to GUI 2D')
        menu_plot_beamprofiles.triggered.connect(self.beamprofile)

        menu_plot_emittance = QtWidgets.QAction('Beam emittance', self)
        # menu_plot_emittance.setShortcut('Ctrl+F')
        # menu_plot_emittance.setStatusTip('Translate Text to GUI 2D')
        menu_plot_emittance.triggered.connect(self.emittance)

        # PLOT #
        # Plot markers
        self.menu_plot_marker = QtWidgets.QAction('marker', self)
        self.menu_plot_marker.setCheckable(True)
        self.menu_plot_marker.setChecked(True)
#               self.connect(self.menu_plot_marker, QtCore.SIGNAL('triggered()'),)
        # Plot geo
        self.menu_plot_geo = QtWidgets.QAction('geometry', self)
        self.menu_plot_geo.setCheckable(True)
        self.menu_plot_geo.setChecked(True)
        # Plot x
        self.menu_plot_x = QtWidgets.QAction('x', self)
        self.menu_plot_x.setCheckable(True)
        self.menu_plot_x.setChecked(True)
        self.menu_plot_x.triggered.connect(self.setxy)
        # Plot y
        self.menu_plot_y = QtWidgets.QAction('y', self)
        self.menu_plot_y.setCheckable(True)
        self.menu_plot_y.setChecked(True)
        self.menu_plot_y.triggered.connect(self.setxy)
        # Background
        self.menu_plot_bg = QtWidgets.QAction('black background', self)
        self.menu_plot_bg.setCheckable(True)
        self.menu_plot_bg.setChecked(False)
        # smoothing
        self.menu_output_smoothing = QtWidgets.QAction("line smoothing", self)
        self.menu_output_smoothing.setStatusTip('draw smoother lines (apply before rendering)')
        self.menu_output_smoothing.setCheckable(True)
        self.menu_output_smoothing.setChecked(True)
        # split view
        self.menu_plot_splitview = QtWidgets.QAction("split view", self)
        self.menu_plot_splitview.setStatusTip('seperate x, y')
        self.menu_plot_splitview.setCheckable(True)
        self.menu_plot_splitview.setChecked(False)

        # INSERT #
        # INPUT
        self.menu_insert_input = QtWidgets.QAction('INPUT[]', self)
        self.menu_insert_input.setShortcut('Ctrl+Shift+I')
        self.menu_insert_input.setStatusTip('insert INPUT[]')
        self.menu_insert_input.triggered.connect(self.InsertINPUT)

        # NAME
        self.menu_insert_name = QtWidgets.QAction('Name', self)
        self.menu_insert_name.setShortcut('Ctrl+Shift+N')
        self.menu_insert_name.setStatusTip('')
        self.menu_insert_name.triggered.connect(self.InsertName)

        # Source
        self.menu_insert_source = QtWidgets.QAction('Source', self)
        self.menu_insert_source.setStatusTip('Use input file')
        self.menu_insert_source.triggered.connect(self.InsertSource)
        # Particle
        self.menu_insert_particle = QtWidgets.QAction('Particle', self)
        self.menu_insert_particle.setStatusTip('Insert a particle')
        self.menu_insert_particle.triggered.connect(self.InsertParticle)
        # Beam
        self.menu_insert_beam = QtWidgets.QAction('Beam', self)
        self.menu_insert_beam.setShortcut('Ctrl+Shift+B')
        self.menu_insert_beam.setStatusTip('Insert a particle beam')
        self.menu_insert_beam.triggered.connect(self.InsertBeam)
        # Beam X
        self.menu_insert_beamx = QtWidgets.QAction('BeamX (3d)', self)
        self.menu_insert_beamx.setStatusTip('Insert a simple 3d beam')
        self.menu_insert_beamx.triggered.connect(self.InsertBeamX)
        # Beam 3D
        self.menu_insert_beam3d = QtWidgets.QAction('Beam3D (3d)', self)
        self.menu_insert_beam3d.setStatusTip('Insert a NICE 3d beam')
        self.menu_insert_beam3d.triggered.connect(self.InsertBeam3d)
        # Beam Gauss
        self.menu_insert_rgauss = QtWidgets.QAction('BeamRandomGauss (3d)', self)
        self.menu_insert_rgauss.setStatusTip('Insert a gaussian beam')
        self.menu_insert_rgauss.triggered.connect(self.InsertRGauss)
        # Gauss Beam
        self.menu_insert_gaussbeam = QtWidgets.QAction('GaussBeam (3d)', self)
        self.menu_insert_gaussbeam.setStatusTip('Insert a gaussian beam (more options)')
        self.menu_insert_gaussbeam.triggered.connect(self.InsertGaussBeam)

        # Slit
        self.menu_insert_slit = QtWidgets.QAction('Slit', self)
        self.menu_insert_slit.setStatusTip('Insert a slit')
        self.menu_insert_slit.setShortcut('Ctrl+Shift+S')
        self.menu_insert_slit.triggered.connect(self.InsertSlit)
        # Aperture
        self.menu_insert_aperture = QtWidgets.QAction('Aperture', self)
        self.menu_insert_aperture.setStatusTip('Insert a aperture')
        # self.menu_insert_aperture.setShortcut('Ctrl+Shift+A')
        self.menu_insert_aperture.triggered.connect(self.InsertAperture)
        # BeamProfile
        self.menu_insert_bpm = QtWidgets.QAction('Beam Profile Monitor', self)
        self.menu_insert_bpm.setStatusTip('Insert a BPM')
        self.menu_insert_bpm.triggered.connect(self.InsertBPM)
        # Modify Emittance
        self.menu_insert_modifyemittance = QtWidgets.QAction('Modify Emittance', self)
        self.menu_insert_modifyemittance.triggered.connect(self.InsertModifyEmittance)
        # Change Parameters
        self.menu_insert_chgparams = QtWidgets.QAction('ChangeBeamParameters', self)
        self.menu_insert_chgparams.triggered.connect(self.InsertChgParams)
        # Foil
        self.menu_insert_foil = QtWidgets.QAction('Foil', self)
        self.menu_insert_foil.triggered.connect(self.InsertFoil)
        # Waist
        self.menu_insert_waist = QtWidgets.QAction('Waist', self)
        self.menu_insert_waist.triggered.connect(self.InsertWaist)
        # general 6x6 Matrix
        self.menu_insert_matrix = QtWidgets.QAction('General 6x6 matrix', self)
        self.menu_insert_matrix.setStatusTip('Insert a general 6x6 transfer matrix')
        self.menu_insert_matrix.triggered.connect(self.InsertMatrix)
        # Drift
        self.menu_insert_drift = QtWidgets.QAction('Drift', self)
        self.menu_insert_drift.setShortcut('Ctrl+Shift+D')
        self.menu_insert_drift.setStatusTip('Insert a drift')
        self.menu_insert_drift.triggered.connect(self.InsertDrift)
        # ESD
        self.menu_insert_esd = QtWidgets.QAction('ESD', self)
        self.menu_insert_esd.setShortcut('Ctrl+Shift+E')
        self.menu_insert_esd.setStatusTip('Insert an ESD')
        self.menu_insert_esd.triggered.connect(self.InsertESD)
        # MSA
        self.menu_insert_hdm = QtWidgets.QAction('MSA', self)
        self.menu_insert_hdm.setShortcut('Ctrl+Shift+M')
        self.menu_insert_hdm.setStatusTip('Insert a homogeneous deflecting magnet')
        self.menu_insert_hdm.triggered.connect(self.InsertHomDeflectingMagnet)
        # MSA_Y
        self.menu_insert_hdmy = QtWidgets.QAction('MSA_Y', self)
        self.menu_insert_hdmy.setStatusTip('Insert a homogeneous deflecting magnet which bends to the Y axis')
        self.menu_insert_hdmy.triggered.connect(self.InsertHomDeflectingMagnetY)
        # Quadrupol
        self.menu_insert_quadrupol = QtWidgets.QAction('Quadrupol', self)
        self.menu_insert_quadrupol.setShortcut('Ctrl+Shift+Q')
        self.menu_insert_quadrupol.setStatusTip('Open Dialog : Insert a quadrupol')
        self.menu_insert_quadrupol.triggered.connect(self.InsertQuadrupol)
        # ThinLens
        self.menu_insert_ThinLens = QtWidgets.QAction('Einzel lens', self)
        self.menu_insert_ThinLens.setShortcut('Ctrl+Shift+L')
        self.menu_insert_ThinLens.setStatusTip('Open Dialog : Insert a thin lense')
        self.menu_insert_ThinLens.triggered.connect(self.InsertThinLens)

        # SO110-EL
        self.menu_insert_so110el = QtWidgets.QAction('AMS SO110-EL', self)
        self.menu_insert_so110el.setStatusTip('Cologne AMS SO110 einzel lens')
        self.menu_insert_so110el.triggered.connect(self.InsertSO110EL)
        # FN-EL
        self.menu_insert_fnel = QtWidgets.QAction('FN EL', self)
        self.menu_insert_fnel.setStatusTip('HVEC-FN 7 einzel lens')
        self.menu_insert_fnel.triggered.connect(self.InsertFNEL)
        # BI-EL
        self.menu_insert_biel = QtWidgets.QAction('AMS BI-EL', self)
        self.menu_insert_biel.setStatusTip('Cologne AMS BI einzel lens')
        self.menu_insert_biel.triggered.connect(self.InsertBIEL)
        # AMSQPT
        self.menu_insert_amsqpt = QtWidgets.QAction('AMS QPT', self)
        self.menu_insert_amsqpt.setStatusTip('Cologne AMS Quadrupole Triplet')
        self.menu_insert_amsqpt.triggered.connect(self.InsertAMSQPT)
        # VBFN
        self.menu_insert_vbfn = QtWidgets.QAction('FN VB', self)
        self.menu_insert_vbfn.setShortcut('Ctrl+Shift+V')
        self.menu_insert_vbfn.setStatusTip('Insert Cologne FN preacceleration segment')
        self.menu_insert_vbfn.triggered.connect(self.InsertVBFN)
        # FNACC
        self.menu_insert_fnacc = QtWidgets.QAction('FN acceleration tube', self)
        # self.menu_insert_fnacc.setShortcut('Ctrl+Shift+F')
        self.menu_insert_fnacc.setStatusTip('Insert HVEC-FN 7 acceleration tube')
        self.menu_insert_fnacc.triggered.connect(self.InsertFNAcc)
        # FNACCNeu
        self.menu_insert_fnaccneu = QtWidgets.QAction('FN acceleration tube', self)
        self.menu_insert_fnaccneu.setShortcut('Ctrl+Shift+F')
        self.menu_insert_fnaccneu.setStatusTip('Insert HVEC-FN 7 acceleration tube')
        self.menu_insert_fnaccneu.triggered.connect(self.InsertFNAccNeu)
        # AMS-Acc
        self.menu_insert_amsacc = QtWidgets.QAction('AMS acceleration tube', self)
        self.menu_insert_amsacc.setStatusTip('Insert Cologne AMS acceleration tube')
        self.menu_insert_amsacc.setShortcut('Ctrl+Shift+A')
        self.menu_insert_amsacc.triggered.connect(self.InsertAMSAcc)

        # TOOLS #
        # dat
        self.menu_output_file = QtWidgets.QAction('-> beam.dat', self)
        self.menu_output_file.setStatusTip('write beams to beam.dat file.')
        self.menu_output_file.setCheckable(True)
        self.menu_output_file.setChecked(False)
        self.menu_output_file.triggered.connect(self.todat)
        # AMS-Spicker
        # self.menu_output_spicker = QtWidgets.QAction("AMS Spicker", self)
        # self.menu_output_spicker.setStatusTip('Cologne AMS Spicker')
        # self.menu_output_spicker.setCheckable(False)
        # self.menu_output_spicker.triggered.connect(self.spicker)
        # plotEmittance
        self.menu_output_emittance = QtWidgets.QAction("Plot Beamprofile", self)
        self.menu_output_emittance.setStatusTip('Plot Beamprofile')
        self.menu_output_emittance.setCheckable(False)
        self.menu_output_emittance.triggered.connect(self.emittance)
        # calculator
        self.menu_output_calc = QtWidgets.QAction("Calculator", self)
        self.menu_output_calc.setStatusTip('calculator')
        self.menu_output_calc.setCheckable(False)
        self.menu_output_calc.triggered.connect(self.calc)

        # Interface
        self.menu_setcom1 = QtWidgets.QAction('COM1', self)
        self.menu_setcom1.setStatusTip('search interface on com1')
        self.menu_setcom1.triggered.connect(self.setcom1)
        self.menu_setcom2 = QtWidgets.QAction('COM2', self)
        self.menu_setcom2.setStatusTip('search interface on com2')
        self.menu_setcom2.triggered.connect(self.setcom2)
        self.menu_setcom3 = QtWidgets.QAction('COM3', self)
        self.menu_setcom3.setStatusTip('search interface on com3')
        self.menu_setcom3.triggered.connect(self.setcom3)
        self.menu_setcom4 = QtWidgets.QAction('COM4', self)
        self.menu_setcom4.setStatusTip('search interface on com4')
        self.menu_setcom4.triggered.connect(self.setcom4)
        self.menu_setcom5 = QtWidgets.QAction('COM5', self)
        self.menu_setcom5.setStatusTip('search interface on com5')
        self.menu_setcom5.triggered.connect(self.setcom5)

        # About
        self.menu_about = QtWidgets.QAction("About", self)
        self.menu_about.triggered.connect(self.About)
        self.menu_licence = QtWidgets.QAction("Licence", self)
        self.menu_licence.triggered.connect(self.Licence)

        # Menubar
        menubar = self.menuBar()

        # FILE
        menu_file = menubar.addMenu('File')
        menu_file.addAction(menu_file_load)

        menu_file.addAction(menu_file_saveas)
        menu_file.addAction(menu_file_save)
        menu_file.addAction(menu_file_exit)

        # TRANSLATE
        menu_translate = menubar.addMenu('Plot')
        # menu_translate.addAction(menu_translate_2d)
        menu_translate.addAction(menu_translate_qt)
        menu_translate.addAction(menu_plot_beamprofiles)
        menu_translate.addAction(menu_plot_emittance)
        # menu_translate.addAction(menu_translate_3d)

        # PLOT
        menu_plot = menubar.addMenu('Options')
        menu_plot.addAction(self.menu_plot_marker)
        menu_plot.addAction(self.menu_plot_geo)
        menu_plot.addAction(self.menu_plot_x)
        menu_plot.addAction(self.menu_plot_y)
        menu_plot.addAction(self.menu_output_smoothing)
        menu_plot.addAction(self.menu_plot_bg)
        menu_plot.addAction(self.menu_plot_splitview)

        # INSERT
        menu_insert = menubar.addMenu('Insert')
        menu_insert.addAction(self.menu_insert_input)
        menu_insert.addAction(self.menu_insert_name)
        menu_insert.addSeparator()
        menu_insert.addAction(self.menu_insert_source)
        menu_insert.addAction(self.menu_insert_particle)
        menu_insert.addAction(self.menu_insert_beam)
        menu_insert.addAction(self.menu_insert_beamx)
        menu_insert.addAction(self.menu_insert_beam3d)
        menu_insert.addAction(self.menu_insert_rgauss)
        menu_insert.addAction(self.menu_insert_gaussbeam)
        menu_insert.addSeparator()
        menu_insert.addAction(self.menu_insert_drift)
        menu_insert.addAction(self.menu_insert_esd)
        menu_insert.addAction(self.menu_insert_hdm)
        menu_insert.addAction(self.menu_insert_hdmy)
        menu_insert.addAction(self.menu_insert_quadrupol)
        menu_insert.addAction(self.menu_insert_ThinLens)
        menu_insert.addAction(self.menu_insert_slit)
        menu_insert.addAction(self.menu_insert_aperture)
        menu_insert.addAction(self.menu_insert_bpm)
        menu_insert.addAction(self.menu_insert_modifyemittance)
        menu_insert.addAction(self.menu_insert_chgparams)
        menu_insert.addAction(self.menu_insert_foil)
        menu_insert.addAction(self.menu_insert_matrix)
        menu_insert.addAction(self.menu_insert_waist)
        menu_insert.addSeparator()
        menu_insert.addAction(self.menu_insert_so110el)
        menu_insert.addAction(self.menu_insert_biel)
        menu_insert.addAction(self.menu_insert_amsqpt)
        menu_insert.addAction(self.menu_insert_amsacc)
        menu_insert.addAction(self.menu_insert_fnel)
        menu_insert.addAction(self.menu_insert_vbfn)
        menu_insert.addAction(self.menu_insert_fnaccneu)

        # TOOLS
        menu_output = menubar.addMenu('Tools')
        menu_output.addAction(self.menu_output_file)
        # menu_output.addAction(self.menu_output_spicker)
        # menu_output.addAction(self.menu_output_calc)
        # menu_output.addAction(self.menu_output_emittance)

        # INTERFACE
        """
        menu_interface = menubar.addMenu('# interface')
        menu_interface.addAction(self.menu_setcom1)
        menu_interface.addAction(self.menu_setcom2)
        menu_interface.addAction(self.menu_setcom3)
        menu_interface.addAction(self.menu_setcom4)
        menu_interface.addAction(self.menu_setcom5)
        """

        # ABOUT
        menu_aboutbar = menubar.addMenu('About')
        menu_aboutbar.addAction(self.menu_about)
        menu_aboutbar.addAction(self.menu_licence)

        # GUI definieren
        self.main_frame = QtWidgets.QWidget()
        self.textedit = myedit(self)
        self.textedit.setLineWrapMode(QtWidgets.QTextEdit.NoWrap)
        self.textedit.setTabStopWidth(30)
        self.setAcceptDrops(True)
        self.highlighter = syntax.PythonHighlighter(self.textedit.document())

        # self.helpwidget = QtGui.QTextEdit()

        hbox1 = QtWidgets.QVBoxLayout()
        hbox1.addWidget(self.textedit)
        # hbox1.addWidget(self.helpwidget)
        self.main_frame.setLayout(hbox1)
        self.setCentralWidget(self.main_frame)

        self.textedit.textChanged.connect(self.setunsaved)
        self.textedit.selectionChanged.connect(self.help)
        # print "started"

        try:
            self.LoadAutosave()
            print("last autosave was restored")
        except:
            print("last savefile was not found")

        # UPDATE
        updatethread = threading.Thread(target=self.update, args=())
        updatethread.start()

    def help(self, selection=""):
        if (not selection):
            selection = self.textedit.textCursor().selectedText()
        if selection in dir(helper):
            helptext = eval("helper.%s" % (selection))
            print("\nHELP:\n------------------------------------------")
            print(helptext)
            print("------------------------------------------")
            self.setToolTip(helptext)
        elif selection in dir(limioptic):
            helptext = selection + " takes:   " + (", ".join(eval("limioptic.%s.__code__.co_varnames" % (selection))))
            print("\nHELP:\n------------------------------------------")
            print(helptext)
            print("------------------------------------------")
            self.setToolTip(helptext)
        else:
            self.setToolTip("double-click on function for help!")

    def closeEvent(self, event):
        print("saving..", end=' ')
        self.SaveNeu(backup_file)
        print("\rsaved to %s" % (backup_file))

    def update(self):
        """ Check uebers Internet ob Updates verfuegbar sind """
        print("checking for updates..", end=' ')
        try:
            a = urllib.request.urlopen("https://raw.githubusercontent.com/alexander-stolz/limioptic/master/LIMIOPTIC/version")
            ver = a.read()
            ver = "".join(map(chr, ver)).strip()
            a.close()
            if (ver == VERSION):
                print("this is the newest version\n")
            else:
                print("there is a newer version online. run \"__update.bat\" to update\n")
                print("Your version:", VERSION, "Latest Version:", ver, "\n")
        except Exception as e:
            print("update check failed\n", e)

#############################

# set com #
    def setcom1(self):
        global PORT
        PORT = "COM1"

    def setcom2(self):
        global PORT
        PORT = "COM2"

    def setcom3(self):
        global PORT
        PORT = "COM3"

    def setcom4(self):
        global PORT
        PORT = "COM4"

    def setcom5(self):
        global PORT
        PORT = "COM5"

##############################
    def setxy(self):
        """ Nur x, nur y, oder beides zeichnen. """
        global xy, plotx, ploty
        xy = 0
        plotx = False
        ploty = False
        if (self.menu_plot_x.isChecked()):
            xy += 1
            plotx = True
        if (self.menu_plot_y.isChecked()):
            xy += 1
            ploty = True

#############################
    # def spicker(self):
    #     """ Der AMS-Spicker """
    #     self.swidget = ams_spicker.Spicker()
    #     self.swidget.show()

#############################
    def calc(self):
        os.system("python function_calculator.py")

#############################
    def emittance(self):
        # plotEmittance.start(limioptic.PROFILEINDEX + 1)
        nBPM = Counter(myapp.textedit.toPlainText().split("\n"))["BeamProfile()"]
        if nBPM:
            os.popen(f"python plotBeamEmittance.pyw {nBPM}")
        else:
            print("ERROR: Add at least one BeamProfile()")
    
    def beamprofile(self):
        # plotEmittance.start(limioptic.PROFILEINDEX + 1)
        nBPM = Counter(myapp.textedit.toPlainText().split("\n"))["BeamProfile()"]
        if nBPM:
            os.popen(f"python plotBeamProfile.pyw {nBPM}")
        else:
            print("ERROR: Add at least one BeamProfile()")

#############################
    def todat(self):
        """ Ausgabe in Datei """
        if (self.menu_output_file.isChecked()):
            # if (RUNNING2D):
            msg = QtWidgets.QMessageBox()
            msg.setText("Please close plotwindow and press Ctrl+F to rerender!\n*.dat files will be produced every time the plot-window opens.")
            msg.exec_()
            # else:
            #     msg = QtWidgets.QMessageBox()
            #     msg.setText("Press Ctrl+G to produce *.dat files!")
            #     msg.exec_()

# save and load #
    def setunsaved(self):
        try:
            self.setWindowTitle("LIMIOPTIC -     changed     - {}".format(self.FileName))
        except:
            self.setWindowTitle("LIMIOPTIC - unsaved")

    def LoadFile(self, filename=None):
        if not filename:
            self.FileName = (QtWidgets.QFileDialog.getOpenFileName(self, "Open file", ".", "*.lim2;;*.lim")[0])
        else:
            self.FileName = filename
        if (self.FileName.endswith("lim")):
            myfile = open(self.FileName, 'r')
            self.textedit.setText(myfile.read())
            myfile.close()

            try:
                global NumberOfInputs, BEZEICHNUNGEN, INPUT
                i = -1
                for line in open(self.FileName.split(".", 1)[0] + ".var"):
                    i += 1
                    BEZEICHNUNGEN[i] = line.split(" = ")[0]
                    INPUT[i] = float(line.split(" = ")[1])
                NumberOfInputs = i + 1
                print("the last values were restored")
            except:
                print("could not load variables")

        elif (self.FileName.endswith("lim2")):
            (text, INPUT, BEZEICHNUNGEN) = pickle.load(open(self.FileName, "rb"))
            self.textedit.setText(text)

        self.setWindowTitle("LIMIOPTIC  -  {}".format(self.FileName))
        print("{} was loaded".format(self.FileName))

    def LoadAutosave(self):
        global NumberOfInputs, BEZEICHNUNGEN, INPUT

        try:
            if backup_file.endswith("lim2"):
                (text, INPUT, BEZEICHNUNGEN) = pickle.load(open(backup_file, "rb"))
                self.textedit.setText(text)
                return 0
            else:
                myfile = open(backup_file + ".lim", "r")
                title  = "LIMIOPTIC  -  _save.lim"
                self.textedit.setText(myfile.read())
                myfile.close()
                self.setWindowTitle(title)
        except:
                print("_save.lim not found!")

        try:
            i = -1
            for line in open(backup_file + ".var", "r"):
                i += 1
                BEZEICHNUNGEN[i] = line.split(" = ")[0]
                INPUT[i] = float(line.split(" = ")[1])
            NumberOfInputs = i + 1
        except:
            print("_save.var not found")

    def SaveAlt(self):
        myfile = open(self.FileName, "w")
        myfile.write(str(self.textedit.toPlainText()))
        myfile.close()
        self.setWindowTitle("LIMIOPTIC")
        print("saved to", self.FileName)
        self.setWindowTitle("LIMIOPTIC  -  {}".format(self.FileName))

        myfile = open(self.FileName.split(".", 1)[0] + ".var", "w")
        for i in range(NumberOfInputs):
                print("{} = {}".format(BEZEICHNUNGEN[i], INPUT[i]), file=myfile)
        myfile.close()

    def SaveNeu(self, filename=None):
        if filename:
            pickle.dump((str(self.textedit.toPlainText()), INPUT, BEZEICHNUNGEN), open(filename, "wb"))
        else:
            pickle.dump((str(self.textedit.toPlainText()), INPUT, BEZEICHNUNGEN), open(self.FileName, "wb"))

    def SaveFileAs(self):
            self.FileName = QtWidgets.QFileDialog.getSaveFileName(self, "Save file", ".", "*.lim2;*.lim")[0]
            if (self.FileName.endswith("lim")):
                self.SaveAlt()
            elif (self.FileName.endswith("lim2")):
                self.SaveNeu()

    def SaveFile(self):
        try:
            if (self.FileName.endswith("lim")):
                self.SaveAlt()
            elif (self.FileName.endswith("lim2")):
                self.SaveNeu()
        except:
            print("Noch kein Filename definiert!")
            self.SaveFileAs()

# plot #
    def plot2d(self):
        global RUNNING2D, RUNNING
        self.SaveNeu(backup_file)

        if not RUNNING:
            print("Sicherungsdatei: ", backup_file)
            self.inputwindow2d = inputcontrol("2d")
            RUNNING2D = RUNNING = True
            self.inputwindow2d.exec_()
            RUNNING2D = RUNNING = False
        elif RUNNING2D:
            self.inputwindow2d.plotwindow.neu()

    def plotqt(self):
        global RUNNINGQT, RUNNING
        self.SaveNeu(backup_file)

        if not RUNNING:
            print("Sicherungsdatei: ", backup_file)
            self.inputwindowqt = inputcontrol("qt")
            RUNNINGQT = RUNNING = True
            self.inputwindowqt.exec_()
            RUNNINGQT = RUNNING = False
        elif RUNNINGQT:
            self.inputwindowqt.plotwindow.update()

    def plot3d(self):
        global RUNNING3D, RUNNING
        self.SaveNeu(backup_file)

        if not RUNNING:
            print("Sicherungsdatei: ", backup_file)
            self.inputwindow3d = inputcontrol("3d")
            RUNNING3D = RUNNING = True
            self.inputwindow3d.exec_()
            RUNNING3D = RUNNING = False
        elif RUNNING3D:
            self.inputwindow3d.plotwindow.neu()

# insert #
    def InsertINPUT(self):
        self.textedit.textCursor().insertText("INPUT[ ]")
        self.textedit.moveCursor(QtGui.QTextCursor.Left)
        self.textedit.moveCursor(QtGui.QTextCursor.Left, QtGui.QTextCursor.KeepAnchor)
        print("define slider number")

    def InsertName(self):
        self.textedit.textCursor().insertText("Name(' ')")
        self.textedit.moveCursor(QtGui.QTextCursor.Left)
        self.textedit.moveCursor(QtGui.QTextCursor.Left)
        self.textedit.moveCursor(QtGui.QTextCursor.Left, QtGui.QTextCursor.KeepAnchor)
        print("define label")

    def InsertParticle(self):
        self.textedit.textCursor().insertText("AddParticle(4, 15, 4, 15, 0, 0)\t\t\t# (x, x\', y, y\', dk, dm)\n")
        self.help("AddParticle")

    def InsertSource(self, _filename=None):
        global SourceObj
        if not _filename:
            _filename = QtWidgets.QFileDialog.getOpenFileName(self, "Open file", ".")[0]
            SourceObj.LoadSource(_filename, filetype=("SRIM" if _filename.endswith("TRANSMIT.txt") else "limioptic"))
        else:
            SourceObj.LoadSource(_filename, filetype=("SRIM" if _filename.endswith("TRANSMIT.txt") else "limioptic"))
        # SourceObj.NormalizeEnergy()
        SourceObj.ShowFits()
        SourceObj.UserInteraction.ChooseFilter()
        # SourceObj.Source = SourceObj.Selection
        # SourceObj.ShowFits()
        print("Source", SourceObj.SourceFile, "loaded")
        self.textedit.textCursor().insertText("Source()\n")
        self.textedit.textCursor().insertText(
            "# ChangeBeamParameters("
            "strag_x={}, "
            "strag_dx={}, "
            "strag_y={}, "
            "strag_dy={}, "
            "strag_k={})\n".format(
                SourceObj.foilparameters["x"],
                SourceObj.foilparameters["x'"],
                SourceObj.foilparameters["y"],
                SourceObj.foilparameters["y'"],
                SourceObj.foilparameters["dk"]))

    def InsertBeam(self):
        self.textedit.textCursor().insertText('############################################\nBeam(4, 15, 4, 15, 0, 0)\t# (xmax, x\'max, ymax, y\'max, dk, dm, delta: 1...360)\n############################################\n\n')
        self.help("Beam")

    def InsertBeamX(self):
        self.textedit.textCursor().insertText('############################################\nBeamX(4,15,4,15,0,0,10)\t# (xmax, x\'max, ymax, y\'max, dk, dm, delta: 1...360)\n############################################\n\n')
        self.help("BeamX")

    def InsertBeam3d(self):
        self.textedit.textCursor().insertText('############################################\nBeam3d(4,15,4,15,0,0,10)\t# (xmax, x\'max, ymax, y\'max, dk, dm, delta_phi: 1...360)\n############################################\n\n')
        self.help("Beam3d")

    def InsertRGauss(self):
        self.textedit.textCursor().insertText('############################################\nBeamRandomGauss(4,15,4,15,0,0,1000)\t# (xmax, x\'max, ymax, y\'max, dk, dm, num)\n############################################\n\n')
        self.help("BeamRandomGauss")

    def InsertGaussBeam(self):
        self.textedit.textCursor().insertText('############################################\nGaussBeam(4, 15, 4, 15)\t# (sigma_x, sigma_x\', sigma_y, sigma_y\', x, x\', y, y\', dk, dm, sigma_k, sigma_m, strag_k, strag_m, number)\n############################################\n\n')
        self.help("GaussBeam")

    def InsertAMSAcc(self):
        self.textedit.textCursor().insertText('AMSAcc(50.e3, 5500.e3, 35.e3, 4)\t# (v_qsnout, v_terminal, v_ext, q)\n\n')
        self.help("AMSAcc")

    def InsertFNAcc(self):
        # self.textedit.textCursor().insertText('AddFNAcc(6000.e3, 100.e3, 5)\t# (v_terminal, v_vorbeschl, q)\n\n')
        self.InsertFNAccNeu()

    def InsertFNAccNeu(self):
        self.textedit.textCursor().insertText('AddFNAccNeu(vt, T0, q, b = 0.57, b1 = -1., b2 = -1., D1 = .088, factor1 = 1., factor2 = 1., beamprofile = False)\n\n')

    def InsertVBFN(self):
        self.textedit.textCursor().insertText('VBFN(extraktion, deltaV, laenge)\t# (v_ext, deltaV, length, [b, b1, b2])\n\n')
        self.help("VBFN")

    def InsertMatrix(self):
        # uebergebe Zeiger auf das TextEdit an den Dialog
        self.dialog = CInsertMatrixDialog(self.textedit)
        self.dialog.exec_()

    def InsertDrift(self):
        self.textedit.textCursor().insertText('Drift(5.)\t\t\t\t\t\t# (length)\n')
        self.help("Drift")

    def InsertSlit(self):
        self.textedit.textCursor().insertText('Slit(0,10,0,10)\t# (x, dx, y, dy)\n')
        self.help("Slit")
    
    def InsertAperture(self):
        self.textedit.textCursor().insertText('Aperture(3)\t# (d)\n')
        self.help("Aperture")

    def InsertBPM(self):
        self.textedit.textCursor().insertText('BeamProfile()\n')
        self.help("BeamProfile")

    def InsertModifyEmittance(self):
        self.textedit.textCursor().insertText('ModifyEmittance(1., 1.)\t# (factor x, factor dx)\n')
        self.help("ModifyEmittance")

    def InsertChgParams(self):
        self.textedit.textCursor().insertText('ChangeBeamParameters(dk=0., dm=0., strag_k=0., strag_m=0.)\t\n')
        self.help("ChangeBeamParameters")

    def InsertFoil(self):
        self.textedit.textCursor().insertText(
            """Foil(
\tdk=0.,\t \t# delta energy in permille
\tstrag_k=0.,\t \t# energy straggling in permille
\tstrag_phi=0.,\t# angular straggling in mrad\t
\tpercentage=.5)\t# percentage of the beam that is affected by energy loss\n""")

    def InsertWaist(self):
        self.textedit.textCursor().insertText('Waist()\n')
        self.help("Waist")

    def InsertESD(self):
        self.textedit.textCursor().insertText('ESD(30., 2., 1.e9)\t# (alpha, r_hor, r_vert, R)\n\n')
        self.help("ESD")

    def InsertHomDeflectingMagnet(self):
        self.textedit.textCursor().insertText('EdgeFocusing(r, beta, K=.45, R)\n')    # Kantenfokussierung
        self.textedit.textCursor().insertText('MSA(r, alpha)\n')    # Magnet
        self.textedit.textCursor().insertText('EdgeFocusing(r, beta, K=.45, R)\n\n')  # Kantenfokussierung
        self.help("MSA")

    def InsertHomDeflectingMagnetY(self):
        self.textedit.textCursor().insertText('EdgeFocusingY(r, beta, K=.45, R)\n')    # Kantenfokussierung
        self.textedit.textCursor().insertText('MSA_Y(r, alpha)\n')    # Magnet
        self.textedit.textCursor().insertText('EdgeFocusingY(r, beta, K=.45, R)\n\n')  # Kantenfokussierung

    def InsertQuadrupol(self):
        self.textedit.textCursor().insertText('QuadrupolRadFoc(k, l, R)\n')   # radial fokussierend
        self.textedit.textCursor().insertText('QuadrupolAxFoc(k, l, R)\n\n')  # axial fokusierend
        self.help("QuadrupolRadFoc")

    def InsertThinLens(self):
        self.textedit.textCursor().insertText('EinzelLens(2.)\t\t\t\t# (f, [R])\n\n')
        self.help("EinzelLens")

    def InsertSO110EL(self):
        self.textedit.textCursor().insertText('AMSSO110EL(vext, vlens)\n\n')

    def InsertBIEL(self):
        self.textedit.textCursor().insertText('AMSBIEL(vext, vlens)\n\n')

    def InsertFNEL(self):
        self.textedit.textCursor().insertText('FNEL(v_ext, v_lens)\n\n')

    def InsertAMSQPT(self):
        self.textedit.textCursor().insertText('AddAMSQPT(gamma2, prozent, astigm, v_terminal, v_ext, q, geo)\n\n')

    def About(self):
        title = "About LIMIOPTIC"
        text = "LIMIOPTIC by Alexander Stolz\nVersion {}\n\n"\
            "Feel free to send me any feedback or suggestions to amstolz@gmail.com.\n\n"\
            "Visit www.limioptic.de for more information.\n\nThanks for using LIMIOPTIC!\n\n".format(VERSION)
        self.dialog = DialogWindow(title, text)

    def Licence(self):
        title = "LIMIOPTIC Licence"
        text = "The software LIMIOPTIC maintained by Alexander Stolz is freely available and distributable. "\
            "However, if you use it for some work whose results are made public, then you have to reference it properly."
        self.dialog = DialogWindow(title, text)


class myedit(QtWidgets.QTextEdit):
    dropped = QtCore.pyqtSignal()

    def __init__(self, type, parent=None):
        super(myedit, self).__init__(parent)
        self.setAcceptDrops(True)

    def dragEnterEvent(self, event):
        if event.mimeData().hasUrls:
            event.accept()
        else:
            event.ignore()

    def dragMoveEvent(self, event):
        if event.mimeData().hasUrls:
            event.setDropAction(QtCore.Qt.CopyAction)
            event.accept()
        else:
            event.ignore()

    def dropEvent(self, event):
        if event.mimeData().hasUrls:
            event.setDropAction(QtCore.Qt.CopyAction)
            event.accept()
            links = []
            for url in event.mimeData().urls():
                links.append(str(url.toLocalFile()))
            self.dropped.emit(links)
            if len(links) == 1:
                if (links[0].endswith(".lim")) or (links[0].endswith(".lim2")):
                    myapp.LoadFile(links[0])
                if links[0].endswith("TRANSMIT.txt") or links[0].endswith(".out"):
                    myapp.InsertSource(links[0])
        else:
            event.ignore()


# Dialoge #
class CInsertParticleDialog(QtWidgets.QDialog):
    """ Wird nicht mehr benoetigt """
    def __init__(self, myarg1):
        QtWidgets.QDialog.__init__(self)
        self.setFixedSize(400, 200)
        self.setWindowTitle('Insert Particle Dialog')

        # Zeiger auf das QTextEdit speichern
        self.parent_textedit = myarg1

        self.insert_syntax_button = QtWidgets.QPushButton('just insert syntax', self)
        self.insert_syntax_button.setGeometry(50, 50, 160, 25)
        self.insert_syntax_button.clicked.connect(self.InsertSyntax)


class DialogWindow(QtWidgets.QDialog):
    def __init__(self, title, text):
        QtWidgets.QDialog.__init__(self)
        self.setWindowFlags(QtCore.Qt.WindowStaysOnTopHint)
        self.setWindowTitle(title)
        self.abouttext = QtWidgets.QTextEdit()
        hbox = QtWidgets.QHBoxLayout()
        hbox.addWidget(self.abouttext)
        self.setLayout(hbox)
        self.abouttext.setText(text)
        self.show()


class CInsertMatrixDialog(QtWidgets.QDialog):
    def __init__(self, myarg1):
        QtWidgets.QDialog.__init__(self)
        self.setFixedSize(200, 100)
        self.setWindowTitle('Insert general 6x6 Matrix Dialog')

        # Zeiger auf das QTextEdit speichern
        self.parent_textedit = myarg1

        self.insert_ThinLens_example_button = QtWidgets.QPushButton('insert unity matrix', self)
        self.insert_ThinLens_example_button.setGeometry(50, 50, 100, 50)
        self.insert_ThinLens_example_button.clicked.connect(self.InsertThinLensExample)

    def InsertThinLensExample(self):
        # nur die Syntax fuer 'AddMatrix' einfuegen
        self.parent_textedit.textCursor().insertText('AddMatrix(n,\n')
        self.parent_textedit.textCursor().insertText('           [1.,0.,0.,0.,0.,0.,\n')
        self.parent_textedit.textCursor().insertText('           0.,1.,0.,0.,0.,0.,\n')
        self.parent_textedit.textCursor().insertText('           0.,0.,1.,0.,0.,0.,\n')
        self.parent_textedit.textCursor().insertText('           0.,0.,0.,1.,0.,0.,\n')
        self.parent_textedit.textCursor().insertText('           0.,0.,0.,0.,1.,0.,\n')
        self.parent_textedit.textCursor().insertText('           0.,0.,0.,0.,0.,1.],\n')
        self.parent_textedit.textCursor().insertText('           length)\n')

################################

VERSION          = open("version", "r").read()
PORT             = "NONE"
INPUT            = []
BEZEICHNUNGEN    = []
OPACITY          = 50
NumberOfInputs   = 8
SCALE3D          = 10.
SUPERSCALE3D     = 10.
SCREENSHOTNUMBER = 0
XCOLOR           = 255, 0, 0, 255
YCOLOR           = 0, 255, 0, 255
LICENCE          = "\nThe software LIMIOPTIC maintained by Alexander Stolz is freely available and distributable. However, if you use it for some work whose results are made public, then you have to reference it properly.\n"

print(LICENCE)

try:
    backup_file = "/".join([os.environ["PROGRAMDATA"], "_save.lim2"])
    _test = open(backup_file + ".test", "w")
    print("writetest", file=_test)
    _test.close()
except Exception as e:
    try:
        backup_file = "/".join([os.path.expanduser("~"), "_save.lim2"])
        _test = open(backup_file + ".test", "w")
        print("writetest", file=_test)
        _test.close()
    except:
        backup_file = "_save.lim2"

for i in range(32):
        INPUT.append(1.)
        BEZEICHNUNGEN.append("")

plotx = True
ploty = True
xy    = 2

RUNNING   = False
RUNNINGQT = False
RUNNING2D = False
RUNNING3D = False


app = QtWidgets.QApplication(sys.argv)

# screen = QtGui.QDesktopWidget().screenGeometry()
screen = QtWidgets.QDesktopWidget().availableGeometry()
# screen = QtGui.QDesktopWidget().desktop()

SourceObj = ImportSource()
myapp = CQtLimioptic()
myapp.setWindowIcon(QtGui.QIcon('logo.png'))

myapp.show()

app.exec_()
del myapp
del app

print("Have a nice day :)")

sys.exit()
