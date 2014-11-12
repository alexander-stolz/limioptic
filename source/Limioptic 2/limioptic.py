#!/usr/bin/env python

"""
Hier werden die Funktionen definiert und an den C++ - Teil uebergeben.
"""

import sys
import ctypes
import array
import math
import random
import vtk
#import pickle
from scipy import optimize
from numpy import linspace

s            = 0.
SOURCE       = []
SOURCEFILES  = []
PROFILEINDEX = -1

geolines = vtk.vtkTable()
geo_s    = vtk.vtkFloatArray()
geo_y    = vtk.vtkFloatArray()
geo_s.SetName("S-Achse")
geo_y.SetName("geometry")

textArray    = []
lastFunction = ""


def AddParticle(x, a, y, b, dk, dm):
    """ Normales Partikel einfuegen """
    global lastFunction
    lastFunction = "AddParticle"
    optic.AddParticle(
        ctypes.c_double(float(x)),
        ctypes.c_double(float(a)), ctypes.c_double(float(y)),
        ctypes.c_double(float(b)), ctypes.c_double(float(dk)),
        ctypes.c_double(float(dm)))
Particle = AddParticle


def AddBeamX(xmax, amax, ymax, bmax, dk, dm, delta):
    """ Einfachen 3d-Strahl einfuegen """
    global lastFunction
    lastFunction = "AddBeamX"
    for degree in xrange(0, 360, int(delta)):
        degtorad = math.radians(degree)

        optic.AddParticle(
            ctypes.c_double(float(xmax * math.cos(degtorad))),
            ctypes.c_double(float(amax * math.sin(degtorad))),
            ctypes.c_double(0.),
            ctypes.c_double(0.),
            ctypes.c_double(float(dk)),
            ctypes.c_double(float(dm)))
        optic.AddParticle(
            ctypes.c_double(0.),
            ctypes.c_double(0.),
            ctypes.c_double(float(ymax * math.cos(degtorad))),
            ctypes.c_double(float(bmax * math.sin(degtorad))),
            ctypes.c_double(float(dk)),
            ctypes.c_double(float(dm)))


def BeamX(xmax, amax, ymax, bmax, dk, dm, num):
    """ Einfachen 3d-Strahl einfuegen """
    global lastFunction
    lastFunction = "BeamX"
    for degree in linspace(0, 360, int(num / 2)):
        degtorad = math.radians(degree)

        optic.AddParticle(
            ctypes.c_double(float(xmax * math.cos(degtorad))),
            ctypes.c_double(float(amax * math.sin(degtorad))),
            ctypes.c_double(0.),
            ctypes.c_double(0.),
            ctypes.c_double(float(dk)),
            ctypes.c_double(float(dm)))
        optic.AddParticle(
            ctypes.c_double(0.),
            ctypes.c_double(0.),
            ctypes.c_double(float(ymax * math.cos(degtorad))),
            ctypes.c_double(float(bmax * math.sin(degtorad))),
            ctypes.c_double(float(dk)),
            ctypes.c_double(float(dm)))


def AddBeam3d(xmax, amax, ymax, bmax, dk, dm, delta):
    """ Komplexen 3d-Strahl einfuegen """
    global lastFunction
    lastFunction = "AddBeam3d"
    for i in xrange(0, 360, int(delta)):
        nu = math.radians(i)
        for j in xrange(0, 360, 10):
            phi = math.radians(j)
            optic.AddParticle(
                ctypes.c_double(xmax * math.cos(nu) * math.cos(phi)),
                ctypes.c_double(amax * math.cos(nu) * math.sin(phi)),
                ctypes.c_double(ymax * math.sin(nu) * math.cos(phi)),
                ctypes.c_double(bmax * math.sin(nu) * math.sin(phi)),
                ctypes.c_double(float(dk)),
                ctypes.c_double(float(dm)))


def Beam3d(xmax, amax, ymax, bmax, dk, dm, num):
    """ Komplexen 3d-Strahl einfuegen """
    global lastFunction
    lastFunction = "Beam3d"
    for i in linspace(0, 360, int(num / 36)):
        nu = math.radians(i)
        for j in xrange(0, 360, 10):
            phi = math.radians(j)
            optic.AddParticle(
                ctypes.c_double(xmax * math.cos(nu) * math.cos(phi)),
                ctypes.c_double(amax * math.cos(nu) * math.sin(phi)),
                ctypes.c_double(ymax * math.sin(nu) * math.cos(phi)),
                ctypes.c_double(bmax * math.sin(nu) * math.sin(phi)),
                ctypes.c_double(float(dk)),
                ctypes.c_double(float(dm)))


def AddGaussBeam(strag_x, strag_a, strag_y, strag_b, x=0., a=0., y=0., b=0., dk=0., dm=0., strag_k=0., strag_m=0., num=250):
    """ Gaussverteilten zufallsbasierten Strahl einfuegen """
    global lastFunction
    lastFunction = "AddGaussBeam"

    optic.AddGaussBeam(
        ctypes.c_double(float(x)),
        ctypes.c_double(float(strag_x)),
        ctypes.c_double(float(a)),
        ctypes.c_double(float(strag_a)),
        ctypes.c_double(float(y)),
        ctypes.c_double(float(strag_y)),
        ctypes.c_double(float(b)),
        ctypes.c_double(float(strag_b)),
        ctypes.c_double(float(dk)),
        ctypes.c_double(float(strag_k)),
        ctypes.c_double(float(dm)),
        ctypes.c_double(float(strag_m)),
        ctypes.c_double(float(num)))


def GaussBeam(strag_x, strag_a, strag_y, strag_b, x=0., a=0., y=0., b=0., dk=0., dm=0., num=250, strag_k=0., strag_m=0., sigma=1.):
    """ Gaussverteilten zufallsbasierten Strahl einfuegen """
    global lastFunction
    lastFunction = "GaussBeam"

    optic.AddGaussBeam(
        ctypes.c_double(float(x)),
        ctypes.c_double(float(strag_x / sigma)),
        ctypes.c_double(float(a)),
        ctypes.c_double(float(strag_a / sigma)),
        ctypes.c_double(float(y)),
        ctypes.c_double(float(strag_y / sigma)),
        ctypes.c_double(float(b)),
        ctypes.c_double(float(strag_b / sigma)),
        ctypes.c_double(float(dk)),
        ctypes.c_double(float(strag_k / sigma)),
        ctypes.c_double(float(dm)),
        ctypes.c_double(float(strag_m / sigma)),
        ctypes.c_double(float(num)))


def AddBeamRandomGauss(xmax, amax, ymax, bmax, dk, dm, number, sigma=1.):
    """ Gaussverteilten zufallsbasierten Strahl einfuegen """
    global lastFunction
    lastFunction = "AddBeamRandomGauss"

    random.seed()
    for i in xrange(number):
        x = random.gauss(0, xmax / sigma)
        y = random.gauss(0, ymax / sigma)
        a = random.gauss(0, amax / sigma)
        b = random.gauss(0, bmax / sigma)

        optic.AddParticle(
            ctypes.c_double(x),
            ctypes.c_double(a),
            ctypes.c_double(y),
            ctypes.c_double(b),
            ctypes.c_double(float(dk)),
            ctypes.c_double(float(dm)))


def BeamRandomGauss(xmax, amax, ymax, bmax, dk, dm, number, sigma=1.):
    """ Gaussverteilten zufallsbasierten Strahl einfuegen """
    global lastFunction
    lastFunction = "AddBeamRandomGauss"

    #random.seed()
    for i in xrange(number):
        x = random.gauss(0, xmax / sigma)
        y = random.gauss(0, ymax / sigma)
        a = random.gauss(0, amax / sigma)
        b = random.gauss(0, bmax / sigma)

        optic.AddParticle(
            ctypes.c_double(x),
            ctypes.c_double(a),
            ctypes.c_double(y),
            ctypes.c_double(b),
            ctypes.c_double(float(dk)),
            ctypes.c_double(float(dm)))


def AddBeam(xmax, amax, ymax, bmax, dk, dm, delta=10):
    """ Einfachen 2d-Strahl einfuegen """
    global lastFunction
    lastFunction = "AddBeam"

    for i in xrange(0, 360, int(delta)):
        degtorad = i * math.pi / 180

        optic.AddParticle(
            ctypes.c_double(float(xmax * math.cos(degtorad))),
            ctypes.c_double(float(amax * math.sin(degtorad))),
            ctypes.c_double(float(ymax * math.cos(degtorad))),
            ctypes.c_double(float(bmax * math.sin(degtorad))),
            ctypes.c_double(float(dk)),
            ctypes.c_double(float(dm)))


def Beam(xmax, amax, ymax, bmax, dk, dm, num=250, x=0, y=0):
    """ Einfachen 2d-Strahl einfuegen """
    global lastFunction
    lastFunction = "AddBeam"

    for i in linspace(0, 360, int(num)):
        degtorad = i * math.pi / 180.

        optic.AddParticle(
            ctypes.c_double(float(x + xmax * math.cos(degtorad))),
            ctypes.c_double(float(amax * math.sin(degtorad))),
            ctypes.c_double(float(y + ymax * math.cos(degtorad))),
            ctypes.c_double(float(bmax * math.sin(degtorad))),
            ctypes.c_double(float(dk)),
            ctypes.c_double(float(dm)))


def AddSource(file=None):
    """ Externe Quelle einfuegen """
    global lastFunction, SOURCE, SOURCEFILES
    lastFunction = "AddSource"

    if (not file is None) and (not file in SOURCEFILES):
        try:
            SOURCEFILES.append(file)
            for line in open(file, "r"):
                SOURCE.append([float(elem) for elem in line.split()])
        except:
            print file, "was not found."

    elif not len(SOURCE) > 0:
        try:
            SOURCE = []
            for line in open("source.dat", "r"):
                SOURCE.append([float(elem) for elem in line.split()])
        except:
            print "source.dat was not found."

    for i in xrange(len(SOURCE)):
        optic.AddParticle(
            ctypes.c_double(SOURCE[i][0]),
            ctypes.c_double(SOURCE[i][1]),
            ctypes.c_double(SOURCE[i][2]),
            ctypes.c_double(SOURCE[i][3]),
            ctypes.c_double(SOURCE[i][4]),
            ctypes.c_double(SOURCE[i][5]))
Source = AddSource


def ClearSource():
    global SOURCE, SOURCEFILES

    SOURCE      = []
    SOURCEFILES = []


def AddMatrix(num, mat, length):
    """ Allgemeinte Matrix einfuegen """
    global lastFunction
    lastFunction = "AddMatrix"

    matrixarray = ctypes.c_double*36
    matrix      = matrixarray()
    for i in xrange(36):
        matrix[i] = ctypes.c_double(float(mat[i]))
    optic.AddMatrix(
        ctypes.c_int(int(num)),
        matrix,
        ctypes.c_double(float(length)))
Matrix = AddMatrix


def AddAMSAcc(v_qsnout, v_gesamt, v_vorbeschl, q):
    """ Beschleunigung Cologne AMS nach HINTERBERGER """
    global lastFunction
    lastFunction = "AddAMSAcc"

    #Widerstand 1-4 Zweig1
    R_a1 = 400.e6 + 10.e3

    #Widerstand 1-4 Zweig2
    R_a2 = 400.e6 + 2.e9

    #Widerstand 1-4
    R_a = float(1 / (1 / R_a1 + 1 / R_a2))

    #Widerstand 5-84
    R_b = 18. * 40.e6 + 62. * 200.e6 + 81. * 200.e6

    #Gesamtwiderstand HE1
    R_gesamt = R_b + R_a

    #Gesamtstrom durch ACC linke Seite
    I_gesamt = float(v_gesamt / R_gesamt)

    #Laenge eines Segmentes
    L = float(25.4e-3)

    #Spannungsabfall 1..4
    V_a  = float(R_a * I_gesamt)
    I_a2 = float(V_a / R_a2)

    # Bias..Anfang Q-Snout
    T0 = float(v_vorbeschl - 4.e3)
    T1 = float(v_vorbeschl + v_qsnout)

    # Q-Snout : Laenge 300.4mm
    T0 = T1
    T1 = T1
    N  = math.sqrt(T1 / T0)
    AddSegment(N, 300.4e-3)

    ### HE1 : Ende Q-Snout..Anfang Terminal
    # 3..4
    T0 = T1
    T1 = float(v_vorbeschl + V_a - 1.e9 * I_a2)
    N  = math.sqrt(T1 / T0)
    AddSegment(N, L)

    # 4..5
    T0 = T1
    T1 = float(v_vorbeschl + V_a)
    N  = math.sqrt(T1 / T0)
    AddSegment(N, L)

    #5..23
    deltaV = 40.e6 * I_gesamt
    for i in xrange(18):
        T0 = T1
        T1 = float(T0 + deltaV)
        N  = math.sqrt(T1 / T0)
        AddSegment(N, L)

    #23..85
    deltaV = 200.e6 * I_gesamt
    for i in xrange(62):
        T0 = T1
        T1 = float(T0 + deltaV)
        N  = math.sqrt(T1 / T0)
        AddSegment(N, L)

    # Space = 76.23mm
    T0 = T1
    T1 = T1
    N  = math.sqrt(T1 / T0)
    AddSegment(N, 76.23e-3)

    ## Space..Terminal
    #85..166
    deltaV = 200.e6 * I_gesamt
    for i in xrange(81):
        T0 = T1
        T1 = float(T0 + deltaV)
        N  = math.sqrt(T1 / T0)
        AddSegment(N, L)

    # Terminal (Laenge 1147 mm) / Stripping
    T0 = T1
    T1 = T1
    N  = math.sqrt(T1 / T0)
    AddSegment(N, 1147.e-3)

    ### HE2 : Ende Terminal..letzte Elektrode
    R_gesamt = float(9.0 * 120.e6 + 156. * 300.e6)

    #Gesamtstrom durch ACC rechte Seite
    I_gesamt = float(v_gesamt / R_gesamt)

    ## Terminal..Space
    #1..9
    deltaV = 120.e6 * I_gesamt
    for i in xrange(9):
        T0 = T1
        T1 = float(T0 + q * deltaV)
        N  = math.sqrt(T1 / T0)
        AddSegment(N, L)

    #9..81
    deltaV = 300.e6 * I_gesamt
    for i in xrange(72):
        T0 = T1
        T1 = float(T0 + q * deltaV)
        N  = math.sqrt(T1 / T0)
        AddSegment(N, L)

    # Space = 76.05 mm
    T0 = T1
    T1 = T1
    N  = math.sqrt(T1 / T0)
    AddSegment(N, 76.05e-3)

    # Space..letzte Elekrode
    #81..165
    deltaV = 300.e6 * I_gesamt
    for i in xrange(84):
        T0 = T1
        T1 = float(T0 + q * deltaV)
        N  = math.sqrt(T1 / T0)
        AddSegment(N, L)

    # Letzte Elektrode..Flansch aussen = 487.8mm
    AddSegment(1, 447.81e-3)
AMSAcc = AddAMSAcc


def AddVBFN(extraktion, deltaV, laenge=.276, b=1.13, b1=-1., b2=-1., segment=0):
    """ Vorbeschleunigung FN """
    global lastFunction
    lastFunction = "AddVBFN"

    if ((b1 == -1.) and (b2 == -1.)):
        b1 = b
        b2 = b

    T0 = float(extraktion)
    T1 = float(extraktion + deltaV)

    if segment == 0:
        AddSegment11(T0, T1, .09, .09, laenge, b1, b2)
    elif segment == 1:
        AddSegment(math.sqrt(T1/T0), laenge)
    elif segment == 2:
        AddSegment1(T0, T1, -1., .09, .09, laenge, b)
    elif segment == 3:
        AddSegment2VBFN(T0, T1, .09, laenge, 0., -(T1-T0)/laenge, 0., -1.)
VBFN = AddVBFN


def AddFNAccNeu(vt, T0, q, b=0.57, b1=-1., b2=-1., D1=.088, factor1=1., factor2=1., beamprofile=False, addwaist=False):
    """ FN Beschleuniger einfach (ohne einzelne Elektroden) """
    global lastFunction
    lastFunction = "AddFNAccNeu"

    if ((b1 == -1.) and (b2 == -1.)):
        b1 = b
        b2 = b
    R1 = 53660.e6
    R2 = 53850.e6
    R3 = 55710.e6
    R4 = 55800.e6
    RLE = R1 + R2
    RHE = R3 + R4

    l1 = 2.438
    l2 = 2.438
    l3 = 2.438
    l4 = 2.438
    space1 = .273
    space2 = .273
    terminal = .609 + .609

    E1 = R1 * vt / RLE / l1
    E2 = R2 * vt / RLE / l2
    E3 = R3 * vt / RHE / l3
    E4 = R4 * vt / RHE / l4

    T1 = T0 + E1 * l1
    T2 = T1 + E2 * l2
    T3 = T2 + q * E3 * l3
    T4 = T3 + q * E4 * l4

    #D1 = .088 #LE1
    D2 = .03  #LE1
    D3 = .03  #LE2
    D4 = .03  #LE2
    D5 = .03  #HE1
    D6 = .03  #HE1
    D7 = .03  #HE2
    D8 = .03  #HE2

    AddSegment11(T0, T1, D1, D2, l1, b1, b2)             #LE1
    AddSegment11(T1, T1, D2, D3, space1, b1, b2)         #Space1
    AddSegment11(T1, T2, D3, D4, l2, b1, b2)             #LE2
    AddSegment11(T2, T2, D4, D5, terminal / 2., b1, b2)  #Terminal

    if (beamprofile):
        print "\n"
        AddBeamProfile()

    if (addwaist):
        AddWaist()

    AddModifyEmittance(factor1, factor2)

    if (beamprofile): AddBeamProfile()

    AddSegment11(T2, T2, D4, D5, terminal / 2., b1, b2)  #Terminal
    AddSegment11(T2, T3, D5, D6, l3, b1, b2)             #HE1
    AddSegment11(T3, T3, D6, D7, space2, b1, b2)         #Space2
    AddSegment11(T3, T4, D7, D8, l4, b1, b2)             #HE2
FNAccNeu = AddFNAccNeu


# Beschleunigung Cologne FN
def AddFNAcc(v_gesamt, v_vorbeschl, q):
    global lastFunction
    lastFunction = "AddFNAcc"

    #Gesamtwiderstand LE
    R_gesamt = 107510.0e6

    #Gesamtstrom durch ACC linke Seite
    I_gesamt = float(-v_gesamt / R_gesamt)

    #Laenge eines Segmentes
    L = float((2.532 + .223) / 100.0)

    #Durchmesser des Elektrodenringes, update: neue Roehren
    D = 5.5 * 2. / 100.0
    q = float(q)

    # Widerstaende in MOhm fuer LE1-2 (erstes Element zw. 1 und 2 usw) Stand 25.10.2011
    R_Segment1 = 0, 0, 0, 580, 570, 0, 600, 610, 600, 600, 0, 600, 590, 550, 600, 600, 590, 590, 580, 590, 600, 600, 600, 590, 600, 600, 600, 590, 600, 600, 610, 590, 600, 0, 610, 610, 600, 600, 600, 600, 600, 600, 590, 610, 600, 0, 600, 600, 600, 610, 600, 600, 600, 600, 600, 610, 600, 600, 610, 600, 610, 600, 610, 600, 600, 610, 600, 600, 610, 590, 610, 610, 560, 620, 610, 620, 610, 610, 610, 630, 600, 620, 620, 610, 640, 610, 620, 610, 630, 610, 600, 630, 630, 620, 610, 600, 0
    R_Segment2 = 570, 550, 580, 580, 580, 550, 580, 580, 570, 580, 570, 580, 570, 560, 570, 580, 570, 580, 580, 580, 610, 570, 580, 580, 570, 560, 570, 580, 570, 570, 610, 580, 560, 570, 590, 560, 610, 0, 540, 560, 580, 570, 550, 570, 560, 600, 570, 570, 580, 600, 570, 580, 590, 570, 580, 570, 570, 570, 570, 560, 580, 540, 520, 550, 550, 550, 540, 530, 540, 550, 550, 540, 550, 550, 560, 550, 570, 570, 560, 540, 560, 580, 550, 550, 540, 570, 540, 560, 560, 550, 560, 580, 570, 580, 550, 580, 0

    # Widerstaende in MOhm fuer HE1-2 (erstes Element zw. 1 und 2 usw) Stand 25.10.2011
    R_Segment3 = 590, 580, 580, 580, 580, 590, 560, 590, 580, 580, 590, 560, 570, 570, 570, 590, 570, 580, 590, 600, 580, 560, 580, 580, 570, 600, 570, 580, 570, 570, 570, 580, 570, 580, 600, 580, 570, 580, 580, 580, 570, 590, 580, 570, 590, 570, 570, 590, 580, 590, 580, 580, 590, 600, 580, 580, 570, 570, 570, 580, 580, 580, 600, 590, 590, 570, 580, 560, 570, 570, 600, 580, 580, 570, 580, 570, 570, 570, 580, 580, 600, 620, 570, 590, 580, 600, 580, 570, 590, 600, 580, 580, 590, 570, 580, 590, 0
    R_Segment4 = 580, 570, 570, 570, 580, 580, 580, 570, 570, 580, 570, 580, 580, 570, 570, 600, 620, 610, 600, 590, 580, 590, 580, 580, 580, 590, 580, 580, 580, 580, 580, 580, 580, 570, 580, 580, 580, 580, 580, 580, 590, 570, 580, 580, 580, 570, 590, 580, 580, 580, 580, 580, 580, 580, 590, 580, 580, 580, 570, 580, 590, 580, 580, 580, 580, 580, 570, 580, 580, 570, 580, 580, 580, 580, 590, 580, 570, 580, 580, 580, 590, 590, 580, 580, 580, 580, 590, 580, 590, 580, 590, 580, 580, 580, 590, 600, 0
    # INFO: Letzter Widerstand (0) dient zur Berechnung von E3

    # LE1 : 4..5
    E1 = 0.
    E2 = R_Segment1[3] * 1.0e6 * I_gesamt / L
    E3 = R_Segment1[4] * 1.0e6 * I_gesamt / L
    T0 = float(v_vorbeschl)
    T1 = float(T0 + float(-R_Segment1[4-1]) * 1.0e6 * I_gesamt)
    AddSegment2(T0, T1, 8.8 / 100.0, L, E1, E2, E3, -1.)

    # LE1 : 5..97
    for i in xrange(5, 9):
        E1 = E2
        E2 = E3
        E3 = R_Segment1[i] * 1.0e6 * I_gesamt / L
        T0 = T1
        T1 = float(T0 + float(-R_Segment1[i-1]) * 1.0e6 * I_gesamt)
        AddSegment2(T0, T1, 7.2 / 100., L, E1, E2, E3, -1.)

    for i in xrange(9, 12):
        E1 = E2
        E2 = E3
        E3 = R_Segment1[i] * 1.e6 * I_gesamt / L
        T0 = T1
        T1 = float(T0 + float(-R_Segment1[i-1]) * 1.e6 * I_gesamt)
        AddSegment2(T0, T1, 6.4 / 100.0, L, E1, E2, E3, -1.)

    for i in xrange(12, 15):
        E1 = E2
        E2 = E3
        E3 = R_Segment1[i] * 1.0e6 * I_gesamt / L
        T0 = T1
        T1 = float(T0 + float(-R_Segment1[i-1]) * 1.0e6 * I_gesamt)
        AddSegment2(T0, T1, 6. / 100., L, E1, E2, E3, -1.)

    for i in xrange(15, 18):
        E1 = E2
        E2 = E3
        E3 = R_Segment1[i] * 1.e6 * I_gesamt / L
        T0 = T1
        T1 = float(T0 + float(-R_Segment1[i-1]) * 1.0e6 * I_gesamt)
        AddSegment2(T0, T1, 5.2 / 100., L, E1, E2, E3, -1.)

    for i in xrange(18, 21):
        E1 = E2
        E2 = E3
        E3 = R_Segment1[i] * 1.e6 * I_gesamt / L
        T0 = T1
        T1 = float(T0 + float(-R_Segment1[i-1]) * 1.e6 * I_gesamt)
        AddSegment2(T0, T1, 4.8 / 100., L, E1, E2, E3, -1.)

    for i in xrange(21, 31):
        E1 = E2
        E2 = E3
        E3 = R_Segment1[i] * 1.e6 * I_gesamt / L
        T0 = T1
        T1 = float(T0 + float(-R_Segment1[i-1]) * 1.e6 * I_gesamt)
        AddSegment2(T0, T1, 4. / 100., L, E1, E2, E3, -1.)

    for i in xrange(31, 76):
        E1 = E2
        E2 = E3
        E3 = R_Segment1[i] * 1.e6 * I_gesamt / L
        T0 = T1
        T1 = float(T0 + float(-R_Segment1[i-1]) * 1.e6 * I_gesamt)
        AddSegment2(T0, T1, 3. / 100., L, E1, E2, E3, -1.)

    for i in xrange(76, 97):
        E1 = E2
        E2 = E3
        E3 = R_Segment1[i] * 1.e6 * I_gesamt / L
        T0 = T1
        T1 = float(T0 + float(-R_Segment1[i-1]) * 1.e6 * I_gesamt)
        AddSegment2(T0, T1, 2. / 100., L, E1, E2, E3, -1.)

    E1 = E2
    E2 = E3
    E3 = 0.
    T0 = T1
    T1 = T1

    # Space = 22.225cm = Drift
    E1 = E2
    E2 = E3
    E3 = R_Segment2[0] * 1.e6 * I_gesamt / L
    T0 = T1
    T1 = T1
    AddSegment2(T0, T1, D, 22.225 / 100.0, E1, E2, E3, -1.0)

    # LE2 : 1..97
    for i in xrange(1, 97):
        E1 = E2
        E2 = E3
        E3 = R_Segment2[i] * 1.e6 * I_gesamt / L
        T0 = T1
        T1 = float(T0 + float(-R_Segment2[i-1]) * 1.e6 * I_gesamt)
        AddSegment2(T0, T1, D, L, E1, E2, E3, -1.0)

    E1 = E2
    E2 = E3
    E3 = 0.
    T0 = T1
    T1 = T1

    # HE1 : R und I
    R_gesamt = 111510.0e6

    #Gesamtstrom durch ACC rechte Seite
    I_gesamt = float(v_gesamt / R_gesamt)

    # Terminal / Stripping = Drift
    E1 = E2
    E2 = E3
    E3 = R_Segment3[0] * 1.e6 * I_gesamt / L
    T0 = T1
    T1 = T1
    AddSegment2(T0, T1, D, 121.92 / 100., E1, E2, E3, -1.)

    # HE1 : 1..97
    for i in xrange(1, 97):
        E1 = E2
        E2 = E3
        E3 = R_Segment3[i] * 1.e6 * I_gesamt / L
        T0 = T1
        T1 = float(T0 + q * float(R_Segment3[i-1]) * 1.e6 * I_gesamt)
        AddSegment2(T0, T1, D, L, E1, E2, E3, q)

    E1 = E2
    E2 = E3
    E3 = 0.
    T0 = T1
    T1 = T1

    # Space = 22.225cm = Drift
    E1 = E2
    E2 = E3
    E3 = R_Segment4[0] * 1.e6 * I_gesamt / L
    T0 = T1
    T1 = T1
    AddSegment2(T0, T1, D, 22.225 / 100., E1, E2, E3, q)

    # HE2 : 1..97
    for i in xrange(1, 97):
        E1 = E2
        E2 = E3
        E3 = R_Segment4[i] * 1.e6 * I_gesamt / L
        T0 = T1
        T1 = float(T0 + q * float(R_Segment4[i-1]) * 1.e6 * I_gesamt)
        AddSegment2(T0, T1, D, L, E1, E2, E3, q)

    E1 = E2
    E2 = E3
    E3 = 0.
    T0 = T1
    T1 = T1

    # Letzte Elektrode..Flansch aussen = Drift
    E1 = E2
    E2 = E3
    E3 = 0.
    T0 = T1
    T1 = T1
    AddSegment2(T0, T1, D, 117.738 / 100., E1, E2, E3, q)
FNAcc = AddFNAcc


def AddFNAcc1(v_gesamt, v_vorbeschl, q):
    global lastFunction
    lastFunction = "AddFNAcc1"

    #Gesamtwiderstand LE
    R_gesamt = 107510.0e6

    #Gesamtstrom durch ACC linke Seite
    I_gesamt = float(-v_gesamt / R_gesamt)

    #Laenge eines Segmentes
    L = float((2.532 + .223) / 100.0)

    #Durchmesser des Elektrodenringes, update: neue Roehren
    D = 5.5 * 2 / 100.
    q = float(q)

    # Widerstaende in MOhm fuer LE1-2 (erstes Element zw. 1 und 2 usw) Stand 25.10.2011
    R_Segment1 = 0, 0, 0, 580, 570, 0, 600, 610, 600, 600, 0, 600, 590, 550, 600, 600, 590, 590, 580, 590, 600, 600, 600, 590, 600, 600, 600, 590, 600, 600, 610, 590, 600, 0, 610, 610, 600, 600, 600, 600, 600, 600, 590, 610, 600, 0, 600, 600, 600, 610, 600, 600, 600, 600, 600, 610, 600, 600, 610, 600, 610, 600, 610, 600, 600, 610, 600, 600, 610, 590, 610, 610, 560, 620, 610, 620, 610, 610, 610, 630, 600, 620, 620, 610, 640, 610, 620, 610, 630, 610, 600, 630, 630, 620, 610, 600, 0
    R_Segment2 = 570, 550, 580, 580, 580, 550, 580, 580, 570, 580, 570, 580, 570, 560, 570, 580, 570, 580, 580, 580, 610, 570, 580, 580, 570, 560, 570, 580, 570, 570, 610, 580, 560, 570, 590, 560, 610, 0, 540, 560, 580, 570, 550, 570, 560, 600, 570, 570, 580, 600, 570, 580, 590, 570, 580, 570, 570, 570, 570, 560, 580, 540, 520, 550, 550, 550, 540, 530, 540, 550, 550, 540, 550, 550, 560, 550, 570, 570, 560, 540, 560, 580, 550, 550, 540, 570, 540, 560, 560, 550, 560, 580, 570, 580, 550, 580, 0

    # Widerstaende in MOhm fuer HE1-2 (erstes Element zw. 1 und 2 usw) Stand 25.10.2011
    R_Segment3 = 590, 580, 580, 580, 580, 590, 560, 590, 580, 580, 590, 560, 570, 570, 570, 590, 570, 580, 590, 600, 580, 560, 580, 580, 570, 600, 570, 580, 570, 570, 570, 580, 570, 580, 600, 580, 570, 580, 580, 580, 570, 590, 580, 570, 590, 570, 570, 590, 580, 590, 580, 580, 590, 600, 580, 580, 570, 570, 570, 580, 580, 580, 600, 590, 590, 570, 580, 560, 570, 570, 600, 580, 580, 570, 580, 570, 570, 570, 580, 580, 600, 620, 570, 590, 580, 600, 580, 570, 590, 600, 580, 580, 590, 570, 580, 590, 0
    R_Segment4 = 580, 570, 570, 570, 580, 580, 580, 570, 570, 580, 570, 580, 580, 570, 570, 600, 620, 610, 600, 590, 580, 590, 580, 580, 580, 590, 580, 580, 580, 580, 580, 580, 580, 570, 580, 580, 580, 580, 580, 580, 590, 570, 580, 580, 580, 570, 590, 580, 580, 580, 580, 580, 580, 580, 590, 580, 580, 580, 570, 580, 590, 580, 580, 580, 580, 580, 570, 580, 580, 570, 580, 580, 580, 580, 590, 580, 570, 580, 580, 580, 590, 590, 580, 580, 580, 580, 590, 580, 590, 580, 590, 580, 580, 580, 590, 600, 0
    # INFO: Letzter Widerstand (0) dient zur Berechnung von E3

    # LE1 : 4..5
    E1 = 0.
    E2 = R_Segment1[3] * 1.e6 * I_gesamt / L
    E3 = R_Segment1[4] * 1.e6 * I_gesamt / L
    T0 = float(v_vorbeschl)
    T1 = float(T0 + float(-R_Segment1[4-1]) * 1.e6 * I_gesamt)

    # //Durch die Grid Lens ist das erste Fringing Field 0
    AddSegment2(T0, T1, 8.8 / 100.0, L, E1, E2, E3, -1.)

    # LE1 : 5..97
    for i in xrange(5, 9):
        E1 = E2
        E2 = E3
        E3 = R_Segment1[i] * 1.e6 * I_gesamt / L
        T0 = T1
        T1 = float(T0 + float(-R_Segment1[i-1]) * 1.e6 * I_gesamt)
        AddSegment2(T0, T1, 7.2 / 100.0, L, E1, E2, E3, -1.)

    for i in xrange(9, 12):
        E1 = E2
        E2 = E3
        E3 = R_Segment1[i]*1.0e6*I_gesamt / L
        T0 = T1
        T1 = float(T0+float(-R_Segment1[i-1])*1.0e6*I_gesamt)
        AddSegment2(T0, T1, 6.4 / 100.0, L, E1, E2, E3, -1.0)

    for i in xrange(12, 15):
        E1 = E2
        E2 = E3
        E3 = R_Segment1[i]*1.0e6*I_gesamt / L
        T0 = T1
        T1 = float(T0+float(-R_Segment1[i-1])*1.0e6*I_gesamt)
        AddSegment2(T0, T1, 6.0 / 100.0, L, E1, E2, E3, -1.0)

    for i in xrange(15, 18):
        E1 = E2
        E2 = E3
        E3 = R_Segment1[i]*1.0e6*I_gesamt / L
        T0 = T1
        T1 = float(T0+float(-R_Segment1[i-1])*1.0e6*I_gesamt)
        AddSegment2(T0, T1, 5.2 / 100.0, L, E1, E2, E3, -1.0)

    for i in xrange(18, 21):
        E1 = E2
        E2 = E3
        E3 = R_Segment1[i]*1.0e6*I_gesamt / L
        T0 = T1
        T1 = float(T0+float(-R_Segment1[i-1])*1.0e6*I_gesamt)
        AddSegment2(T0, T1, 4.8 / 100.0, L, E1, E2, E3, -1.0)

    for i in xrange(21, 31):
        E1 = E2
        E2 = E3
        E3 = R_Segment1[i]*1.0e6*I_gesamt / L
        T0 = T1
        T1 = float(T0+float(-R_Segment1[i-1])*1.0e6*I_gesamt)
        AddSegment2(T0, T1, 4.0 / 100.0, L, E1, E2, E3, -1.0)

    for i in xrange(31, 76):
        E1 = E2
        E2 = E3
        E3 = R_Segment1[i]*1.0e6*I_gesamt / L
        T0 = T1
        T1 = float(T0+float(-R_Segment1[i-1])*1.0e6*I_gesamt)
        AddSegment2(T0, T1, 3.6 / 100.0, L, E1, E2, E3, -1.0)

    for i in xrange(76, 97):
        E1 = E2
        E2 = E3
        E3 = R_Segment1[i]*1.0e6*I_gesamt / L
        T0 = T1
        T1 = float(T0+float(-R_Segment1[i-1])*1.0e6*I_gesamt)
        AddSegment2(T0, T1, 2.8 / 100.0, L, E1, E2, E3, -1.0)

    E1 = E2
    E2 = E3
    E3 = 0.0
    T0 = T1
    T1 = T1

    # Space = 22.225cm = Drift
    E1 = E2
    E2 = E3
    E3 = R_Segment2[0] * 1.e6 * I_gesamt / L
    T0 = T1
    T1 = T1
    AddSegment2(T0, T1, D, 22.225 / 100.0, E1, E2, E3, -1.0)

    # LE2 : 1..97
    for i in xrange(1, 97):
        E1 = E2
        E2 = E3
        E3 = R_Segment2[i] * 1.e6 * I_gesamt / L
        T0 = T1
        T1 = float(T0+float(-R_Segment2[i-1]) * 1.0e6 * I_gesamt)
        AddSegment2(T0, T1, D, L, E1, E2, E3, -1.)

    E1 = E2
    E2 = E3
    E3 = 0.0
    T0 = T1
    T1 = T1

    # HE1 : R und I
    R_gesamt = 111510.0e6

    #Gesamtstrom durch ACC rechte Seite
    I_gesamt = float(v_gesamt / R_gesamt)

    # Terminal / Stripping = Drift
    E1 = E2
    E2 = E3
    E3 = R_Segment3[0] * 1.e6 * I_gesamt / L
    T0 = T1
    T1 = T1
    N  = math.sqrt(T1 / T0)
    AddSegment(N, 121.92 / 100.0)

    # HE1 : 1..97
    for i in xrange(1, 97):
        E1 = E2
        E2 = E3
        E3 = R_Segment3[i]*1.0e6*I_gesamt / L
        T0 = T1
        T1 = float(T0 + q * float(R_Segment3[i-1]) * 1.e6 * I_gesamt)
        N  = math.sqrt(T1 / T0)
        AddSegment(N, L)

    E1 = E2
    E2 = E3
    E3 = 0.
    T0 = T1
    T1 = T1

    # Space = 22.225cm = Drift
    E1 = E2
    E2 = E3
    E3 = R_Segment4[0]*1.0e6*I_gesamt / L
    T0 = T1
    T1 = T1
    N  = math.sqrt(T1 / T0)
    AddSegment(N, 22.225 / 100.0)

    # HE2 : 1..97
    for i in xrange(1, 97):
        E1 = E2
        E2 = E3
        E3 = R_Segment4[i]*1.0e6*I_gesamt / L
        T0 = T1
        T1 = float(T0 + q * float(R_Segment4[i-1]) * 1.e6 * I_gesamt)
        N  = math.sqrt(T1 / T0)
        AddSegment(N, L)

    E1 = E2
    E2 = E3
    E3 = 0.
    T0 = T1
    T1 = T1

    # Letzte Elektrode..Flansch aussen = Drift
    E1 = E2
    E2 = E3
    E3 = 0.
    T0 = T1
    T1 = T1
    AddSegment2(T0, T1, D, 117.738 / 100.0, E1, E2, E3, q)


def AddFNAcc2(v_gesamt, v_vorbeschl, q):
    global lastFunction
    lastFunction = "AddFNAcc2"

    #Gesamtwiderstand LE
    R_gesamt = 107510.0e6

    #Gesamtstrom durch ACC linke Seite
    I_gesamt = float(-v_gesamt / R_gesamt)

    #Laenge eines Segmentes
    L = float((2.532+.223)/100.0)

    #Durchmesser des Elektrodenringes, update: neue Roehren
    D = 6.35 * 2 / 100
    q = float(q)

    # Widerstaende in MOhm fuer LE1-2 (erstes Element zw. 1 und 2 usw) Stand 25.10.2011
    R_Segment1 = 0, 0, 0, 580, 570, 0, 600, 610, 600, 600, 0, 600, 590, 550, 600, 600, 590, 590, 580, 590, 600, 600, 600, 590, 600, 600, 600, 590, 600, 600, 610, 590, 600, 0, 610, 610, 600, 600, 600, 600, 600, 600, 590, 610, 600, 0, 600, 600, 600, 610, 600, 600, 600, 600, 600, 610, 600, 600, 610, 600, 610, 600, 610, 600, 600, 610, 600, 600, 610, 590, 610, 610, 560, 620, 610, 620, 610, 610, 610, 630, 600, 620, 620, 610, 640, 610, 620, 610, 630, 610, 600, 630, 630, 620, 610, 600, 0
    R_Segment2 = 570, 550, 580, 580, 580, 550, 580, 580, 570, 580, 570, 580, 570, 560, 570, 580, 570, 580, 580, 580, 610, 570, 580, 580, 570, 560, 570, 580, 570, 570, 610, 580, 560, 570, 590, 560, 610, 0, 540, 560, 580, 570, 550, 570, 560, 600, 570, 570, 580, 600, 570, 580, 590, 570, 580, 570, 570, 570, 570, 560, 580, 540, 520, 550, 550, 550, 540, 530, 540, 550, 550, 540, 550, 550, 560, 550, 570, 570, 560, 540, 560, 580, 550, 550, 540, 570, 540, 560, 560, 550, 560, 580, 570, 580, 550, 580, 0

    # Widerstaende in MOhm fuer HE1-2 (erstes Element zw. 1 und 2 usw) Stand 25.10.2011
    R_Segment3 = 590, 580, 580, 580, 580, 590, 560, 590, 580, 580, 590, 560, 570, 570, 570, 590, 570, 580, 590, 600, 580, 560, 580, 580, 570, 600, 570, 580, 570, 570, 570, 580, 570, 580, 600, 580, 570, 580, 580, 580, 570, 590, 580, 570, 590, 570, 570, 590, 580, 590, 580, 580, 590, 600, 580, 580, 570, 570, 570, 580, 580, 580, 600, 590, 590, 570, 580, 560, 570, 570, 600, 580, 580, 570, 580, 570, 570, 570, 580, 580, 600, 620, 570, 590, 580, 600, 580, 570, 590, 600, 580, 580, 590, 570, 580, 590, 0
    R_Segment4 = 580, 570, 570, 570, 580, 580, 580, 570, 570, 580, 570, 580, 580, 570, 570, 600, 620, 610, 600, 590, 580, 590, 580, 580, 580, 590, 580, 580, 580, 580, 580, 580, 580, 570, 580, 580, 580, 580, 580, 580, 590, 570, 580, 580, 580, 570, 590, 580, 580, 580, 580, 580, 580, 580, 590, 580, 580, 580, 570, 580, 590, 580, 580, 580, 580, 580, 570, 580, 580, 570, 580, 580, 580, 580, 590, 580, 570, 580, 580, 580, 590, 590, 580, 580, 580, 580, 590, 580, 590, 580, 590, 580, 580, 580, 590, 600, 0
    # INFO: Letzter Widerstand (0) dient zur Berechnung von E3

    # LE1 : 4..5
    E1 = 0.
    E2 = R_Segment1[3] * 1.e6 * I_gesamt / L
    E3 = R_Segment1[4] * 1.e6 * I_gesamt / L
    T0 = float(v_vorbeschl)
    T1 = float(T0+float(-R_Segment1[4-1])*1.0e6*I_gesamt)
    N = math.sqrt(T1 / T0)
    AddSegment(N, L)

    # LE1 : 5..97
    for i in xrange(5, 9):
        E1 = E2
        E2 = E3
        E3 = R_Segment1[i] * 1.e6 * I_gesamt / L
        T0 = T1
        T1 = float(T0 + float(-R_Segment1[i-1]) * 1.e6 * I_gesamt)
        N = math.sqrt(T1 / T0)
        AddSegment(N, L)

    for i in xrange(9, 12):
        E1 = E2
        E2 = E3
        E3 = R_Segment1[i] * 1.e6 * I_gesamt / L
        T0 = T1
        T1 = float(T0 + float(-R_Segment1[i-1]) * 1.e6 * I_gesamt)
        N = math.sqrt(T1 / T0)
        AddSegment(N, L)

    for i in xrange(12, 15):
        E1 = E2
        E2 = E3
        E3 = R_Segment1[i] * 1.e6 * I_gesamt / L
        T0 = T1
        T1 = float(T0 + float(-R_Segment1[i-1]) * 1.e6 * I_gesamt)
        N = math.sqrt(T1 / T0)
        AddSegment(N, L)

    for i in xrange(15, 18):
        E1 = E2
        E2 = E3
        E3 = R_Segment1[i] * 1.e6 * I_gesamt / L
        T0 = T1
        T1 = float(T0 + float(-R_Segment1[i-1]) * 1.e6 * I_gesamt)
        N = math.sqrt(T1 / T0)
        AddSegment(N, L)

    for i in xrange(18, 21):
        E1 = E2
        E2 = E3
        E3 = R_Segment1[i] * 1.e6 * I_gesamt / L
        T0 = T1
        T1 = float(T0 + float(-R_Segment1[i-1]) * 1.e6 * I_gesamt)
        N = math.sqrt(T1 / T0)
        AddSegment(N, L)

    for i in xrange(21, 31):
        E1 = E2
        E2 = E3
        E3 = R_Segment1[i]*1.0e6*I_gesamt / L
        T0 = T1
        T1 = float(T0 + float(-R_Segment1[i-1]) * 1.e6 * I_gesamt)
        N  = math.sqrt(T1 / T0)
        AddSegment(N, L)

    for i in xrange(31, 76):
        E1 = E2
        E2 = E3
        E3 = R_Segment1[i]*1.0e6*I_gesamt / L
        T0 = T1
        T1 = float(T0 + float(-R_Segment1[i-1]) * 1.e6 * I_gesamt)
        N  = math.sqrt(T1 / T0)
        AddSegment(N, L)

    for i in xrange(76, 97):
        E1 = E2
        E2 = E3
        E3 = R_Segment1[i]*1.0e6*I_gesamt / L
        T0 = T1
        T1 = float(T0 + float(-R_Segment1[i-1]) * 1.e6 * I_gesamt)
        N  = math.sqrt(T1 / T0)
        AddSegment(N, L)

    E1 = E2
    E2 = E3
    E3 = 0.
    T0 = T1
    T1 = T1

    # Space = 22.225cm = Drift
    E1 = E2
    E2 = E3
    E3 = R_Segment2[0] * 1.0e6 * I_gesamt / L
    T0 = T1
    T1 = T1
    N  = math.sqrt(T1 / T0)
    AddSegment(N, 22.225 / 100.0)

    ### LE2 : 1..97
    for i in xrange(1, 97):
        E1 = E2
        E2 = E3
        E3 = R_Segment2[i] * 1.0e6 * I_gesamt / L
        T0 = T1
        T1 = float(T0 + float(-R_Segment2[i-1]) * 1.e6 * I_gesamt)
        N = math.sqrt(T1 / T0)
        AddSegment(N, L)

    E1 = E2
    E2 = E3
    E3 = 0.
    T0 = T1
    T1 = T1

    ### HE1 : R und I
    R_gesamt = 111510.0e6

    #Gesamtstrom durch ACC rechte Seite
    I_gesamt = float(v_gesamt / R_gesamt)

    # Terminal / Stripping = Drift
    E1 = E2
    E2 = E3
    E3 = R_Segment3[0] * 1.0e6 * I_gesamt / L
    T0 = T1
    T1 = T1
    N = math.sqrt(T1 / T0)
    AddSegment(N, 121.92 / 100.0)

    # HE1 : 1..97
    for i in xrange(1, 97):
        E1 = E2
        E2 = E3
        E3 = R_Segment3[i] * 1.0e6 * I_gesamt / L
        T0 = T1
        T1 = float(T0 + q * float(R_Segment3[i-1]) * 1.0e6 * I_gesamt)
        N = math.sqrt(T1 / T0)
        AddSegment(N, L)
        #AddSegment2(T0, T1, D, L, E1, E2, E3, q)

    E1 = E2
    E2 = E3
    E3 = 0.
    T0 = T1
    T1 = T1

    # Space = 22.225cm = Drift
    E1 = E2
    E2 = E3
    E3 = R_Segment4[0] * 1.0e6 * I_gesamt / L
    T0 = T1
    T1 = T1
    N = math.sqrt(T1 / T0)
    AddSegment(N, 22.225 / 100.)
    #AddSegment2(T0, T1, D, 22.225 / 100.0, E1, E2, E3, q)

    # HE2 : 1..97
    for i in xrange(1, 97):
        E1 = E2
        E2 = E3
        E3 = R_Segment4[i] * 1.0e6 * I_gesamt / L
        T0 = T1
        T1 = float(T0 + q * float(R_Segment4[i-1]) * 1.0e6 * I_gesamt)
        N = math.sqrt(T1 / T0)
        AddSegment(N, L)
        #AddSegment2(T0, T1, D, L, E1, E2, E3, q)

    E1 = E2
    E2 = E3
    E3 = 0.
    T0 = T1
    T1 = T1

    # Letzte Elektrode..Flansch aussen = Drift
    E1 = E2
    E2 = E3
    E3 = 0.0
    T0 = T1
    T1 = T1
    #AddSegment(1, 117.738 / 100.0)
    AddSegment2(T0, T1, D, 117.738 / 100., E1, E2, E3, q)


# ACC Segment einfuegen nach Hinterberger
def AddSegment(N, L):
    global s, geo_x, geo_y
    if (s > -0.5):
        geo_s.InsertNextValue(s)
        geo_s.InsertNextValue(s+L)
        geo_y.InsertNextValue(55)
        geo_y.InsertNextValue(55)
        s = s + L

    try:
        m11  = (3. - N) / 2.
        m12  = 2. * L / (1. + N)
        m21  = -3. / 8. / N / N / L * (N * N - 1.) * (N - 1.)
        m22  = (3. * N - 1.) / (2. * N * N)
        dkdk = 1. / N / N
        dmdm = 1.

        AddMatrix(1, [
            m11, m12, 0., 0., 0., 0.,
            m21, m22, 0., 0., 0., 0.,
            0., 0., m11, m12, 0., 0.,
            0., 0., m21, m22, 0., 0.,
            0., 0., 0., 0., dkdk, 0.,
            0., 0., 0., 0., 0., dmdm], L)
    except:
        print "Division durch Null! (AddSegment)"


# Korrektur Apertur
def AddSegment1(T0, T1, q, D1, D2, L, b=.57):
    global s, geo_x, geo_y
    if (s > -0.5):
        geo_s.InsertNextValue(s)
        geo_s.InsertNextValue(s)
        geo_s.InsertNextValue(s+L)
        geo_s.InsertNextValue(s+L)
        geo_y.InsertNextValue(D1*1000./2.)
        geo_y.InsertNextValue(55)
        geo_y.InsertNextValue(55)
        geo_y.InsertNextValue(D2*1000./2.)
        s = s + L

    T0 = float(T0)
    T1 = float(T1)
    q  = float(q)
    L  = float(L)
    D1 = float(D1)
    D2 = float(D2)

    E = (T1-T0) / q / L
    N = math.sqrt(T1/T0)

    f1  = -q * (0. - E) / (2. * T0 * (1. + math.sqrt(1. - q * (0. - E) * b * D1 / T0)))  # 1/f1
    f2  = -q * (E - 0.) / (2. * T1 * (1. + math.sqrt(1. - q * (E - 0.) * b * D2 / T1)))  # 1/f2
    LL  = 2. * L / (1. + N)
    m11 = 1. - LL * f1
    m12 = LL
    m21 = -f2 - (1./N - f2*LL)*f1
    m22 = 1/N - f2*LL

    dkdk = 1./N/N
    dmdm = 1.

    AddMatrix(1, [
        m11, m12, 0., 0., 0., 0.,
        m21, m22, 0., 0., 0., 0.,
        0., 0., m11, m12, 0., 0.,
        0., 0., m21, m22, 0., 0.,
        0., 0., 0., 0., dkdk, 0.,
        0., 0., 0., 0., 0., dmdm], L)


# Korrektur Apertur
def AddSegment11(T0, T1, D1, D2, L, b1=.57, b2=.57, f1=True, f2=True):
    global s, geo_x, geo_y
    if (s > -0.5):
        geo_s.InsertNextValue(s)
        geo_s.InsertNextValue(s)
        geo_s.InsertNextValue(s+L)
        geo_s.InsertNextValue(s+L)
        geo_y.InsertNextValue(D1*1000./2.)
        geo_y.InsertNextValue(55)
        geo_y.InsertNextValue(55)
        geo_y.InsertNextValue(D2*1000./2.)
        s = s + L

    T0 = float(T0)
    T1 = float(T1)
    L  = float(L)
    D1 = float(D1)
    D2 = float(D2)
    deltaT = T1-T0

    try:
        N = math.sqrt(T1/T0)
    except:
        print "T0 must not be zero!"

    if f1:
        f1  = + deltaT / (L * 2. * T0 * (1. + math.sqrt(1. + deltaT * b1 * D1 / T0 / L)))  # 1/f1
    else:
        f1 = 0.
    if f2:
        f2  = - deltaT / (L * 2. * T1 * (1. + math.sqrt(1. - deltaT * b2 * D2 / T1 / L)))  # 1/f2
    else:
        f2 = 0.

    LL  = 2. * L / (1. + N)
    m11 = 1. - LL * f1
    m12 = LL
    m21 = -f2 - (1./N - f2*LL)*f1
    m22 = 1/N - f2*LL

    dkdk = 1./N/N
    dmdm = 1.

    AddMatrix(1, [
        m11, m12, 0., 0., 0., 0.,
        m21, m22, 0., 0., 0., 0.,
        0., 0., m11, m12, 0., 0.,
        0., 0., m21, m22, 0., 0.,
        0., 0., 0., 0., dkdk, 0.,
        0., 0., 0., 0., 0., dmdm], L)


# ACC Segment einfuegen nach Elkind und T. Joy. Nur Eintritt beruecksichtigen.
def AddSegment2(T0, T1, D, L, E1, E2, E3, q):
    global s, geo_x, geo_y
    if (s > -0.5):
        geo_s.InsertNextValue(s)
        geo_y.InsertNextValue(D*1000./2.)
        geo_s.InsertNextValue(s)
        geo_y.InsertNextValue(55)
        geo_s.InsertNextValue(s+L)
        geo_y.InsertNextValue(55)
        geo_s.InsertNextValue(s+L)
        geo_y.InsertNextValue(D*1000./2.)
        s = s + L

    N = float(math.sqrt(T1/T0))
    D = float(D)*100.

    f1 = -q*(E1-E2)/2.0/T0/(1.0+math.sqrt(1.0-q*(E1-E2)/T0*0.57*D))  # 1/f1
    LL = 2.0*L/(1.0+N)

    # Nur der Eintritt und die Beschleunigung werden beruecksichtigt:
    m11 = 1.0-LL*f1
    m12 = LL
    m21 = -f1/N
    m22 = 1/N

    dkdk = 1.0/N/N
    dmdm = 1.0

    AddMatrix(1, [
        m11, m12, 0., 0., 0., 0.,
        m21, m22, 0., 0., 0., 0.,
        0., 0., m11, m12, 0., 0.,
        0., 0., m21, m22, 0., 0.,
        0., 0., 0., 0., dkdk, 0.,
        0., 0., 0., 0., 0., dmdm], L)


def AddSegment2VBFN(T0, T1, D, L, E1, E2, E3, q):
    global s, geo_x, geo_y
    if (s > -0.5):
        geo_s.InsertNextValue(s)
        geo_y.InsertNextValue(D*1000./2.)
        geo_s.InsertNextValue(s)
        geo_y.InsertNextValue(55)
        geo_s.InsertNextValue(s+L)
        geo_y.InsertNextValue(55)
        geo_s.InsertNextValue(s+L)
        geo_y.InsertNextValue(D*1000./2.)
        s = s + L

    N = float(math.sqrt(T1/T0))
    D = float(D) * 100.

    f1 = -q * (E1-E2) / (2. * T0 * (1. + math.sqrt(1. - q * (E1-E2) * 1.13 * D / T0)))  # 1/f1
    f2 = -q * (E2-E3) / (2. * T1 * (1. + math.sqrt(1. - q * (E2-E3) * 1.13 * D / T1)))  # 1/f2
    LL = 2. * L / (1. + N)

    m11 = 1.0-LL*f1
    m12 = LL
    m21 = -f2-(1.0/N-f2*LL)*f1
    m22 = 1/N-f2*LL

    dkdk = 1.0/N/N
    dmdm = 1.0

    AddMatrix(1, [
        m11, m12, 0., 0., 0., 0.,
        m21, m22, 0., 0., 0., 0.,
        0., 0., m11, m12, 0., 0.,
        0., 0., m21, m22, 0., 0.,
        0., 0., 0., 0., dkdk, 0.,
        0., 0., 0., 0., 0., dmdm], L)


def AddDrift(num, gamma2, length):
    global s, geo_s, geo_y
    global lastFunction
    lastFunction = "AddDrift"

    if (s > -0.5):
        geo_s.InsertNextValue(s)
        geo_s.InsertNextValue(s + length)
        geo_y.InsertNextValue(55)
        geo_y.InsertNextValue(55)
        s = s + length

    optic.AddDrift(
        ctypes.c_int(int(num)),
        ctypes.c_double(float(gamma2)),
        ctypes.c_double(float(length)))


def Drift(length):
    global s, geo_s, geo_y
    global lastFunction
    lastFunction = "AddDrift"

    num = 1
    gamma2 = 1.

    if (s > -0.5):
        geo_s.InsertNextValue(s)
        geo_s.InsertNextValue(s + length)
        geo_y.InsertNextValue(55)
        geo_y.InsertNextValue(55)
        s = s + length

    optic.AddDrift(
        ctypes.c_int(int(num)),
        ctypes.c_double(float(gamma2)),
        ctypes.c_double(float(length)))


def AddBeamProfile():
    global PROFILEINDEX
    PROFILEINDEX += 1
    optic.AddBeamProfile(PROFILEINDEX)
BeamProfile = AddBeamProfile


def AddWaist():
    optic.AddWaist()
Waist = AddWaist


def ErrFkt(_param_, beamline, INPUT, source):
    optic.Clear()
    ExecText(beamline, INPUT=INPUT, SOURCE1=SOURCE, _param_=_param_)
    optic.CalculateTrajectories()

    return float(optic.GetSpotSize())


def AddModifyEmittance(factor1, factor2):
    """ Aendern der Emittanz um einen Faktor (z.B. beim Strippingprozess) """
    optic.AddModifyEmittance(
        ctypes.c_double(float(factor1)),
        ctypes.c_double(float(factor2)))
ModifyEmittance = AddModifyEmittance


def ChangeBeamParameters(dk=0., dm=0., strag_k=0., strag_m=0., strag_x=0., strag_y=0., strag_dx=0., strag_dy=0.):
    """ Z. B. fuer Folie in der Beamline """
    optic.ChangeBeamParameters(
        ctypes.c_double(float(dk)),
        ctypes.c_double(float(dm)),
        ctypes.c_double(float(strag_k)),
        ctypes.c_double(float(strag_m)),
        ctypes.c_double(float(strag_x)),
        ctypes.c_double(float(strag_y)),
        ctypes.c_double(float(strag_dx)),
        ctypes.c_double(float(strag_dy)),
        ctypes.c_double(1.))


def AddFoil(dk=0., strag_k=0., strag_phi=0.):
    """ Fuer Degraderfolie """
    optic.ChangeBeamParameters(
        ctypes.c_double(float(dk)),
        ctypes.c_double(0.),
        ctypes.c_double(float(strag_k)),
        ctypes.c_double(0.),
        ctypes.c_double(0.),
        ctypes.c_double(0.),
        ctypes.c_double(float(strag_phi)),
        ctypes.c_double(float(strag_phi)),
        ctypes.c_double(.5))
Foil = AddFoil

def ChangeBeamParameters2(dk=0., dm=0., strag_k_over_E=0., strag_m=0.):
    """ Z. B. fuer Folie in der Beamline """
    optic.ChangeBeamParameters2(
        ctypes.c_double(float(dk/1000.)),
        ctypes.c_double(float(dm/1000.)),
        ctypes.c_double(float(strag_k_over_E)),
        ctypes.c_double(float(strag_m)))


def AddSlit(x, dx, y, dy):
    global lastFunction
    lastFunction = "AddSlit"

    global s, geo_s, geo_y
    if (s > -0.5):
        geo_s.InsertNextValue(s)
        geo_y.InsertNextValue(x+dx/2.)
        geo_s.InsertNextValue(s+.05)
        geo_y.InsertNextValue(x+dx/2.)
        geo_s.InsertNextValue(s+.05)
        geo_y.InsertNextValue(x-dx/2.)
        geo_s.InsertNextValue(s)
        geo_y.InsertNextValue(x-dx/2.)
        geo_s.InsertNextValue(s)
        geo_y.InsertNextValue(-20.)
    optic.AddSlit(
        ctypes.c_double(float(x)),
        ctypes.c_double(float(dx)),
        ctypes.c_double(float(y)),
        ctypes.c_double(float(dy)))
Slit = AddSlit


def AddThinLens(arg1, arg2, geo=25):
    global s, geo_s, geo_y
    global lastFunction
    lastFunction = "AddThinLens"

    if (s > -0.5):
        geo_s.InsertNextValue(s)
        geo_y.InsertNextValue(geo)

    optic.AddThinLens(
        ctypes.c_int(1),
        ctypes.c_double(float(arg1)),
        ctypes.c_double(float(arg2)),
        ctypes.c_double(0.))
ThinLens = AddThinLens


def AddEinzelLens(f, geo=25.):
    AddThinLens(f, f, geo)
EinzelLens = AddEinzelLens


def AddAMSSO110EL(phi1, v_el):
    global lastFunction
    lastFunction = "AddAMSSO110EL"

    global s, geo_s, geo_y
    if (s > -0.5):
        geo_s.InsertNextValue(s)
        geo_y.InsertNextValue(50)

    # phi2/phi1
    x = (phi1 - v_el) / float(phi1)

    #f = .001/(0.022648131734415336-0.13770935811064514*x+0.44804113697363235*x**2-0.9575118663530009*x**3+1.3523795047227753*x**4-1.2071207469218252*x**5+0.6161254239513653*x**6-0.13696307518260592*x**7)#reihenentwicklung fuer x=.3 ... .8

    #reihenentwicklung fuer x=.2 ... .7
    f = .001 / (
        + 0.024838874085440145
        - 0.1721716204656418 * x
        + 0.6764569821500281 * x**2
        - 1.7844240368333488 * x**3
        + 3.118746647284872  * x**4
        - 3.434443885979744  * x**5
        + 2.152129217984039  * x**6
        - 0.5841534279352003 * x**7)

    optic.AddThinLens(
        ctypes.c_int(1),
        ctypes.c_double(float(f)),
        ctypes.c_double(float(f)),
        ctypes.c_double(0.))
AMSSO110EL = AddAMSSO110EL


def AddAMSBIEL(phi1, v_el):
    global lastFunction
    lastFunction = "AddAMSBIEL"

    global s, geo_s, geo_y
    if (s > -0.5):
        geo_s.InsertNextValue(s)
        geo_y.InsertNextValue(56)

    # phi2/phi1
    x = (phi1 - v_el) / phi1

    #reihenentwicklung fuer x=.2 ... .7
    f = .001 / (
        + 0.01993277008153473
        - 0.13793678662472572 * x
        + 0.5408115698616647  * x**2
        - 1.4238381952266255  * x**3
        + 2.4846514325496316  * x**4
        - 2.7328746285908734  * x**5
        + 1.7109251000894101  * x**6
        - 0.4640648127510327  * x**7)

    optic.AddThinLens(
        ctypes.c_int(1),
        ctypes.c_double(float(f)),
        ctypes.c_double(float(f)),
        ctypes.c_double(0.))
AMSBIEL = AddAMSBIEL


def AddFNEL(phi1, v_el):
    global lastFunction
    lastFunction = "AddFNEL"

    global s, geo_s, geo_y
    if (s > -0.5):
        geo_s.InsertNextValue(s)
        geo_y.InsertNextValue(63)


    # phi2/phi1
    x = (phi1 - v_el) / phi1

    #reihenentwicklung fuer x=.2 ... .7
    f = 1. / (
        + 11.542244476708195
        - 61.57930337802449  * x
        + 183.54364546860486 * x**2
        - 385.71960791043387 * x**3
        + 562.1897241491755  * x**4
        - 529.9264768291281  * x**5
        + 288.3485380347237  * x**6
        - 68.5022805970122   * x**7)

    """
    x = (phi1 + v_el) / phi1
    # reihenentwicklung fuer x=1..2
    f = 1. / (
        + 6.76764
        - 22.3617   * x
        + 30.8469   * x**2
        - 24.4108   * x**3
        + 12.5201   * x**4
        - 4.05526   * x**5
        + 0.754707  * x**6
        - 0.0615853 * x**7)
    """
    optic.AddThinLens(
        ctypes.c_int(1),
        ctypes.c_double(float(f)),
        ctypes.c_double(float(f)),
        ctypes.c_double(0.))
FNEL = AddFNEL


def AddFNSputterEL(phi1, v_el):
    global lastFunction
    lastFunction = "AddFNSputterEL"

    global s, geo_s, geo_y
    if (s > -0.5):
        geo_s.InsertNextValue(s)
        geo_y.InsertNextValue(31)

    # phi2/phi1
    x = (phi1 - v_el) / phi1

    if (x < .44) or (x > .57):
        print "WARNUNG: Sputterlinse nicht im optimalen Bereich!\nx = ", x, "\n"

    # reihenentwicklung fuer x=.44...57 nach Messung
    f = 1. / (
        - 2.14464643
        + 38.80407837 * x
        - 58.46349099 * x**2)

    optic.AddThinLens(
        ctypes.c_int(1),
        ctypes.c_double(float(f)),
        ctypes.c_double(float(f)),
        ctypes.c_double(0.))
FNSputterEL = AddFNSputterEL


def AddFNSputterEL_SIMION(phi1, v_el):
    global lastFunction
    lastFunction = "AddFNSputterEL_SIMION"


def AddESD(num, arg1, arg2, arg3, arg4, arg5, geo, korrektur=None, V_ist=1., V_soll=1.):
    global lastFunction
    lastFunction = "AddESD"

    if not korrektur:
        korrektur = float(V_ist) / float(V_soll)

    global s, geo_s, geo_y
    if (s > -0.5):
        geo_s.InsertNextValue(s)
        geo_s.InsertNextValue(s+num*arg2*arg3)
        geo_y.InsertNextValue(geo)
        geo_y.InsertNextValue(geo)
        s = s + num * arg2 * arg3

    optic.AddESD(
        ctypes.c_int(int(num)),
        ctypes.c_double(float(arg1)),
        ctypes.c_double(float(arg2)),
        ctypes.c_double(float(arg3)),
        ctypes.c_double(float(arg4)),
        ctypes.c_double(float(arg5)),
        ctypes.c_double(float(korrektur)))


def AddEdgeFocusing(r, beta, K, geo):
    global lastFunction
    lastFunction = "AddEdgeFocusing"

    g = 2. * geo / 1.e3

    betaeff = beta - g / r * K * (1. + math.sin(beta) * math.sin(beta)) / (math.cos(beta))

    optic.AddEdgeFocusing(
        ctypes.c_int(1),
        ctypes.c_double(float(r)),
        ctypes.c_double(float(beta)),
        ctypes.c_double(float(betaeff)))


def AddMSA(num, arg1, arg2, arg3, geo, korrektur=None, B_ist=1., B_soll=1.):
    global lastFunction
    lastFunction = "AddMSA"

    if not korrektur:
        korrektur = float(B_ist) / float(B_soll)

    global s, geo_s, geo_y
    if (s > -0.5):
        geo_s.InsertNextValue(s)
        geo_s.InsertNextValue(s+num*arg2*arg3)
        geo_y.InsertNextValue(geo)
        geo_y.InsertNextValue(geo)
        s = s + num * arg2 * arg3

    optic.AddHomDeflectingMagnet(
        ctypes.c_int(int(num)),
        ctypes.c_double(float(arg1)),
        ctypes.c_double(float(arg2)),
        ctypes.c_double(float(arg3)),
        ctypes.c_double(float(korrektur)))


def ESD(alpha, arg3, arg4, geo=25., korrektur=None, V_ist=1., V_soll=1.):

    global lastFunction
    lastFunction = "AddESD"

    alpha = math.radians(alpha)
    beta0 = 0.
    num = 10

    if not korrektur:
        korrektur = float(V_ist) / float(V_soll)

    global s, geo_s, geo_y
    if (s > -0.5):
        geo_s.InsertNextValue(s)
        geo_s.InsertNextValue(s + alpha * arg3)
        geo_y.InsertNextValue(geo)
        geo_y.InsertNextValue(geo)
        s = s + alpha * arg3

    optic.AddESD(
        ctypes.c_int(int(num)),
        ctypes.c_double(float(1)),
        ctypes.c_double(float(alpha)),
        ctypes.c_double(float(arg3)),
        ctypes.c_double(float(arg4)),
        ctypes.c_double(float(beta0)),
        ctypes.c_double(float(korrektur)))


def EdgeFocusing(r, beta, K=.45, geo=30.):

    global lastFunction
    lastFunction = "AddEdgeFocusing"

    g = 2. * geo / 1.e3
    beta = math.radians(beta)

    betaeff = beta - g / r * K * (1. + math.sin(beta) * math.sin(beta)) / (math.cos(beta))

    optic.AddEdgeFocusing(
        ctypes.c_int(1),
        ctypes.c_double(float(r)),
        ctypes.c_double(float(beta)),
        ctypes.c_double(float(betaeff)))


def MSA(r, alpha, geo=30., korrektur=None, B_ist=1., B_soll=1.):

    global lastFunction
    lastFunction = "AddMSA"

    alpha = math.radians(alpha)
    num = 10

    if not korrektur:
        korrektur = float(B_ist) / float(B_soll)

    global s, geo_s, geo_y
    if (s > -0.5):
        geo_s.InsertNextValue(s)
        geo_s.InsertNextValue(s + r * alpha)
        geo_y.InsertNextValue(geo)
        geo_y.InsertNextValue(geo)
        s = s + r * alpha

    optic.AddHomDeflectingMagnet(
        ctypes.c_int(int(num)),
        ctypes.c_double(float(1)),
        ctypes.c_double(float(r)),
        ctypes.c_double(float(alpha)),
        ctypes.c_double(float(korrektur)))


def AddInhomMSA(num, rho, phi, n, geo):
    global s, geo_s, geo_y

    if (s > -0.5):
        geo_s.InsertNextValue(s)
        geo_s.InsertNextValue(s+num*rho*phi)
        geo_y.InsertNextValue(geo)
        geo_y.InsertNextValue(geo)
        s = s + num * rho * phi

    optic.AddInhomDeflectingMagnet(
        ctypes.c_int(int(num)),
        ctypes.c_double(float(rho)),
        ctypes.c_double(float(phi)),
        ctypes.c_double(float(n)))


def AddQuadrupolRadFoc(num, arg1, arg2, arg3, geo):
    global lastFunction, s, geo_s, geo_y
    lastFunction = "AddQuadrupolRadFoc"

    if (s > -0.5):
        geo_s.InsertNextValue(s)
        geo_s.InsertNextValue(s+num*arg3)
        geo_y.InsertNextValue(geo)
        geo_y.InsertNextValue(geo)
        s = s + num * arg3
    optic.AddQuadrupolRadFoc(
        ctypes.c_int(int(num)),
        ctypes.c_double(float(arg1)),
        ctypes.c_double(float(arg2)),
        ctypes.c_double(float(arg3)))


def AddQuadrupolAxFoc(num, arg1, arg2, arg3, geo):
    global lastFunction, s, geo_s, geo_y
    lastFunction = "AddQuadrupolAxFoc"

    if (s > -0.5):
        geo_s.InsertNextValue(s)
        geo_s.InsertNextValue(s+num*arg3)
        geo_y.InsertNextValue(geo)
        geo_y.InsertNextValue(geo)
        s = s + num * arg3

    optic.AddQuadrupolAxFoc(
        ctypes.c_int(int(num)),
        ctypes.c_double(float(arg1)),
        ctypes.c_double(float(arg2)),
        ctypes.c_double(float(arg3)))


def QuadrupolRadFoc(arg2, arg3, geo=25.):
    global lastFunction, s, geo_s, geo_y
    lastFunction = "AddQuadrupolRadFoc"

    num = 10

    if (s > -0.5):
        geo_s.InsertNextValue(s)
        geo_s.InsertNextValue(s + arg3)
        geo_y.InsertNextValue(geo)
        geo_y.InsertNextValue(geo)
        s = s + arg3
    optic.AddQuadrupolRadFoc(
        ctypes.c_int(int(num)),
        ctypes.c_double(float(1.)),
        ctypes.c_double(float(arg2)),
        ctypes.c_double(float(arg3)))


def QuadrupolAxFoc(arg2, arg3, geo=25.):

    global lastFunction, s, geo_s, geo_y
    lastFunction = "AddQuadrupolAxFoc"

    num = 10

    if (s > -0.5):
        geo_s.InsertNextValue(s)
        geo_s.InsertNextValue(s + arg3)
        geo_y.InsertNextValue(geo)
        geo_y.InsertNextValue(geo)
        s = s + arg3

    optic.AddQuadrupolAxFoc(
        ctypes.c_int(int(num)),
        ctypes.c_double(float(1.)),
        ctypes.c_double(float(arg2)),
        ctypes.c_double(float(arg3)))


def AddAMSQPT_XYX(gamma2, prozent, astigm, v_terminal, v_ext, q, geo):
    global lastFunction
    lastFunction = "AddAMSQPT_XYX"

    global s, geo_s, geo_y
    if (s > -0.5):
        geo_s.InsertNextValue(s)
        geo_y.InsertNextValue(geo)
        geo_s.InsertNextValue(s+.25)
        geo_y.InsertNextValue(geo)
        s = s+.25
        geo_s.InsertNextValue(s)
        geo_y.InsertNextValue(55)
        geo_s.InsertNextValue(s+.076)
        geo_y.InsertNextValue(55)
        s = s+.076
        geo_s.InsertNextValue(s)
        geo_y.InsertNextValue(geo)
        geo_s.InsertNextValue(s+.5)
        geo_y.InsertNextValue(geo)
        s = s+.5
        geo_s.InsertNextValue(s)
        geo_y.InsertNextValue(55)
        geo_s.InsertNextValue(s+.076)
        geo_y.InsertNextValue(55)
        s = s+.076
        geo_s.InsertNextValue(s)
        geo_y.InsertNextValue(geo)
        geo_s.InsertNextValue(s+.25)
        geo_y.InsertNextValue(geo)
        s = s+.25

    #empirisch bestimmt
    vx = 30.e3 * (prozent + astigm) / 100.
    vy = 30.e3 * (prozent - astigm) / 100.
    T  = v_ext + (q + 1) * float(v_terminal)

    kx  = abs(vx) / float(geo * 1.e-3)**2 * q / T
    k   = (vx + vy) / 2. / float(geo * 1.e-3)**2 * q / T

    optic.AddQuadrupolRadFoc(
        ctypes.c_int(10),
        ctypes.c_double(float(gamma2)),
        ctypes.c_double(float(k)),
        ctypes.c_double(float(.025)))
    optic.AddDrift(
        ctypes.c_int(1),
        ctypes.c_double(float(gamma2)),
        ctypes.c_double(.076))
    optic.AddQuadrupolAxFoc(
        ctypes.c_int(10),
        ctypes.c_double(float(gamma2)),
        ctypes.c_double(float(kx)),
        ctypes.c_double(float(.05)))
    optic.AddDrift(
        ctypes.c_int(1),
        ctypes.c_double(float(gamma2)),
        ctypes.c_double(.076))
    optic.AddQuadrupolRadFoc(
        ctypes.c_int(10),
        ctypes.c_double(float(gamma2)),
        ctypes.c_double(float(k)),
        ctypes.c_double(float(.025)))


def AddAMSQPT_YXY(gamma2, prozent, astigm, v_terminal, v_ext, q, geo):
    global lastFunction
    lastFunction = "AddAMSQPT_YXY"

    global s, geo_s, geo_y
    if (s > -0.5):
        geo_s.InsertNextValue(s)
        geo_y.InsertNextValue(geo)
        geo_s.InsertNextValue(s+.25)
        geo_y.InsertNextValue(geo)
        s = s+.25
        geo_s.InsertNextValue(s)
        geo_y.InsertNextValue(55)
        geo_s.InsertNextValue(s+.076)
        geo_y.InsertNextValue(55)
        s = s+.076
        geo_s.InsertNextValue(s)
        geo_y.InsertNextValue(geo)
        geo_s.InsertNextValue(s+.5)
        geo_y.InsertNextValue(geo)
        s = s+.5
        geo_s.InsertNextValue(s)
        geo_y.InsertNextValue(55)
        geo_s.InsertNextValue(s+.076)
        geo_y.InsertNextValue(55)
        s = s+.076
        geo_s.InsertNextValue(s)
        geo_y.InsertNextValue(geo)
        geo_s.InsertNextValue(s+.25)
        geo_y.InsertNextValue(geo)
        s = s+.25

    vx  = 30.e3 * (prozent + astigm) / 100.
    vy  = 30.e3 * (prozent - astigm) / 100.
    T   = v_ext + (q + 1) * float(v_terminal)

    kx  = abs(vx)/float(geo*1.e-3)**2*q/T
    k   = (vx+vy)/2./float(geo*1.e-3)**2*q/T

    optic.AddQuadrupolAxFoc(
        ctypes.c_int(10),
        ctypes.c_double(float(gamma2)),
        ctypes.c_double(float(k)),
        ctypes.c_double(float(.025)))
    optic.AddDrift(
        ctypes.c_int(1),
        ctypes.c_double(float(gamma2)),
        ctypes.c_double(.076))
    optic.AddQuadrupolRadFoc(
        ctypes.c_int(10),
        ctypes.c_double(float(gamma2)),
        ctypes.c_double(float(kx)),
        ctypes.c_double(float(.05)))
    optic.AddDrift(
        ctypes.c_int(1),
        ctypes.c_double(float(gamma2)),
        ctypes.c_double(.076))
    optic.AddQuadrupolAxFoc(
        ctypes.c_int(10),
        ctypes.c_double(float(gamma2)),
        ctypes.c_double(float(k)),
        ctypes.c_double(float(.025)))


def AddAMSQPT(gamma2, prozent, astigm, v_terminal, v_ext, q, geo):
    global lastFunction
    lastFunction = "AddAMSQPT"

    AddAMSQPT_XYX(gamma2, prozent, astigm, v_terminal, v_ext, q, geo)


def AddGeo(s1, x1, s2, x2):
    """ zum zeichnen eigener geometrie """
    global lastFunction
    lastFunction = "AddGeo"

    global s, geo_s, geo_y
    if (s > -0.5):
        geo_s.InsertNextValue(s+s1)
        geo_s.InsertNextValue(s+s2)
        geo_y.InsertNextValue(x1)
        geo_y.InsertNextValue(x2)


def AddQPT_XYX(gamma2, vx, vy, v_terminal, v_ext, q, geo):
    global s, geo_s, geo_y
    vx = float(vx)
    vy = float(vy)
    T  = v_ext + (q + 1) * float(v_terminal)

    kx = abs(vx * 1.e3) / float(geo * 1.e-3)**2 * q / T
    ky = abs(vy * 1.e3) / float(geo * 1.e-3)**2 * q / T

    optic.AddQuadrupolAxFoc(
        ctypes.c_int(10),
        ctypes.c_double(float(gamma2)),
        ctypes.c_double(float(ky)),
        ctypes.c_double(float(.025)))
    optic.AddDrift(
        ctypes.c_int(1),
        ctypes.c_double(float(gamma2)),
        ctypes.c_double(.076))
    optic.AddQuadrupolRadFoc(
        ctypes.c_int(10),
        ctypes.c_double(float(gamma2)),
        ctypes.c_double(float(kx)),
        ctypes.c_double(float(.05)))
    optic.AddDrift(
        ctypes.c_int(1),
        ctypes.c_double(float(gamma2)),
        ctypes.c_double(.076))
    optic.AddQuadrupolAxFoc(
        ctypes.c_int(10),
        ctypes.c_double(float(gamma2)),
        ctypes.c_double(float(ky)),
        ctypes.c_double(float(.025)))


###########################################################


def GetTrajectories():
    particlenum  = int(optic.GetParticleNum())
    particlesize = int(optic.GetParticleSize())
    trajsize     = int(optic.GetTrajectoriesSize())

    # erstelle ein Array vom C-Type double*
    trajectories = (ctypes.c_double * trajsize)()

    optic.GetTrajectories(trajectories)

    traj = array.array('d', trajectories)
    return (traj, particlenum, particlesize)


def GetTrajectory(iparticle, iproperty):
    traj = (ctypes.c_double * optic.GetTrajectorySize())()
    optic.GetTrajectory(ctypes.c_int(iparticle), ctypes.c_int(iproperty), traj)
    traj = array.array('d', traj)
    return traj


def PrintTrajectories(traj):
    (trajectories, particlenum, particlesize) = GetTrajectories()
    for i in xrange(0, particlenum * particlesize, particlesize):
        for j in xrange(0, len(trajectories), particlenum * particlesize):
            for k in xrange(0, particlesize, 1):
                print "%f" % (trajectories[j+i+k]),
            print "\n",
        print "\n",


def ExecText(text, INPUT=None, SOURCE1=None, _param_=None):
    global textArray, lastFunction, PROFILEINDEX

    textArray = []
    lastFunction = "unknown"
    PROFILEINDEX = -1

    global SOURCE
    if len(SOURCE1) > 0:
        SOURCE = SOURCE1
    exec(text)


def Limioptic(filename):
    optic.Clear()
    execfile(filename)
    optic.CalculateTrajectories()
    traj = GetTrajectories()
    PrintTrajectories(traj)


def Name(text):
    textArray.append([s, text])


try:
    winver = sys.winver
except:
    winver = None

if winver:
    try:
        optic = ctypes.CDLL("./liblimioptic-win.so")
    except:
        print "\n\nFirst start?\nTrying to compile C++ code. ",
        from os import system
        print "\nThis is a windows OS (right?). So this will only work if MinGW is installed!\n(Make sure that X:\\..\\MinGW\\bin is in your PATH)\n"
        system("PAUSE")
        print "compiling climioptic.cpp"
        print "------------------------"
        system("g++ -std=c++0x -Wall -fPIC -O -c climioptic.cpp")
        print "\n"

        print "compiling limioptic.cpp"
        print "------------------------"
        system("g++ -Wall -fPIC -O -c limioptic.cpp")
        print "\n"

        print "building library"
        print "------------------------"
        system("g++ -fPIC -shared -Wl,-soname,liblimioptic-win.so -O -o liblimioptic-win.so limioptic.o climioptic.o")
        print "completed."

        optic = ctypes.CDLL("./liblimioptic-win.so")

else:
    try:
        optic = ctypes.CDLL("./liblimioptic-linux.so")
    except:
        print "\n\nFirst start?\nTrying to compile C++ code. ",
        from os import system

        print "compiling climioptic.cpp"
        print "------------------------"
        system("g++ -std=c++0x -Wall -fPIC -O -c climioptic.cpp")
        print "\n"

        print "compiling limioptic.cpp"
        print "------------------------"
        system("g++ -Wall -fPIC -O -c limioptic.cpp")
        print "\n"

        print "building library"
        print "------------------------"
        system("g++ -fPIC -shared -Wl,-soname,liblimioptic-linux.so -O -o liblimioptic-linux.so limioptic.o climioptic.o")
        print "completed."

        optic = ctypes.CDLL("./liblimioptic-linux.so")


optic.GetSpotSize.restype = ctypes.c_double
optic.GetSigmaX.restype = ctypes.c_double
optic.GetSigmaY.restype = ctypes.c_double

if __name__ == "__main__":
    if (len(sys.argv) != 2):
        print("Usage : %s filename" % (sys.argv[0]))
    else:
        Limioptic(sys.argv[1])

#TestOptFkt