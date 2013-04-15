#!/usr/bin/env python

"""
Hier werden die Funktionen definiert und an den C++ - Teil uebergeben.
"""

import sys
import ctypes
import re
import array
import math
import random
import vtk
#import serial
#import time
#import thread

###########################################################

def AddParticle(p1,p2,p3,p4,p5,p6):
    """ Normales Partikel einfuegen """
    global lastFunction
    lastFunction = "AddParticle"
    optic.AddParticle(ctypes.c_double(float(p1)),
        ctypes.c_double(float(p2)),ctypes.c_double(float(p3)),
        ctypes.c_double(float(p4)),ctypes.c_double(float(p5)),
        ctypes.c_double(float(p6)))
        

def AddBeamX(xmax,amax,ymax,bmax,dk,dm,delta):
    """ Einfachen 3d-Strahl einfuegen """
    global lastFunction
    lastFunction = "AddBeamX"
    for i in range(0,360,int(delta)):
        degtorad=i*math.pi/180
        
        optic.AddParticle(ctypes.c_double(float(xmax*math.cos(degtorad))),
            ctypes.c_double(float(amax*math.sin(degtorad))),
            ctypes.c_double(0.),
            ctypes.c_double(0.),
            ctypes.c_double(float(dk)),
            ctypes.c_double(float(dm)))         
        optic.AddParticle(ctypes.c_double(0.),
            ctypes.c_double(0.),
            ctypes.c_double(float(ymax*math.cos(degtorad))),
            ctypes.c_double(float(bmax*math.sin(degtorad))),
            ctypes.c_double(float(dk)),
            ctypes.c_double(float(dm)))
            

def AddBeam3d(xmax,amax,ymax,bmax,dk,dm,delta):
    """ Komplexen 3d-Strahl einfuegen """
    global lastFunction
    lastFunction = "AddBeam3d"
    for i in range(0,360,int(delta)):
        nu=i*math.pi/180
        for j in range(0,360,10):
            phi = j*math.pi/180
            optic.AddParticle(ctypes.c_double(xmax*math.cos(nu)*math.cos(phi)),
                ctypes.c_double(amax*math.cos(nu)*math.sin(phi)),
                ctypes.c_double(ymax*math.sin(nu)*math.cos(phi)),
                ctypes.c_double(bmax*math.sin(nu)*math.sin(phi)),
                ctypes.c_double(float(dk)),
                ctypes.c_double(float(dm)))                                 

            
def AddBeamRandomGauss(xmax,amax,ymax,bmax,dk,dm,number):
    """ Gaussverteilten zufallsbasierten Strahl einfuegen """
    global lastFunction
    lastFunction = "AddBeamRandomGauss"
    random.seed()
    for i in range(0,number,1):
        try:    
            #phix = random.gauss(math.pi/4.,math.pi/2.)
            #phiy = random.gauss(math.pi/4.,math.pi/2.)
            x   =   random.gauss(0,xmax/3.)
            y   =   random.gauss(0,ymax/3.)
            a   =   random.gauss(0,amax/3.)
            b   =   random.gauss(0,bmax/3.)
            #a = 0.
            #b = 0.
            #if (abs(x)<=1.):   a   = math.sin(math.acos(x))
            #if (abs(y)<=1.):   b   = math.sin(math.acos(y))
            
            optic.AddParticle(ctypes.c_double(x),
                ctypes.c_double(a),
                ctypes.c_double(y),
                ctypes.c_double(b),
                ctypes.c_double(float(dk)),
                ctypes.c_double(float(dm)))
        except:
            #number = number + 1
            print "except in AddBeamRandomGauss()"
            
                
def AddBeam(xmax,amax,ymax,bmax,dk,dm,delta):
    """ Einfachen 2d-Strahl einfuegen """
    global lastFunction
    lastFunction = "AddBeam"
    for i in range(0,360,int(delta)):
        degtorad=i*math.pi/180
        
        optic.AddParticle(ctypes.c_double(float(xmax*math.cos(degtorad))),
            ctypes.c_double(float(amax*math.sin(degtorad))),
            ctypes.c_double(float(ymax*math.cos(degtorad))),
            ctypes.c_double(float(bmax*math.sin(degtorad))),
            ctypes.c_double(float(dk)),
            ctypes.c_double(float(dm)))     

            
def AddGaussBeam(xmax,dx,dy,dk,dm,delta):
    """ Wird nicht mehr verwendet """
    global lastFunction
    lastFunction = "AddGaussBeam"
    q=0.0
    for i in range(0,rng,1):
        q=q+math.exp(-0.5*i*i/(rng*20.))
        if (q>=1.0):        
            optic.AddParticle(ctypes.c_double(float(i)/rng*xmax+dx),
                ctypes.c_double(0.0),
                ctypes.c_double(float(i)/rng*xmax+dy),
                ctypes.c_double(0.0),
                ctypes.c_double(float(dk)),
                ctypes.c_double(float(dm)))
            optic.AddParticle(ctypes.c_double(-float(i)/rng*xmax+dx),
                ctypes.c_double(0.0),
                ctypes.c_double(-float(i)/rng*xmax+dy),
                ctypes.c_double(0.0),
                ctypes.c_double(float(dk)),
                ctypes.c_double(float(dm)))
            q=0.0
            
def AddSource():
    """ Externe Quelle einfuegen """
    global lastFunction
    lastFunction = "AddSource"
    for i in range(0,len(SOURCE)):
        optic.AddParticle(ctypes.c_double(SOURCE[i][0]),
            ctypes.c_double(SOURCE[i][1]),
            ctypes.c_double(SOURCE[i][2]),
            ctypes.c_double(SOURCE[i][3]),
            ctypes.c_double(SOURCE[i][4]),
            ctypes.c_double(SOURCE[i][5]))  
            
    #print "Anzahl Partikel =", len(SOURCE)
            
            
def AddMatrix(num,mat,length):
    """ Allgemeinte Matrix einfuegen """
    global lastFunction
    lastFunction = "AddMatrix"
    matrixarray=ctypes.c_double*36
    matrix=matrixarray()
    for i in range(0,36):
        matrix[i]=ctypes.c_double(float(mat[i]))
    optic.AddMatrix(ctypes.c_int(int(num)),matrix,ctypes.c_double(float(length)))

    
def AddAMSAcc(v_qsnout, v_gesamt, v_vorbeschl, q):
    """ Beschleunigung Cologne AMS nach HINTERBERGER """
    global lastFunction
    lastFunction = "AddAMSAcc"
    R_a1 = 400.0e6+10.0e3 #Widerstand 1-4 Zweig1
    R_a2 = 400.0e6+2.0e9 #Widerstand 1-4 Zweig2
    R_a = float(1/(1/R_a1+1/R_a2)) #Widerstand 1-4
    R_b = 18.0*40.0e6+62.0*200.0e6+81.0*200.0e6 #Widerstand 5-84
    R_gesamt = R_b+R_a #Gesamtwiderstand HE1
    I_gesamt = float(v_gesamt/R_gesamt) #Gesamtstrom durch ACC linke Seite
    L=float(25.4e-3) #Laenge eines Segmentes
    V_a = float(R_a * I_gesamt) #Spannungsabfall 1..4
    I_a2 = float(V_a / R_a2)
#   V_qsnout = float(I_a2 * (400.0e6+100.0e9))
#   a = s
    
#   print V_a, V_a-1.e9*I_a2, V_a-2.e9*I_a2, V_a-(2.e9+200.e6)*I_a2, V_a-(2.e9+400.e6)*I_a2
    
    # Bias..Anfang Q-Snout
    T0 = float(v_vorbeschl-4.e3)
    T1 = float(v_vorbeschl+v_qsnout)
#   N = math.sqrt(T1/T0)
#   AddSegment(N,223.290615e-3)
#   print "Spannung00 = %1.0f" % (T1)

    # Q-Snout : Laenge 300.4mm
    T0 = T1
    T1 = T1
    N = math.sqrt(T1/T0)
    AddSegment(N,300.4e-3)
    # print "Q-Snout", "Spannung = %1.0f "% (T1-v_vorbeschl), "Potential = %1.0f" % (T1), "E4 = %1.0f"%(float(v_vorbeschl+V_a-1.e9*I_a2))

# HE1 : Ende Q-Snout..Anfang Terminal
# 3..4
    T0 = T1
    T1 = float(v_vorbeschl+V_a-1.e9*I_a2)
    N = math.sqrt(T1/T0)
    AddSegment(N,L)
#   print "Spannung04 = %1.0f" % (T1)
    
# 4..5
    T0 = T1
    T1 = float(v_vorbeschl+V_a)
    N = math.sqrt(T1/T0)
    AddSegment(N,L)
#   print "Spannung05 = %1.0f" % (T1)
    
    deltaV = 40.0e6*I_gesamt #5..23
    for i in range(0,18):
        T0 = T1
        T1 = float(T0+deltaV)
        N = math.sqrt(T1/T0)
        AddSegment(N,L)
        #print "Spannung0%1.0f = %1.0f" % (i+6, T1)
    
    deltaV = 200.0e6*I_gesamt #23..85
    for i in range(0,62):
        T0 = T1
        T1 = float(T0+deltaV)
        N = math.sqrt(T1/T0)
        AddSegment(N,L)
#   print "Spannung85 = %1.0f" % (T1)
    
# Space = 76.23mm
    T0 = T1
    T1 = T1
    N = math.sqrt(T1/T0)
    AddSegment(N,76.23e-3)
#   print "Spannung85 = %1.0f" % (T1)
    
# Space..Terminal
    deltaV = 200.0e6*I_gesamt #85..166
    for i in range(0,81):
        T0 = T1
        T1 = float(T0+deltaV)
        N = math.sqrt(T1/T0)
        AddSegment(N,L)
        
#   print "Terminalspannung = %1.0f" % (T1)
        
# Terminal / Stripping
    T0 = T1
    T1 = T1
    N = math.sqrt(T1/T0)
    AddSegment(N,1147.0e-3) #Terminallaenge in mm

# HE2 : Ende Terminal..letzte Elektrode
    R_gesamt = float(9.0*120.0e6+156.0*300.0e6)
    I_gesamt = float(v_gesamt/R_gesamt) #Gesamtstrom durch ACC rechte Seite

# Terminal..Space
    deltaV = 120.0e6*I_gesamt #1..9
    for i in range(0,9):
        T0 = T1
        T1 = float(T0+q*deltaV)
        N = math.sqrt(T1/T0)
        AddSegment(N,L)

    deltaV = 300.0e6*I_gesamt #9..81
    for i in range(0,72):
        T0 = T1
        T1 = float(T0+q*deltaV)
        N = math.sqrt(T1/T0)
        AddSegment(N,L)

# Space = 76.05 mm 
    T0 = T1
    T1 = T1
    N = math.sqrt(T1/T0)
    AddSegment(N,76.05e-3)

# Space..letzte Elekrode
    deltaV = 300.0e6*I_gesamt #81..165
    for i in range(0,84):
        T0 = T1
        T1 = float(T0+q*deltaV)
        N = math.sqrt(T1/T0)
        AddSegment(N,L)
        
# Letzte Elektrode..Flansch aussen = 487.8mm
    AddSegment(1,447.81e-3)
#   print "Tanklaenge:", s-a
        

    

def AddVBFN(extraktion, deltaV, laenge = .276, b = 1.13, b1 = -1., b2 = -1.):
    """ Vorbeschleunigung FN """
    global lastFunction
    lastFunction = "AddVBFN"
    if ((b1 == -1.) and (b2 == -1.)):
        b1 = b
        b2 = b
    #E1 = -extraktion/0.95
    #E2 = -deltaV/float(laenge)
    #E3 = 0.0
    T0 = float(extraktion)
    T1 = float(extraktion + deltaV)
    #AddSegment(float(math.sqrt(T1/T0)),.28)#
    #AddSegment2VBFN(T0,T1,9./100.0,float(laenge),E1,E2,E3,-1.0)
    AddSegment11(T0, T1, .09, .09, laenge, b1, b2)

    
def AddFNAccNeu(vt, T0, q, b = 0.57, b1 = -1., b2 = -1., D1 = .088, factor1 = 1., factor2 = 1., beamprofile = False):
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
    #print E1,E2,E3,E4
    

    T1 = T0 + E1 * l1
    T2 = T1 + E2 * l2
    T3 = T2 + q * E3 * l3
    T4 = T3 + q * E4 * l4
    #print T1, T2, T3, T4

    #D1 = .088 #LE1
    D2 = .03  #LE1
    D3 = .03  #LE2
    D4 = .03  #LE2
    D5 = .03  #HE1
    D6 = .03  #HE1
    D7 = .03  #HE2
    D8 = .03  #HE2

    AddSegment11(T0,T1,D1,D2,l1,b1,b2)          #LE1
    AddSegment11(T1,T1,D2,D3,space1,b1,b2)      #Space1
    AddSegment11(T1,T2,D3,D4,l2,b1,b2)          #LE2
    AddSegment11(T2,T2,D4,D5,terminal/2.,b1,b2) #Terminal
    if (beamprofile): 
           print "\n"
           AddBeamProfile()
    AddModifyEmittance(factor1, factor2)
    if (beamprofile): AddBeamProfile()
    AddSegment11(T2,T2,D4,D5,terminal/2.,b1,b2) #Terminal
    AddSegment11(T2,T3,D5,D6,l3,b1,b2)          #HE1
    AddSegment11(T3,T3,D6,D7,space2,b1,b2)      #Space2
    AddSegment11(T3,T4,D7,D8,l4,b1,b2)          #HE2

    
# Beschleunigung Cologne FN
def AddFNAcc(v_gesamt, v_vorbeschl, q):
    global lastFunction
    lastFunction = "AddFNAcc"
    
    R_gesamt = 107510.0e6 #Gesamtwiderstand LE
    I_gesamt = float(-v_gesamt/R_gesamt) #Gesamtstrom durch ACC linke Seite
    
    L=float((2.532+.223)/100.0) #Laenge eines Segmentes
    D = 5.5*2./100.0    #Durchmesser des Elektrodenringes, update: neue Roehren
    q=float(q)  
    
    # Widerstaende in MOhm fuer LE1-2 (erstes Element zw. 1 und 2 usw) Stand 25.10.2011
    R_Segment1 = 0,0,0,580,570,0,600,610,600,600,0,600,590,550,600,600,590,590,580,590,600,600,600,590,600,600,600,590,600,600,610,590,600,0,610,610,600,600,600,600,600,600,590,610,600,0,600,600,600,610,600,600,600,600,600,610,600,600,610,600,610,600,610,600,600,610,600,600,610,590,610,610,560,620,610,620,610,610,610,630,600,620,620,610,640,610,620,610,630,610,600,630,630,620,610,600,0
    R_Segment2 = 570,550,580,580,580,550,580,580,570,580,570,580,570,560,570,580,570,580,580,580,610,570,580,580,570,560,570,580,570,570,610,580,560,570,590,560,610,0,540,560,580,570,550,570,560,600,570,570,580,600,570,580,590,570,580,570,570,570,570,560,580,540,520,550,550,550,540,530,540,550,550,540,550,550,560,550,570,570,560,540,560,580,550,550,540,570,540,560,560,550,560,580,570,580,550,580,0
    
    # Widerstaende in MOhm fuer HE1-2 (erstes Element zw. 1 und 2 usw) Stand 25.10.2011
    R_Segment3 = 590,580,580,580,580,590,560,590,580,580,590,560,570,570,570,590,570,580,590,600,580,560,580,580,570,600,570,580,570,570,570,580,570,580,600,580,570,580,580,580,570,590,580,570,590,570,570,590,580,590,580,580,590,600,580,580,570,570,570,580,580,580,600,590,590,570,580,560,570,570,600,580,580,570,580,570,570,570,580,580,600,620,570,590,580,600,580,570,590,600,580,580,590,570,580,590,0
    R_Segment4 = 580,570,570,570,580,580,580,570,570,580,570,580,580,570,570,600,620,610,600,590,580,590,580,580,580,590,580,580,580,580,580,580,580,570,580,580,580,580,580,580,590,570,580,580,580,570,590,580,580,580,580,580,580,580,590,580,580,580,570,580,590,580,580,580,580,580,570,580,580,570,580,580,580,580,590,580,570,580,580,580,590,590,580,580,580,580,590,580,590,580,590,580,580,580,590,600,0
    # INFO: Letzter Widerstand (0) dient zur Berechnung von E3


# LE1 : 4..5

    E1=0.
    E2=R_Segment1[3]*1.0e6*I_gesamt/L
    E3=R_Segment1[4]*1.0e6*I_gesamt/L
    T0 = float(v_vorbeschl)
    T1 = float(T0+float(-R_Segment1[4-1])*1.0e6*I_gesamt)
    #N = math.sqrt(T1/T0)#
    #AddSegment(N,L)#
    AddSegment2(T0,T1,8.8/100.0,L,E1,E2,E3,-1.0)
    
# LE1 : 5..97

    for i in range(5,9):
        E1=E2
        E2=E3
        E3=R_Segment1[i]*1.0e6*I_gesamt/L
        T0 = T1
        T1 = float(T0+float(-R_Segment1[i-1])*1.0e6*I_gesamt)
        #N = math.sqrt(T1/T0)#
        #AddSegment(N,L)#
        AddSegment2(T0,T1,7.2/100.0,L,E1,E2,E3,-1.0)#


    for i in range(9,12):
        E1=E2
        E2=E3
        E3=R_Segment1[i]*1.0e6*I_gesamt/L
        T0 = T1
        T1 = float(T0+float(-R_Segment1[i-1])*1.0e6*I_gesamt)
        #N = math.sqrt(T1/T0)#
        #AddSegment(N,L)#
        AddSegment2(T0,T1,6.4/100.0,L,E1,E2,E3,-1.0)#

    for i in range(12,15):
        E1=E2
        E2=E3
        E3=R_Segment1[i]*1.0e6*I_gesamt/L
        T0 = T1
        T1 = float(T0+float(-R_Segment1[i-1])*1.0e6*I_gesamt)
        #N = math.sqrt(T1/T0)#
        #AddSegment(N,L)#
        AddSegment2(T0,T1,6.0/100.0,L,E1,E2,E3,-1.0)#
        
    for i in range(15,18):
        E1=E2
        E2=E3
        E3=R_Segment1[i]*1.0e6*I_gesamt/L
        T0 = T1
        T1 = float(T0+float(-R_Segment1[i-1])*1.0e6*I_gesamt)
        #N = math.sqrt(T1/T0)#
        #AddSegment(N,L)#
        AddSegment2(T0,T1,5.2/100.0,L,E1,E2,E3,-1.0)#
        
    for i in range(18,21):
        E1=E2
        E2=E3
        E3=R_Segment1[i]*1.0e6*I_gesamt/L
        T0 = T1
        T1 = float(T0+float(-R_Segment1[i-1])*1.0e6*I_gesamt)
        #N = math.sqrt(T1/T0)#
        #AddSegment(N,L)#
        AddSegment2(T0,T1,4.8/100.0,L,E1,E2,E3,-1.0)#

    for i in range(21,31):
        E1=E2
        E2=E3
        E3=R_Segment1[i]*1.0e6*I_gesamt/L
        T0 = T1
        T1 = float(T0+float(-R_Segment1[i-1])*1.0e6*I_gesamt)
        #N = math.sqrt(T1/T0)#
        #AddSegment(N,L)#
        AddSegment2(T0,T1,4.0/100.0,L,E1,E2,E3,-1.0)#       

    for i in range(31,76):
        E1=E2
        E2=E3
        E3=R_Segment1[i]*1.0e6*I_gesamt/L
        T0 = T1
        T1 = float(T0+float(-R_Segment1[i-1])*1.0e6*I_gesamt)
        #N = math.sqrt(T1/T0)#
        #AddSegment(N,L)#
        AddSegment2(T0,T1,3.6/100.0,L,E1,E2,E3,-1.0)#
        
    for i in range(76,97):
        E1=E2
        E2=E3
        E3=R_Segment1[i]*1.0e6*I_gesamt/L
        T0 = T1
        T1 = float(T0+float(-R_Segment1[i-1])*1.0e6*I_gesamt)
        #N = math.sqrt(T1/T0)#
        #AddSegment(N,L)#
        AddSegment2(T0,T1,2.8/100.0,L,E1,E2,E3,-1.0)#
        
    E1=E2
    E2=E3
    E3=0.0
    T0 = T1
    T1 = T1
    
# Space = 22.225cm = Drift

    E1=E2
    E2=E3
    E3=R_Segment2[0]*1.0e6*I_gesamt/L
    T0 = T1
    T1 = T1
    #N = math.sqrt(T1/T0)#
    #AddSegment(N,22.225/100.0)#
    AddSegment2(T0,T1,D,22.225/100.0,E1,E2,E3,-1.0)
    
# LE2 : 1..97

    for i in range(1,97):
        E1=E2
        E2=E3
        E3=R_Segment2[i]*1.0e6*I_gesamt/L
        T0 = T1
        T1 = float(T0+float(-R_Segment2[i-1])*1.0e6*I_gesamt)
        #N = math.sqrt(T1/T0)#
        #AddSegment(N,L)#
        AddSegment2(T0,T1,D,L,E1,E2,E3,-1.0)#
    
    E1=E2
    E2=E3
    E3=0.0
    T0 = T1
    T1 = T1
    
# HE1 : R und I

    R_gesamt = 111510.0e6
    I_gesamt = float(v_gesamt/R_gesamt) #Gesamtstrom durch ACC rechte Seite
            
# Terminal / Stripping = Drift

    E1=E2
    E2=E3
    E3=R_Segment3[0]*1.0e6*I_gesamt/L
    T0 = T1
    T1 = T1
    #N = math.sqrt(T1/T0)#
    #AddSegment(N,121.92/100.0)#
    AddSegment2(T0,T1,D,121.92/100.0,E1,E2,E3,-1.0)#

# HE1 : 1..97

    for i in range(1,97):
        E1=E2
        E2=E3
        E3=R_Segment3[i]*1.0e6*I_gesamt/L
        T0 = T1
        T1 = float(T0+q*float(R_Segment3[i-1])*1.0e6*I_gesamt)
        #N = math.sqrt(T1/T0)#
        #AddSegment(N,L)#
        AddSegment2(T0,T1,D,L,E1,E2,E3,q)#

    E1=E2
    E2=E3
    E3=0.0
    T0 = T1
    T1 = T1

# Space = 22.225cm = Drift

    E1=E2
    E2=E3
    E3=R_Segment4[0]*1.0e6*I_gesamt/L
    T0 = T1
    T1 = T1
    #N = math.sqrt(T1/T0)#
    #AddSegment(N,22.225/100.0)#
    AddSegment2(T0,T1,D,22.225/100.0,E1,E2,E3,q)#
    
# HE2 : 1..97

    for i in range(1,97):
        E1=E2
        E2=E3
        E3=R_Segment4[i]*1.0e6*I_gesamt/L
        T0 = T1
        T1 = float(T0+q*float(R_Segment4[i-1])*1.0e6*I_gesamt)
        #N = math.sqrt(T1/T0)#
        #AddSegment(N,L)#
        AddSegment2(T0,T1,D,L,E1,E2,E3,q)#

    E1=E2
    E2=E3
    E3=0.0  
    T0 = T1
    T1 = T1

# Letzte Elektrode..Flansch aussen = Drift
    
    E1=E2
    E2=E3
    E3=0.0
    T0=T1
    T1=T1
    #AddSegment(1,117.738/100.0)#
    AddSegment2(T0,T1,D,117.738/100.0,E1,E2,E3,q)#

    
def AddFNAcc1(v_gesamt, v_vorbeschl, q):
    global lastFunction
    lastFunction = "AddFNAcc1"
    
    R_gesamt = 107510.0e6 #Gesamtwiderstand LE
    I_gesamt = float(-v_gesamt/R_gesamt) #Gesamtstrom durch ACC linke Seite
    
    L=float((2.532+.223)/100.0) #Laenge eines Segmentes
    D = 5.5*2/100.        #Durchmesser des Elektrodenringes, update: neue Roehren
    q=float(q)  
    
    # Widerstaende in MOhm fuer LE1-2 (erstes Element zw. 1 und 2 usw) Stand 25.10.2011
    R_Segment1 = 0,0,0,580,570,0,600,610,600,600,0,600,590,550,600,600,590,590,580,590,600,600,600,590,600,600,600,590,600,600,610,590,600,0,610,610,600,600,600,600,600,600,590,610,600,0,600,600,600,610,600,600,600,600,600,610,600,600,610,600,610,600,610,600,600,610,600,600,610,590,610,610,560,620,610,620,610,610,610,630,600,620,620,610,640,610,620,610,630,610,600,630,630,620,610,600,0
    R_Segment2 = 570,550,580,580,580,550,580,580,570,580,570,580,570,560,570,580,570,580,580,580,610,570,580,580,570,560,570,580,570,570,610,580,560,570,590,560,610,0,540,560,580,570,550,570,560,600,570,570,580,600,570,580,590,570,580,570,570,570,570,560,580,540,520,550,550,550,540,530,540,550,550,540,550,550,560,550,570,570,560,540,560,580,550,550,540,570,540,560,560,550,560,580,570,580,550,580,0
    
    # Widerstaende in MOhm fuer HE1-2 (erstes Element zw. 1 und 2 usw) Stand 25.10.2011
    R_Segment3 = 590,580,580,580,580,590,560,590,580,580,590,560,570,570,570,590,570,580,590,600,580,560,580,580,570,600,570,580,570,570,570,580,570,580,600,580,570,580,580,580,570,590,580,570,590,570,570,590,580,590,580,580,590,600,580,580,570,570,570,580,580,580,600,590,590,570,580,560,570,570,600,580,580,570,580,570,570,570,580,580,600,620,570,590,580,600,580,570,590,600,580,580,590,570,580,590,0
    R_Segment4 = 580,570,570,570,580,580,580,570,570,580,570,580,580,570,570,600,620,610,600,590,580,590,580,580,580,590,580,580,580,580,580,580,580,570,580,580,580,580,580,580,590,570,580,580,580,570,590,580,580,580,580,580,580,580,590,580,580,580,570,580,590,580,580,580,580,580,570,580,580,570,580,580,580,580,590,580,570,580,580,580,590,590,580,580,580,580,590,580,590,580,590,580,580,580,590,600,0
    # INFO: Letzter Widerstand (0) dient zur Berechnung von E3


# LE1 : 4..5

    E1=0.
    E2=R_Segment1[3]*1.0e6*I_gesamt/L
    E3=R_Segment1[4]*1.0e6*I_gesamt/L
    T0 = float(v_vorbeschl)
    T1 = float(T0+float(-R_Segment1[4-1])*1.0e6*I_gesamt)
    #N = math.sqrt(T1/T0)#
    #AddSegment(N,L)#
    AddSegment2(T0,T1,8.8/100.0,L,E1,E2,E3,-1.0)# //Durch die Grid Lens ist das erste Fringing Field 0
    
# LE1 : 5..97

    for i in range(5,9):
        E1=E2
        E2=E3
        E3=R_Segment1[i]*1.0e6*I_gesamt/L
        T0 = T1
        T1 = float(T0+float(-R_Segment1[i-1])*1.0e6*I_gesamt)
        #N = math.sqrt(T1/T0)#
        #AddSegment(N,L)#
        AddSegment2(T0,T1,7.2/100.0,L,E1,E2,E3,-1.0)#


    for i in range(9,12):
        E1=E2
        E2=E3
        E3=R_Segment1[i]*1.0e6*I_gesamt/L
        T0 = T1
        T1 = float(T0+float(-R_Segment1[i-1])*1.0e6*I_gesamt)
        #N = math.sqrt(T1/T0)#
        #AddSegment(N,L)#
        AddSegment2(T0,T1,6.4/100.0,L,E1,E2,E3,-1.0)#

    for i in range(12,15):
        E1=E2
        E2=E3
        E3=R_Segment1[i]*1.0e6*I_gesamt/L
        T0 = T1
        T1 = float(T0+float(-R_Segment1[i-1])*1.0e6*I_gesamt)
        #N = math.sqrt(T1/T0)#
        #AddSegment(N,L)#
        AddSegment2(T0,T1,6.0/100.0,L,E1,E2,E3,-1.0)#
        
    for i in range(15,18):
        E1=E2
        E2=E3
        E3=R_Segment1[i]*1.0e6*I_gesamt/L
        T0 = T1
        T1 = float(T0+float(-R_Segment1[i-1])*1.0e6*I_gesamt)
        #N = math.sqrt(T1/T0)#
        #AddSegment(N,L)#
        AddSegment2(T0,T1,5.2/100.0,L,E1,E2,E3,-1.0)#
        
    for i in range(18,21):
        E1=E2
        E2=E3
        E3=R_Segment1[i]*1.0e6*I_gesamt/L
        T0 = T1
        T1 = float(T0+float(-R_Segment1[i-1])*1.0e6*I_gesamt)
        #N = math.sqrt(T1/T0)#
        #AddSegment(N,L)#
        AddSegment2(T0,T1,4.8/100.0,L,E1,E2,E3,-1.0)#

    for i in range(21,31):
        E1=E2
        E2=E3
        E3=R_Segment1[i]*1.0e6*I_gesamt/L
        T0 = T1
        T1 = float(T0+float(-R_Segment1[i-1])*1.0e6*I_gesamt)
        #N = math.sqrt(T1/T0)#
        #AddSegment(N,L)#
        AddSegment2(T0,T1,4.0/100.0,L,E1,E2,E3,-1.0)#       

    for i in range(31,76):
        E1=E2
        E2=E3
        E3=R_Segment1[i]*1.0e6*I_gesamt/L
        T0 = T1
        T1 = float(T0+float(-R_Segment1[i-1])*1.0e6*I_gesamt)
        #N = math.sqrt(T1/T0)#
        #AddSegment(N,L)#
        AddSegment2(T0,T1,3.6/100.0,L,E1,E2,E3,-1.0)#
        
    for i in range(76,97):
        E1=E2
        E2=E3
        E3=R_Segment1[i]*1.0e6*I_gesamt/L
        T0 = T1
        T1 = float(T0+float(-R_Segment1[i-1])*1.0e6*I_gesamt)
        #N = math.sqrt(T1/T0)#
        #AddSegment(N,L)#
        AddSegment2(T0,T1,2.8/100.0,L,E1,E2,E3,-1.0)#
        
    E1=E2
    E2=E3
    E3=0.0
    T0 = T1
    T1 = T1
    
# Space = 22.225cm = Drift

    E1=E2
    E2=E3
    E3=R_Segment2[0]*1.0e6*I_gesamt/L
    T0 = T1
    T1 = T1
    #N = math.sqrt(T1/T0)#
    #AddSegment(N,22.225/100.0)#
    AddSegment2(T0,T1,D,22.225/100.0,E1,E2,E3,-1.0)
    
# LE2 : 1..97

    for i in range(1,97):
        E1=E2
        E2=E3
        E3=R_Segment2[i]*1.0e6*I_gesamt/L
        T0 = T1
        T1 = float(T0+float(-R_Segment2[i-1])*1.0e6*I_gesamt)
        #N = math.sqrt(T1/T0)#
        #AddSegment(N,L)#
        AddSegment2(T0,T1,D,L,E1,E2,E3,-1.0)#
    
    E1=E2
    E2=E3
    E3=0.0
    T0 = T1
    T1 = T1
    
# HE1 : R und I

    R_gesamt = 111510.0e6
    I_gesamt = float(v_gesamt/R_gesamt) #Gesamtstrom durch ACC rechte Seite
            
# Terminal / Stripping = Drift

    E1=E2
    E2=E3
    E3=R_Segment3[0]*1.0e6*I_gesamt/L
    T0 = T1
    T1 = T1
    N = math.sqrt(T1/T0)#
    AddSegment(N,121.92/100.0)#
    #AddSegment2(T0,T1,D,121.92/100.0,E1,E2,E3,-1.0)#

# HE1 : 1..97

    for i in range(1,97):
        E1=E2
        E2=E3
        E3=R_Segment3[i]*1.0e6*I_gesamt/L
        T0 = T1
        T1 = float(T0+q*float(R_Segment3[i-1])*1.0e6*I_gesamt)
        N = math.sqrt(T1/T0)#
        AddSegment(N,L)#
        #AddSegment2(T0,T1,D,L,E1,E2,E3,q)#

    E1=E2
    E2=E3
    E3=0.0
    T0 = T1
    T1 = T1

# Space = 22.225cm = Drift

    E1=E2
    E2=E3
    E3=R_Segment4[0]*1.0e6*I_gesamt/L
    T0 = T1
    T1 = T1
    N = math.sqrt(T1/T0)#
    AddSegment(N,22.225/100.0)#
    #AddSegment2(T0,T1,D,22.225/100.0,E1,E2,E3,q)#
    
# HE2 : 1..97

    for i in range(1,97):
        E1=E2
        E2=E3
        E3=R_Segment4[i]*1.0e6*I_gesamt/L
        T0 = T1
        T1 = float(T0+q*float(R_Segment4[i-1])*1.0e6*I_gesamt)
        N = math.sqrt(T1/T0)#
        AddSegment(N,L)#
        #AddSegment2(T0,T1,D,L,E1,E2,E3,q)#

    E1=E2
    E2=E3
    E3=0.0  
    T0 = T1
    T1 = T1

# Letzte Elektrode..Flansch aussen = Drift
    
    E1=E2
    E2=E3
    E3=0.0
    T0=T1
    T1=T1
    #AddSegment(1,117.738/100.0)#
    AddSegment2(T0,T1,D,117.738/100.0,E1,E2,E3,q)#
    
    

def AddFNAcc2(v_gesamt, v_vorbeschl, q):
    global lastFunction
    lastFunction = "AddFNAcc2"
    
    R_gesamt = 107510.0e6 #Gesamtwiderstand LE
    I_gesamt = float(-v_gesamt/R_gesamt) #Gesamtstrom durch ACC linke Seite
    
    L=float((2.532+.223)/100.0) #Laenge eines Segmentes
    D = 6.35*2/100        #Durchmesser des Elektrodenringes, update: neue Roehren
    q=float(q)  
    
    # Widerstaende in MOhm fuer LE1-2 (erstes Element zw. 1 und 2 usw) Stand 25.10.2011
    R_Segment1 = 0,0,0,580,570,0,600,610,600,600,0,600,590,550,600,600,590,590,580,590,600,600,600,590,600,600,600,590,600,600,610,590,600,0,610,610,600,600,600,600,600,600,590,610,600,0,600,600,600,610,600,600,600,600,600,610,600,600,610,600,610,600,610,600,600,610,600,600,610,590,610,610,560,620,610,620,610,610,610,630,600,620,620,610,640,610,620,610,630,610,600,630,630,620,610,600,0
    R_Segment2 = 570,550,580,580,580,550,580,580,570,580,570,580,570,560,570,580,570,580,580,580,610,570,580,580,570,560,570,580,570,570,610,580,560,570,590,560,610,0,540,560,580,570,550,570,560,600,570,570,580,600,570,580,590,570,580,570,570,570,570,560,580,540,520,550,550,550,540,530,540,550,550,540,550,550,560,550,570,570,560,540,560,580,550,550,540,570,540,560,560,550,560,580,570,580,550,580,0
    
    # Widerstaende in MOhm fuer HE1-2 (erstes Element zw. 1 und 2 usw) Stand 25.10.2011
    R_Segment3 = 590,580,580,580,580,590,560,590,580,580,590,560,570,570,570,590,570,580,590,600,580,560,580,580,570,600,570,580,570,570,570,580,570,580,600,580,570,580,580,580,570,590,580,570,590,570,570,590,580,590,580,580,590,600,580,580,570,570,570,580,580,580,600,590,590,570,580,560,570,570,600,580,580,570,580,570,570,570,580,580,600,620,570,590,580,600,580,570,590,600,580,580,590,570,580,590,0
    R_Segment4 = 580,570,570,570,580,580,580,570,570,580,570,580,580,570,570,600,620,610,600,590,580,590,580,580,580,590,580,580,580,580,580,580,580,570,580,580,580,580,580,580,590,570,580,580,580,570,590,580,580,580,580,580,580,580,590,580,580,580,570,580,590,580,580,580,580,580,570,580,580,570,580,580,580,580,590,580,570,580,580,580,590,590,580,580,580,580,590,580,590,580,590,580,580,580,590,600,0
    # INFO: Letzter Widerstand (0) dient zur Berechnung von E3


# LE1 : 4..5

    E1=0.
    E2=R_Segment1[3]*1.0e6*I_gesamt/L
    E3=R_Segment1[4]*1.0e6*I_gesamt/L
    T0 = float(v_vorbeschl)
    T1 = float(T0+float(-R_Segment1[4-1])*1.0e6*I_gesamt)
    N = math.sqrt(T1/T0)#
    AddSegment(N,L)#
    #AddSegment2(T0,T1,8.8/100.0,L,E1,E2,E3,-1.0)# //Durch die Grid Lens ist das erste Fringing Field 0
    
# LE1 : 5..97

    for i in range(5,9):
        E1=E2
        E2=E3
        E3=R_Segment1[i]*1.0e6*I_gesamt/L
        T0 = T1
        T1 = float(T0+float(-R_Segment1[i-1])*1.0e6*I_gesamt)
        N = math.sqrt(T1/T0)#
        AddSegment(N,L)#
        #AddSegment2(T0,T1,7.2/100.0,L,E1,E2,E3,-1.0)#


    for i in range(9,12):
        E1=E2
        E2=E3
        E3=R_Segment1[i]*1.0e6*I_gesamt/L
        T0 = T1
        T1 = float(T0+float(-R_Segment1[i-1])*1.0e6*I_gesamt)
        N = math.sqrt(T1/T0)#
        AddSegment(N,L)#
        #AddSegment2(T0,T1,6.4/100.0,L,E1,E2,E3,-1.0)#

    for i in range(12,15):
        E1=E2
        E2=E3
        E3=R_Segment1[i]*1.0e6*I_gesamt/L
        T0 = T1
        T1 = float(T0+float(-R_Segment1[i-1])*1.0e6*I_gesamt)
        N = math.sqrt(T1/T0)#
        AddSegment(N,L)#
        #AddSegment2(T0,T1,6.0/100.0,L,E1,E2,E3,-1.0)#
        
    for i in range(15,18):
        E1=E2
        E2=E3
        E3=R_Segment1[i]*1.0e6*I_gesamt/L
        T0 = T1
        T1 = float(T0+float(-R_Segment1[i-1])*1.0e6*I_gesamt)
        N = math.sqrt(T1/T0)#
        AddSegment(N,L)#
        #AddSegment2(T0,T1,5.2/100.0,L,E1,E2,E3,-1.0)#
        
    for i in range(18,21):
        E1=E2
        E2=E3
        E3=R_Segment1[i]*1.0e6*I_gesamt/L
        T0 = T1
        T1 = float(T0+float(-R_Segment1[i-1])*1.0e6*I_gesamt)
        N = math.sqrt(T1/T0)#
        AddSegment(N,L)#
        #AddSegment2(T0,T1,4.8/100.0,L,E1,E2,E3,-1.0)#

    for i in range(21,31):
        E1=E2
        E2=E3
        E3=R_Segment1[i]*1.0e6*I_gesamt/L
        T0 = T1
        T1 = float(T0+float(-R_Segment1[i-1])*1.0e6*I_gesamt)
        N = math.sqrt(T1/T0)#
        AddSegment(N,L)#
        #AddSegment2(T0,T1,4.0/100.0,L,E1,E2,E3,-1.0)#      

    for i in range(31,76):
        E1=E2
        E2=E3
        E3=R_Segment1[i]*1.0e6*I_gesamt/L
        T0 = T1
        T1 = float(T0+float(-R_Segment1[i-1])*1.0e6*I_gesamt)
        N = math.sqrt(T1/T0)#
        AddSegment(N,L)#
        #AddSegment2(T0,T1,3.6/100.0,L,E1,E2,E3,-1.0)#
        
    for i in range(76,97):
        E1=E2
        E2=E3
        E3=R_Segment1[i]*1.0e6*I_gesamt/L
        T0 = T1
        T1 = float(T0+float(-R_Segment1[i-1])*1.0e6*I_gesamt)
        N = math.sqrt(T1/T0)#
        AddSegment(N,L)#
        #AddSegment2(T0,T1,2.8/100.0,L,E1,E2,E3,-1.0)#
        
    E1=E2
    E2=E3
    E3=0.0
    T0 = T1
    T1 = T1
    
# Space = 22.225cm = Drift

    E1=E2
    E2=E3
    E3=R_Segment2[0]*1.0e6*I_gesamt/L
    T0 = T1
    T1 = T1
    N = math.sqrt(T1/T0)#
    AddSegment(N,22.225/100.0)#
    #AddSegment2(T0,T1,D,22.225/100.0,E1,E2,E3,-1.0)
    
# LE2 : 1..97

    for i in range(1,97):
        E1=E2
        E2=E3
        E3=R_Segment2[i]*1.0e6*I_gesamt/L
        T0 = T1
        T1 = float(T0+float(-R_Segment2[i-1])*1.0e6*I_gesamt)
        N = math.sqrt(T1/T0)#
        AddSegment(N,L)#
        #AddSegment2(T0,T1,D,L,E1,E2,E3,-1.0)#
    
    E1=E2
    E2=E3
    E3=0.0
    T0 = T1
    T1 = T1
    
# HE1 : R und I

    R_gesamt = 111510.0e6
    I_gesamt = float(v_gesamt/R_gesamt) #Gesamtstrom durch ACC rechte Seite
            
# Terminal / Stripping = Drift

    E1=E2
    E2=E3
    E3=R_Segment3[0]*1.0e6*I_gesamt/L
    T0 = T1
    T1 = T1
    N = math.sqrt(T1/T0)#
    AddSegment(N,121.92/100.0)#
    #AddSegment2(T0,T1,D,121.92/100.0,E1,E2,E3,-1.0)#

# HE1 : 1..97

    for i in range(1,97):
        E1=E2
        E2=E3
        E3=R_Segment3[i]*1.0e6*I_gesamt/L
        T0 = T1
        T1 = float(T0+q*float(R_Segment3[i-1])*1.0e6*I_gesamt)
        N = math.sqrt(T1/T0)#
        AddSegment(N,L)#
        #AddSegment2(T0,T1,D,L,E1,E2,E3,q)#

    E1=E2
    E2=E3
    E3=0.0
    T0 = T1
    T1 = T1

# Space = 22.225cm = Drift

    E1=E2
    E2=E3
    E3=R_Segment4[0]*1.0e6*I_gesamt/L
    T0 = T1
    T1 = T1
    N = math.sqrt(T1/T0)#
    AddSegment(N,22.225/100.0)#
    #AddSegment2(T0,T1,D,22.225/100.0,E1,E2,E3,q)#
    
# HE2 : 1..97

    for i in range(1,97):
        E1=E2
        E2=E3
        E3=R_Segment4[i]*1.0e6*I_gesamt/L
        T0 = T1
        T1 = float(T0+q*float(R_Segment4[i-1])*1.0e6*I_gesamt)
        N = math.sqrt(T1/T0)#
        AddSegment(N,L)#
        #AddSegment2(T0,T1,D,L,E1,E2,E3,q)#

    E1=E2
    E2=E3
    E3=0.0  
    T0 = T1
    T1 = T1

# Letzte Elektrode..Flansch aussen = Drift
    
    E1=E2
    E2=E3
    E3=0.0
    T0=T1
    T1=T1
    #AddSegment(1,117.738/100.0)#
    AddSegment2(T0,T1,D,117.738/100.0,E1,E2,E3,q)#



# ACC Segment einfuegen nach Hinterberger
def AddSegment(N,L):
    global s, geo_x, geo_y
    if (s > -0.5):
        geo_s.InsertNextValue(s)
        geo_s.InsertNextValue(s+L)
        geo_y.InsertNextValue(55)   
        geo_y.InsertNextValue(55)   
        s = s+L
    
    try:
        m11=(3.0-N)/2.0
        m12=2.0*L/(1.0+N)
        m21=-3.0/8.0/N/N/L*(N*N-1.0)*(N-1.0)
        m22=(3.0*N-1.0)/(2.0*N*N)
        dkdk= 1.0/N/N
        dmdm= 1.0
        AddMatrix(1,[m11,m12,0.0,0.0,0.0,0.0,m21,m22,0.0,0.0,0.0,0.0,0.0,0.0,m11,m12,0.0,0.0,0.0,0.0,m21,m22,0.0,0.0,0.0,0.0,0.0,0.0,dkdk,0.0,0.0,0.0,0.0,0.0,0.0,dmdm],L)
    except:
        print "Division durch Null! (AddSegment)"

    
# Korrektur Apertur
def AddSegment1(T0, T1, q, D1, D2, L, b = .57):
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
    q = float(q)
    L = float(L)
    D1 = float(D1)
    D2 = float(D2)
    
    E = (T1-T0) / q / L
    N = math.sqrt(T1/T0)
    
    f1 = -q * (0. - E) / (2. * T0 * (1. + math.sqrt(1. - q * (0. - E) * b * D1 / T0)))  # 1/f1
    f2 = -q * (E - 0.) / (2. * T1 * (1. + math.sqrt(1. - q * (E - 0.) * b * D2 / T1)))  # 1/f2
    LL = 2. * L / (1. + N)
    m11= 1. - LL * f1
    m12= LL
    m21= -f2 - (1./N - f2*LL)*f1
    m22= 1/N - f2*LL

    dkdk= 1./N/N
    dmdm= 1.

    AddMatrix(1, [m11,m12,0.,0.,0.,0.,m21,m22,0.,0.,0.,0.,0.,0.,m11,m12,0.,0.,0.,0.,m21,m22,0.,0.,0.,0.,0.,0.,dkdk,0.,0.,0.,0.,0.,0.,dmdm], L)

# Korrektur Apertur
def AddSegment11(T0, T1, D1, D2, L, b1 = .57, b2 = .57):
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
    deltaT = T1-T0
    L = float(L)
    D1 = float(D1)
    D2 = float(D2)
    
    N = math.sqrt(T1/T0)
    
    f1 = + deltaT / (L * 2. * T0 * (1. + math.sqrt(1. + deltaT * b1 * D1 / T0 / L)))  # 1/f1
    f2 = - deltaT / (L * 2. * T1 * (1. + math.sqrt(1. - deltaT * b2 * D2 / T1 / L)))  # 1/f2
    LL = 2. * L / (1. + N)
    m11= 1. - LL * f1
    m12= LL
    m21= -f2 - (1./N - f2*LL)*f1
    m22= 1/N - f2*LL

    dkdk= 1./N/N
    dmdm= 1.

    AddMatrix(1, [m11,m12,0.,0.,0.,0.,m21,m22,0.,0.,0.,0.,0.,0.,m11,m12,0.,0.,0.,0.,m21,m22,0.,0.,0.,0.,0.,0.,dkdk,0.,0.,0.,0.,0.,0.,dmdm], L)

    
# ACC Segment einfuegen nach Elkind und T. Joy. Nur Eintritt beruecksichtigen.
def AddSegment2(T0,T1,D,L,E1,E2,E3,q):
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
        s = s+L
    
    N=float(math.sqrt(T1/T0))
    deltaT = float(T1-T0)
    D = float(D)*100.
    
    f1 = -q*(E1-E2)/2.0/T0/(1.0+math.sqrt(1.0-q*(E1-E2)/T0*0.57*D))  # 1/f1
    #f2 = -q*(E2-E3)/2.0/T1/(1.0+math.sqrt(1.0-q*(E2-E3)/T1*0.57*D))  # 1/f2
    LL = 2.0*L/(1.0+N)

    # Nur der Eintritt und die Beschleunigung werden beruecksichtigt:
    m11=1.0-LL*f1
    m12=LL
    m21=-f1/N
    m22=1/N
    
    dkdk= 1.0/N/N
    dmdm= 1.0
        
    AddMatrix(1,[m11,m12,0.0,0.0,0.0,0.0,m21,m22,0.0,0.0,0.0,0.0,0.0,0.0,m11,m12,0.0,0.0,0.0,0.0,m21,m22,0.0,0.0,0.0,0.0,0.0,0.0,dkdk,0.0,0.0,0.0,0.0,0.0,0.0,dmdm],L)

    
def AddSegment2VBFN(T0,T1,D,L,E1,E2,E3,q):
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
    
    N=float(math.sqrt(T1/T0))
    deltaT = float(T1-T0)
    D = float(D)*100.
    
    f1 = -q * (E1-E2) / (2. * T0 * (1. + math.sqrt(1. - q * (E1-E2) * 1.13 * D / T0)))  # 1/f1
    f2 = -q * (E2-E3) / (2. * T1 * (1. + math.sqrt(1. - q * (E2-E3) * 1.13 * D / T1)))  # 1/f2
    LL = 2. * L / (1. + N)
    

    m11= 1.0-LL*f1
    m12= LL
    m21= -f2-(1.0/N-f2*LL)*f1
    m22= 1/N-f2*LL
    
    
    dkdk= 1.0/N/N
    dmdm= 1.0
        
    AddMatrix(1,[m11,m12,0.0,0.0,0.0,0.0,m21,m22,0.0,0.0,0.0,0.0,0.0,0.0,m11,m12,0.0,0.0,0.0,0.0,m21,m22,0.0,0.0,0.0,0.0,0.0,0.0,dkdk,0.0,0.0,0.0,0.0,0.0,0.0,dmdm],L)  

    
def AddDrift(num,gamma2,length):
    global s, geo_s, geo_y
    global lastFunction
    lastFunction = "AddDrift"
    if (s > -0.5):
        geo_s.InsertNextValue(s)
        geo_s.InsertNextValue(s+length)
        geo_y.InsertNextValue(55)
        geo_y.InsertNextValue(55)
        s = s+length
    optic.AddDrift(ctypes.c_int(int(num)),ctypes.c_double(float(gamma2)), ctypes.c_double(float(length))) 
    

def AddBeamProfile():
    optic.AddBeamProfile()
    
def AddModifyEmittance(factor1, factor2):
    """
    Aendern der Emittanz um einen Faktor (z.B. beim Strippingprozess).
    """
    optic.AddModifyEmittance(ctypes.c_double(float(factor1)),ctypes.c_double(float(factor2)))
    
def AddSlit(x,dx,y,dy):
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
    optic.AddSlit(ctypes.c_double(float(x)),ctypes.c_double(float(dx)), ctypes.c_double(float(y)), ctypes.c_double(float(dy))) 
  
def AddThinLens(arg1,arg2,geo):
    global s, geo_s, geo_y
    global lastFunction
    lastFunction = "AddThinLens"
    if (s > -0.5):
        geo_s.InsertNextValue(s)
        geo_y.InsertNextValue(geo)
    optic.AddThinLens(ctypes.c_int(1),ctypes.c_double(float(arg1)), ctypes.c_double(float(arg2)),ctypes.c_double(0.))   

def AddAMSSO110EL(phi1,v_el):
    global lastFunction
    lastFunction = "AddAMSSO110EL"
    
    global s, geo_s, geo_y
    if (s > -0.5):
        geo_s.InsertNextValue(s)
        geo_y.InsertNextValue(50)
        
    x = (phi1-v_el)/float(phi1) # phi2/phi1
    #f =  .001/(0.022648131734415336-0.13770935811064514*x+0.44804113697363235*x**2-0.9575118663530009*x**3+1.3523795047227753*x**4-1.2071207469218252*x**5+0.6161254239513653*x**6-0.13696307518260592*x**7)#reihenentwicklung fuer x=.3 ... .8
    f =  .001/(0.024838874085440145-0.1721716204656418*x+0.6764569821500281*x**2-1.7844240368333488*x**3+3.118746647284872*x**4-3.434443885979744*x**5+2.152129217984039*x**6-0.5841534279352003*x**7)#reihenentwicklung fuer x=.2 ... .7
    optic.AddThinLens(ctypes.c_int(1),ctypes.c_double(float(f)), ctypes.c_double(float(f)),ctypes.c_double(0.)) 

def AddAMSBIEL(phi1,v_el):
    global lastFunction
    lastFunction = "AddAMSBIEL"
    
    global s, geo_s, geo_y
    if (s > -0.5):
        geo_s.InsertNextValue(s)
        geo_y.InsertNextValue(56)
        
    x = (phi1-v_el)/phi1 # phi2/phi1
    f =  .001/(0.01993277008153473-0.13793678662472572*x+0.5408115698616647*x**2-1.4238381952266255*x**3+2.4846514325496316*x**4-2.7328746285908734*x**5+1.7109251000894101*x**6-0.4640648127510327*x**7)#reihenentwicklung fuer x=.2 ... .7
    optic.AddThinLens(ctypes.c_int(1),ctypes.c_double(float(f)), ctypes.c_double(float(f)),ctypes.c_double(0.)) 
    
def AddFNEL(phi1,v_el):
    global lastFunction
    lastFunction = "AddFNEL"
    
    global s, geo_s, geo_y
    if (s > -0.5):
        geo_s.InsertNextValue(s)
        geo_y.InsertNextValue(63)
        
    x = (phi1-v_el)/phi1 # phi2/phi1
    f =  1./(11.542244476708195-61.57930337802449*x+183.54364546860486*x**2-385.71960791043387*x**3+562.1897241491755*x**4-529.9264768291281*x**5+288.3485380347237*x**6-68.5022805970122*x**7)#reihenentwicklung fuer x=.2 ... .7
    optic.AddThinLens(ctypes.c_int(1),ctypes.c_double(float(f)), ctypes.c_double(float(f)),ctypes.c_double(0.)) 
    
def AddESD(num,arg1,arg2,arg3,arg4,arg5,geo, korrektur = None, V_ist = 1., V_soll = 1.):
    global lastFunction
    lastFunction = "AddESD"
    
    if korrektur == None: korrektur = float(V_ist) / float(V_soll)
    global s, geo_s, geo_y
    if (s > -0.5):
        geo_s.InsertNextValue(s)
        geo_s.InsertNextValue(s+num*arg2*arg3)
        geo_y.InsertNextValue(geo)
        geo_y.InsertNextValue(geo)
        s = s+num*arg2*arg3
    optic.AddESD(ctypes.c_int(int(num)),ctypes.c_double(float(arg1)),ctypes.c_double(float(arg2)),ctypes.c_double(float(arg3)),ctypes.c_double(float(arg4)),ctypes.c_double(float(arg5)), ctypes.c_double(float(korrektur)))

def AddEdgeFocusing(r,beta,K,geo):
    global lastFunction
    lastFunction = "AddEdgeFocusing"
    
    g = 2.*geo/1.e3
    betaeff = beta - g/r*K*(1.+math.sin(beta)*math.sin(beta))/(math.cos(beta))
    optic.AddEdgeFocusing(ctypes.c_int(1),ctypes.c_double(float(r)), ctypes.c_double(float(beta)),ctypes.c_double(float(betaeff)))

def AddMSA(num,arg1,arg2,arg3,geo, korrektur = None, B_ist = 1., B_soll = 1.):
    global lastFunction
    lastFunction = "AddMSA"
    
    if korrektur == None: korrektur = float(B_ist) / float(B_soll)
    global s, geo_s, geo_y
    if (s > -0.5):
        geo_s.InsertNextValue(s)
        geo_s.InsertNextValue(s+num*arg2*arg3)
        geo_y.InsertNextValue(geo)
        geo_y.InsertNextValue(geo)
        s = s+num*arg2*arg3
    optic.AddHomDeflectingMagnet(ctypes.c_int(int(num)),ctypes.c_double(float(arg1)), ctypes.c_double(float(arg2)),ctypes.c_double(float(arg3)), ctypes.c_double(float(korrektur)))
    
def AddInhomMSA(num,rho,phi,n,geo):
    global s, geo_s, geo_y
    if (s > -0.5):
        geo_s.InsertNextValue(s)
        geo_s.InsertNextValue(s+num*rho*phi)
        geo_y.InsertNextValue(geo)
        geo_y.InsertNextValue(geo)
        s = s+num*rho*phi
    optic.AddInhomDeflectingMagnet(ctypes.c_int(int(num)),ctypes.c_double(float(rho)), ctypes.c_double(float(phi)),ctypes.c_double(float(n)))

def AddQuadrupolRadFoc(num,arg1,arg2,arg3,geo):
    global lastFunction
    lastFunction = "AddQuadrupolRadFoc"
    global s, geo_s, geo_y
    if (s > -0.5):
        geo_s.InsertNextValue(s)
        geo_s.InsertNextValue(s+num*arg3)
        geo_y.InsertNextValue(geo)
        geo_y.InsertNextValue(geo)
        s = s+num*arg3
    optic.AddQuadrupolRadFoc(ctypes.c_int(int(num)),ctypes.c_double(float(arg1)), ctypes.c_double(float(arg2)),ctypes.c_double(float(arg3)))

def AddQuadrupolAxFoc(num,arg1,arg2,arg3,geo):
    global lastFunction
    lastFunction = "AddQuadrupolAxFoc"
    
    global s, geo_s, geo_y
    if (s > -0.5):
        geo_s.InsertNextValue(s)
        geo_s.InsertNextValue(s+num*arg3)
        geo_y.InsertNextValue(geo)
        geo_y.InsertNextValue(geo)
        s = s+num*arg3
    optic.AddQuadrupolAxFoc(ctypes.c_int(int(num)),ctypes.c_double(float(arg1)), ctypes.c_double(float(arg2)),ctypes.c_double(float(arg3)))

def AddAMSQPT_YXY(gamma2,prozent,astigm,v_terminal,v_ext,q,geo):
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
        
    vx  =   30.e3*(prozent+astigm)/100. #empirisch bestimmt
    vy  =   30.e3*(prozent-astigm)/100. #empirisch bestimmt
    T       = v_ext+(q+1)*float(v_terminal)
    
    kx  = abs(vx)/float(geo*1.e-3)**2*q/T
#   ky  = abs(vy)/float(geo*1.e-3)**2*q/T
    k       = (vx+vy)/2./float(geo*1.e-3)**2*q/T
    
    #print "kx=\t",kx,"ky=\t",k
    
    optic.AddQuadrupolAxFoc(ctypes.c_int(10),ctypes.c_double(float(gamma2)), ctypes.c_double(float(k)),ctypes.c_double(float(.025)))
    optic.AddDrift(ctypes.c_int(1),ctypes.c_double(float(gamma2)), ctypes.c_double(.076)) 
    optic.AddQuadrupolRadFoc(ctypes.c_int(10),ctypes.c_double(float(gamma2)), ctypes.c_double(float(kx)),ctypes.c_double(float(.05)))
    optic.AddDrift(ctypes.c_int(1),ctypes.c_double(float(gamma2)), ctypes.c_double(.076)) 
    optic.AddQuadrupolAxFoc(ctypes.c_int(10),ctypes.c_double(float(gamma2)), ctypes.c_double(float(k)),ctypes.c_double(float(.025)))
    
def AddAMSQPT_XYX(gamma2,prozent,astigm,v_terminal,v_ext,q,geo):
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
        
    vx  =   30.e3*(prozent+astigm)/100. #empirisch bestimmt
    vy  =   30.e3*(prozent-astigm)/100. #empirisch bestimmt
    T       = v_ext+(q+1)*float(v_terminal)
    
    kx  = abs(vx)/float(geo*1.e-3)**2*q/T
#   ky  = abs(vy)/float(geo*1.e-3)**2*q/T
    k       = (vx+vy)/2./float(geo*1.e-3)**2*q/T
    
    #print "kx=\t",kx,"ky=\t",k
    
    optic.AddQuadrupolRadFoc(ctypes.c_int(10),ctypes.c_double(float(gamma2)), ctypes.c_double(float(k)),ctypes.c_double(float(.025)))
    optic.AddDrift(ctypes.c_int(1),ctypes.c_double(float(gamma2)), ctypes.c_double(.076)) 
    optic.AddQuadrupolAxFoc(ctypes.c_int(10),ctypes.c_double(float(gamma2)), ctypes.c_double(float(kx)),ctypes.c_double(float(.05)))
    optic.AddDrift(ctypes.c_int(1),ctypes.c_double(float(gamma2)), ctypes.c_double(.076)) 
    optic.AddQuadrupolRadFoc(ctypes.c_int(10),ctypes.c_double(float(gamma2)), ctypes.c_double(float(k)),ctypes.c_double(float(.025)))

def AddAMSQPT(gamma2,prozent,astigm,v_terminal,v_ext,q,geo):
    global lastFunction
    lastFunction = "AddAMSQPT"
    
    AddAMSQPT_XYX(gamma2,prozent,astigm,v_terminal,v_ext,q,geo)

    
#zum zeichnen eigener geometrie
def AddGeo(s1,x1,s2,x2):
    global lastFunction
    lastFunction = "AddGeo"
    
    global s, geo_s, geo_y
    if (s > -0.5):
        geo_s.InsertNextValue(s+s1)
        geo_s.InsertNextValue(s+s2)
        geo_y.InsertNextValue(x1)
        geo_y.InsertNextValue(x2)
        
    
def AddQPT_XYX(gamma2,vx,vy,v_terminal,v_ext,q,geo):
    global s, geo_s, geo_y
    vx  =   float(vx)
    vy  =   float(vy)
    T       = v_ext+(q+1)*float(v_terminal)
    
    kx  = abs(vx*1.e3)/float(geo*1.e-3)**2*q/T
    ky  = abs(vy*1.e3)/float(geo*1.e-3)**2*q/T
    print vx,vy
    print kx,ky
    
    optic.AddQuadrupolAxFoc(ctypes.c_int(10),ctypes.c_double(float(gamma2)), ctypes.c_double(float(ky)),ctypes.c_double(float(.025)))
    optic.AddDrift(ctypes.c_int(1),ctypes.c_double(float(gamma2)), ctypes.c_double(.076)) 
    optic.AddQuadrupolRadFoc(ctypes.c_int(10),ctypes.c_double(float(gamma2)), ctypes.c_double(float(kx)),ctypes.c_double(float(.05)))
    optic.AddDrift(ctypes.c_int(1),ctypes.c_double(float(gamma2)), ctypes.c_double(.076)) 
    optic.AddQuadrupolAxFoc(ctypes.c_int(10),ctypes.c_double(float(gamma2)), ctypes.c_double(float(ky)),ctypes.c_double(float(.025)))
    
###########################################################
def GetTrajectories():
    particlenum=int(optic.GetParticleNum())
    particlesize=int(optic.GetParticleSize())
    trajsize=int(optic.GetTrajectoriesSize())
    trajectories=(ctypes.c_double*trajsize)() # erstelle ein Array vom C-Type double*
    optic.GetTrajectories(trajectories)
    traj=array.array('d',trajectories)
    return (traj,particlenum,particlesize)

###########################################################
def GetTrajectory(iparticle,iproperty):
    traj=(ctypes.c_double*optic.GetTrajectorySize())()
    optic.GetTrajectory(ctypes.c_int(iparticle),ctypes.c_int(iproperty),traj)
    traj=array.array('d',traj)
    return traj

###########################################################
def PrintTrajectories(traj):
    (trajectories,particlenum,particlesize)=GetTrajectories()
    for i in range(0,particlenum*particlesize,particlesize):
        for j in range(0,len(trajectories),particlenum*particlesize):
            for k in range(0,particlesize,1):
                print "%f" % (trajectories[j+i+k]),
            print "\n",
        print "\n",

###########################################################

#def DoTheLoop(text):
#   ser = serial.Serial('/dev/ttyACM0', 9600)
#   INPUT1 = float(ser.readline())
#   LOOP = int(ser.readline())
#   while LOOP == 1:    
#      INPUT1 = float(ser.readline())
#      LOOP = int(ser.readline())
#      exec(text)
#      print(INPUT1)
#      print(LOOP)
#      time.sleep(1)

def ExecText(text, INPUT, SOURCE1):
    global textArray, lastFunction
    textArray = []
    lastFunction = "unknown"
    
    global SOURCE
    SOURCE = SOURCE1
    exec(text)
    

###########################################################
def Limioptic(filename):
    optic.Clear()
    execfile(filename)
    optic.CalculateTrajectories()
    traj=GetTrajectories()
    PrintTrajectories(traj)

    
    
def Name(text):
    textArray.append([s, text])
    
    
    
######################################################################################################################

"""
Momentan nicht verwendete Funktionen!
"""

def AddAMSAcc1(v_qsnout, v_gesamt, v_vorbeschl, q):
    R_a1 = 400.0e6+10.0e3 #Widerstand 1-4 Zweig1
    R_a2 = 400.0e6+2.0e9 #Widerstand 1-4 Zweig2
    R_a = float(1/(1/R_a1+1/R_a2)) #Widerstand 1-4
    R_b = 18.0*40.0e6+62.0*200.0e6+81.0*200.0e6 #Widerstand 5-84
    R_gesamt = R_b+R_a #Gesamtwiderstand HE1
    I_gesamt = float(v_gesamt/R_gesamt) #Gesamtstrom durch ACC linke Seite
    L=float(25.4e-3) #Laenge eines Segmentes
    V_a = float(R_a * I_gesamt) #Spannungsabfall 1..4
    I_a2 = float(V_a / R_a2)
    a = s
    
    # Q-Snout : Laenge 300.4mm
    T0 = float(v_vorbeschl+v_qsnout)
    T1 = T0
    AddSegment1(T0,T1,q,44.e-3,95.e-3,300.4e-3)

# HE1 : Ende Q-Snout..Anfang Terminal
# 3..4
    T0 = T1
    T1 = float(v_vorbeschl+V_a-1.e9*I_a2)
    AddSegment1(T0,T1,q,95.e-3,0.,L)
#   print "Spannung04 = %1.0f" % (T1)
    
# 4..5
    T0 = T1
    T1 = float(v_vorbeschl+V_a)
    N = math.sqrt(T1/T0)
    AddSegment(N,L)
#   print "Spannung05 = %1.0f" % (T1)
    
    deltaV = 40.0e6*I_gesamt #5..23
    for i in range(0,18):
        T0 = T1
        T1 = float(T0+deltaV)
        N = math.sqrt(T1/T0)
        AddSegment(N,L)
        #print "Spannung0%1.0f = %1.0f" % (i+6, T1)
    
    deltaV = 200.0e6*I_gesamt #23..85
    for i in range(0,62):
        T0 = T1
        T1 = float(T0+deltaV)
        N = math.sqrt(T1/T0)
        AddSegment(N,L)
#   print "Spannung85 = %1.0f" % (T1)
    
# Space = 76.23mm
    T0 = T1
    T1 = T1
    N = math.sqrt(T1/T0)
    AddSegment(N,76.23e-3)
#   print "Spannung85 = %1.0f" % (T1)
    
# Space..Terminal
    deltaV = 200.0e6*I_gesamt #85..166
    for i in range(0,81):
        T0 = T1
        T1 = float(T0+deltaV)
        N = math.sqrt(T1/T0)
        AddSegment(N,L)
        
#   print "Terminalspannung = %1.0f" % (T1)
        
# Terminal / Stripping
    T0 = T1
    T1 = T1
    N = math.sqrt(T1/T0)
    AddSegment(N,1147.0e-3) #Terminallaenge in mm

# HE2 : Ende Terminal..letzte Elektrode
    R_gesamt = float(9.0*120.0e6+156.0*300.0e6)
    I_gesamt = float(v_gesamt/R_gesamt) #Gesamtstrom durch ACC rechte Seite

# Terminal..Space
    deltaV = 120.0e6*I_gesamt #1..9
    for i in range(0,9):
        T0 = T1
        T1 = float(T0+q*deltaV)
        N = math.sqrt(T1/T0)
        AddSegment(N,L)

    deltaV = 300.0e6*I_gesamt #9..81
    for i in range(0,72):
        T0 = T1
        T1 = float(T0+q*deltaV)
        N = math.sqrt(T1/T0)
        AddSegment(N,L)

# Space = 76.05 mm 
    T0 = T1
    T1 = T1
    N = math.sqrt(T1/T0)
    AddSegment(N,76.05e-3)

# Space..letzte Elekrode
    deltaV = 300.0e6*I_gesamt #81..165
    for i in range(0,84):
        T0 = T1
        T1 = float(T0+q*deltaV)
        N = math.sqrt(T1/T0)
        AddSegment(N,L)
        
# Letzte Elektrode..Flansch aussen = 487.8mm
    AddSegment(1,447.81e-3)
#   print "Tanklaenge:", s-a



    
#Alex:  Beschleunigung Cologne AMS nach Elkind und T. Joy
def AddAMSAcc2(v_qsnout, v_gesamt, v_vorbeschl, q):
    R_a1 = 400.0e6+10.0e3 #Widerstand 1-4 Zweig1
    R_a2 = 400.0e6+200.0e9 #Widerstand 1-4 Zweig2
    R_a = float(1/(1/R_a1+1/R_a2)) #Widerstand 1-4
    R_b = 18.0*40.0e6+62.0*200.0e6+81.0*200.0e6 #Widerstand 5-84
    R_gesamt = R_b+R_a #Gesamtwiderstand HE1
    I_gesamt = float(-v_gesamt/R_gesamt) #Gesamtstrom durch ACC linke Seite
    L=float(25.38e-3) #Laenge eines Segmentes
    V_a = float(-R_a * I_gesamt) #Spannungsabfall 1..4
    I_a2 = float(V_a / R_a2)
#   V_qsnout = float(I_a2 * (400.0e6+100.0e9))
    
    
# Bias..Anfang Q-Snout
    T0 = float(v_vorbeschl-4.e3)
    T1 = float(T0+v_qsnout+4.e3)
    E1 = 4.e3/646.4255e-3
    E2 = -(v_qsnout+4.e3)/(223.290615e-3) #zu ergaenzen: bias-schlitze
    E3 = 0. # kein feld in der qsnout
    AddSegment2(T0,T1,10.,223.290615e-3,E1,E2,E3,-1.) #oeffnung gross

# Q-Snout : Laenge 300.4mm
    T0 = T1
    T1 = T1
    E1 = E2
    E2 = E3
    E3 = (v_qsnout-V_a)/L
    AddSegment2(T0,T1,50.e-3,300.4e-3,E1,E2,E3,-1.)


# HE1 : Ende Q-Snout..Anfang Terminal
    T0 = T1
    T1 = float(v_vorbeschl+V_a)
    E1 = E2
    E2 = E3
    E3 = 40.0e6*I_gesamt/L
    AddSegment2(T0,T1,475.e-3,L,E1,E2,E3,-1.)  #3..4
    
    deltaV = 40.0e6*I_gesamt #4..22
    for i in range(0,18):
        T0 = T1
        T1 = float(T0-deltaV) # (-) da q=-1
        E1 = E2
        E2 = E3
        E3 = deltaV/L
        AddSegment2(T0,T1,58.e-3,L,E1,E2,E3,-1.)
    
    deltaV = 200.0e6*I_gesamt #22..84
    for i in range(0,62):
        T0 = T1
        T1 = float(T0-deltaV)
        E1 = E2
        E2 = E3
        E3 = deltaV/L
        AddSegment2(T0,T1,58.e-3,L,E1,E2,E3,-1.)

# Space = 76.23mm
    T0 = T1
    T1 = T1
    E1 = E2
    E2 = 0. #elektroden verbunden
    E3 = 200.0e6*I_gesamt/L
    AddSegment2(T0,T1,58.e-3,76.23e-3,E1,E2,E3,-1.)
    
# Space..Terminal
    deltaV = 200.0e6*I_gesamt #84..165
    for i in range(0,81):
        T0 = T1
        T1 = float(T0-deltaV)
        E1 = E2
        E2 = E3
        E3 = deltaV/L
        AddSegment2(T0,T1,58.e-3,L,E1,E2,E3,-1.)
        
        
# HE2
    R_gesamt = float(9.0*120.0e6+156.0*300.0e6)
    I_gesamt = float(v_gesamt/R_gesamt) #Gesamtstrom durch ACC rechte Seite
        
# Terminal / Stripping
    T0 = T1
    T1 = T1
    E1 = E2
    E2 = 0.
    E3 = 120.0e6*I_gesamt/L
    AddSegment2(T0,T1,58.e-3,1147.0e-3,E1,E2,E3,-1.)

# Terminal..Space
    deltaV = 120.0e6*I_gesamt #1..9
    for i in range(0,9):
        T0 = T1
        T1 = float(T0+q*deltaV)
        E1 = E2
        E2 = E3
        E3 = deltaV/L
        AddSegment2(T0,T1,58.e-3,L,E1,E2,E3,q)

    deltaV = 300.0e6*I_gesamt #9..81
    for i in range(0,72):
        T0 = T1
        T1 = float(T0+q*deltaV)
        E1 = E2
        E2 = E3
        E3 = deltaV/L
        AddSegment2(T0,T1,58.e-3,L,E1,E2,E3,q)

# Space = 76.05 mm 
    T0 = T1
    T1 = T1
    E1 = E2
    E2 = 0.
    E3 = 300.0e6*I_gesamt/L
    AddSegment2(T0,T1,58.e-3,76.05e-3,E1,E2,E3,q)

# Space..letzte Elekrode
    deltaV = 300.0e6*I_gesamt #81..165
    for i in range(0,84):
        T0 = T1
        T1 = float(T0+q*deltaV)
        E1 = E2
        E2 = E3
        E3 = deltaV/L
        AddSegment2(T0,T1,58.e-3,L,E1,E2,E3,q)
        
# Letzte Elektrode..Flansch aussen = 487.8mm
    E1 = E2
    E2 = 0.
    E3 = 0.
    AddSegment2(T0,T1,58.e-3,487.75e-3,E1,E2,E3,q)


######################################################################################################################













    


s = 0.

SOURCE = []

geolines = vtk.vtkTable()
geo_s = vtk.vtkFloatArray()
geo_s.SetName("S-Achse")
geo_y = vtk.vtkFloatArray()
geo_y.SetName("geometry")
textArray = []
lastFunction = ""

try:  
    optic=ctypes.CDLL("./liblimioptic.so")
except:
    try: 
        optic=ctypes.CDLL("./liblimioptic-win.so")
    except:    
        print "\n\nFirst start?\nTrying to compile C++ code. ", 
        from os import system
        print "\nIf this is a windows OS: This will only work if MinGW is installed!\n"
        system("PAUSE")
        print "compiling climioptic.cpp"
        print "------------------------"
        system("g++ -Wall -fPIC -c climioptic.cpp")
        print "\n"
        
        print "compiling limioptic.cpp"
        print "------------------------"
        system("g++ -Wall -fPIC -c limioptic.cpp")
        print "\n"
        
        print "building library"
        print "------------------------"
        system("g++ -fPIC -shared -Wl,-soname,liblimioptic.so -o liblimioptic.so limioptic.o climioptic.o")
        print "completed."
        
        optic=ctypes.CDLL("./liblimioptic.so")
        
if (__name__=="__main__"):
    if (len(sys.argv)!=2):
        print("Usage : %s filename" % (sys.argv[0]))
    else:
        Limioptic(sys.argv[1])

###########################################################

