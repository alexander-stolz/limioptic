from pylab import *
import random
from PyQt4 import QtGui as gui
from PyQt4 import QtCore as core
import sys
import time

class ImportSource():
    def __init__(self):
        self.mittel = [0]*6
        self.sigma = [0]*4
        self.Source = []
        self.Selection = []
        self.UserInteraction = UserInteraction(self)

    def LoadSource(self,File):
        self.SourceFile = File
        self.Source = []
        i = 0
        for line in open(self.SourceFile):        
            try:
                temp = [float(line.split()[0]),float(line.split()[1]),float(line.split()[2]),float(line.split()[3]),float(line.split()[4]),float(line.split()[5])]
                self.Source.append(temp)
                i += 1
                for j in range(0,6):
                    self.mittel[j]+=temp[j]
            except:
                try:
                    temp = [float(line.split()[0]),float(line.split()[1]),float(line.split()[2]),float(line.split()[3]),0.,0.]
                    self.Source.append(temp)
                    i += 1
                    for j in range(0,6):
                        self.mittel[j]+=temp[j]
                except:
                    print "something went wrong in line", i+1
                    #print line.split()
                
        for j in range(0,6):
            self.mittel[j] = self.mittel[j]/len(self.Source)
        print len(self.Source),"particles read."
        print "average energy deviation is",self.mittel[4],"permille."
        print "average mass deviation is",self.mittel[5],"permille."

    def NormalizeEnergy(self):
        for i in range(0,len(self.Source)):
            self.Source[i][4] -= self.mittel[4]

    def NormalizeMass(self):
        for i in range(0,len(self.Source)):
            self.Source[i][5] -= self.mittel[5]
        
    def ShowFits(self):
        art = "x","x'","y","y'"
        ax = [0]*4
        fig, ((ax[0],ax[1]),(ax[2],ax[3])) = subplots(nrows=2,ncols=2)
        for z in range(0,4):
            sigma = 0.
            for j in range(0,len(self.Source)):
                sigma += (self.Source[j][z]-self.mittel[z])**2
            self.sigma[z] = sqrt(sigma/len(self.Source))

            n, bins, patches = ax[z].hist([row[z] for row in self.Source],bins=arange(self.mittel[z]-3*self.sigma[z],self.mittel[z]+3*self.sigma[z],self.sigma[z]/10.),normed=True)

            y = normpdf(arange(self.mittel[z]-3*self.sigma[z],self.mittel[z]+3*self.sigma[z],self.sigma[z]/100.),self.mittel[z],self.sigma[z])
            ax[z].plot(arange(self.mittel[z]-3*self.sigma[z],self.mittel[z]+3*self.sigma[z],self.sigma[z]/100.),y)
            ax[z].grid(True)
            #ax[z].set_ylabel("Probability")
            ax[z].set_title(art[z])
            print "mu {}\t=".format(art[z]),self.mittel[z],"\nsigma\t=",self.sigma[z]

        fig.canvas.set_window_title("Source Beam - {}".format(self.SourceFile))
        tight_layout()
        show()

    def SelectRunaways(self,f):
        self.Selection = []
        for i in range(0,len(self.Source)):
            x = self.Source[i][0]
            xx = self.Source[i][1]
            y = self.Source[i][2]
            yy = self.Source[i][3]
            phi1 = arctan(xx/x)
            phi2 = arctan(yy/y)
            if(((f*self.sigma[0]*cos(phi1))**2+(f*self.sigma[1]*sin(phi1))**2 < x**2+xx**2) or ((f*self.sigma[2]*cos(phi2))**2+(f*self.sigma[3]*sin(phi2))**2 < y**2+yy**2)):
                self.Selection.append(self.Source[i])
                #print phi1, self.sigma[0],(self.sigma[0]*cos(phi1))**2+(self.sigma[1]*sin(phi1))**2 , x**2+xx**2
        print "Filter applied. Before:",len(self.Source),"Particles, after:",len(self.Selection),"Particles."
        
    def SelectRandom(self,num):
        self.Selection = []
        random.seed()
        for i in range(0,num):
            self.Selection.append(self.Source[random.randint(0,len(self.Source)-1)])
        print "Filter applied. Before:",len(self.Source),"Particles, after:",len(self.Selection),"Particles."

    def SelectBorder(self,f,step):
        self.Selection = []
        for i in range(0,360,int(step)):
            radi = radians(i)
            for j in range (0,360,int(step)):
                radj = radians(j)
                self.Selection.append([f*self.sigma[0]*cos(radi)*cos(radj),f*self.sigma[1]*cos(radi)*sin(radj),f*self.sigma[2]*sin(radi)*cos(radj),f*self.sigma[3]*sin(radi)*sin(radj),self.Source[random.choice(range(0,len(self.Source)))][4],0.])
        print "Filter applied. Before:",len(self.Source),"Particles, after:",len(self.Selection),"Particles."


class UserInteraction(gui.QDialog):
    def __init__(self,parent):
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
        for i in range(0,len(self.button)):
            vbox.addWidget(self.button[i])
        self.connect(self.button[0],core.SIGNAL("clicked()"),self.filter1)
        self.connect(self.button[1],core.SIGNAL("clicked()"),self.filter2)
        self.connect(self.button[2],core.SIGNAL("clicked()"),self.filter3)
        self.connect(self.button[3],core.SIGNAL("clicked()"),self.filter4)
        self.setLayout(vbox)
        self.adjustSize()
        self.show()

    def filter1(self):
        self.parent.SelectRunaways(3.)
        self.parent.Source = self.parent.Selection
        self.parent.ShowFits()
        self.close()
    def filter2(self):
        self.parent.SelectRandom(1000)
        self.parent.Source = self.parent.Selection
        self.parent.ShowFits()
        self.close()
    def filter3(self):
        self.parent.SelectBorder(3.,10)
        self.parent.Source = self.parent.Selection
        self.parent.ShowFits()
        self.close()
    def filter4(self):
        self.parent.Selection = self.parent.Source
        self.close()
    

if __name__ == "__main__":
    source = "C:\\Users\\Alexander\\Dropbox\\uni\\aaaa\\pavel_beam_x_v_x_y_v_y.dat"
    #source = "C:\\Users\\astolz\\Dropbox\\uni\\aaaa\\pp_to_alex.out"
    app = gui.QApplication(sys.argv)
    myapp = ImportSource(source)
    myapp.LoadSource()
    myapp.NormalizeEnergy()
    myapp.ShowFits()
    myapp.UserInteraction.ChooseFilter()
    app.exec_()
    #myapp.Source = myapp.Selection
    #myapp.ShowFits()
    sys.exit()

