import sys
from PyQt4 import QtCore as core
from PyQt4 import QtGui as gui
import math

class Spicker(gui.QDialog):
    def __init__(self):
        gui.QDialog.__init__(self)
        vbox = gui.QVBoxLayout()

        inputbox = gui.QGridLayout()
        self.q = gui.QSpinBox()
        self.q.setPrefix("q= ")
        self.q.setSuffix(" e")
        self.q.setSingleStep(1)
        self.q.setMinimum(1)
        self.q.setValue(4)
        self.vext = gui.QDoubleSpinBox()
        self.vext.setPrefix("Ext= ")
        self.vext.setSuffix(" kV")
        self.vext.setValue(35)
        self.term = gui.QDoubleSpinBox()
        self.term.setPrefix("Terminal= ")
        self.term.setSuffix(" MV")
        self.term.setDecimals(5)
        self.term.setSingleStep(.5)
        self.term.setValue(5.5)
        self.a = gui.QDoubleSpinBox()
        self.a.setPrefix("a= ")
        self.a.setSuffix(" mm")
        self.a.setValue(25.17)
        self.t = gui.QDoubleSpinBox()
        self.t.setDecimals(3)
        self.t.setRange(-100,100)
        self.t.setPrefix("Korrektur= ")
        self.t.setSuffix(" %")
        self.t.setValue(5.686)

        mbox = gui.QHBoxLayout()
        self.min = gui.QSpinBox()
        self.min.setSingleStep(1)
        self.min.setPrefix("in= ")
        self.min.setSuffix(" u")
        self.min.setRange(1,1000)
        self.min.setValue(14)
        self.mout = gui.QSpinBox()
        self.mout.setSingleStep(1)
        self.mout.setPrefix("out= ")
        self.mout.setSuffix(" u")
        self.mout.setRange(1,1000)
        self.mout.setValue(14)
        mbox.addWidget(self.min)
        mbox.addWidget(self.mout)

        inputbox.addWidget(self.q,0,0)
        inputbox.addWidget(self.vext,0,1)
        inputbox.addWidget(self.term,0,2)
        inputbox.addWidget(self.a,1,0)
        inputbox.addWidget(self.t,1,1)
        inputbox.addLayout(mbox,1,2)

        vbox.addLayout(inputbox)

        self.GIC = gui.QRadioButton("GIC")
        self.TOF = gui.QRadioButton("TOF")
        self.SIG = gui.QRadioButton("SiGI")
        self.AFD = gui.QRadioButton("AFD")
        detbox = gui.QGridLayout()
        detbox.addWidget(gui.QLabel("Ausgang: "),0,0)
        detbox.addWidget(self.GIC,0,1)
        detbox.addWidget(self.TOF,0,2)
        detbox.addWidget(self.SIG,0,3)
        detbox.addWidget(self.AFD,0,4)
        vbox.addLayout(detbox)

        #vbox.addWidget(space)

        qplayout = gui.QGridLayout()
        labels = ["<b>QPT:<\b>","ACC","HES","HEE","DSW"]
        self.label = []
        for i in range(0,5):
            self.label.append(gui.QLabel(labels[i]))
            qplayout.addWidget(self.label[i],i,0)
        #self.percentlable = gui.QLabel("[%]")
        #qplayout.addWidget(self.percentlable,0,1)
        self.percent = []
        for i in range(0,4):
            self.percent.append(gui.QDoubleSpinBox())
            self.percent[i].setDecimals(3)
            self.percent[i].setRange(-150,150)
            self.percent[i].setPrefix("[%]= ")
            self.percent[i].setSuffix(" %")
            qplayout.addWidget(self.percent[i],i+1,1)
        #self.astlabel = gui.QLabel("[Ast]")
        #qplayout.addWidget(self.astlabel,0,2)
        self.ast = []
        for i in range(0,4):
            self.ast.append(gui.QDoubleSpinBox())
            self.ast[i].setDecimals(3)
            self.ast[i].setPrefix("[Ast]= ")
            self.ast[i].setSuffix(" %")
            self.ast[i].setRange(-100,100)
            qplayout.addWidget(self.ast[i],i+1,2)
        qplayout.setColumnStretch(0,0)
        qplayout.setColumnStretch(1,100)
        qplayout.setColumnStretch(2,100)
        vbox.addLayout(qplayout)

        #vbox.addWidget(space)

        msalayout = gui.QGridLayout()
        mlabels   = ["<b>MSA:<\b>","BI","HEM","DSW","AFD"]
        #msalayout.addWidget(gui.QLabel("Strom"),0,1)
        #msalayout.addWidget(gui.QLabel("Feld"),0,2)
        self.mlabel = []
        for i in range(0,5):
            self.mlabel.append(gui.QLabel(mlabels[i]))
            msalayout.addWidget(self.mlabel[i],i,0)
        self.mstrom = []
        self.mfeld  = []
        for i in range(1,5):
            self.mstrom.append(gui.QDoubleSpinBox())
            self.mstrom[i-1].setDecimals(3)
            self.mstrom[i-1].setPrefix("Strom= ")
            self.mstrom[i-1].setSuffix(" A")
            self.mstrom[i-1].setRange(-1000,1000)
            msalayout.addWidget(self.mstrom[i-1],i,1)

            self.mfeld.append(gui.QDoubleSpinBox())
            self.mfeld[i-1].setDecimals(6)
            self.mfeld[i-1].setPrefix("Feld= ")
            self.mfeld[i-1].setSuffix(" T")
            self.mfeld[i-1].setRange(-100,100)
            msalayout.addWidget(self.mfeld[i-1],i,2)
        msalayout.setColumnStretch(1,100)
        msalayout.setColumnStretch(2,100)
        vbox.addLayout(msalayout)


        esdlayout = gui.QGridLayout()
        elabels = ["<b>ESD:<\b>","BI","HEE"]

        self.elabel = []
        for i in range(0,3):
            self.elabel.append(gui.QLabel(elabels[i]))
            esdlayout.addWidget(self.elabel[i],0,i)
        self.espannung = []
        for i in range(1,3):
            self.espannung.append(gui.QDoubleSpinBox())
            self.espannung[i-1].setDecimals(3)
            self.espannung[i-1].setPrefix("Spannung= ")
            self.espannung[i-1].setSuffix(" kV")
            esdlayout.addWidget(self.espannung[i-1],1,i)
        esdlayout.setColumnStretch(1,100)
        esdlayout.setColumnStretch(2,100)
        vbox.addLayout(esdlayout)


        self.setLayout(vbox)


        self.connect(self.q,    core.SIGNAL("valueChanged(int)"),self.rechne)
        self.connect(self.min,  core.SIGNAL("valueChanged(int)"),self.rechne)
        self.connect(self.mout, core.SIGNAL("valueChanged(int)"),self.rechne)
        self.connect(self.vext, core.SIGNAL("valueChanged(double)"),self.rechne)
        self.connect(self.term, core.SIGNAL("valueChanged(double)"),self.rechne)
        self.connect(self.a,    core.SIGNAL("valueChanged(double)"),self.rechne)
        self.connect(self.t,    core.SIGNAL("valueChanged(double)"),self.rechne)
        self.connect(self.term, core.SIGNAL("valueChanged(double)"),self.rechne)
        self.connect(self.GIC,  core.SIGNAL("clicked()"),self.rechne)
        self.connect(self.TOF,  core.SIGNAL("clicked()"),self.rechne)
        self.connect(self.SIG,  core.SIGNAL("clicked()"),self.rechne)
        self.connect(self.AFD,  core.SIGNAL("clicked()"),self.rechne)

        self.setWindowTitle("Cologne AMS Spicker - Limioptic 2 Widget")

        self.rechne()

    def rechne(self):
        q = self.q.value()
        T = self.vext.value()*1000.+self.term.value()*(float(self.mout.value())/self.min.value()+q)*1.e6*(self.t.value()/100.+1.)
        a = self.a.value()/1000.
        xi = 30000.*q/T/a**2
        chi = math.sqrt(2.*self.mout.value()*T)/q*math.sqrt(1.660538921e-27 / 1.602176565e-19)

        for  i in range(0,4):
            self.percent[i].setValue(ky[i]/xi*100.)
            self.ast[i].setValue((kx[i]-ky[i])/xi*100.)

        if (self.GIC.isChecked()): radius[2] = 1.545
        if (self.TOF.isChecked()): radius[2] = 1.545*30./20.
        if (self.SIG.isChecked()): radius[2] = 1.545*30./(-15.)
        if (self.AFD.isChecked()): radius[2] = -1.545

        for i in range(1,4):
            feld = chi/radius[i]
            self.mfeld[i].setValue(feld)
            self.mstrom[i].setValue(self.magnetstrom(i,feld))
        feld = math.sqrt(2.*self.min.value()*self.vext.value()*1000.)*math.sqrt(1.660538921e-27 / 1.602176565e-19)/radius[0]
        self.mfeld[0].setValue(feld)
        self.mstrom[0].setValue(self.magnetstrom(0,feld))

        self.espannung[0].setValue(self.vext.value()*.05325/.469)
        self.espannung[1].setValue(T/q*.032/2.816/1000.)

    def magnetstrom(self,m,x):
        if (m == 0): return -1.65809+234.238*x-336.824*x**2+1237.62*x**3-2011.29*x**4+1467.22*x**5-378.914*x**6
        if (m == 1): return -0.680739+187.943*x-7.24334*x**2-30.7645*x**3+94.3876*x**4-88.3585*x**5+28.931*x**6
        if (m == 2):
            if (x >= 0): return .357713+127.112*x+58.7234*x**2-192.05*x**3+306.316*x**4-230.378*x**5+66.9716*x**6
            if (x <  0): return -(.357713-127.112*x+58.7234*x**2+192.05*x**3+306.316*x**4+230.378*x**5+66.9716*x**6)
        if (m == 3): return 0.370571+163.385*x+4.48565*x**2



radius = [.403,2.112,1.611,.455] #msas
steigung = [203.56,183.38,135.95,0.] #AFD fehlt noch
kx = [4.32,3.74,2.33,3.65]
ky = [4.72,3.91,3.05,3.67]

if __name__ == "__main__":

    app = gui.QApplication(sys.argv)
    myapp = Spicker()

    myapp.show()

    sys.exit(app.exec_())
