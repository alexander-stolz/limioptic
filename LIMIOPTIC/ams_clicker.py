import sys
from PyQt4 import QtCore as core
from PyQt4 import QtGui as gui
import time
import thread
import win32com.client as comclt

class clicker(gui.QDialog):
    def __init__(self):
        gui.QDialog.__init__(self)
        hbox = gui.QVBoxLayout()
        
        self.a = gui.QSpinBox()
        self.a.setValue(6)
        self.a.setPrefix("erhoehe alle ")
        self.a.setSuffix(" min")
        self.b = gui.QPushButton()
        self.b.setText("start")
        
        hbox.addWidget(self.a)
        hbox.addWidget(self.b)

        self.setLayout(hbox)
        
        self.connect(self.b, core.SIGNAL("clicked()"),self.los)

        self.setWindowTitle("Cologne AMS Clicker - Limioptic 2 Widget")

        

    def los(self):
        global RUN
        RUN = not RUN
        if (RUN):
            thread.start_new_thread(self.los2,(0,))
            self.b.setText("STOP")
        else:
            self.b.setText("Start")
    def los2(self,idd):
        while RUN:
            time.sleep(60*self.a.value())
            wsh.SendKeys("{RIGHT}")
            
            

RUN = False
wsh= comclt.Dispatch("WScript.Shell")

if __name__ == "__main__":
    
    app = gui.QApplication(sys.argv)
    myapp = clicker()

    myapp.show()

    sys.exit(app.exec_())
