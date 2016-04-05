# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 12:31:44 2015

@author: LH242250
"""

# Main program for launching STM_Data_Analysis
import pyqtv2
from PyQt4.QtGui import *
from PyQt4.QtCore import *
import sys
import citswidget

def main(args):
    app=QApplication(args)
    mw=citswidget.CitsWidget(None)
    mw.show()
    app.exec_()
    
if __name__ == "__main__":
    main(sys.argv)