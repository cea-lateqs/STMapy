# -*- coding: utf-8 -*-

# Main program for launching Scampy
import PyQt5.QtWidgets as QtWidgets
import sys
from scampy import citswidget


def main(args):
    app = QtWidgets.QApplication(args)
    mw = citswidget.CitsWidget(None)
    mw.show()
    app.exec_()


if __name__ == "__main__":
    main(sys.argv)
