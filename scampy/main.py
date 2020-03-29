# -*- coding: utf-8 -*-

# Main program for launching Scampy
import PyQt5.QtWidgets as QtWidgets
import sys
from scampy import citswidget


def main():
    app = QtWidgets.QApplication(sys.argv)
    mw = citswidget.CitsWidget(None)
    mw.show()
    app.exec_()


if __name__ == "__main__":
    main()
