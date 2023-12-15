# -*- coding: utf-8 -*-

# Main program for launching Stmapy
import PyQt5.QtWidgets as QtWidgets
import sys
import logging
from stmapy import citswidget


def main():
    logging.basicConfig(format="%(levelname)s:%(message)s", level=logging.DEBUG)
    logging.getLogger('matplotlib.font_manager').disabled = True
    app = QtWidgets.QApplication(sys.argv)
    mw = citswidget.CitsWidget(None)
    mw.show()
    app.exec_()


if __name__ == "__main__":
    main()
