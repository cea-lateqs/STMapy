# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'ui_citswidget.ui'
#
# Created: Thu Aug 13 09:50:39 2015
#      by: PyQt4 UI code generator 4.9.6
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

try:
    _encoding = QtGui.QApplication.UnicodeUTF8
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig)

class Ui_CitsWidget(object):
    def setupUi(self, CitsWidget):
        CitsWidget.setObjectName(_fromUtf8("CitsWidget"))
        CitsWidget.resize(698, 535)
        self.dockWidgetContents = QtGui.QWidget()
        self.dockWidgetContents.setObjectName(_fromUtf8("dockWidgetContents"))
        self.gridLayout_2 = QtGui.QGridLayout(self.dockWidgetContents)
        self.gridLayout_2.setObjectName(_fromUtf8("gridLayout_2"))
        self.spec_layout = QtGui.QVBoxLayout()
        self.spec_layout.setSizeConstraint(QtGui.QLayout.SetDefaultConstraint)
        self.spec_layout.setContentsMargins(0, -1, -1, -1)
        self.spec_layout.setObjectName(_fromUtf8("spec_layout"))
        self.m_specWidget = MatplotlibWidget(self.dockWidgetContents)
        self.m_specWidget.setObjectName(_fromUtf8("m_specWidget"))
        self.spec_layout.addWidget(self.m_specWidget)
        self.gridLayout_2.addLayout(self.spec_layout, 3, 4, 1, 1)
        self.horizontalLayout = QtGui.QHBoxLayout()
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        self.m_openButton = QtGui.QPushButton(self.dockWidgetContents)
        self.m_openButton.setObjectName(_fromUtf8("m_openButton"))
        self.horizontalLayout.addWidget(self.m_openButton)
        self.label = QtGui.QLabel(self.dockWidgetContents)
        self.label.setObjectName(_fromUtf8("label"))
        self.horizontalLayout.addWidget(self.label)
        self.m_voltageBox = QtGui.QDoubleSpinBox(self.dockWidgetContents)
        self.m_voltageBox.setDecimals(0)
        self.m_voltageBox.setObjectName(_fromUtf8("m_voltageBox"))
        self.horizontalLayout.addWidget(self.m_voltageBox)
        self.gridLayout_2.addLayout(self.horizontalLayout, 0, 3, 1, 1)
        self.m_mapWidget = MatplotlibWidget(self.dockWidgetContents)
        self.m_mapWidget.setObjectName(_fromUtf8("m_mapWidget"))
        self.gridLayout_2.addWidget(self.m_mapWidget, 3, 3, 1, 1)
        self.horizontalLayout_2 = QtGui.QHBoxLayout()
        self.horizontalLayout_2.setSizeConstraint(QtGui.QLayout.SetDefaultConstraint)
        self.horizontalLayout_2.setObjectName(_fromUtf8("horizontalLayout_2"))
        self.m_fwdButton = QtGui.QRadioButton(self.dockWidgetContents)
        self.m_fwdButton.setObjectName(_fromUtf8("m_fwdButton"))
        self.horizontalLayout_2.addWidget(self.m_fwdButton)
        self.m_bwdButton = QtGui.QRadioButton(self.dockWidgetContents)
        self.m_bwdButton.setObjectName(_fromUtf8("m_bwdButton"))
        self.horizontalLayout_2.addWidget(self.m_bwdButton)
        self.gridLayout_2.addLayout(self.horizontalLayout_2, 1, 3, 1, 1)
        CitsWidget.setWidget(self.dockWidgetContents)

        self.retranslateUi(CitsWidget)
        QtCore.QMetaObject.connectSlotsByName(CitsWidget)

    def retranslateUi(self, CitsWidget):
        CitsWidget.setWindowTitle(_translate("CitsWidget", "DockWidget", None))
        self.m_openButton.setText(_translate("CitsWidget", "Open CITS", None))
        self.label.setText(_translate("CitsWidget", "Voltage", None))
        self.m_fwdButton.setText(_translate("CitsWidget", "Forward", None))
        self.m_bwdButton.setText(_translate("CitsWidget", "Backward", None))

from matplotlibwidget import MatplotlibWidget
