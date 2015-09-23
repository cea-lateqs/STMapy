# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'ui_citswidget.ui'
#
# Created: Tue Sep 22 18:27:16 2015
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
        CitsWidget.resize(941, 588)
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
        self.gridLayout_2.addLayout(self.spec_layout, 4, 4, 1, 1)
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
        self.m_clearSpec = QtGui.QPushButton(self.dockWidgetContents)
        self.m_clearSpec.setObjectName(_fromUtf8("m_clearSpec"))
        self.gridLayout_2.addWidget(self.m_clearSpec, 1, 4, 1, 1)
        self.m_avgSpec = QtGui.QPushButton(self.dockWidgetContents)
        self.m_avgSpec.setObjectName(_fromUtf8("m_avgSpec"))
        self.gridLayout_2.addWidget(self.m_avgSpec, 0, 4, 1, 1)
        self.map_layout = QtGui.QVBoxLayout()
        self.map_layout.setContentsMargins(0, 0, -1, -1)
        self.map_layout.setObjectName(_fromUtf8("map_layout"))
        self.m_mapWidget = MatplotlibWidget(self.dockWidgetContents)
        self.m_mapWidget.setMinimumSize(QtCore.QSize(440, 440))
        self.m_mapWidget.setObjectName(_fromUtf8("m_mapWidget"))
        self.map_layout.addWidget(self.m_mapWidget)
        self.gridLayout_2.addLayout(self.map_layout, 4, 3, 1, 1)
        self.m_scaleVoltage = QtGui.QCheckBox(self.dockWidgetContents)
        self.m_scaleVoltage.setChecked(True)
        self.m_scaleVoltage.setObjectName(_fromUtf8("m_scaleVoltage"))
        self.gridLayout_2.addWidget(self.m_scaleVoltage, 2, 4, 1, 1)
        self.m_avgBox = QtGui.QCheckBox(self.dockWidgetContents)
        self.m_avgBox.setObjectName(_fromUtf8("m_avgBox"))
        self.gridLayout_2.addWidget(self.m_avgBox, 2, 3, 1, 1)
        self.m_derivBox = QtGui.QCheckBox(self.dockWidgetContents)
        self.m_derivBox.setObjectName(_fromUtf8("m_derivBox"))
        self.gridLayout_2.addWidget(self.m_derivBox, 3, 3, 1, 1)
        self.m_vLineBox = QtGui.QCheckBox(self.dockWidgetContents)
        self.m_vLineBox.setObjectName(_fromUtf8("m_vLineBox"))
        self.gridLayout_2.addWidget(self.m_vLineBox, 3, 4, 1, 1)
        CitsWidget.setWidget(self.dockWidgetContents)

        self.retranslateUi(CitsWidget)
        QtCore.QMetaObject.connectSlotsByName(CitsWidget)

    def retranslateUi(self, CitsWidget):
        CitsWidget.setWindowTitle(_translate("CitsWidget", "STM Data Analysis", None))
        self.m_openButton.setText(_translate("CitsWidget", "Open CITS", None))
        self.label.setText(_translate("CitsWidget", "Voltage", None))
        self.m_fwdButton.setText(_translate("CitsWidget", "Forward", None))
        self.m_bwdButton.setText(_translate("CitsWidget", "Backward", None))
        self.m_clearSpec.setText(_translate("CitsWidget", "Clear spectra", None))
        self.m_avgSpec.setText(_translate("CitsWidget", "Average on map", None))
        self.m_scaleVoltage.setText(_translate("CitsWidget", "Scale in Volts", None))
        self.m_avgBox.setText(_translate("CitsWidget", "Average spectra by selection", None))
        self.m_derivBox.setText(_translate("CitsWidget", "Plot derivative", None))
        self.m_vLineBox.setText(_translate("CitsWidget", "Voltage guideline", None))

from matplotlibwidget import MatplotlibWidget
