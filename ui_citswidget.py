# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'ui_citswidget.ui'
#
# Created by: PyQt5 UI code generator 5.10.1
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_CitsWidget(object):
    def setupUi(self, CitsWidget):
        CitsWidget.setObjectName("CitsWidget")
        CitsWidget.resize(1392, 1281)
        self.centralwidget = QtWidgets.QWidget(CitsWidget)
        self.centralwidget.setObjectName("centralwidget")
        self.gridLayout_2 = QtWidgets.QGridLayout(self.centralwidget)
        self.gridLayout_2.setObjectName("gridLayout_2")
        self.verticalLayout_2 = QtWidgets.QVBoxLayout()
        self.verticalLayout_2.setContentsMargins(-1, -1, -1, 0)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.m_openButton = QtWidgets.QPushButton(self.centralwidget)
        self.m_openButton.setObjectName("m_openButton")
        self.horizontalLayout.addWidget(self.m_openButton)
        spacerItem = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout.addItem(spacerItem)
        self.m_topoButton = QtWidgets.QPushButton(self.centralwidget)
        self.m_topoButton.setObjectName("m_topoButton")
        self.horizontalLayout.addWidget(self.m_topoButton)
        spacerItem1 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout.addItem(spacerItem1)
        self.m_wholeLengthCutButton = QtWidgets.QPushButton(self.centralwidget)
        self.m_wholeLengthCutButton.setObjectName("m_wholeLengthCutButton")
        self.horizontalLayout.addWidget(self.m_wholeLengthCutButton)
        self.verticalLayout_2.addLayout(self.horizontalLayout)
        self.groupBox_3 = QtWidgets.QGroupBox(self.centralwidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.groupBox_3.sizePolicy().hasHeightForWidth())
        self.groupBox_3.setSizePolicy(sizePolicy)
        self.groupBox_3.setObjectName("groupBox_3")
        self.verticalLayout_5 = QtWidgets.QVBoxLayout(self.groupBox_3)
        self.verticalLayout_5.setObjectName("verticalLayout_5")
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_2.setSizeConstraint(QtWidgets.QLayout.SetDefaultConstraint)
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.label_7 = QtWidgets.QLabel(self.groupBox_3)
        self.label_7.setObjectName("label_7")
        self.horizontalLayout_2.addWidget(self.label_7)
        self.m_normButton = QtWidgets.QPushButton(self.groupBox_3)
        self.m_normButton.setObjectName("m_normButton")
        self.horizontalLayout_2.addWidget(self.m_normButton)
        self.m_channelBox = QtWidgets.QComboBox(self.groupBox_3)
        self.m_channelBox.setObjectName("m_channelBox")
        self.horizontalLayout_2.addWidget(self.m_channelBox)
        self.verticalLayout_5.addLayout(self.horizontalLayout_2)
        self.horizontalLayout_4 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_4.setObjectName("horizontalLayout_4")
        self.label = QtWidgets.QLabel(self.groupBox_3)
        self.label.setObjectName("label")
        self.horizontalLayout_4.addWidget(self.label)
        self.m_voltageBox = QtWidgets.QSpinBox(self.groupBox_3)
        self.m_voltageBox.setMaximum(100)
        self.m_voltageBox.setObjectName("m_voltageBox")
        self.horizontalLayout_4.addWidget(self.m_voltageBox)
        self.verticalLayout_5.addLayout(self.horizontalLayout_4)
        self.verticalLayout_2.addWidget(self.groupBox_3)
        self.groupBox = QtWidgets.QGroupBox(self.centralwidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.groupBox.sizePolicy().hasHeightForWidth())
        self.groupBox.setSizePolicy(sizePolicy)
        self.groupBox.setObjectName("groupBox")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.groupBox)
        self.verticalLayout.setObjectName("verticalLayout")
        self.horizontalLayout_13 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_13.setObjectName("horizontalLayout_13")
        self.label_16 = QtWidgets.QLabel(self.groupBox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_16.sizePolicy().hasHeightForWidth())
        self.label_16.setSizePolicy(sizePolicy)
        self.label_16.setObjectName("label_16")
        self.horizontalLayout_13.addWidget(self.label_16)
        self.line_4 = QtWidgets.QFrame(self.groupBox)
        self.line_4.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_4.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_4.setObjectName("line_4")
        self.horizontalLayout_13.addWidget(self.line_4)
        self.verticalLayout.addLayout(self.horizontalLayout_13)
        self.m_avgBox = QtWidgets.QCheckBox(self.groupBox)
        self.m_avgBox.setObjectName("m_avgBox")
        self.verticalLayout.addWidget(self.m_avgBox)
        self.horizontalLayout_7 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_7.setObjectName("horizontalLayout_7")
        self.label_2 = QtWidgets.QLabel(self.groupBox)
        self.label_2.setObjectName("label_2")
        self.horizontalLayout_7.addWidget(self.label_2)
        self.m_rcAvgBox = QtWidgets.QSpinBox(self.groupBox)
        self.m_rcAvgBox.setMaximum(100)
        self.m_rcAvgBox.setProperty("value", 2)
        self.m_rcAvgBox.setObjectName("m_rcAvgBox")
        self.horizontalLayout_7.addWidget(self.m_rcAvgBox)
        self.verticalLayout.addLayout(self.horizontalLayout_7)
        self.m_avgCheckBox = QtWidgets.QCheckBox(self.groupBox)
        self.m_avgCheckBox.setObjectName("m_avgCheckBox")
        self.verticalLayout.addWidget(self.m_avgCheckBox)
        self.m_avgWidget = QtWidgets.QWidget(self.groupBox)
        self.m_avgWidget.setObjectName("m_avgWidget")
        self.gridLayout_3 = QtWidgets.QGridLayout(self.m_avgWidget)
        self.gridLayout_3.setObjectName("gridLayout_3")
        self.m_belowLine = QtWidgets.QLineEdit(self.m_avgWidget)
        self.m_belowLine.setReadOnly(True)
        self.m_belowLine.setObjectName("m_belowLine")
        self.gridLayout_3.addWidget(self.m_belowLine, 3, 0, 1, 1)
        self.m_avgVButton = QtWidgets.QPushButton(self.m_avgWidget)
        self.m_avgVButton.setObjectName("m_avgVButton")
        self.gridLayout_3.addWidget(self.m_avgVButton, 4, 2, 1, 1)
        self.label_4 = QtWidgets.QLabel(self.m_avgWidget)
        self.label_4.setObjectName("label_4")
        self.gridLayout_3.addWidget(self.label_4, 0, 2, 1, 1)
        self.m_aboveBar = QtWidgets.QProgressBar(self.m_avgWidget)
        self.m_aboveBar.setProperty("value", 50)
        self.m_aboveBar.setTextVisible(False)
        self.m_aboveBar.setOrientation(QtCore.Qt.Vertical)
        self.m_aboveBar.setInvertedAppearance(True)
        self.m_aboveBar.setObjectName("m_aboveBar")
        self.gridLayout_3.addWidget(self.m_aboveBar, 2, 2, 1, 1)
        self.m_belowBar = QtWidgets.QProgressBar(self.m_avgWidget)
        self.m_belowBar.setProperty("value", 50)
        self.m_belowBar.setTextVisible(False)
        self.m_belowBar.setOrientation(QtCore.Qt.Vertical)
        self.m_belowBar.setObjectName("m_belowBar")
        self.gridLayout_3.addWidget(self.m_belowBar, 2, 0, 1, 1)
        self.label_3 = QtWidgets.QLabel(self.m_avgWidget)
        self.label_3.setObjectName("label_3")
        self.gridLayout_3.addWidget(self.label_3, 0, 0, 1, 1)
        self.label_6 = QtWidgets.QLabel(self.m_avgWidget)
        self.label_6.setObjectName("label_6")
        self.gridLayout_3.addWidget(self.label_6, 1, 3, 1, 1)
        self.m_aboveBox = QtWidgets.QSpinBox(self.m_avgWidget)
        self.m_aboveBox.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.m_aboveBox.setMaximum(100)
        self.m_aboveBox.setProperty("value", 50)
        self.m_aboveBox.setObjectName("m_aboveBox")
        self.gridLayout_3.addWidget(self.m_aboveBox, 1, 2, 1, 1)
        self.label_5 = QtWidgets.QLabel(self.m_avgWidget)
        self.label_5.setObjectName("label_5")
        self.gridLayout_3.addWidget(self.label_5, 1, 1, 1, 1)
        self.m_belowBox = QtWidgets.QSpinBox(self.m_avgWidget)
        self.m_belowBox.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.m_belowBox.setMaximum(100)
        self.m_belowBox.setProperty("value", 50)
        self.m_belowBox.setObjectName("m_belowBox")
        self.gridLayout_3.addWidget(self.m_belowBox, 1, 0, 1, 1)
        self.m_aboveLine = QtWidgets.QLineEdit(self.m_avgWidget)
        self.m_aboveLine.setReadOnly(True)
        self.m_aboveLine.setObjectName("m_aboveLine")
        self.gridLayout_3.addWidget(self.m_aboveLine, 3, 2, 1, 1)
        self.verticalLayout.addWidget(self.m_avgWidget)
        self.horizontalLayout_14 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_14.setContentsMargins(-1, 0, -1, -1)
        self.horizontalLayout_14.setObjectName("horizontalLayout_14")
        self.label_17 = QtWidgets.QLabel(self.groupBox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_17.sizePolicy().hasHeightForWidth())
        self.label_17.setSizePolicy(sizePolicy)
        self.label_17.setObjectName("label_17")
        self.horizontalLayout_14.addWidget(self.label_17)
        self.line_5 = QtWidgets.QFrame(self.groupBox)
        self.line_5.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_5.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_5.setObjectName("line_5")
        self.horizontalLayout_14.addWidget(self.line_5)
        self.verticalLayout.addLayout(self.horizontalLayout_14)
        self.m_avgSpec = QtWidgets.QPushButton(self.groupBox)
        self.m_avgSpec.setObjectName("m_avgSpec")
        self.verticalLayout.addWidget(self.m_avgSpec)
        self.horizontalLayout_15 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_15.setContentsMargins(-1, 10, -1, -1)
        self.horizontalLayout_15.setObjectName("horizontalLayout_15")
        self.label_14 = QtWidgets.QLabel(self.groupBox)
        self.label_14.setObjectName("label_14")
        self.horizontalLayout_15.addWidget(self.label_14)
        self.m_CitsAvgBox = QtWidgets.QSpinBox(self.groupBox)
        self.m_CitsAvgBox.setMinimum(1)
        self.m_CitsAvgBox.setMaximum(10000000)
        self.m_CitsAvgBox.setObjectName("m_CitsAvgBox")
        self.horizontalLayout_15.addWidget(self.m_CitsAvgBox)
        self.m_averageCitsButton = QtWidgets.QPushButton(self.groupBox)
        self.m_averageCitsButton.setObjectName("m_averageCitsButton")
        self.horizontalLayout_15.addWidget(self.m_averageCitsButton)
        self.verticalLayout.addLayout(self.horizontalLayout_15)
        self.verticalLayout_2.addWidget(self.groupBox)
        self.gridLayout_2.addLayout(self.verticalLayout_2, 0, 0, 1, 1)
        self.verticalLayout_3 = QtWidgets.QVBoxLayout()
        self.verticalLayout_3.setContentsMargins(-1, -1, -1, 0)
        self.verticalLayout_3.setObjectName("verticalLayout_3")
        self.groupBox_4 = QtWidgets.QGroupBox(self.centralwidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.groupBox_4.sizePolicy().hasHeightForWidth())
        self.groupBox_4.setSizePolicy(sizePolicy)
        self.groupBox_4.setObjectName("groupBox_4")
        self.verticalLayout_6 = QtWidgets.QVBoxLayout(self.groupBox_4)
        self.verticalLayout_6.setObjectName("verticalLayout_6")
        self.horizontalLayout_6 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_6.setObjectName("horizontalLayout_6")
        self.m_waterfallButton = QtWidgets.QRadioButton(self.groupBox_4)
        self.m_waterfallButton.setChecked(False)
        self.m_waterfallButton.setObjectName("m_waterfallButton")
        self.horizontalLayout_6.addWidget(self.m_waterfallButton)
        self.m_plot2DButton = QtWidgets.QRadioButton(self.groupBox_4)
        self.m_plot2DButton.setChecked(True)
        self.m_plot2DButton.setObjectName("m_plot2DButton")
        self.horizontalLayout_6.addWidget(self.m_plot2DButton)
        self.verticalLayout_6.addLayout(self.horizontalLayout_6)
        self.m_viewSelectedBox = QtWidgets.QCheckBox(self.groupBox_4)
        self.m_viewSelectedBox.setObjectName("m_viewSelectedBox")
        self.verticalLayout_6.addWidget(self.m_viewSelectedBox)
        self.verticalLayout_3.addWidget(self.groupBox_4)
        self.groupBox_2 = QtWidgets.QGroupBox(self.centralwidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.groupBox_2.sizePolicy().hasHeightForWidth())
        self.groupBox_2.setSizePolicy(sizePolicy)
        self.groupBox_2.setObjectName("groupBox_2")
        self.verticalLayout_4 = QtWidgets.QVBoxLayout(self.groupBox_2)
        self.verticalLayout_4.setObjectName("verticalLayout_4")
        self.m_clearSpec = QtWidgets.QPushButton(self.groupBox_2)
        self.m_clearSpec.setObjectName("m_clearSpec")
        self.verticalLayout_4.addWidget(self.m_clearSpec)
        self.horizontalLayout_5 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_5.setObjectName("horizontalLayout_5")
        self.label_9 = QtWidgets.QLabel(self.groupBox_2)
        self.label_9.setObjectName("label_9")
        self.horizontalLayout_5.addWidget(self.label_9)
        self.m_shiftXBox = QtWidgets.QLineEdit(self.groupBox_2)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.m_shiftXBox.sizePolicy().hasHeightForWidth())
        self.m_shiftXBox.setSizePolicy(sizePolicy)
        self.m_shiftXBox.setObjectName("m_shiftXBox")
        self.horizontalLayout_5.addWidget(self.m_shiftXBox)
        self.label_13 = QtWidgets.QLabel(self.groupBox_2)
        self.label_13.setObjectName("label_13")
        self.horizontalLayout_5.addWidget(self.label_13)
        self.m_shiftYBox = QtWidgets.QLineEdit(self.groupBox_2)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.m_shiftYBox.sizePolicy().hasHeightForWidth())
        self.m_shiftYBox.setSizePolicy(sizePolicy)
        self.m_shiftYBox.setObjectName("m_shiftYBox")
        self.horizontalLayout_5.addWidget(self.m_shiftYBox)
        self.verticalLayout_4.addLayout(self.horizontalLayout_5)
        self.m_logBox = QtWidgets.QCheckBox(self.groupBox_2)
        self.m_logBox.setObjectName("m_logBox")
        self.verticalLayout_4.addWidget(self.m_logBox)
        self.horizontalLayout_11 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_11.setObjectName("horizontalLayout_11")
        self.m_derivBox = QtWidgets.QCheckBox(self.groupBox_2)
        self.m_derivBox.setObjectName("m_derivBox")
        self.horizontalLayout_11.addWidget(self.m_derivBox)
        self.m_derivWindowWidget = QtWidgets.QWidget(self.groupBox_2)
        self.m_derivWindowWidget.setObjectName("m_derivWindowWidget")
        self.horizontalLayout_12 = QtWidgets.QHBoxLayout(self.m_derivWindowWidget)
        self.horizontalLayout_12.setObjectName("horizontalLayout_12")
        self.label_15 = QtWidgets.QLabel(self.m_derivWindowWidget)
        self.label_15.setObjectName("label_15")
        self.horizontalLayout_12.addWidget(self.label_15)
        self.m_derivNBox = QtWidgets.QSpinBox(self.m_derivWindowWidget)
        self.m_derivNBox.setEnabled(False)
        self.m_derivNBox.setMinimum(1)
        self.m_derivNBox.setMaximum(501)
        self.m_derivNBox.setSingleStep(2)
        self.m_derivNBox.setProperty("value", 11)
        self.m_derivNBox.setObjectName("m_derivNBox")
        self.horizontalLayout_12.addWidget(self.m_derivNBox)
        self.horizontalLayout_11.addWidget(self.m_derivWindowWidget)
        self.verticalLayout_4.addLayout(self.horizontalLayout_11)
        self.m_fitWidget = QtWidgets.QWidget(self.groupBox_2)
        self.m_fitWidget.setObjectName("m_fitWidget")
        self.gridLayout_4 = QtWidgets.QGridLayout(self.m_fitWidget)
        self.gridLayout_4.setObjectName("gridLayout_4")
        self.m_fitLowerBox = QtWidgets.QLineEdit(self.m_fitWidget)
        self.m_fitLowerBox.setEnabled(True)
        self.m_fitLowerBox.setObjectName("m_fitLowerBox")
        self.gridLayout_4.addWidget(self.m_fitLowerBox, 3, 0, 1, 1)
        self.m_fitUpperBox = QtWidgets.QLineEdit(self.m_fitWidget)
        self.m_fitUpperBox.setEnabled(True)
        self.m_fitUpperBox.setObjectName("m_fitUpperBox")
        self.gridLayout_4.addWidget(self.m_fitUpperBox, 3, 1, 1, 1)
        self.m_fitLowerLabel = QtWidgets.QLabel(self.m_fitWidget)
        self.m_fitLowerLabel.setObjectName("m_fitLowerLabel")
        self.gridLayout_4.addWidget(self.m_fitLowerLabel, 2, 0, 1, 1)
        self.m_fitUpperLabel = QtWidgets.QLabel(self.m_fitWidget)
        self.m_fitUpperLabel.setObjectName("m_fitUpperLabel")
        self.gridLayout_4.addWidget(self.m_fitUpperLabel, 2, 1, 1, 1)
        self.m_fitSpec = QtWidgets.QPushButton(self.m_fitWidget)
        self.m_fitSpec.setObjectName("m_fitSpec")
        self.gridLayout_4.addWidget(self.m_fitSpec, 0, 0, 1, 1)
        self.m_fitCustomCheckbox = QtWidgets.QCheckBox(self.m_fitWidget)
        self.m_fitCustomCheckbox.setObjectName("m_fitCustomCheckbox")
        self.gridLayout_4.addWidget(self.m_fitCustomCheckbox, 0, 1, 1, 1)
        self.verticalLayout_4.addWidget(self.m_fitWidget)
        self.verticalLayout_3.addWidget(self.groupBox_2)
        self.groupBox_5 = QtWidgets.QGroupBox(self.centralwidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.groupBox_5.sizePolicy().hasHeightForWidth())
        self.groupBox_5.setSizePolicy(sizePolicy)
        self.groupBox_5.setObjectName("groupBox_5")
        self.verticalLayout_7 = QtWidgets.QVBoxLayout(self.groupBox_5)
        self.verticalLayout_7.setObjectName("verticalLayout_7")
        self.horizontalLayout_8 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_8.setObjectName("horizontalLayout_8")
        self.m_scaleVoltage = QtWidgets.QCheckBox(self.groupBox_5)
        self.m_scaleVoltage.setChecked(True)
        self.m_scaleVoltage.setObjectName("m_scaleVoltage")
        self.horizontalLayout_8.addWidget(self.m_scaleVoltage)
        self.m_scaleMetric = QtWidgets.QCheckBox(self.groupBox_5)
        self.m_scaleMetric.setChecked(True)
        self.m_scaleMetric.setObjectName("m_scaleMetric")
        self.horizontalLayout_8.addWidget(self.m_scaleMetric)
        self.verticalLayout_7.addLayout(self.horizontalLayout_8)
        self.horizontalLayout_10 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_10.setSizeConstraint(QtWidgets.QLayout.SetFixedSize)
        self.horizontalLayout_10.setSpacing(6)
        self.horizontalLayout_10.setObjectName("horizontalLayout_10")
        self.m_vLineBox = QtWidgets.QCheckBox(self.groupBox_5)
        self.m_vLineBox.setObjectName("m_vLineBox")
        self.horizontalLayout_10.addWidget(self.m_vLineBox)
        self.m_legendBox = QtWidgets.QCheckBox(self.groupBox_5)
        self.m_legendBox.setObjectName("m_legendBox")
        self.horizontalLayout_10.addWidget(self.m_legendBox)
        self.verticalLayout_7.addLayout(self.horizontalLayout_10)
        self.m_cbarCheckBox = QtWidgets.QCheckBox(self.groupBox_5)
        self.m_cbarCheckBox.setObjectName("m_cbarCheckBox")
        self.verticalLayout_7.addWidget(self.m_cbarCheckBox)
        self.m_cbarWidget = QtWidgets.QWidget(self.groupBox_5)
        self.m_cbarWidget.setObjectName("m_cbarWidget")
        self.gridLayout = QtWidgets.QGridLayout(self.m_cbarWidget)
        self.gridLayout.setObjectName("gridLayout")
        self.m_colorBarBox = QtWidgets.QComboBox(self.m_cbarWidget)
        self.m_colorBarBox.setObjectName("m_colorBarBox")
        self.gridLayout.addWidget(self.m_colorBarBox, 0, 1, 1, 1)
        self.label_10 = QtWidgets.QLabel(self.m_cbarWidget)
        self.label_10.setObjectName("label_10")
        self.gridLayout.addWidget(self.label_10, 0, 0, 1, 1)
        self.m_cbarLowerBox = QtWidgets.QLineEdit(self.m_cbarWidget)
        self.m_cbarLowerBox.setEnabled(False)
        self.m_cbarLowerBox.setObjectName("m_cbarLowerBox")
        self.gridLayout.addWidget(self.m_cbarLowerBox, 3, 0, 1, 1)
        self.m_cbarCustomCheckbox = QtWidgets.QCheckBox(self.m_cbarWidget)
        self.m_cbarCustomCheckbox.setObjectName("m_cbarCustomCheckbox")
        self.gridLayout.addWidget(self.m_cbarCustomCheckbox, 1, 0, 1, 1)
        self.label_12 = QtWidgets.QLabel(self.m_cbarWidget)
        self.label_12.setObjectName("label_12")
        self.gridLayout.addWidget(self.label_12, 2, 1, 1, 1)
        self.m_cbarUpperBox = QtWidgets.QLineEdit(self.m_cbarWidget)
        self.m_cbarUpperBox.setEnabled(False)
        self.m_cbarUpperBox.setObjectName("m_cbarUpperBox")
        self.gridLayout.addWidget(self.m_cbarUpperBox, 3, 1, 1, 1)
        self.label_11 = QtWidgets.QLabel(self.m_cbarWidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_11.sizePolicy().hasHeightForWidth())
        self.label_11.setSizePolicy(sizePolicy)
        self.label_11.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignVCenter)
        self.label_11.setObjectName("label_11")
        self.gridLayout.addWidget(self.label_11, 2, 0, 1, 1)
        self.verticalLayout_7.addWidget(self.m_cbarWidget)
        self.verticalLayout_3.addWidget(self.groupBox_5)
        self.gridLayout_2.addLayout(self.verticalLayout_3, 0, 1, 1, 1)
        self.map_layout = QtWidgets.QVBoxLayout()
        self.map_layout.setContentsMargins(0, 0, -1, -1)
        self.map_layout.setObjectName("map_layout")
        self.m_mapWidget = MatplotlibWidget(self.centralwidget)
        self.m_mapWidget.setMinimumSize(QtCore.QSize(440, 440))
        self.m_mapWidget.setObjectName("m_mapWidget")
        self.map_layout.addWidget(self.m_mapWidget)
        self.gridLayout_2.addLayout(self.map_layout, 1, 0, 1, 1)
        self.spec_layout = QtWidgets.QVBoxLayout()
        self.spec_layout.setSizeConstraint(QtWidgets.QLayout.SetDefaultConstraint)
        self.spec_layout.setContentsMargins(0, -1, -1, -1)
        self.spec_layout.setObjectName("spec_layout")
        self.m_specWidget = MatplotlibWidget(self.centralwidget)
        self.m_specWidget.setObjectName("m_specWidget")
        self.spec_layout.addWidget(self.m_specWidget)
        self.gridLayout_2.addLayout(self.spec_layout, 1, 1, 1, 1)
        CitsWidget.setCentralWidget(self.centralwidget)
        self.m_menuBar = QtWidgets.QMenuBar(CitsWidget)
        self.m_menuBar.setGeometry(QtCore.QRect(0, 0, 1392, 22))
        self.m_menuBar.setObjectName("m_menuBar")
        CitsWidget.setMenuBar(self.m_menuBar)
        self.m_statusBar = QtWidgets.QStatusBar(CitsWidget)
        self.m_statusBar.setObjectName("m_statusBar")
        CitsWidget.setStatusBar(self.m_statusBar)

        self.retranslateUi(CitsWidget)
        QtCore.QMetaObject.connectSlotsByName(CitsWidget)

    def retranslateUi(self, CitsWidget):
        _translate = QtCore.QCoreApplication.translate
        CitsWidget.setWindowTitle(_translate("CitsWidget", "STM Data Analysis"))
        self.m_openButton.setText(_translate("CitsWidget", "Open CITS"))
        self.m_topoButton.setText(_translate("CitsWidget", "Draw topo"))
        self.m_wholeLengthCutButton.setText(_translate("CitsWidget", "Whole length cut"))
        self.groupBox_3.setTitle(_translate("CitsWidget", "CITS parameters"))
        self.label_7.setText(_translate("CitsWidget", "Displayed channel"))
        self.m_normButton.setText(_translate("CitsWidget", "Normalize current channel"))
        self.label.setText(_translate("CitsWidget", "V/Z Index"))
        self.groupBox.setTitle(_translate("CitsWidget", "Averaging"))
        self.label_16.setText(_translate("CitsWidget", "Multiple spectra"))
        self.m_avgBox.setText(_translate("CitsWidget", "Select spectra with left-click (uncheck to get the result)"))
        self.label_2.setText(_translate("CitsWidget", "Number of pixels for right-click average"))
        self.m_avgCheckBox.setText(_translate("CitsWidget", "Average with respect to values"))
        self.m_avgVButton.setText(_translate("CitsWidget", "Start the averaging"))
        self.label_4.setText(_translate("CitsWidget", "Maximum value for the above averaging"))
        self.label_3.setText(_translate("CitsWidget", "Minimum value for the below averaging"))
        self.label_6.setText(_translate("CitsWidget", "% of colormap"))
        self.label_5.setText(_translate("CitsWidget", "% of colormap"))
        self.label_17.setText(_translate("CitsWidget", "Whole CITS"))
        self.m_avgSpec.setText(_translate("CitsWidget", "Plot average spectrum"))
        self.label_14.setText(_translate("CitsWidget", "Number of spectra to average along X direction"))
        self.m_averageCitsButton.setText(_translate("CitsWidget", "Average CITS"))
        self.groupBox_4.setTitle(_translate("CitsWidget", "Cuts"))
        self.m_waterfallButton.setText(_translate("CitsWidget", "Waterfall"))
        self.m_plot2DButton.setText(_translate("CitsWidget", "2D plot"))
        self.m_viewSelectedBox.setText(_translate("CitsWidget", "View selected spectra"))
        self.groupBox_2.setTitle(_translate("CitsWidget", "Spectra plot"))
        self.m_clearSpec.setText(_translate("CitsWidget", "Clear spectra"))
        self.label_9.setText(_translate("CitsWidget", "Shift plot along X"))
        self.m_shiftXBox.setText(_translate("CitsWidget", "0.0"))
        self.label_13.setText(_translate("CitsWidget", "Shift plot along Y"))
        self.m_shiftYBox.setText(_translate("CitsWidget", "0.0"))
        self.m_logBox.setText(_translate("CitsWidget", "Plot log of the spectrum"))
        self.m_derivBox.setText(_translate("CitsWidget", "Plot derivative"))
        self.label_15.setText(_translate("CitsWidget", "Window length (in number of pts) :"))
        self.m_fitLowerBox.setText(_translate("CitsWidget", "0.0"))
        self.m_fitUpperBox.setText(_translate("CitsWidget", "0.0"))
        self.m_fitLowerLabel.setText(_translate("CitsWidget", "Lower limit"))
        self.m_fitUpperLabel.setText(_translate("CitsWidget", "Upper limit"))
        self.m_fitSpec.setText(_translate("CitsWidget", "Fit spectra"))
        self.m_fitCustomCheckbox.setText(_translate("CitsWidget", "Use a custom range"))
        self.groupBox_5.setTitle(_translate("CitsWidget", "Display parameters"))
        self.m_scaleVoltage.setText(_translate("CitsWidget", "Scale in Volts"))
        self.m_scaleMetric.setText(_translate("CitsWidget", "Scale in metric units"))
        self.m_vLineBox.setText(_translate("CitsWidget", "Voltage index guideline"))
        self.m_legendBox.setText(_translate("CitsWidget", "Deactivate legend"))
        self.m_cbarCheckBox.setText(_translate("CitsWidget", "Colorbar settings"))
        self.label_10.setText(_translate("CitsWidget", "Colorbar to use"))
        self.m_cbarLowerBox.setPlaceholderText(_translate("CitsWidget", "0"))
        self.m_cbarCustomCheckbox.setText(_translate("CitsWidget", "Use custom limits"))
        self.label_12.setText(_translate("CitsWidget", "Custom upper limit"))
        self.m_cbarUpperBox.setPlaceholderText(_translate("CitsWidget", "0"))
        self.label_11.setText(_translate("CitsWidget", "Custom lower limit"))

from matplotlibwidget import MatplotlibWidget
