"""
This file defines the optical tab of ErwinJr, for simulation and optimization
of 1D waveguiding
"""

from functools import partial
import numpy as np
from PyQt5.QtCore import pyqtSignal, Qt
from PyQt5.QtWidgets import (QWidget, QVBoxLayout, QHBoxLayout, QGridLayout,
                             QLabel, QComboBox, QGroupBox, QDoubleSpinBox,
                             QSpinBox,
                             QPushButton, QTableWidget, QTableWidgetItem,
                             QTextEdit)
from OptStrata import OptStratum, rIdx, Alloy, Dopable
from EJcanvas import EJcanvas
from customQTClass import mtrlComboBox
mtrlList = list(rIdx.keys()) + list(Alloy.keys())
facetList = ('cleaved', 'perfect AR', 'perfect HR', 'custom')


class OpticalTab(QWidget):
    """The Optical Tab of ErwinJr. This is designed to be a GUI wrapper of
    the class Stratum
    Member variable

    --- GUI widget ---
    1st column (settingBox):
        wlBox
        mtrlBox
        ridgeBox

    2nd column (strataBox):
        insertButton  deleteButton
        layerTable
        solveButton   OptimizeButton

    3rd column (figureBox):
        optCanvas
    """
    dirty = pyqtSignal()

    def __init__(self, stratum=None, parent=None):
        super(OpticalTab, self).__init__(parent)
        self.stratum = stratum if stratum else OptStratum(3.0)
        # TODO: put gain as part of ridx?
        self.alphas = {}
        self.periods = {}
        self.ridgeLength = 3.0
        self.ridgeLoss = 0.0
        self.Lps = {}
        self.facet1 = facetList[0]
        self.facet2 = facetList[0]
        opticalLayout = QHBoxLayout()
        opticalLayout.setSpacing(0)
        opticalLayout.addLayout(self._settingBox(350))
        opticalLayout.addLayout(self._strataBox(320))
        opticalLayout.addLayout(self._figureBox())
        self.setLayout(opticalLayout)
        self.dirty.connect(self.update)

    def _settingBox(self, width):
        settingBox = QVBoxLayout()

        settingBox.addWidget(QLabel('<b><center>Wavelength</center></b>'))
        self.wlBox = QDoubleSpinBox()
        self.wlBox.setDecimals(1)
        self.wlBox.setValue(self.stratum.wl)
        self.wlBox.setSuffix(' μm')
        self.wlBox.setRange(0.0, 100.0)
        self.wlBox.setMaximumWidth(width)
        self.wlBox.valueChanged[float].connect(self.input_wl)
        settingBox.addWidget(self.wlBox)

        mtrlGroupBox = QGroupBox("Custom Material")
        mtrlGroupBox.setMaximumWidth(width)
        mtrlLayout = QGridLayout()
        mtrlLayout.addWidget(QLabel('<center><b>material Name</b></center>'),
                             0, 0, 1, 2)
        self.mtrlsBox = QComboBox()
        self.mtrlsBox.addItems(self.stratum.custom.keys())
        self.mtrlsBox.currentIndexChanged[str].connect(self.set_MtrlBox)
        mtrlLayout.addWidget(self.mtrlsBox, 1, 0, 1, 2)
        mtrlLayout.addWidget(QLabel('<center><b>n<sub>eff</sub></b></center>'),
                             2, 0)
        self.rIdxBox = QDoubleSpinBox()
        self.rIdxBox.valueChanged[float].connect(self.input_rIdx)
        mtrlLayout.addWidget(self.rIdxBox, 3, 0)
        mtrlLayout.addWidget(QLabel('<center><b>gain/loss α</b></center>'),
                             2, 1)
        self.alphaBox = QDoubleSpinBox()
        self.alphaBox.setSuffix(" /cm")
        self.alphaBox.valueChanged[float].connect(self.input_alpha)
        mtrlLayout.addWidget(self.alphaBox, 3, 1)
        mtrlLayout.addWidget(QLabel('<center><b>Period Length</b></center>'),
                             4, 0)
        self.periodBox = QDoubleSpinBox()
        self.periodBox.setSuffix(" Å")
        if self.stratum.custom:
            mtrl = self.stratum.custom.keys()[1]
            self.mtrlsBox.setCurrentText(mtrl)
            self.rIdxBox.setValue(self.stratum.custom[mtrl])
            self.alphaBox.setValue(self.alphas[mtrl])
            self.periodBox.setValue(self.Lps[mtrl])
        else:
            self.rIdxBox.setValue(1.0)
            self.rIdxBox.setEnabled(False)
            self.alphaBox.setValue(0.0)
            self.alphaBox.setEnabled(False)
            self.periodBox.setValue(0.0)
            self.periodBox.setEnabled(False)
        mtrlLayout.addWidget(self.periodBox, 5, 0)
        mtrlLayout.addWidget(QLabel('<center><b>Periods</b></center>'), 4, 1)
        self.repeatBox = QSpinBox()
        self.repeatBox.setValue(1)
        self.repeatBox.valueChanged[int].connect(self.input_repeats)
        mtrlLayout.addWidget(self.repeatBox, 5, 1)
        mtrlLayout.setHorizontalSpacing(1)
        mtrlLayout.setVerticalSpacing(3)
        mtrlGroupBox.setLayout(mtrlLayout)
        settingBox.addWidget(mtrlGroupBox)

        ridgeBox = QGroupBox("Ridge Geometry")
        ridgeBox.setMaximumWidth(width)
        ridgeLayout = QGridLayout()
        ridgeLayout.addWidget(QLabel("<center><b>Fecet1</b></center>"), 0, 0)
        self.facetBox1 = QComboBox()
        self.facetBox1.addItems(facetList)
        self.facetBox1.setCurrentText(self.facet1)
        self.facetBox1.currentIndexChanged[int].connect(self.input_facet1)
        ridgeLayout.addWidget(self.facetBox1, 1, 0)
        ridgeLayout.addWidget(QLabel("<center><b>Reflectivity</b></center>"),
                              2, 0, 1, 2)
        self.refBox1 = QDoubleSpinBox()
        self.refBox1.setDecimals(1)
        self.refBox1.setRange(0.0, 100.0)
        self.refBox1.setSuffix(" %")
        self.refBox1.setEnabled(False)
        self.refBox1.setValue(0.0)
        ridgeLayout.addWidget(self.refBox1, 3, 0)

        ridgeLayout.addWidget(QLabel("<center><b>Fecet2</b></center>"), 0, 1)
        self.facetBox2 = QComboBox()
        self.facetBox2.addItems(facetList)
        self.facetBox2.setCurrentText(self.facet2)
        self.facetBox2.currentIndexChanged[int].connect(self.input_facet2)
        ridgeLayout.addWidget(self.facetBox2, 1, 1)
        self.refBox2 = QDoubleSpinBox()
        self.refBox2.setDecimals(1)
        self.refBox2.setRange(0.0, 100.0)
        self.refBox2.setSuffix(" %")
        self.refBox2.setEnabled(False)
        self.refBox2.setValue(0.0)
        ridgeLayout.addWidget(self.refBox2, 3, 1)

        ridgeLayout.addWidget(QLabel("<center><b>Ridge Length</b></center>"),
                              4, 0)
        self.ridgeLengthBox = QDoubleSpinBox()
        self.ridgeLengthBox.setDecimals(1)
        self.ridgeLengthBox.setRange(0.0, 20.0)
        self.ridgeLengthBox.setSingleStep(1)
        self.ridgeLengthBox.setValue(self.ridgeLength)
        self.ridgeLengthBox.setSuffix(" mm")
        self.ridgeLengthBox.valueChanged[float].connect(self.input_ridgeL)
        ridgeLayout.addWidget(self.ridgeLengthBox, 5, 0)
        ridgeLayout.addWidget(QLabel("<center><b>Mirror Loss: </b></center>"),
                              4, 1)
        self.mirrorLoss = QLabel("<center>0.0 %</center>")
        ridgeLayout.addWidget(self.mirrorLoss, 5, 1)
        ridgeLayout.setHorizontalSpacing(1)
        ridgeLayout.setVerticalSpacing(3)
        ridgeBox.setLayout(ridgeLayout)
        settingBox.addWidget(ridgeBox)

        settingBox.addStretch()
        return settingBox
        # _settingBox end

    def _strataBox(self, width):
        strataBox = QGridLayout()
        strataBox.setSpacing(5)
        self.insertButton = QPushButton("Insert Strata")
        self.insertButton.clicked.connect(self.insert_strata)
        self.deleteButton = QPushButton("Delete Strata")
        self.deleteButton.clicked.connect(self.delete_strata)
        strataBox.addWidget(self.insertButton, 0, 0)
        strataBox.addWidget(self.deleteButton, 0, 1)
        self.strataTable = QTableWidget()
        self.strataTable.setMaximumWidth(width)
        self.strataTable.setMinimumWidth(width)
        self.strataTable.setMinimumHeight(300)
        self.strataTable.setMaximumHeight(400)
        self.strataTable.setSelectionBehavior(QTableWidget.SelectRows)
        self.strataTable.setSelectionMode(QTableWidget.SingleSelection)
        self.strataTable.itemChanged.connect(self.strataTable_item)
        self.strataTable.itemSelectionChanged.connect(self.strataTable_select)
        self.strataTable_refresh()
        strataBox.addWidget(self.strataTable, 1, 0, 1, 2)
        self.solveButton = QPushButton("Solve Mode")
        self.solveButton.clicked.connect(self.solve)
        self.optimizeButton = QPushButton("Optimize Strata")
        self.optimizeButton.clicked.connect(self.optimizeStrata)
        strataBox.addWidget(self.solveButton, 2, 0)
        strataBox.addWidget(self.optimizeButton, 2, 1)
        strataBox.addWidget(QLabel("Performance"), 3, 0, 1, 2)
        self.resultBox = QTextEdit()
        self.resultBox.setMinimumWidth(width)
        self.resultBox.setMaximumWidth(width)
        self.resultBox.setMinimumHeight(120)
        strataBox.addWidget(self.resultBox, 4, 0, 1, 2)

        strataBox.setRowStretch(5, 0)
        return strataBox

    def _figureBox(self):
        figureBox = QVBoxLayout()
        self.optCanvas = EJcanvas(xlabel='Position (μm)',
                                  ylabel='Refractive index')
        figureBox.addWidget(self.optCanvas)
        return figureBox

    def input_wl(self, wl):
        self.stratum.setWl(wl)

    def set_MtrlBox(self, mtrl):
        if mtrl not in self.stratum.custom:
            self.stratum.custom[mtrl] = 1.0
        if mtrl not in self.alphas:
            self.alphas[mtrl] = 0.0
        self.rIdxBox.setValue(self.stratum.custom[mtrl])
        self.alphaBox.setValue(self.alphas[mtrl])

    def input_rIdx(self, value):
        self.stratum.custom[self.mtrlsBox.currentText()] = value

    def input_alpha(self, value):
        self.alphas[self.mtrlsBox.currentText()] = value

    def input_repeats(self, value):
        pass

    def input_facet1(self, idx):
        self.facet1 = facetList[idx]
        self.updateLoss()

    def input_facet2(self, idx):
        self.facet2 = facetList[idx]
        self.updateLoss()

    def input_ridgeL(self, value):
        self.ridgeLength = value
        self.updateLoss()

    def updateLoss(self):
        pass

    def strataTable_select(self):
        pass

    def strataTable_item(self, item):
        pass

    def insert_strata(self):
        """SLOT connected to self.insertButton.clicked()"""
        row = self.strataTable.currentRow()
        if row >= len(self.stratum.materials) or row < 0:
            # Add new lines in the last strata
            row = len(self.stratum.materials) - 1
        elif row == 0:
            row = 1
        # len(self.stratum.materials) is at least 2 for top and substrate
        # s.t. row >= 1 and <= len(self.stratum.materials) - 1
        self.stratum.insert(row)
        self.dirty.emit()

    def delete_strata(self):
        """SLOT connected to self.deleteButton.clicked()"""
        row = self.strataTable.currentRow()
        if row >= len(self.stratum.materials)-1 or row < 1:
            return
        # s.t. len(self.stratum.materials) > 2 and Ls is not empty
        self.stratum.delete(row)
        self.dirty.emit()
        print(self.stratum)

    def strataTable_refresh(self):
        self.stratum.updateIndices()
        self.strataTable.clear()
        self.strataTable.setColumnCount(5)
        self.strataTable.setHorizontalHeaderLabels([
            "Material", "x", "Width",
            "Doping", "n"])
        self.strataTable.setRowCount(len(self.stratum.materials))
        self.strataTable.setVerticalHeaderLabels(
            str(n+1) for n in range(len(self.stratum.materials)))

        for q, mtrl in enumerate(self.stratum.materials):
            # material name
            mtrlName = mtrlComboBox()
            mtrlName.addItems(mtrlList)
            mtrlName.addItems(self.stratum.custom.keys())
            mtrlName.setCurrentText(mtrl)
            mtrlName.currentTextChanged.connect(
                partial(self.strataTable_mtrlChanged, q))
            self.strataTable.setCellWidget(q, 0, mtrlName)

            # Mole Frac
            if mtrl not in Alloy:
                moleFrac = QTableWidgetItem('N/A')
                moleFrac.setFlags(Qt.NoItemFlags)
            else:
                moleFrac = QTableWidgetItem("%.2f" % self.stratum.moleFracs[q])
            # moleFrac.setTextAlignment(Qt.AlignCenter)
            self.strataTable.setItem(q, 1, moleFrac)

            # Thickness
            if q == 0 or q == len(self.stratum.materials)-1:
                thickness = QTableWidgetItem('1.0')
                thickness.setFlags(Qt.NoItemFlags)
            else:
                thickness = QTableWidgetItem("%.2f" % self.stratum.Ls[q-1])
                if mtrl in self.stratum.custom:
                    thickness.setFlags(Qt.NoItemFlags)
            self.strataTable.setItem(q, 2, thickness)

            # doping
            doping = QTableWidgetItem("%.1f" % self.stratum.dopings[q])
            if mtrl not in Dopable:
                doping.setFlags(Qt.NoItemFlags)
            self.strataTable.setItem(q, 3, doping)

            # refractive index
            if q == 0:
                ridx = self.stratum.index0
            elif q == len(self.stratum.materials)-1:
                ridx = self.stratum.indexs
            else:
                ridx = self.stratum.indices[q-1]
            ridx = QTableWidgetItem("%.2f + %.2fi" % (ridx.real, ridx.imag))
            ridx.setFlags(Qt.NoItemFlags)
            self.strataTable.setItem(q, 4, ridx)

        self.strataTable.resizeColumnsToContents()

    def strataTable_mtrlChanged(self, row, selection):
        """SLOT as partial(self.strataTable_mtrlChanged, q) connected to
        mtrlName.currentTextChanged(int)"""
        self.stratum.materials[row] = selection
        self.strataTable.selectRow(row)
        self.dirty.emit()

    def update(self):
        self.stratum.updateIndices()
        self.strataTable_refresh()

    def solve(self):
        self.xs = np.linspace(-1, sum(self.stratum.Ls)+1, 1000)
        self.beta = self.stratum.boundModeTM()
        self.Ey, _, _ = self.stratum.populateMode(self.beta, self.xs)
        self.nx = self.stratum.populateIndices(self.xs)
        return self.Ey, self.nx

    def optimizeStrata(self):
        pass

# vim: ts=4 sw=4 sts=4 expandtab
