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
                             QTextEdit, QMessageBox)
from OptStrata import OptStrata, rIdx, Alloy, Dopable
from EJcanvas import EJcanvas
from EJcanvas import config as canvasConfig
from customQTClass import mtrlComboBox
from versionAndName import ejError
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
        resultBox

    3rd column (figureBox):
        optCanvas
    """
    dirty = pyqtSignal()

    def __init__(self, stratum=None, parent=None):
        super(OpticalTab, self).__init__(parent)
        self.stratum = stratum if stratum else OptStrata(3.0)
        # TODO: put gain as part of ridx?
        self.periods = {}
        self.ridgeLength = 3.0
        self.ridgeLoss = 0.0
        self.Lps = {}
        self.facet1 = facetList[0]
        self.facet2 = facetList[0]
        opticalLayout = QHBoxLayout()
        opticalLayout.setSpacing(0)
        opticalLayout.addLayout(self._settingBox(350))
        opticalLayout.addLayout(self._strataBox(375))
        figureBox = QVBoxLayout()
        self.optCanvas = EJcanvas(xlabel='Position $x$ (μm)',
                                  ylabel='Refractive index $n$')
        self.ridxAxis = self.optCanvas.axes
        self.modeAxis = self.ridxAxis.twinx()
        self.modeAxis.set_frame_on(False)
        self.modeAxis.get_yaxis().set_visible(False)
        # self.modeAxis = self.ridxAxis
        # self.modeAxis.set_ylabel("TM Mode $E_y$ (a.u.)", color='C0', 
        #                          fontsize=canvasConfig["fontsize"])
        # self.modeAxis.tick_params(axis='y', labelcolor='C0')
        figureBox.addWidget(self.optCanvas)
        opticalLayout.addLayout(figureBox)

        self.beta = None
        self.select = None

        self.setLayout(opticalLayout)
        self.dirty.connect(self.update_xs)
        self.update_xs()
        self.update_canvas()

    def _settingBox(self, width):
        """Return a Qt Layout object containning all setting widgets"""
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
        self.mtrlsBox.addItems(self.stratum.cstmIndx.keys())
        self.mtrlsBox.currentIndexChanged[str].connect(self.set_MtrlBox)
        mtrlLayout.addWidget(self.mtrlsBox, 1, 0, 1, 2)
        mtrlLayout.addWidget(QLabel(
            '<center><b>RI n<sub>eff</sub></b></center>'), 2, 0)
        self.rIdxRealBox = QDoubleSpinBox()
        mtrlLayout.addWidget(self.rIdxRealBox, 3, 0)
        mtrlLayout.addWidget(QLabel('<center><b>gain/loss k</b></center>'),
                             2, 1)
        self.rIdxImagBox = QDoubleSpinBox()
        mtrlLayout.addWidget(self.rIdxImagBox, 3, 1)
        mtrlLayout.addWidget(QLabel('<center><b>Period Length</b></center>'),
                             4, 0)
        self.periodBox = QDoubleSpinBox()
        self.periodBox.setSuffix(" Å")
        self.periodBox.setMaximum(9999.99)
        self.repeatBox = QSpinBox()
        self.repeatBox.setValue(1)
        if self.stratum.cstmIndx:
            self.mtrlsBox.setCurrentIndex(0)
            mtrl = self.mtrlsBox.currentText()
            mtrlIdx = self.stratum.cstmIndx[mtrl]
            self.rIdxRealBox.setValue(mtrlIdx.real)
            self.rIdxImagBox.setValue(mtrlIdx.imag)
            if (mtrl in self.stratum.cstmPrd and
                    self.stratum.cstmPrd[mtrl][0] > 0):
                self.periodBox.setValue(self.stratum.cstmPrd[mtrl][0])
                self.repeatBox.setValue(self.stratum.cstmPrd[mtrl][1])
            else:
                self.periodBox.setValue(0.0)
                self.repeatBox.setValue(0)
                self.repeatBox.setEnabled(False)
        else:
            self.rIdxRealBox.setValue(1.0)
            self.rIdxRealBox.setEnabled(False)
            self.rIdxImagBox.setValue(0.0)
            self.rIdxImagBox.setEnabled(False)
            self.periodBox.setValue(0.0)
            self.periodBox.setEnabled(False)
        self.rIdxRealBox.valueChanged[float].connect(self.input_rIdx)
        self.rIdxImagBox.valueChanged[float].connect(self.input_alpha)
        self.periodBox.valueChanged[float].connect(self.input_period)
        self.repeatBox.valueChanged[int].connect(self.input_repeats)
        mtrlLayout.addWidget(self.periodBox, 5, 0)
        mtrlLayout.addWidget(QLabel('<center><b>Periods</b></center>'), 4, 1)
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
        self.refBox1.valueChanged[float].connect(self.input_ref1)
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
        self.refBox2.valueChanged[float].connect(self.input_ref2)
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
        """Return a Qt Layout object containning all strata table widgets"""
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
        self.strataTable.setMinimumHeight(375)
        self.strataTable.setMaximumHeight(500)
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

    def input_wl(self, wl):
        self.stratum.setWl(wl)
        self.dirty.emit()

    def set_MtrlBox(self, mtrl):
        if mtrl not in self.stratum.cstmIndx:
            self.stratum.cstmIndx[mtrl] = 1.0
        self.rIdxRealBox.setValue(self.stratum.cstmIndx[mtrl].real)
        self.rIdxImagBox.setValue(self.stratum.cstmIndx[mtrl].imag)

    def input_rIdx(self, value):
        self.stratum.cstmIndx[self.mtrlsBox.currentText()] = value

    def input_alpha(self, value):
        pass

    def input_period(self, value):
        mtrl = self.mtrlsBox.currentText()
        if mtrl not in self.stratum.cstmPrd:
            self.stratum.cstmPrd[mtrl] = [value, 1]
        else:
            self.stratum.cstmPrd[self.mtrlsBox.currentText()][0] = value
        if value == 0:
            self.repeatBox.setEnabled(False)
        else:
            self.repeatBox.setEnabled(True)
            self.update_customLength()

    def input_repeats(self, value):
        mtrl = self.mtrlsBox.currentText()
        if mtrl not in self.stratum.cstmPrd:
            self.stratum.cstmPrd[mtrl] = [0, value]
        else:
            self.stratum.cstmPrd[self.mtrlsBox.currentText()][1] = value
        self.update_customLength()

    def update_customLength(self):
        length = self.periodBox.value() * self.repeatBox.value() / 1E4  # to um
        for n, mtrl in enumerate(self.stratum.materials):
            if mtrl == self.mtrlsBox.currentText():
                self.stratum.Ls[n] = length
        self.dirty.emit()

    def input_facet1(self, idx):
        self.facet1 = facetList[idx]
        self.updateLoss()

    def input_facet2(self, idx):
        self.facet2 = facetList[idx]
        self.updateLoss()

    def input_ref1(self, ref):
        pass

    def input_ref2(self, ref):
        pass

    def input_ridgeL(self, value):
        self.ridgeLength = value
        self.updateLoss()

    def updateLoss(self):
        pass

    def strataTable_select(self):
        row = self.strataTable.currentRow()
        if row < len(self.stratum.materials):
            self.select = row
        else:
            self.select = None
            self.stratTable.clearSelection()
        self.update_canvas()

    def strataTable_item(self, item):
        row = item.row()
        column = item.column()
        if column in (1, 2, 3):
            try:
                value = float(item.text())
            except ValueError:
                # invalid input
                QMessageBox.warning(self, ejError,
                                    "This value should be a number")
                self.strataTable_refresh()
                return
            if column == 1:
                # molefrac
                self.stratum.moleFracs[row] = value
            if column == 2:
                # width
                self.stratum.Ls[row] = value
            if column == 3:
                # doping
                self.stratum.dopings[row] = value

        self.stratum.updateIndices()
        self.dirty.emit()

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
        self.strataTable.selectRow(row)

    def delete_strata(self):
        """SLOT connected to self.deleteButton.clicked()"""
        row = self.strataTable.currentRow()
        if row >= len(self.stratum.materials)-1 or row < 1:
            return
        # s.t. len(self.stratum.materials) > 2 and Ls is not empty
        self.stratum.delete(row)
        self.dirty.emit()
        self.strataTable.selectRow(row)

    def strataTable_refresh(self):
        self.strataTable.blockSignals(True)
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
            mtrlName.addItems(self.stratum.cstmIndx.keys())
            mtrlName.setCurrentText(mtrl)
            mtrlName.currentTextChanged.connect(
                partial(self.strataTable_mtrlChanged, q))
            self.strataTable.setCellWidget(q, 0, mtrlName)

            # Mole Frac
            if mtrl not in Alloy:
                moleFrac = QTableWidgetItem('N/A')
                moleFrac.setFlags(Qt.ItemIsSelectable)
            else:
                moleFrac = QTableWidgetItem("%.2f" % self.stratum.moleFracs[q])
            # moleFrac.setTextAlignment(Qt.AlignCenter)
            self.strataTable.setItem(q, 1, moleFrac)

            # Thickness
            if q == 0 or q == len(self.stratum.materials)-1:
                thickness = QTableWidgetItem('1.0' if q==0 else '2.0')
                thickness.setFlags(Qt.ItemIsSelectable)
            else:
                thickness = QTableWidgetItem("%.2f" % self.stratum.Ls[q])
                if mtrl in self.stratum.cstmIndx:
                    thickness.setFlags(Qt.ItemIsSelectable)
            self.strataTable.setItem(q, 2, thickness)

            # doping
            if mtrl in Dopable:
                doping = QTableWidgetItem("%.1f" % self.stratum.dopings[q])
            else:
                doping = QTableWidgetItem("N/A")
                doping.setFlags(Qt.ItemIsSelectable)
            self.strataTable.setItem(q, 3, doping)

            # refractive index
            if q == 0:
                ridx = self.stratum.index0
            elif q == len(self.stratum.materials)-1:
                ridx = self.stratum.indexs
            else:
                ridx = self.stratum.indices[q-1]
            ridx = QTableWidgetItem("%.2f + %.2fi" % (ridx.real, ridx.imag))
            ridx.setFlags(Qt.ItemIsSelectable)
            self.strataTable.setItem(q, 4, ridx)

        self.strataTable.resizeColumnsToContents()
        self.strataTable.blockSignals(False)

    def strataTable_mtrlChanged(self, row, selection):
        """SLOT as partial(self.strataTable_mtrlChanged, q) connected to
        mtrlName.currentTextChanged(int)"""
        self.stratum.materials[row] = selection
        self.strataTable.selectRow(row)
        self.dirty.emit()

    def update_xs(self):
        self.beta = None
        self.stratum.updateIndices()
        self.strataTable_refresh()
        self.xs = np.linspace(-1, sum(self.stratum.Ls[1:]), 1000)

    def solve(self):
        try:
            self.beta = self.stratum.boundModeTM()
        except (TimeoutError, ValueError):
            QMessageBox.warning(self, ejError, "Failed to solve for modes")
            return
        self.Ey, _, _ = self.stratum.populateMode(self.beta, self.xs)
        self.update_canvas()
        self.resultBox.setText("β = %.3f+%.3gi" % (
            self.beta.real, self.beta.imag))
        return self.Ey

    def update_canvas(self):
        self.optCanvas.clear()
        nx = self.stratum.populateIndices(self.xs).real
        self.ridxAxis.plot(self.xs, nx, 'k', lw=1)
        # plot select strata
        if self.select is not None:
            if self.select == 0:
                idx = self.xs < 0
            elif self.select == len(self.stratum.Ls) - 1:
                idx = self.xs > sum(self.stratum.Ls[1:-1])
            else:
                lsum = sum(self.stratum.Ls[1:self.select])
                idx = (self.xs >= lsum) & (
                    self.xs < lsum+self.stratum.Ls[self.select])
            self.ridxAxis.plot(self.xs[idx], nx[idx], 'b', lw=1.5)
        if self.beta is not None:
            self.modeAxis.plot(self.xs, np.abs(self.Ey)**2, color='C0')
        # self.optCanvas.figure.tight_layout()
        self.optCanvas.draw()

    def optimizeStrata(self):
        pass

# vim: ts=4 sw=4 sts=4 expandtab
