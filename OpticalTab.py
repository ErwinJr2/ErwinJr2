"""
This file defines the optical tab of ErwinJr, for simulation and optimization
of 1D waveguiding
"""

from functools import partial
import sys
import numpy as np
from numpy import log, pi
from PyQt5.QtCore import pyqtSignal, pyqtSlot, Qt
from PyQt5.QtGui import QPalette
from PyQt5.QtWidgets import (QWidget, QVBoxLayout, QHBoxLayout, QGridLayout,
                             QLabel, QComboBox, QGroupBox, QDoubleSpinBox,
                             QSpinBox,
                             QPushButton, QTableWidget, QTableWidgetItem,
                             QTextEdit, QMessageBox)
from OptStrata import OptStrata, rIdx, Alloy, Dopable
from EJcanvas import EJcanvas
# from EJcanvas import config as canvasConfig
from customQTClass import mtrlComboBox
from versionAndName import ejError
mtrlList = list(rIdx.keys()) + list(Alloy.keys())
facetList = ('cleaved', 'perfect AR', 'perfect HR', 'custom')


class OpticalTab(QWidget):
    """The Optical Tab of ErwinJr. This is designed to be a GUI wrapper of
    the class Stratum
    Member variable
        stratum: A class to describe the physics of 1D optical slab waveguide

        customized SIGNAL
            dirty: Show if there's any new changes to qclayers that need to
                   be saved

        --- GUI widget ---
        1st column (settingBox):
            wlBox
            mtrlsBox
            rIdxRealBox     rIdxImagBox
            periodBox       repeatBox
            fieldBox        gainCoeffBox
            facetBox[0]     facetBox[1]
            refBox[0]       refBox[1]
            ridgeLengthBox  mirrorLoss(QLabel)

        2nd column (strataBox):
            insertButton    deleteButton
            strataTable
            solveButton     OptimizeButton
            infoBox

        3rd column (figureBox):
            optCanvas

    Member method
        setupActive
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
        self.facet = [facetList[0], facetList[0]]
        self.beta = None
        self.select = None
        self.redActive = False

        self.setAutoFillBackground(True)
        self.setBackgroundRole(QPalette.Window)

        if sys.platform.startswith('win'):
            settingBoxWidth = 200
            strataBoxWidth = 350
        elif sys.platform.startswith('darwin'):
            settingBoxWidth = 350
            strataBoxWidth = 370
        elif sys.platform.startswith('linux'):
            settingBoxWidth = 350
            strataBoxWidth = 350
        else:
            settingBoxWidth = 350
            strataBoxWidth = 375

        opticalLayout = QHBoxLayout()
        opticalLayout.setSpacing(0)
        opticalLayout.addLayout(self._settingBox(settingBoxWidth))
        opticalLayout.addLayout(self._strataBox(strataBoxWidth))
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
        settingBox.addWidget(self.wlBox)

        mtrlGroupBox = QGroupBox("Custom Material")
        mtrlGroupBox.setMaximumWidth(width)
        mtrlLayout = QGridLayout()
        mtrlLayout.addWidget(QLabel('<center>Material Name</center>'),
                             0, 0, 1, 2)
        self.mtrlsBox = QComboBox()
        self.mtrlsBox.addItems(self.stratum.cstmIndx.keys())
        mtrlLayout.addWidget(self.mtrlsBox, 1, 0, 1, 2)
        mtrlLayout.addWidget(QLabel(
            '<center>index n<sub>eff</sub></center>'), 2, 0)
        self.rIdxRealBox = QDoubleSpinBox()
        self.rIdxRealBox.setSingleStep(0.1)
        mtrlLayout.addWidget(self.rIdxRealBox, 3, 0)
        mtrlLayout.addWidget(QLabel('<center>passive loss k</center>'),
                             2, 1)
        self.rIdxImagBox = QDoubleSpinBox()
        self.rIdxImagBox.setSingleStep(0.1)
        mtrlLayout.addWidget(self.rIdxImagBox, 3, 1)
        mtrlLayout.addWidget(QLabel('<center>Period Length</center>'),
                             4, 0)
        self.periodBox = QDoubleSpinBox()
        self.periodBox.setSuffix(" Å")
        self.periodBox.setMaximum(19999.99)
        self.repeatBox = QSpinBox()
        mtrlLayout.addWidget(self.periodBox, 5, 0)
        mtrlLayout.addWidget(QLabel('<center>Periods</center>'), 4, 1)
        mtrlLayout.addWidget(self.repeatBox, 5, 1)
        mtrlLayout.addWidget(QLabel('<center>Active property</center>'),
                             6, 0, 1, 2)
        mtrlLayout.addWidget(QLabel('<center>Bias Field</center>'), 7, 0)
        self.fieldBox = QDoubleSpinBox()
        self.fieldBox.setSuffix(' kV/cm')
        self.fieldBox.setEnabled(False)
        self.fieldBox.setDecimals(1)
        self.fieldBox.setMaximum(500)
        mtrlLayout.addWidget(self.fieldBox, 8, 0)
        mtrlLayout.addWidget(QLabel('<center>Gain Coef.</center>'), 7, 1)
        self.gainCoeffBox = QDoubleSpinBox()
        self.gainCoeffBox.setSuffix(' cm/kA')
        self.gainCoeffBox.setEnabled(False)
        self.gainCoeffBox.setDecimals(1)
        mtrlLayout.addWidget(self.gainCoeffBox, 8, 1)
        # TODO
        self.updateMtrl()
        mtrlLayout.setHorizontalSpacing(1)
        mtrlLayout.setVerticalSpacing(3)
        mtrlGroupBox.setLayout(mtrlLayout)
        settingBox.addWidget(mtrlGroupBox)

        self.wlBox.valueChanged[float].connect(self.input_wl)
        self.mtrlsBox.currentIndexChanged[int].connect(self.updateMtrl)
        self.rIdxRealBox.valueChanged[float].connect(self.input_rIdx)
        self.rIdxImagBox.valueChanged[float].connect(self.input_alpha)
        self.periodBox.valueChanged[float].connect(self.input_period)
        self.repeatBox.valueChanged[int].connect(self.input_repeats)
        self.gainCoeffBox.valueChanged[float].connect(self.input_gainCoeff)

        ridgeBox = QGroupBox("Ridge Geometry")
        ridgeBox.setMaximumWidth(width)
        ridgeLayout = QGridLayout()
        self.facetBox = [None]*2
        self.refBox = [None]*2
        for n in (0, 1):
            ridgeLayout.addWidget(QLabel(
                "<center>Fecet%d</center>" % (n+1)), 0, n)
            self.facetBox[n] = QComboBox()
            self.facetBox[n].addItems(facetList)
            self.facetBox[n].setCurrentText(self.facet[n])
            ridgeLayout.addWidget(self.facetBox[n], 1, n)
            self.refBox[n] = QDoubleSpinBox()
            self.refBox[n].setDecimals(1)
            self.refBox[n].setRange(0.0, 100.0)
            self.refBox[n].setSuffix(" %")
            self.refBox[n].setEnabled(self.facet[n] == 'custom')
            self.refBox[n].setValue(self.facetRefct(n))
            ridgeLayout.addWidget(self.refBox[n], 3, n)
            self.facetBox[n].currentIndexChanged[int].connect(
                partial(self.input_facet, n))
            self.refBox[n].valueChanged[float].connect(
                partial(self.input_ref, n))
        ridgeLayout.addWidget(QLabel(
            "<center>Reflectivity</center>"), 2, 0, 1, 2)

        ridgeLayout.addWidget(QLabel("<center>Ridge Length</center>"),
                              4, 0)
        self.ridgeLengthBox = QDoubleSpinBox()
        self.ridgeLengthBox.setDecimals(1)
        self.ridgeLengthBox.setRange(0.0, 20.0)
        self.ridgeLengthBox.setSingleStep(1)
        self.ridgeLengthBox.setValue(self.ridgeLength)
        self.ridgeLengthBox.setSuffix(" mm")
        self.ridgeLengthBox.valueChanged[float].connect(self.input_ridgeL)
        ridgeLayout.addWidget(self.ridgeLengthBox, 5, 0)
        ridgeLayout.addWidget(QLabel("<center>Mirror Loss: </center>"),
                              4, 1)
        self.mirrorLoss = QLabel("<center>0.0 cm<sup>-1</sup></center>")
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
        self.infoBox = QTextEdit()
        self.infoBox.setReadOnly(True)
        self.infoBox.setMinimumWidth(width)
        self.infoBox.setMaximumWidth(width)
        self.infoBox.setMinimumHeight(120)
        strataBox.addWidget(self.infoBox, 4, 0, 1, 2)

        strataBox.setRowStretch(5, 0)
        return strataBox

    @pyqtSlot(float)
    def input_wl(self, wl):
        """SLOT connected to self.wlBox.valueChanged[float]"""
        self.stratum.setWl(wl)
        self.dirty.emit()

    @pyqtSlot(int)
    def updateMtrl(self, idx=0):
        """SLOT connect to mtrlsBox.currentIndexChanged[int]"""
        if self.stratum.cstmIndx:
            self.mtrlsBox.setCurrentIndex(idx)
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
            if mtrl in self.stratum.cstmGain:
                self.gainCoeffBox.setValue(self.stratum.cstmGain[mtrl])
        else:
            self.rIdxRealBox.setValue(1.0)
            self.rIdxRealBox.setEnabled(False)
            self.rIdxImagBox.setValue(0.0)
            self.rIdxImagBox.setEnabled(False)
            self.periodBox.setValue(0.0)
            self.periodBox.setEnabled(False)
            return

    @pyqtSlot(float)
    def input_rIdx(self, value):
        """SLOT connected to self.rIdxRealBox.valueChanged[float]"""
        mtrl = self.mtrlsBox.currentText()
        alpha = self.stratum.cstmIndx[mtrl].imag
        self.stratum.cstmIndx[mtrl] = value + 1j * alpha
        self.dirty.emit()

    @pyqtSlot(float)
    def input_alpha(self, value):
        """SLOT connected to self.rIdxImagBox.valueChanged[float]"""
        mtrl = self.mtrlsBox.currentText()
        neff = self.stratum.cstmIndx[mtrl].real
        self.stratum.cstmIndx[mtrl] = neff + 1j * value
        self.dirty.emit()

    @pyqtSlot(float)
    def input_period(self, value):
        """SLOT connected to self.periodBox.valueChanged[float]"""
        mtrl = self.mtrlsBox.currentText()
        if mtrl not in self.stratum.cstmPrd:
            self.stratum.cstmPrd[mtrl] = [value, 1]
        else:
            self.stratum.cstmPrd[mtrl][0] = value
        if value == 0:
            self.repeatBox.setEnabled(False)
        else:
            self.repeatBox.setEnabled(True)
            self.update_customLength()

    @pyqtSlot(int)
    def input_repeats(self, value):
        """SLOT connected to self.repeatBox.valueChanged[int]"""
        mtrl = self.mtrlsBox.currentText()
        if mtrl not in self.stratum.cstmPrd:
            self.stratum.cstmPrd[mtrl] = [1, value]
        else:
            self.stratum.cstmPrd[mtrl][1] = value
        self.update_customLength()

    @pyqtSlot(float)
    def input_gainCoeff(self, value):
        """SLOT connected to self.gainCoeffBox.valueChanged[float]"""
        mtrl = self.mtrlsBox.currentText()
        self.stratum.cstmGain[mtrl] = value

    def setupActive(self, wl, EField, gaincoeff, Lp):
        """Interface to get parameters from quantum tab"""
        self.wlBox.setValue(wl)
        # set up active core
        self.mtrlsBox.setCurrentText("Active Core")
        self.fieldBox.setValue(EField)
        self.periodBox.setValue(Lp)
        self.gainCoeffBox.setValue(gaincoeff)

    def update_customLength(self):
        length = self.periodBox.value() * self.repeatBox.value() / 1E4  # to um
        for n, mtrl in enumerate(self.stratum.materials):
            if mtrl == self.mtrlsBox.currentText():
                self.stratum.Ls[n] = length
        self.dirty.emit()

    def input_facet(self, n, idx):
        """SLOT as partial(self.input_facet, n) connected to
        facetBox[n].currentIndexChanged(int)"""
        self.facet[n] = facetList[idx]
        self.refBox[n].setValue(self.facetRefct(n)*100)
        self.refBox[n].setEnabled(self.facet[n] == 'custom')
        self.update_Loss()

    def input_ref(self, n, ref):
        """SLOT as partial(self.input_facet, n) connected to
        facetBox[n].currentIndexChanged(int)"""
        self.update_Loss()

    def facetRefct(self, n):
        """Return the reflectivity of facet n"""
        if self.facet[n] == 'cleaved':
            if self.beta is None:
                return -1
            R = abs((self.beta.real - 1)/(self.beta.real + 1))**2
            self.refBox[n].setValue(100 * R)
            return R
        if self.facet[n] == 'perfect AR':
            return 1E-9
        if self.facet[n] == 'perfect HR':
            return 1
        if self.facet[n] == 'custom':
            return self.refBox[n].value()/100

    @pyqtSlot(float)
    def input_ridgeL(self, value):
        """SLOT connected to ridgeLengthBox.valueChanged[float]"""
        self.ridgeLength = value
        self.update_Loss()

    def update_Loss(self):
        """Update the mirror loss, should be called whenever the facet
        settings are changed"""
        # print(self.facetRefct(0), self.facetRefct(1))
        perRunLoss = self.facetRefct(1) * self.facetRefct(0)
        self.alpham = -log(perRunLoss)/(2*self.ridgeLength/10)  # to cm-1
        self.mirrorLoss.setText(
            "<center>%.1f cm<sup>-1</sup></center>" % self.alpham)
        self.dirty.emit()

    @pyqtSlot()
    def strataTable_select(self):
        """SLOT connected to SIGNAL strataTable.itemSelectionChanged"""
        row = self.strataTable.currentRow()
        if row < len(self.stratum.materials):
            self.select = row
        else:
            self.select = None
            self.stratTable.clearSelection()
        self.update_canvas()

    @pyqtSlot()
    def strataTable_item(self, item):
        """SLOT connected to strataTable.itemChanged"""
        row = item.row()
        column = item.column()
        if column in (1, 2, 3):
            try:
                value = float(item.text())
                if value < 0:
                    raise ValueError
            except ValueError:
                # invalid input
                QMessageBox.warning(self, ejError, "This value should be "
                                    "a non-negative number")
                self.strataTable_refresh()
                return
            if column == 1:
                # molefrac
                if value > 1:
                    QMessageBox.warning(
                        self, ejError, "Mole Fraction must be between 0 and 1")
                    self.strataTable.setItem(row, column, QTableWidgetItem(
                        "%.2f" % self.stratum.moleFracs[row]))
                    return
                self.stratum.moleFracs[row] = value
            if column == 2:
                # width
                self.stratum.Ls[row] = value
            if column == 3:
                # doping
                self.stratum.dopings[row] = value

        self.stratum.updateIndices()
        self.dirty.emit()

    @pyqtSlot()
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

    @pyqtSlot()
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
        """Update strataTable content, should be called whenever the strata
        structure is changed"""
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
                thickness = QTableWidgetItem('1.0' if q == 0 else '2.0')
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

    @pyqtSlot()
    def update_xs(self):
        """Update position vector for plotting and confinement factor integral,
        also as SLOT to self.dirty SGINAL"""
        self.beta = None
        self.stratum.updateIndices()
        self.strataTable_refresh()
        self.xs = np.linspace(-1, sum(self.stratum.Ls[1:]), 5000)

    @pyqtSlot()
    def solve(self):
        """SLOT connected to self.solveButton.clicked"""
        try:
            self.beta = self.stratum.boundModeTM()
        except (TimeoutError, ValueError):
            QMessageBox.warning(self, ejError, "Failed to solve for modes")
            return
        self.Ey, _, _ = self.stratum.populateMode(self.beta, self.xs)
        # TODO
        nx = self.stratum.populateIndices(self.xs).real
        self.confinement = 0
        for ar in self.stratum.populateMtrl(self.xs):
            self.confinement += np.trapz(
                (nx*np.abs(self.Ey)**2)[ar], self.xs[ar])
        self.confinement = self.beta.real * self.confinement / np.trapz(
            (nx * np.abs(self.Ey))**2, self.xs)
        self.alphaw = 4*pi/(self.stratum.wl/1E4) * self.beta.imag  # cm^-1
        self.update_canvas()
        self.update_Loss()
        self.update_info()
        return self.Ey

    def update_info(self):
        """Update information in the info box"""
        info = ""
        if self.beta is not None:
            info += "Effective refractive index:\n"
            info += "  β = %.3f + (%.3g)i\n" % (self.beta.real, self.beta.imag)
            info += "Waveguide loss:\n"
            info += "  α<sub>w</sub> = %.3f cm<sup>-1</sup>\n" % self.alphaw
            info += "Confinement factor:\n"
            info += "  Γ = %.1f%%\n" % (self.confinement * 100)
            try:
                self.jth = (self.alpham + self.alphaw)/(
                    self.stratum.cstmGain["Active Core"] * self.confinement)
                info += "Threshold current:\n"
                info += "  J<sub>th</sub> = %.1f kA/cm<sup>-2</sup>" % self.jth
            except (AttributeError, ZeroDivisionError, KeyError):
                info += "\nUse the quantum tab to define Active region"
        self.infoBox.setHtml(info.replace("\n", "<br>"))

    def update_canvas(self):
        """Update figure in optCanvas"""
        self.optCanvas.clear()
        nx = self.stratum.populateIndices(self.xs).real
        self.ridxAxis.plot(self.xs, nx, 'k', lw=1)
        if self.redActive:
            for ar in self.stratum.populateIndices(self.xs):
                self.ridxAxis.plot(self.xs[ar], nx[ar], 'r', lw=2)
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

    @pyqtSlot()
    def optimizeStrata(self):
        # TODO
        QMessageBox.warning(self, ejError,
                            "Optimization is not yet implemented")

# vim: ts=4 sw=4 sts=4 expandtab
