#!/usr/bin/env python
# -*- coding:utf-8 -*-

# ===========================================================================
# ErwinJr2 is a simulation program for quantum semiconductor lasers.
# Copyright (C) 2012 Kale J. Franz, PhD
# Copyright (C) 2018 Ming Lyu (CareF)
#
# A portion of this code is Copyright (c) 2011, California Institute of
# Technology ("Caltech"). U.S. Government sponsorship acknowledged.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
# ===========================================================================

# TODO:
# In plot controls, add "show one period"
# save and load pickle for qclayers
# Reverse layers
# Bug: when delete all layers
# export excel file for growth sheet

import sys
import traceback
import numpy as np
from numpy import pi, sqrt
from functools import partial, wraps

from QCLayers import QCLayers, qcMaterial, h, c0, e0
from EJcanvas import EJcanvas, EJplotControl

from PyQt5.QtCore import pyqtSignal, pyqtSlot, Qt
from PyQt5.QtGui import QPalette, QColor
from PyQt5.QtWidgets import (QWidget, QVBoxLayout, QLabel, QComboBox,
                             QSpinBox, QDoubleSpinBox, QGroupBox,
                             QCheckBox, QTextEdit, QSizePolicy,
                             QGridLayout, QHBoxLayout, QPushButton,
                             QTableWidget, QTableWidgetItem,
                             QApplication, QMessageBox)

from Material import AParm
ejError = "ErwinJr2 - Error"
ejWarning = "ErwinJr2 - Warning"

class QuantumTab(QWidget):
    """The Quantum Tab of ErwinJr. This is designed to be a GUI wrapper of
    the class QCLayers
    Member viable (! label public interface):
        !qclayers: A class to describe the physics of QC layers
        mtrlList: for Material column in the layerTable
        colors: colors for different wavefunctions

        customized SIGNAL
            dirty: Show if there's any new changes to qclayers that need to
                   be saved

        --- state selection ---
        stateHolder: the list containing the index of seleted states. the
                    index is defined in qclayers
        pairSelected: Boolean flag for if a pair of states is selected
        --- plot flags ---
        solveType: 'basis' or 'whole', decide different kinds of solving
        !plotVX, !plotVL, !plotLH, !plotSO: show X point of conduction band 
                                            (L point, LH band, SO band)

        --- GUI widget ---
        1st column (settingBox):
            inputSubstrateBox
            inputEFieldBox
            inputxResBox
            inputEresBox
            inputRepeatsBox
            inputARInjectorCheck
            inputInjectorARCheck
            LpFirstSpinbox          LpLastSpinbox

        2nd column (layerBox):
            insertLayerAboveButton  deleteLayerButton
            layerTable

        3rd column (solveBox):
            solveBasisButton
            solveWholeButton
            descBox
            mtrlTable
            offsetLabel
            netStrainLabel
            LOPhononLabel
            wellSelectButton
            zoominButton            zoomOutButton
            panButton               clearWFsButton
            pairSelectButton
            FoMButton
            pairSelectString

        4th column (figureBox):
            !quantumCanvas**

    Member method (! label public interface):
        !reload
        !rotate_layer
        !view_VXBand(VLBand, LHBand, SOBand)
        !export_quantumCanvas
        !export_band_data
        !set_temperature
        (private is omitted)
        """
    dirty = pyqtSignal()

    def __init__(self, qclayers=None, parent=None):
        super(QuantumTab, self).__init__(parent)
        self.qclayers = qclayers if qclayers else QCLayers()

        # colors for different wavefunctions
        self.colors = ((0.584, 0.450, 0.701), (0.431, 0.486, 0.745),
                       (0.576, 0.694, 0.517), (0.682, 0.780, 0.321),
                       (0.501, 0.501, 0.509), (0.854, 0.741, 0.247),
                       (0.874, 0.607, 0.290), (0.823, 0.341, 0.278),
                       (0.725, 0.321, 0.623), (0.411, 0.741, 0.270),
                       (0.078, 0.078, 0.078), (0.431, 0.803, 0.870),
                       (0.223, 0.321, 0.643))
        # colors for different materials in tables
        self.mtrlcolors = tuple([QColor(*rgb) for rgb in (
            (255, 255, 255), (230, 230, 240), (230, 240, 230), 
            (240, 230, 230), (230, 240, 240), (240, 230, 240), 
            (240, 240, 230), (230, 230, 230))])

        self.updating = False
        self.solveType = None
        self.plotVX = False
        self.plotVL = False
        self.plotLH = False
        self.plotSO = False

        self.stateHolder = []
        self.pairSelected = False
        # plotType can be mode, wf or DoS (TODO)
        #  self.plotType = "wf"
        self.plotType = "mode"
        self.fillplot = 0.3  # alpha of fill; False for not fill
        self.fillplot = False

        # Platform dependent settings, eg. layerout size settings
        if sys.platform.startswith('win'):
            settingBoxWidth = 150
            layerBoxWidth = 400
            solveBoxWidth = 220
        elif sys.platform.startswith('darwin'):
            settingBoxWidth = 90
            layerBoxWidth = 250
            solveBoxWidth = 190
        elif sys.platform.startswith('linux'):
            settingBoxWidth = 150
            layerBoxWidth = 365
            solveBoxWidth = 240
        else:
            QMessageBox.warning(self, ejWarning, 
                                'Platform %s not tested.' % sys.platform)
            settingBoxWidth = 150
            layerBoxWidth = 400
            solveBoxWidth = 190


        quantumLayout = QHBoxLayout()
        settingBox = self._generateSettingBox(settingBoxWidth)
        layerBox = self._generateLayerBox(layerBoxWidth)
        figureBox, plotControlGrid = self._generateFigureBox()
        solveBox = self._generateSolveBox(plotControlGrid, solveBoxWidth)
        quantumLayout.addLayout(settingBox)
        quantumLayout.addLayout(layerBox)
        quantumLayout.addLayout(solveBox)
        quantumLayout.addLayout(figureBox)
        self.setLayout(quantumLayout)
        self.setAutoFillBackground(True)
        self.setBackgroundRole(QPalette.Window)

        self.reload()
    # __init__ end

    def _generateSettingBox(self, width):
        """ Return a Qt Layout object containning all setting parameters """
        settingBox = QVBoxLayout()

        settingBox.addWidget(QLabel(
            "<center><b>Substrate</b></center>"))
        self.inputSubstrateBox = QComboBox()
        self.inputSubstrateBox.addItems(qcMaterial.keys())
        self.inputSubstrateBox.currentIndexChanged[str].connect(
            self.input_substrate)
        settingBox.addWidget(self.inputSubstrateBox)

        settingBox.addWidget(QLabel(
            '<center><b><i>E<sub>field</sub></i></b></center>'))
        self.inputEFieldBox = QDoubleSpinBox()
        self.inputEFieldBox.setDecimals(1)
        self.inputEFieldBox.setSuffix(' kV/cm')
        self.inputEFieldBox.setRange(-250.0, 250.0)
        self.inputEFieldBox.valueChanged[float].connect(self.input_EField)
        settingBox.addWidget(self.inputEFieldBox)

        settingBox.addWidget(QLabel(
            '<center><b>Position<br>Resolution</b></center>'))
        self.inputxResBox = QComboBox()
        self.xresTuple = (1.0, 0.5, 0.25, 0.1, 0.05)
        self.inputxResBox.addItems([u'%.2f \u212B'%res for res in
                                       self.xresTuple])
        self.inputxResBox.currentIndexChanged[int].connect(self.input_xres)
        settingBox.addWidget(self.inputxResBox)

        settingBox.addWidget(QLabel(
            '<center><b>Energy<br>Resolution</b></center>'))
        self.inputEresBox = QDoubleSpinBox()
        self.inputEresBox.setDecimals(2)
        self.inputEresBox.setRange(0.0, 10.0)
        self.inputEresBox.setSingleStep(0.1)
        self.inputEresBox.setSuffix(' meV')
        self.inputEresBox.valueChanged[float].connect(self.input_Eres)
        settingBox.addWidget(self.inputEresBox)

        settingBox.addWidget(QLabel(
            '<center><b>Structure Repeats</b></center>'))
        self.inputRepeatsBox = QSpinBox()
        self.inputRepeatsBox.setRange(1, 5)
        self.inputRepeatsBox.valueChanged[int].connect(self.input_repeats)
        settingBox.addWidget(self.inputRepeatsBox)

        # Basis solver devider setting
        basisGroupBox = QGroupBox("Basis Divisions")
        self.inputARInjectorCheck = QCheckBox("AR->Injector")
        self.inputInjectorARCheck = QCheckBox("Injector->AR")
        self.inputARInjectorCheck.setChecked(True)
        self.inputInjectorARCheck.setChecked(True)
        basisLayout = QVBoxLayout()
        basisLayout.addWidget(self.inputARInjectorCheck)
        basisLayout.addWidget(self.inputInjectorARCheck)
        basisGroupBox.setLayout(basisLayout)
        self.inputARInjectorCheck.stateChanged.connect(self.input_basis)
        self.inputInjectorARCheck.stateChanged.connect(self.input_basis)
        settingBox.addWidget(basisGroupBox)

        # Period information groupbox
        LpLayoutGroupBox = QGroupBox("Period Info")
        self.LpFirstSpinbox = QSpinBox()
        self.LpFirstSpinbox.setValue(1)
        self.LpFirstSpinbox.setRange(1, 1)
        self.LpFirstSpinbox.valueChanged.connect(self.update_Lp_box)
        self.LpLastSpinbox = QSpinBox()
        self.LpLastSpinbox.setRange(1, 1)
        self.LpLastSpinbox.valueChanged.connect(self.update_Lp_box)
        self.LpStringBox = QTextEdit('')
        self.LpStringBox.setReadOnly(True)
        self.LpStringBox.setSizePolicy(
            QSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed))
        self.LpStringBox.setMaximumHeight(95)
        self.LpStringBox.setMaximumWidth(width)
        LpLayout = QGridLayout()
        LpLayout.addWidget(QLabel('<b>first</b>'), 0, 0)
        LpLayout.addWidget(QLabel('<b>last</b>'), 0, 1)
        LpLayout.addWidget(self.LpFirstSpinbox, 1, 0)
        LpLayout.addWidget(self.LpLastSpinbox, 1, 1)
        LpLayout.addWidget(self.LpStringBox, 2, 0, 1, 2)
        LpLayoutGroupBox.setLayout(LpLayout)
        settingBox.addWidget(LpLayoutGroupBox)

        settingBox.addStretch()
        return settingBox
        # _generateSettingBox end

    def _generateLayerBox(self, width):
        """ Return a Qt Layout object containning all layer parameters """
        layerBox = QGridLayout()
        self.insertLayerAboveButton = QPushButton("Insert Layer Above")
        self.insertLayerAboveButton.clicked.connect(self.insert_layerAbove)
        layerBox.addWidget(self.insertLayerAboveButton, 0, 0)
        self.deleteLayerButton = QPushButton("Delete Layer")
        self.deleteLayerButton.clicked.connect(self.delete_layer)
        layerBox.addWidget(self.deleteLayerButton, 0, 1)

        # set up layerTable
        self.layerTable = QTableWidget()
        self.layerTable.setSelectionBehavior(QTableWidget.SelectRows)
        self.layerTable.setSelectionMode(QTableWidget.SingleSelection)
        self.layerTable.setMaximumWidth(width)
        self.layerTable.setSizePolicy(
            QSizePolicy(QSizePolicy.Fixed, QSizePolicy.Expanding))
        self.layerTable.itemChanged.connect(self.layerTable_itemChanged)
        self.layerTable.itemSelectionChanged.connect(
            self.layerTable_itemSelectionChanged)
        layerBox.addWidget(self.layerTable, 1, 0, 1, 2)

        return layerBox
        # _generateLayerBox end

    def _generateFigureBox(self):
        """Return:
            A Qt Layout containning a canvas for plotting; 
            A QGiridLayout containning plot control;"""
        self.quantumCanvas = EJcanvas(xlabel=u'Position (Å)',
                                      ylabel='Energy (eV)', parent=self)
        self.plotControl = EJplotControl(self.quantumCanvas, self)
        figureBox = QVBoxLayout()
        figureBox.addWidget(self.quantumCanvas)

        self.zoominButton = QPushButton("Zoom")
        self.plotControl.set_action('zoom', self.zoominButton)
        self.zoomOutButton = QPushButton("Reset")
        self.plotControl.set_action('home', self.zoomOutButton)
        self.panButton = QPushButton("Pan")  # to move
        self.plotControl.set_action('pan', self.panButton)
        self.wellSelectButton = QPushButton("Layer Select")
        self.plotControl.set_custom('wellselect', self.wellSelectButton,
                                    self.well_select)
        self.wellSelectButton.clicked.connect(
            partial(self.plotControl.custom, 'wellselect'))
        self.clearWFsButton = QPushButton("Clear")
        self.clearWFsButton.clicked.connect(self.clear_WFs)
        plotControlGrid = QGridLayout()
        plotControlGrid.addWidget(self.wellSelectButton, 0, 0, 1, 2)
        plotControlGrid.addWidget(self.zoominButton, 1, 0, 1, 1)
        plotControlGrid.addWidget(self.zoomOutButton, 1, 1, 1, 1)
        plotControlGrid.addWidget(self.panButton, 2, 0, 1, 1)
        plotControlGrid.addWidget(self.clearWFsButton, 2, 1, 1, 1)

        return figureBox, plotControlGrid
        # _generateFigureBox end

    def _generateSolveBox(self, plotControlGrid, width):
        """ Return a Qt Layout containning material information,
        eigensolve control, states properties calculation and plot control"""
        solveBox = QVBoxLayout()
        self.solveBasisButton = QPushButton("Solve Basis")
        self.solveBasisButton.clicked.connect(self.solve_basis)
        solveBox.addWidget(self.solveBasisButton)
        self.solveWholeButton = QPushButton("Solve Whole")
        self.solveWholeButton.clicked.connect(self.solve_whole)
        solveBox.addWidget(self.solveWholeButton)

        # set up description box
        self.descBox = QTextEdit('')
        self.descBox.setReadOnly(False)
        self.descBox.setSizePolicy(QSizePolicy(
            QSizePolicy.Fixed, QSizePolicy.Minimum))
        self.descBox.setMinimumHeight(60)
        self.descBox.setMaximumWidth(width)
        self.descBox.textChanged.connect(self.input_description)
        descLayout = QVBoxLayout()
        descLayout.addWidget(self.descBox)
        descLayoutGroupBox = QGroupBox("Description")
        descLayoutGroupBox.setLayout(descLayout)
        solveBox.addWidget(descLayoutGroupBox)

        # set up material composition inputs
        self.mtrlTable = QTableWidget()
        #  self.mtrlTable.setSelectionBehavior(QTableWidget.SelectItems)
        self.mtrlTable.setSelectionMode(QTableWidget.NoSelection)
        self.mtrlTable.setMaximumWidth(width)
        self.mtrlTable.setMinimumWidth(100)
        self.mtrlTable.setSizePolicy(
            QSizePolicy(QSizePolicy.Fixed, QSizePolicy.Minimum))
        self.mtrlTable.itemChanged.connect(self.mtrlTable_itemChanged)
        self.offsetLabel = QLabel('')
        self.netStrainLabel = QLabel('')
        self.LOPhononLabel = QLabel('')
        self.addMtrlButton = QPushButton("Add")
        self.addMtrlButton.clicked.connect(self.add_mtrl)
        self.delMtrlButton = QPushButton("Del")
        self.delMtrlButton.clicked.connect(self.del_mtrl)
        mtrlGrid = QGridLayout()
        mtrlGrid.addWidget(self.addMtrlButton, 0, 0)
        mtrlGrid.addWidget(self.delMtrlButton, 0, 1)
        mtrlGrid.addWidget(self.mtrlTable, 1, 0, 1, 2)
        mtrlGrid.addWidget(self.offsetLabel, 2, 0, 1, 2)
        mtrlGrid.addWidget(self.netStrainLabel, 3, 0, 1, 2)
        mtrlGrid.addWidget(self.LOPhononLabel, 4, 0, 1, 2)
        mtrlGroupBox = QGroupBox("Materials")
        mtrlGroupBox.setLayout(mtrlGrid)
        solveBox.addWidget(mtrlGroupBox)

        # set up plot control inputs
        plotControlGroupBox = QGroupBox("Plot Controls")
        plotControlGroupBox.setLayout(plotControlGrid)
        solveBox.addWidget(plotControlGroupBox)

        # set up Calculate controls
        self.pairSelectButton = QPushButton("Pair Select")
        self.pairSelectButton.setEnabled(False)
        self.plotControl.set_custom('pairselect', self.pairSelectButton,
                                    self.state_pick)
        self.pairSelectButton.clicked.connect(
            partial(self.plotControl.custom, 'pairselect'))
        self.FoMButton = QPushButton("FoM")
        self.FoMButton.setEnabled(False)
        self.FoMButton.clicked.connect(self.updateFoM)
        self.pairSelectString = QTextEdit('')
        self.pairSelectString.setReadOnly(True)
        self.pairSelectString.setMaximumWidth(width)
        self.pairSelectString.setSizePolicy(QSizePolicy(
            QSizePolicy.Fixed, QSizePolicy.Fixed))
        calculateControlGrid = QGridLayout()
        calculateControlGrid.addWidget(self.pairSelectButton, 0, 0, 1, 2)
        calculateControlGrid.addWidget(self.FoMButton, 1, 0, 1, 1)
        calculateControlGrid.addWidget(self.pairSelectString, 2, 0, 1, 2)
        calculateControlGroupBox = QGroupBox("Calculate")
        calculateControlGroupBox.setLayout(calculateControlGrid)
        solveBox.addWidget(calculateControlGroupBox)

        solveBox.addStretch()
        return solveBox
        # _generateSolveBox end

    def settingslot(fn):
        """ A decorator to ask slots to skip reaction when doing massive 
        updating"""
        @wraps(fn)
        def wrapper(self, *args, **kwargs):
            if self.updating:
                return
            else:
                return fn(self, *args, **kwargs)
        return wrapper

    def reload(self):
        """ Reload everything from qclayers, typically when it's changed
        externally or loaded"""
        # TODO
        self.updating = True
        self._update_mtrlList()
        self.update_Lp_limits()
        self.update_Lp_box()
        self.mtrlTable_refresh()
        self.layerTable_refresh()
        self.layerTable.selectRow(1)
        self.qclayers.populate_x()
        self.update_quantumCanvas()
        self.updating = False

    def _update_mtrlList(self):
        """Update self.mtrlList according to self.qclayers material
        information. mtrlList is used for mtrl column in layerTable"""
        self.mtrlList = []
        for n, mtrl in enumerate(self.qclayers.materials):
            name = AParm[mtrl]['name']
            name = name.replace("1-x", str(1-self.qclayers.moleFracs[n]))
            name = name.replace("x", str(self.qclayers.moleFracs[n]))
            name += "#%d"%(n+1) + name
            self.mtrlList.append(name)

# ===========================================================================
# Input Controls
# ===========================================================================
    @pyqtSlot('QString')
    @settingslot
    def input_substrate(self, substrateType):
        """ SLOT connected to inputSubstrateBox.currentIndexChanged(QString)
        update substrate chosen """
        if substrateType in ('InP', 'GaAs', 'GaSb'):
            self.qclayers.set_substrate(substrateType)
        else:
            QMessageBox.information(
                self, ejError,
                '%s substrates have not yet been implemented.'%substrateType)
            self.inputSubstrateBox.setCurrentIndex(
                self.inputSubstrateBox.findText(self.qclayers.substrate))
            return
        self.layerTable_refresh()
        self.qclayers.populate_x()
        self.update_quantumCanvas()
        self.dirty.emit()

    @pyqtSlot(float)
    @settingslot
    def input_EField(self, EField):
        """ SLOT connected to inputEFieldBox.valueChanged(double)
        update external E field in unit kV/cm """
        self.qclayers.EField = EField
        self.qclayers.populate_x()
        self.update_quantumCanvas()
        self.dirty.emit()

    @pyqtSlot(int)
    @settingslot
    def input_xres(self, indx):
        """ SLOT connected to inputxResBox.currentIndexChanged(int)
        update position resolution (xres), in angstrom """
        self.qclayers.xres = self.xresTuple[indx]
        self.qclayers.populate_x()
        self.update_quantumCanvas()
        self.dirty.emit()

    @pyqtSlot(float)
    @settingslot
    def input_Eres(self, ERes):
        """ SLOT connected to inputEresBox.valueChanged
        Update initial energy resolution for eigensolver. Set this too small
        may result in loss of some eigenvalue. """
        self.qclayers.Eres = ERes
        self.dirty.emit()

    @pyqtSlot(int)
    @settingslot
    def input_repeats(self, repeat):
        """SLOT connected to SINGAL self.inputRepeatsBox.valueChanged(int)
        update number of repeats for the whole structure."""
        self.qclayers.repeats = repeat
        self.qclayers.populate_x()
        self.update_quantumCanvas()
        self.dirty.emit()

    @pyqtSlot()
    def input_basis(self, state):
        """SLOT connected to self.inputARInjectorCheck.stateChanged(int) and
        self.inputInjectorARCheck.stateChanged(int)
        update dividers info
        """
        self.qclayers.basisARInjector = self.inputARInjectorCheck.isChecked()
        self.qclayers.basisInjectorAR = self.inputInjectorARCheck.isChecked()

    def update_Lp_limits(self):
        """Update Lp select range in the Period Info box (GUI)
        This should be called whenever number of layers is changed
        """
        self.LpFirstSpinbox.setRange(1, len(self.qclayers.layerWidths)-1)
        self.LpFirstSpinbox.setValue(1)
        self.LpLastSpinbox.setRange(1, len(self.qclayers.layerWidths))
        self.LpLastSpinbox.setValue(len(self.qclayers.layerWidths))

    @pyqtSlot()
    @settingslot
    def update_Lp_box(self):
        """Update Lp box in the Period Info box (GUI):
            Lp:total length
            nD: average doping (cm-3)
            ns: 2D carrier density in 1E11 cm-2
        This needs layerWidths and layerDopings information of self.qclayers
        SLOT connected to LpFirstSpinbox/LpLastSpinBox.valueChanged(int)
        """
        LpFirst = self.LpFirstSpinbox.value()
        LpLast = self.LpLastSpinbox.value() + 1
        # +1 because range is not inclusive of last value
        # total length of the layers (1 period)
        Lp = sum(self.qclayers.layerWidths[LpFirst:LpLast])
        Lp_string = u"Lp: %.1f \u212B<br>" % Lp
        # average doping of the layers
        if Lp == 0:
            Lp_string += (u"n<sub>D</sub>: NA\u00D710<sup>17</sup>"
                          u"cm<sup>-3</sup><br>")
        else:
            nD = sum(self.qclayers.layerDopings[LpFirst:LpLast] * 
                     self.qclayers.layerWidths[LpFirst:LpLast]) / Lp
            Lp_string += (u"n<sub>D</sub>: %6.3f\u00D710<sup>17</sup>"
                          u"cm<sup>-3</sup><br>") % nD
        # 2D carrier density in 1E11cm-2
        ns = sum(self.qclayers.layerDopings[LpFirst:LpLast] *
                 self.qclayers.layerWidths[LpFirst:LpLast]) * 1e-2
        Lp_string += (u"n<sub>s</sub>: %6.3f\u00D710<sup>11</sup>"
                      u"cm<sup>-2</sup") % ns
        self.LpStringBox.setText(Lp_string)

    @pyqtSlot()
    @settingslot
    def input_description(self):
        """ SLOT connected to self.descBox.textChanged()
        Change description string for the design"""
        self.qclayers.description = self.descBox.toPlainText()
        self.dirty.emit()

    @settingslot
    def set_temperature(self, T):
        self.qclayers.set_temperature(T)
        self.qclayers.populate_x()
        self.update_quantumCanvas()

#===========================================================================
# Layer Table Control
#===========================================================================
    def layerTable_refresh(self):
        """Refresh layer table, called every time after data update"""
        # Block itemChanged SIGNAL while refreshing
        self.layerTable.blockSignals(True)
        self.layerTable.clear()
        self.layerTable.setColumnCount(5)
        # An extra blanck line for adding new layers
        self.layerTable.setRowCount(len(self.qclayers.layerWidths) + 1)
        self.layerTable.setHorizontalHeaderLabels(
            ('Width', 'ML', 'Mtrl', 'AR', 'Doping'))
        vertLabels = [str(n) for n in range(len(self.qclayers.layerWidths))]
        self.layerTable.setVerticalHeaderLabels(vertLabels)

        gray2 = QColor(230, 230, 230)  # for unchangable background

        for q, layerWidths in enumerate(self.qclayers.layerWidths):
            color = self.mtrlcolors[
                self.qclayers.layerMaterialIdxs[q]%len(self.mtrlcolors)]
            # Width Setup
            width = QTableWidgetItem("%5.1f" % layerWidths)
            width.setTextAlignment(Qt.AlignCenter)
            width.setBackground(color)
            self.layerTable.setItem(q, 0, width)

            # "ML" number of monolayer Setup
            mlThickness = self.qclayers.mtrlAlloys[
                self.qclayers.layerMaterialIdxs[q]].a_perp
            numML = QTableWidgetItem("%5.1f" % (layerWidths / mlThickness))
            numML.setTextAlignment(Qt.AlignCenter)
            numML.setBackground(color)
            self.layerTable.setItem(q, 1, numML)

            # Material Setup
            mtrlWidget = QComboBox()
            mtrlWidget.addItems(self.mtrlList)
            mtrlWidget.setCurrentIndex(self.qclayers.layerMaterialIdxs[q])
            mtrlWidget.currentIndexChanged.connect(
                partial(self.layerTable_materialChanged, q))
            self.layerTable.setCellWidget(q, 2, mtrlWidget)

            # Active Region Layer Setup
            ARitem = QTableWidgetItem()
            ARitem.setCheckState(
                Qt.Checked if self.qclayers.layerARs[q] == 1
                else Qt.Unchecked)
            ARitem.setBackground(color)
            self.layerTable.setItem(q, 3, ARitem)

            # Layer Doping Setup
            doping = QTableWidgetItem(str(self.qclayers.layerDopings[q]))
            doping.setTextAlignment(Qt.AlignCenter)
            doping.setBackground(color)
            self.layerTable.setItem(q, 4, doping)

        self.layerTable.resizeColumnsToContents()
        self.layerTable.blockSignals(False)

    @pyqtSlot()
    def insert_layerAbove(self):
        """SLOT connected to self.insertLayerAboveButton.clicked()"""
        row = self.layerTable.currentRow()
        N = len(self.qclayers.materials)
        if row == -1:
            return
        elif row >= N: 
            # Add new lines in the last layer
            row = N
            AR = self.qclayers.layerARs[row-1] and self.qclayers.layerARs[0]
            doping = self.qclayers.layerDopings[row-1]
            mtrlIdx = (self.qclayers.layerMaterialIdxs[row-1] + 1)%N
        else:
            AR = self.qclayers.layerARs[row] and self.qclayers.layerARs[row-1]
            doping = self.qclayers.layerDopings[row]
            mtrlIdx = (self.qclayers.layerMaterialIdxs[row-1] - 1)%N

        self.qclayers.add_layer(row, 0.0, mtrlIndx, AR, doping)
        self.update_Lp_limits()
        self.update_Lp_box()
        self.layerTable_refresh()
        self.layerTable.selectRow(row)
        self.dirty.emit()

    @pyqtSlot()
    def delete_layer(self):
        """SLOT connected to self.deleteLayerButton.clicked()"""
        row = self.layerTable.currentRow()
        if row == -1 or row >= len(self.qclayers.layerWidths):
            return
        # don't delete last layer
        if len(self.qclayers.layerWidths) == 1:
            self.qclayers.layerWidths[0] = 0.0
            return

        self.qclayers.del_layer(row)
        self.update_Lp_limits()
        self.update_Lp_box()
        self.layerTable_refresh()
        self.layerTable.selectRow(row)
        self.dirty.emit()

    @pyqtSlot(QTableWidgetItem)
    def layerTable_itemChanged(self, item):
        """SLOT connected to layerTable.itemChanged(QTableWidgetItem*)
        Update layer profile after user input"""
        row = item.row()
        if row > len(self.qclayers.layerWidths):
            raise ValueError("Bad layer width input row number")
        column = item.column()
        if column in (0, 1, 4):
            try: 
                value = float(item.text())
            except ValueError:
                # invalid input
                QMessageBox.warning(self, "ErwinJr - Error",
                                    "This value should be a number")
                self.layerTable_refresh()
                return
        if column in (0, 1):
            if column == 0:  # column == 0 for Widths column
                new_width = value
            else:  # column == 1 for "ML" number of monolayers
                new_width = value * self.qclayers.mtrlAlloys[ 
                    self.layerMaterialIdxs[row]].a_perp

            if row == len(self.qclayers.layerWidths):
                # add row at end of list
                AR = self.qclayers.layerARs[row-1]
                doping = self.qclayers.layerDopings[row-1]
                mtrlIdx = (self.qclayers.layerMaterialIdxs[row] + 1)%N
                self.qclayers.add_layer(row, new_width, mtrlIndx, AR, doping)
                row += 1  # used so that last (blank) row is again selected
                self.update_Lp_limits()
                self.update_Lp_box()
            else:  # change Width of selected row in-place
                self.qclayers.layerWidths[row] = new_width
                self.update_Lp_box()

        elif column == 2:
            # column == 2 for item change in mtrl column, should be 
            # controled by layerTable_materialChanged
            raise Exception("Should not be here")

        elif column == 3:
            # column == 3 for item change in AR column
            if row == len(self.qclayers.layerWidths):
                # don't do anything if row is last row
                return
            self.qclayers.layerARs[row] = (item.checkState() == Qt.Checked)

        elif column == 4:
            # column == 4 for item change in Doping column
            doping = value
            if row == len(self.qclayers.layerWidths):
                AR = self.qclayers.layerARs[row-1]
                mtrlIdx = (self.qclayers.layerMaterialIdxs[row] + 1)%N
                self.qclayers.add_layer(row, 0.0, mtrlIndx, AR, doping)
                row += 1  # used so that last (blank) row is again selected
                self.update_Lp_limits()
                self.update_Lp_box()
            else: 
                self.qclayers.layerDopings[row] = doping

        else:
            raise Exception("Should not be here")

        self.layerTable_refresh() 
        self.layerTable.setCurrentCell(row, column)
        self.update_quantumCanvas()
        self.dirty.emit()

    @pyqtSlot()
    @settingslot
    def layerTable_itemSelectionChanged(self):
        """SLOT connected to layerTable.itemSelectionChanged()"""
        # This is the primary call to update_quantumCanvas
        self.qclayers.layerSelected = self.layerTable.currentRow()
        self.qclayers.populate_select()
        self.update_quantumCanvas()

    def layerTable_materialChanged(self, row, selection):
        """SLOT as partial(self.layerTable_materialChanged, q)) connected to
        materialWidget.currentIndexChanged(int) """
        self.qclayers.layerMaterialIdxs[row] = selection
        self.qclayers.populate_x()
        self.layerTable.selectRow(row)
        self.dirty.emit()

    @pyqtSlot()
    def rotate_layer(self):
        """Move last layer to first layer, SLOT open to public"""
        self.qclayers.rotate_layer()
        self.layerTable_refresh()
        self.layerTable.setCurrentCell(1, 0)
        self.update_quantumCanvas()
        self.dirty.emit()

    @pyqtSlot()
    def copy_structure(self):
        clipboard = QApplication.clipboard()
        string = ''
        for width in self.qclayers.layerWidth:
            string += '%.1f\n' % width
        clipboard.setText(string)

#=========================================================================
# mtrlTable Control
#=========================================================================
    def mtrlTable_refresh(self):
        """Set up material table, both format and consistency with qclayers"""
        self.mtrlTable.blockSignals(True)
        self.mtrlTable.clear()
        self.mtrlTable.setColumnCount(3)
        self.mtrlTable.setRowCount(len(self.qclayers.materials))
        self.mtrlTable.setHorizontalHeaderLabels(["#", "mtrl", "x"])
        # TODO: material name support

        possibleMtrl = tuple([AParm[m]['name'] for m in
                              qcMaterial[self.qclayers.substrate]])
        for n, mtrl in enumerate(self.qclayers.materials): 
            color = self.mtrlcolors[n%len(self.mtrlcolors)]
            name = QTableWidgetItem(str(n+1))
            name.setTextAlignment(Qt.AlignCenter)
            name.setBackground(color)
            name.setFlags(Qt.NoItemFlags)
            self.mtrlTable.setItem(n, 0, name)

            # Choose from available materials, according to substrate 
            mtrlItem = QComboBox()
            mtrlItem.addItems(possibleMtrl)
            mtrlItem.setCurrentText(AParm[mtrl]['name'])
            mtrlItem.currentIndexChanged.connect(
                partial(self.mtrlTable_mtrlChanged, n))
            self.mtrlTable.setCellWidget(n, 1, mtrlItem)

            # Set mole fraction for materials
            # TODO: number only with QItemDelegate
            moleFrac = QTableWidgetItem(str(self.qclayers.moleFracs[n]))
            moleFrac.setTextAlignment(Qt.AlignCenter)
            moleFrac.setBackground(color)
            self.mtrlTable.setItem(n, 2, moleFrac)

        self.mtrlTable.resizeColumnsToContents()
        self.update_mtrl_info()
        self.mtrlTable.blockSignals(False)

    def mtrlTable_mtrlChanged(self, row, selection):
        """SLOT as partial(self.mtrlTable_mrtlChanged, q)) connected to
        mtrlItem.currentIndexChanged(int) """
        self.qclayers.materials[row] = qcMaterial[
            self.qclayers.substrate][selection]
        self.update_mtrl_info()
        self.dirty.emit()

    @pyqtSlot(QTableWidgetItem)
    def mtrlTable_itemChanged(self, item):
        """SLOT connected to mtrlTable.itemChanged(QTableWidgetItem*)
        Update material definition after user input"""
        column = item.column()
        row = item.row()
        if column == 0: 
            # change "#" or name column, TODO
            return
        elif column == 1:
            # change "mtrl" or material column, should be handled by
            # mtrlTable_mtrlChanged
            return
        elif column == 3:
            # change "x" or mole fraction
            try:
                mf = float(item.text())
                assert(mf <= 1)
                self.qclayers.set_mtrl(row, moleFrac=int(100*mf)/100)
            except ValueError, AssertionError:
                # Input is not a number or is larger than 1
                QMessageBox.warning(self, ejError, "invalid mole Fraction")
            item.setText(str(self.qclayers.moleFracs[row]))
            self.update_mtrl_info()
        else:
            # Should never be here
            raise ValueError

    @pyqtSlot()
    def add_mtrl(self): 
        """SLOT connected to self.addMtrlButton.clicked()"""
        pass

    @pyqtSlot()
    def del_mtrl(self): 
        """SLOT connected to self.delMtrlButton.clicked()"""
        pass

    def update_mtrl_info(self):
        """Update labels below mtrlTable"""
        # TODO: why averaging?
        self.offsetLabel.setText(
            '<center>ΔE<sub>c</sub>: <b>%6.0f meV </b></center>' % (
            self.qclayers.offset() * 1000))
        self.netStrainLabel.setText(
            "<center>Net Strain: <b>%6.3f%%</b></center>" % 
            self.qclayers.netStrain())
        self.LOPhononLabel.setText(
            "<center>E<sub>LO</sub>: <b>%4.1f meV</b></center>" % 
            self.qclayers.avghwLO())

# =========================================================================
# Quantum Tab Plotting and Plot Control
# =========================================================================
    def update_quantumCanvas(self):
        #  print "update "+inspect.stack()[1][3]
        #  print inspect.stack()[2][3]
        self.quantumCanvas.clear()
        self.quantumCanvas.axes.plot(self.qclayers.xPoints,
                                     self.qclayers.xVc, 'k', linewidth=1)

        # plot Conduction Band L-Valley/X-Valley, Light Hole Valence Bnad and
        # Spin-Orbit coupling Valence Band
        for bandFlag, xv, conf in (
                (self.plotVL, self.qclayers.xVL, 'g--'),
                (self.plotVX, self.qclayers.xVX, 'm-.'),
                (self.plotLH, self.qclayers.xVLH, 'k'),
                (self.plotSO, self.qclayers.xVSO, 'r--')):
            if bandFlag:
                self.quantumCanvas.axes.plot(self.qclayers.xPoints, xv,
                                             conf, linewidth=1)

        # highlight selected layer & make AR layers bold
        self.quantumCanvas.axes.plot(
            self.qclayers.xPoints, 
            self.qclayers.xLayerAR(), 
            'k', linewidth=1.5)
        if self.qclayers.layerSelected:
            self.quantumCanvas.axes.plot(
                self.qclayers.xPoints,
                self.qclayers.xLayerSelected(), 
                'b', linewidth=1.5 if self.qclayers.layerARs[
                    self.qclayers.layerSelected] else 1)

        if hasattr(self.qclayers, 'eigenE'):
            self.curveWF = []
            if self.plotType == "mode":
                y = self.qclayers.psis**2
            elif self.plotType == "wf":
                y = self.qclayers.psis
            for n in range(len(self.qclayers.eigenE)):
                curve, = self.quantumCanvas.axes.plot(
                    self.qclayers.xPoints,
                    y[n, :] + self.qclayers.eigenE[n],
                    color=self.colors[n % len(self.colors)])
                if self.fillplot:
                    self.quantumCanvas.axes.fill_between(
                        self.qclayers.xPoints,
                        y[:, n] + self.qclayers.eigenE[n],
                        self.qclayers.eigenE[n],
                        facecolor=self.colors[n % len(self.colors)],
                        alpha=self.fillplot)
                self.curveWF.append(curve)
            for n in self.stateHolder:
                curve, = self.quantumCanvas.axes.plot(
                    self.qclayers.xPoints,
                    y[:, n] + self.qclayers.eigenE[n],
                    'k', linewidth=2)

        self.quantumCanvas.draw()

    def well_select(self, event):
        """callback registered in plotControl when it's in wellselect mode.
        It's mpl_connect to button_release_event of quantumCanvas """
        if event.button == 1:  # left button clicked
            x = event.xdata
            xLayerNum = np.argmin((self.qclayers.xPoints - x)**2)
            layerNum = self.qclayers.xLayerNums[xLayerNum]
            self.layerTable.selectRow(layerNum)

    def clear_WFs(self):
        if hasattr(self.qclayers, 'eigenE'):
            delattr(self.qclayers, 'eigenE')
        self.pairSelectButton.setEnabled(False)
        self.update_quantumCanvas()


# ===========================================================================
# Export Functions
# ===========================================================================
    def export_quantumCanvas(self, filename=None):
        self.plotControl.save_figure(
            "ErwinJr - Export Band Structure Image", filename, u'png')

    def export_band_data(self, fname):
        np.savetxt(fname.split('.')[0] + '_CB' + '.csv',
                   np.column_stack([self.qclayers.xPoints, self.qclayers.xVc]),
                   delimiter=',')

        if hasattr(self.qclayers, 'psis'):
            # otherwise band structure hasn't been solved yet
            # TODO: make it consistant with plotting
            xyPsiPsiEig = np.zeros(self.qclayers.psis.shape)
            for q in range(len(self.qclayers.eigenE)):
                xyPsiPsiEig[:, q] = (self.qclayers.psis[:, q] +
                                     self.qclayers.eigenE[q])
            np.savetxt(
                fname.split('.')[0] + '_States' + '.csv',
                np.column_stack([self.qclayers.xPoints, xyPsiPsiEig]),
                delimiter=',')


# ===========================================================================
# View Band Items
# ===========================================================================
    def view_Band(self, band):
        # TODO: combine following slot to this one
        self.bandFlags[band] = not self.bandFlags[band]
        self.update_quantumCanvas()

    def view_VXBand(self):
        if self.plotVX:
            self.plotVX = False
        else:
            self.plotVX = True
        self.update_quantumCanvas()

    def view_VLBand(self):
        if self.plotVL:
            self.plotVL = False
        else:
            self.plotVL = True
        self.update_quantumCanvas()

    def view_LHBand(self):
        if self.plotLH:
            self.plotLH = False
        else:
            self.plotLH = True
        self.update_quantumCanvas()

    def view_SOBand(self):
        if self.plotSO:
            self.plotSO = False
        else:
            self.plotSO = True
        self.update_quantumCanvas()


# ===========================================================================
# Calculations
# ===========================================================================
    def Calculating(self, is_doing):
        """UI repaint for doing calculating """
        for button in (self.solveWholeButton, self.solveBasisButton,
                       self.pairSelectButton):
            button.setEnabled(not is_doing)
            button.repaint()
        if self.pairSelected:
            self.FoMButton.setEnabled(not is_doing)
            self.FoMButton.repaint()

    @pyqtSlot()
    def solve_whole(self):  # solves whole structure
        """SLOT connected to solveWholeButton.clicked()
        Whole solver """
        if hasattr(self.qclayers, "eigenE"):
            self.clear_WFs()
        self.pairSelected = False
        self.Calculating(True)

        try:
            self.stateHolder = []
            self.qclayers.solve_whole()
            self.solveType = 'whole'
            self.update_quantumCanvas()
            self.pairSelectButton.setEnabled(True)
        except (IndexError, TypeError):
            QMessageBox.warning(self, 'ErwinJr - Error',
                                traceback.format_exc())

        self.Calculating(False)

    @pyqtSlot()
    def solve_basis(self):  # solves structure with basis
        """SLOT connected to solveBasisButton.clicked()
        Basis solver """
        if hasattr(self.qclayers, "eigenE"):
            self.clear_WFs()
        self.pairSelected = False
        self.Calculating(True)

        try:
            self.stateHolder = []
            self.qclayers.solve_basis()
            self.solveType = 'basis'
            self.update_quantumCanvas()
            self.pairSelectButton.setEnabled(True)
        except (ValueError, IndexError):
            QMessageBox.warning(self, "ErwinJr - Error",
                                traceback.format_exc())

        self.Calculating(False)

    def state_pick(self, event):
        """Callback registered in plotControl when it's in pairselect mode.
        It's mpl_connect to button_release_event of quantumCanvas """
        if not hasattr(self.qclayers, 'eigenE'):
            # Not yet solved
            return

        if event.button == 1:  # left button clicked: select a state
            if len(self.stateHolder) >= 2:
                # start new pair selection
                self.stateHolder = []
                self.pairSelected = False

            # TODO: replace follwoing by matplotlib lines picker
            x = event.xdata
            y = event.ydata
            xData = np.tile(self.qclayers.xPoints,
                            (self.qclayers.psis.shape[1], 1)).T
            if self.plotType in ("mode", "wf"):
                yData = self.qclayers.eigenE + (self.qclayers.psis
                                                if self.plotType == "mode"
                                                else self.qclayers.psis)

            width, height = self.quantumCanvas.axes.bbox.size
            xmin, xmax = self.quantumCanvas.axes.get_xlim()
            ymin, ymax = self.quantumCanvas.axes.get_ylim()
            xScale = (xmax - xmin) / width
            yScale = (ymax - ymin) / height

            r = np.nanmin(sqrt(((xData - x) / xScale)**2 +
                               ((yData - y) / yScale)**2), axis=0)
            ss = np.nanargmin(r)
            self.stateHolder.append(ss)
            #  self.curveWF[ss].set_color('black')
            #  self.curveWF[ss].set_linewidth(2)
        elif event.button == 3:  # right button clicked: remove last selection
            if len(self.stateHolder) == 2:
                self.pairSelected = False
            self.stateHolder.pop()
        else:
            return
        self.update_quantumCanvas()

        if len(self.stateHolder) == 1:
            self.pairString = (u"selected: %d, ..<br>" % self.stateHolder[0])
        elif len(self.stateHolder) == 2:
            self.pairSelected = True
            # TODO: put these enablement to a functions
            self.FoMButton.setEnabled(True)
            E_i = self.qclayers.EigenE[self.stateHolder[0]]
            E_j = self.qclayers.EigenE[self.stateHolder[1]]
            if E_i > E_j:
                upper = self.stateHolder[0]
                lower = self.stateHolder[1]
            else:
                upper = self.stateHolder[1]
                lower = self.stateHolder[0]

            self.eDiff = 1000 * (E_i - E_j)
            self.wavelength = h * c0 / (e0 * np.abs(E_i - E_j)) * 1e6

            if self.solveType is 'basis':
                couplingEnergy = self.qclayers.coupling_energy(
                    self.dCL, upper, lower)
                self.transitionBroadening = self.qclayers.broadening_energy(
                    upper, lower)
                self.opticalDipole = self.qclayers.dipole(upper, lower)
                self.tauUpperLower = 1 / self.qclayers.lo_transition_rate(
                    upper, lower)
                self.pairString = (
                    u"selected: %d, %d<br>"
                    u"energy diff: <b>%6.1f meV</b> (%6.1f \u00B5m)<br>"
                    u"coupling: %6.1f meV<br>broadening: %6.1f meV<br>"
                    u"dipole: <b>%6.1f \u212B</b>"
                    u"<br>LO scattering: <b>%6.3g ps</b><br>") % (
                        self.stateHolder[0],
                        self.stateHolder[1],
                        self.eDiff, self.wavelength,
                        couplingEnergy,
                        self.transitionBroadening,
                        self.opticalDipole,
                        self.tauUpperLower)

            elif self.solveType is 'whole':
                self.opticalDipole = self.qclayers.dipole(upper, lower)
                self.tauUpperLower = 1 / self.qclayers.lo_transition_rate(
                    upper, lower)
                self.transitionBroadening = 0.1 * self.eDiff  # TODO
                self.pairString = (
                    u"selected: %d, %d<br>"
                    u"energy diff: <b>%6.1f meV</b> (%6.1f \u00B5m)<br>"
                    u"dipole: %6.1f \u212B<br>" u"LO scattering: %6.3g ps<br>"
                ) % (self.stateHolder[0], self.stateHolder[1], self.eDiff,
                     self.wavelength, self.opticalDipole, self.tauUpperLower)
            else:
                self.FoMButton.setEnabled(False)

        self.pairSelectString.clear()
        self.pairSelectString.setText(self.pairString)

    def updateFoM(self):
        """ SLOT connected to FoMButton.clicked()
        Calculate Figure of merit.  """
        if len(self.stateHolder) < 2:
            return
        self.Calculating(True)
        self.FoMButton.setEnabled(False)
        self.FoMButton.repaint()

        upper = self.stateHolder[1]
        lower = self.stateHolder[0]
        if upper < lower:
            upper, lower = lower, upper

        self.tauLower = self.qclayers.lo_life_time(lower)
        self.tauUpper = self.qclayers.lo_life_time(upper)
        self.FoM = self.opticalDipole**2 * self.tauUpper \
            * (1 - self.tauLower / self.tauUpperLower)
        # tauUpperLower is the inverse of transition rate (lifetime)
        self.alphaISB = self.qclayers.alphaISB(upper, lower)

        self.FoMString = (
            u"<i>\u03C4<sub>upper</sub></i> : %6.3f ps<br>"
            u"<i>\u03C4<sub>lower</sub></i> : %6.3f ps"
            u"<br>FoM: <b>%6.0f ps \u212B<sup>2</sup></b>"
            u"<br><i>\u03B1<sub>ISB</sub></i> : %.3f cm<sup>2</sup>") % (
                self.tauUpper, self.tauLower, self.FoM, self.alphaISB)
        self.pairSelectString.setText(self.pairString + self.FoMString)

        self.Calculating(False)
        self.FoMButton.setEnabled(True)

# vim: ts=4 sw=4 sts=4 expandtab
