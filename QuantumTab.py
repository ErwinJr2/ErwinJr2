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

import sys
import traceback
import numpy as np
from numpy import pi, sqrt
from functools import partial, wraps

from QCLayers import QCLayers, qcMaterial
from QCLayers import h, c0, e0
from EJcanvas import EJcanvas, EJplotControl

from PyQt5.QtCore import pyqtSignal, QString, Qt
from PyQt5.QtGui import QPalette, QColor
from PyQt5.QtWidgets import (QWidget, QVBoxLayout, QLabel, QComboBox,
                             QSpinBox, QDoubleSpinBox, QGroupBox,
                             QCheckBox, QTextEdit, QSizePolicy,
                             QGridLayout, QHBoxLayout, QPushButton,
                             QTableWidget, QTableWidgetItem,
                             QApplication, QMessageBox)

# ===========================================================================
# Debug options
# ===========================================================================
DEBUG = 2
if DEBUG >= 3:
    import pickle
    #  import inspect


class QuantumTab(QWidget):
    """The Quantum Tab of ErwinJr. This is designed to be a GUI wrapper of
    the class QCLayers
    Member viable (! label public interface):
        !qclayers: A class to describe the physics of QC layers
        materialList: for Material column in the layerTable
        colors: colors for different wavefunctions

        customized SIGNAL: 'dirty'

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
            inputHorzResBox
            inputVertResBox
            inputRepeatsBox
            inputARInjectorCheck
            inputInjectorARCheck
            LpFirstSpinbox          LpLastSpinbox

        2nd column (layerBox):
            insertLayerAboveButton  deleteLayerButton
            OptimizeFoMButton       OptimizeDipoleButton
            layerTable

        3rd column (solveBox):
            solveBasisButton
            solveWholeButton
            DescriptionBox
            mtrl_well               mtrl_barr
            MoleFracWellBox[n]      MoleFracBarrBox[n]     offsetLabel[n]
            strainDescription
            LOPhononDescription
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
        !layerNum
        !bump_first_layer
        !copy_structure
        !view_VXBand(VLBand, LHBand, SOBand)
        !export_quantumCanvas
        !export_band_data
        !set_temperature
        (private is omitted)
        """
    #  if pyqt5:
    dirty = pyqtSignal()

    def __init__(self, parent=None):
        super(QuantumTab, self).__init__(parent)
        self.qclayers = QCLayers()

        # colors for different wavefunctions
        self.colors = ((0.584, 0.450, 0.701), (0.431, 0.486, 0.745),
                       (0.576, 0.694, 0.517), (0.682, 0.780, 0.321),
                       (0.501, 0.501, 0.509), (0.854, 0.741, 0.247),
                       (0.874, 0.607, 0.290), (0.823, 0.341, 0.278),
                       (0.725, 0.321, 0.623), (0.411, 0.741, 0.270),
                       (0.078, 0.078, 0.078), (0.431, 0.803, 0.870),
                       (0.223, 0.321, 0.643))

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
        if sys.platform == 'win32':
            layerTableSize = 400
            DescriptionBoxWidth = 220
            LpStringBoxWidth = 150
        elif sys.platform == 'darwin':
            layerTableSize = 250
            DescriptionBoxWidth = 190
            LpStringBoxWidth = 90
        elif sys.platform == 'linux2':
            layerTableSize = 365
            DescriptionBoxWidth = 240
            LpStringBoxWidth = 150
        else:
            QMessageBox.warning(self, 'ErwinJr - Warning',
                                'Platform %s not tested.' % sys.platform)
            layerTableSize = 400
            DescriptionBoxWidth = 190
            LpStringBoxWidth = 150
        pairSelectStringWidth = DescriptionBoxWidth

        self.updating = True
        # ####################################################
        # settingBox, containing all setting parameter
        # ####################################################
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
        self.inputEFieldBox.setRange(0.0, 250.0)
        self.inputEFieldBox.valueChanged.connect(self.input_EField)
        settingBox.addWidget(self.inputEFieldBox)

        settingBox.addWidget(QLabel(
            '<center><b>Position<br>Resolution</b></center>'))
        self.inputHorzResBox = QComboBox()
        self.HorzRes = (1.0, 0.5, 0.25, 0.1, 0.05)
        self.inputHorzResBox.addItems([u'%.2f \u212B'%res for res in
                                       self.HorzRes])
        self.inputHorzResBox.currentIndexChanged.connect(self.input_horzRes)
        settingBox.addWidget(self.inputHorzResBox)

        settingBox.addWidget(QLabel(
            '<center><b>Energy<br>Resolution</b></center>'))
        self.inputVertResBox = QDoubleSpinBox()
        self.inputVertResBox.setDecimals(2)
        self.inputVertResBox.setValue(0.5)
        self.inputVertResBox.setRange(0.0, 10.0)
        self.inputVertResBox.setSingleStep(0.1)
        self.inputVertResBox.setSuffix(' meV')
        self.inputVertResBox.valueChanged.connect(self.input_vertRes)
        # TODO: check the SLOT if it can make use of the input
        settingBox.addWidget(self.inputVertResBox)

        settingBox.addWidget(QLabel(
            '<center><b>Structure Repeats</b></center>'))
        self.inputRepeatsBox = QSpinBox()
        self.inputRepeatsBox.setValue(1)
        self.inputRepeatsBox.setRange(1, 5)
        self.inputRepeatsBox.valueChanged.connect(self.input_repeats)
        # TODO: check the SLOT if it can make use of the input
        settingBox.addWidget(self.inputRepeatsBox)

        # Basis solver devider setting
        basis_groupBox = QGroupBox("Basis Divisions")
        self.inputARInjectorCheck = QCheckBox("AR->Injector")
        self.inputInjectorARCheck = QCheckBox("Injector->AR")
        self.inputARInjectorCheck.setChecked(True)
        self.inputInjectorARCheck.setChecked(True)
        basisLayout = QVBoxLayout()
        basisLayout.addWidget(self.inputARInjectorCheck)
        basisLayout.addWidget(self.inputInjectorARCheck)
        basis_groupBox.setLayout(basisLayout)
        self.inputARInjectorCheck.stateChanged.connect(self.input_basis)
        self.inputInjectorARCheck.stateChanged.connect(self.input_basis)
        settingBox.addWidget(basis_groupBox)

        # Period information groupbox
        LpLayout_groupBox = QGroupBox("Period Info")
        self.LpFirstSpinbox = QSpinBox()
        self.LpFirstSpinbox.setValue(1)
        self.LpFirstSpinbox.setRange(1, 1)
        self.LpFirstSpinbox.valueChanged.connect(self.update_inputBoxes)
        # TODO: change this slot to prevent unnecessary data updates
        self.LpLastSpinbox = QSpinBox()
        self.LpLastSpinbox.setValue(1)
        self.LpLastSpinbox.setRange(1, 1)
        self.LpLastSpinbox.valueChanged.connect(self.update_inputBoxes)
        # TODO: change this slot to prevent unnecessary data updates
        self.LpStringBox = QTextEdit('')
        self.LpStringBox.setReadOnly(True)
        self.LpStringBox.setSizePolicy(
            QSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed))
        self.LpStringBox.setMaximumHeight(95)
        self.LpStringBox.setMaximumWidth(LpStringBoxWidth)
        LpLayout = QGridLayout()
        LpLayout.addWidget(QLabel('<b>first</b>'), 0, 0)
        LpLayout.addWidget(QLabel('<b>last</b>'), 0, 1)
        LpLayout.addWidget(self.LpFirstSpinbox, 1, 0)
        LpLayout.addWidget(self.LpLastSpinbox, 1, 1)
        LpLayout.addWidget(self.LpStringBox, 2, 0, 1, 2)
        LpLayout_groupBox.setLayout(LpLayout)
        settingBox.addWidget(LpLayout_groupBox)

        settingBox.addStretch()
        # vBox1: settingBox end

        # ######################################################
        # layerBox, containing layer table and related buttons
        # ######################################################
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
        self.layerTable.setMaximumWidth(layerTableSize)
        self.layerTable.setMinimumWidth(layerTableSize)
        self.layerTable.itemChanged.connect(self.layerTable_itemChanged)
        self.layerTable.itemSelectionChanged.connect(
            self.layerTable_itemSelectionChanged)
        layerBox.addWidget(self.layerTable, 1, 0, 1, 2)
        # vBox2: layerBox end

        # ######################################################
        # figureBox, containing the band and wavefunc figure
        # ######################################################
        # set up quantumCanvas for band structure plot
        self.quantumCanvas = EJcanvas(xlabel=u'Position (Å)',
                                      ylabel='Energy (eV)', parent=self)
        self.plotControl = EJplotControl(self.quantumCanvas, self)
        figureBox = QVBoxLayout()
        figureBox.addWidget(self.quantumCanvas)
        # vBox4: figureBox end

        #######################################################
        # solveBox, containing all eigensolver and calculation control
        ########################################################
        solveBox = QVBoxLayout()
        self.solveBasisButton = QPushButton("Solve Basis")
        self.solveBasisButton.clicked.connect(self.solve_basis)
        solveBox.addWidget(self.solveBasisButton)
        self.solveWholeButton = QPushButton("Solve Whole")
        self.solveWholeButton.clicked.connect(self.solve_whole)
        solveBox.addWidget(self.solveWholeButton)

        # set up description box
        self.DescriptionBox = QTextEdit('')
        self.DescriptionBox.setReadOnly(False)
        self.DescriptionBox.setSizePolicy(QSizePolicy(
            QSizePolicy.Fixed, QSizePolicy.Fixed))
        self.DescriptionBox.setMaximumHeight(60)
        self.DescriptionBox.setMaximumWidth(DescriptionBoxWidth)
        self.DescriptionBox.textChanged.connect(self.input_description)
        DescLayout = QVBoxLayout()
        DescLayout.addWidget(self.DescriptionBox)
        DescLayout_groupBox = QGroupBox("Description")
        DescLayout_groupBox.setLayout(DescLayout)
        solveBox.addWidget(DescLayout_groupBox)

        # set up material composition inputs
        self.mtrl_well = QLabel()
        self.mtrl_barr = QLabel()
        self.input_substrate('InP')
        # TODO: set these according to substrate dict
        #  self.MoleFracWellBox = []
        #  self.MoleFracBarrBox = []
        #  self.offsetLabel = []
        #  for n in range(self.numMaterials // 2):
        #      self.MoleFracWellBox.append(QDoubleSpinBox())
        #      self.MoleFracWellBox[n].setDecimals(3)
        #      self.MoleFracWellBox[n].setValue(0.53)
        #      self.MoleFracWellBox[n].setRange(0.0, 1.0)
        #      self.MoleFracWellBox[n].setSingleStep(0.001)
        #      self.MoleFracWellBox[n].editingFinished.connect(
        #          partial(self.input_moleFrac, 2 * n))
        #      self.MoleFracBarrBox.append(QDoubleSpinBox())
        #      self.MoleFracBarrBox[n].setDecimals(3)
        #      self.MoleFracBarrBox[n].setValue(0.52)
        #      self.MoleFracBarrBox[n].setRange(0.0, 1.0)
        #      self.MoleFracBarrBox[n].setSingleStep(0.001)
        #      self.MoleFracBarrBox[n].editingFinished.connect(
        #          partial(self.input_moleFrac, 2 * n + 1))
        #      self.offsetLabel.append(QLabel(''))
        self.MoleFracWellBox = QDoubleSpinBox()
        self.MoleFracWellBox.setDecimals(3)
        self.MoleFracWellBox.setValue(0.53)
        self.MoleFracWellBox.setRange(0.0, 1.0)
        self.MoleFracWellBox.setSingleStep(0.001)
        self.MoleFracWellBox.editingFinished.connect(
            partial(self.input_moleFrac, 0))
        self.MoleFracBarrBox = QDoubleSpinBox()
        self.MoleFracBarrBox.setDecimals(3)
        self.MoleFracBarrBox.setValue(0.52)
        self.MoleFracBarrBox.setRange(0.0, 1.0)
        self.MoleFracBarrBox.setSingleStep(0.001)
        self.MoleFracBarrBox.editingFinished.connect(
            partial(self.input_moleFrac, 1))
        self.offsetLabel = QLabel('')
        self.strainDescription = QLabel('')
        self.LOPhononDescription = QLabel('')
        mtrl_grid = QGridLayout()
        #  mtrl_grid.addWidget(QLabel(
        #      '<center><b>Mole Fractions</b></center>'), 0, 0, 1, 3)
        mtrl_grid.addWidget(QLabel(u'<center><b>#</b></center>'), 1, 0)
        mtrl_grid.addWidget(self.mtrl_well, 1, 1)
        mtrl_grid.addWidget(self.mtrl_barr, 1, 2)
        self.mtrl_indxBox = QComboBox()
        self.mtrl_indxBox.addItems(['1', '2', '3', '4'])
        #  self.mtrl_indxBox.setSizeAdjustPolicy(0)
        self.mtrl_indxBox.setSizePolicy(
            QSizePolicy(QSizePolicy.Maximum, QSizePolicy.Fixed))
        mtrl_grid.addWidget(self.mtrl_indxBox, 2, 0)
        mtrl_grid.addWidget(self.MoleFracWellBox, 2, 1)
        mtrl_grid.addWidget(self.MoleFracBarrBox, 2, 2)
        #  for n in range(self.numMaterials // 2):
        #      mtrl_grid.addWidget(QLabel(
        #          '<center><b>#%d</b></center>' % (n + 1)), 2 + n, 0)
        #      mtrl_grid.addWidget(self.MoleFracWellBox[n], 2 + n, 1)
        #      mtrl_grid.addWidget(self.MoleFracBarrBox[n], 2 + n, 2)
        #      mtrl_grid.addWidget(self.offsetLabel[n], 2 + n, 3)
        mtrl_grid.addWidget(QLabel('<center>(well)</center>'), 3, 1)
        mtrl_grid.addWidget(QLabel('<center>(barrier)</center>'), 3, 2)
        mtrl_grid.addWidget(self.offsetLabel, 4, 0, 1, 3)
        mtrl_grid.addWidget(self.strainDescription, 5, 0, 1, 3)
        mtrl_grid.addWidget(self.LOPhononDescription, 6, 0, 1, 3)
        mtrl_groupBox = QGroupBox("Mole Fractions")
        mtrl_groupBox.setLayout(mtrl_grid)
        solveBox.addWidget(mtrl_groupBox)

        # set up plot control inputs
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
        plotControl_grid = QGridLayout()
        plotControl_grid.addWidget(self.wellSelectButton, 0, 0, 1, 2)
        plotControl_grid.addWidget(self.zoominButton, 1, 0, 1, 1)
        plotControl_grid.addWidget(self.zoomOutButton, 1, 1, 1, 1)
        plotControl_grid.addWidget(self.panButton, 2, 0, 1, 1)
        plotControl_grid.addWidget(self.clearWFsButton, 2, 1, 1, 1)
        plotControl_groupBox = QGroupBox("Plot Controls")
        plotControl_groupBox.setLayout(plotControl_grid)
        solveBox.addWidget(plotControl_groupBox)

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
        # SLOT connected to SIGNAL of tranferOpticalParametersButton should
        # be defined out side this tab, in the main window
        self.pairSelectString = QTextEdit('')
        self.pairSelectString.setReadOnly(True)
        self.pairSelectString.setMaximumWidth(pairSelectStringWidth)
        self.pairSelectString.setSizePolicy(QSizePolicy(
            QSizePolicy.Fixed, QSizePolicy.Fixed))
        calculateControl_grid = QGridLayout()
        calculateControl_grid.addWidget(self.pairSelectButton, 0, 0, 1, 2)
        calculateControl_grid.addWidget(self.FoMButton, 1, 0, 1, 1)
        calculateControl_grid.addWidget(self.pairSelectString, 2, 0, 1, 2)
        calculateControl_groupBox = QGroupBox("Calculate")
        calculateControl_groupBox.setLayout(calculateControl_grid)
        solveBox.addWidget(calculateControl_groupBox)

        solveBox.addStretch()
        # vBox3: solveBox end

        quantumLayout = QHBoxLayout()
        quantumLayout.addLayout(settingBox)
        quantumLayout.addLayout(layerBox)
        quantumLayout.addLayout(solveBox)
        quantumLayout.addLayout(figureBox)

        self.setLayout(quantumLayout)
        self.setAutoFillBackground(True)
        self.setBackgroundRole(QPalette.Window)

        self.layerTable_refresh()
        self.layerTable.selectRow(1)
        self.updating = False

        self.update_inputBoxes()
        self.layerTable.setFocus()
    # __init__ end

    def settingslot(fn):
        @wraps(fn)
        def wrapper(self, *args, **kwargs):
            if self.updating:
                return
            else:
                #  print args
                #  print kwargs
                return fn(self, *args, **kwargs)
        return wrapper

    def reload(self):
        self.update_Lp_limits()
        self.update_inputBoxes()
        self.layerTable_refresh()
        self.updating = True
        self.layerTable.selectRow(1)
        self.layerTable.setFocus()
        self.updating = False
        self.qclayers.populate_x()
        self.update_quantumCanvas()

# ===========================================================================
# Input Controls
# ===========================================================================
    def update_inputBoxes(self):
        """ Update all input boxes except for the layerTable.
        SLOT connected to LpFirstSpinbox/LpLastSpinBox.valueChanged(int)
        """
        self.updating = True
        try:
            self.inputSubstrateBox.setCurrentText(self.qclayers.substrate)
        except Exception:
            QMessageBox.warning(self, "ErwinJr - Warning",
                                "Substrate data wrong.\n" +
                                traceback.format_exc())

        self.qclayers.update_alloys()
        self.qclayers.update_strain()
        self.qclayers.populate_x()
        mtrl_indx = int(unicode(self.mtrl_indxBox.currentText())) - 1
        self.MoleFracWellBox.setValue(
            self.qclayers.moleFrac[2 * mtrl_indx])
        self.MoleFracBarrBox.setValue(
            self.qclayers.moleFrac[2 * mtrl_indx + 1])
        self.offsetLabel.setText(
            u'<center>ΔE<sub>c</sub>: <b>%6.0f meV </b></center>' % (
                (self.qclayers.EcG[2 * mtrl_indx + 1] -
                 self.qclayers.EcG[2 * mtrl_indx]) * 1000))
        #  for n in range(self.numMaterials // 2):
        #      self.MoleFracWellBox[n].setValue(self.qclayers.moleFrac[2 * n])
        #      self.MoleFracBarrBox[n].setValue(self.qclayers.moleFrac[2*n+1])
        #      self.offsetLabel[n].setText("%6.0f meV" % (
        #          (self.qclayers.EcG[2 * n + 1] -
        #           self.qclayers.EcG[2 * n]) * 1000))

        self.DescriptionBox.setText(self.qclayers.description)
        strainString = ("<center>Net Strain: <b>%6.3f%%</b></center>" %
                        self.qclayers.netStrain)
        self.strainDescription.setText(strainString)
        hwLOString = ("<center>E<sub>LO</sub>:"
                      "<b>%4.1f ~ %4.1f meV</b></center>") % (
                          self.qclayers.hwLO[2 * mtrl_indx] * 1000,
                          self.qclayers.hwLO[2 * mtrl_indx + 1] * 1000)
        self.LOPhononDescription.setText(hwLOString)

        self.inputVertResBox.setValue(self.qclayers.vertRes)
        self.inputEFieldBox.setValue(self.qclayers.EField)
        try:
            n = self.HorzRes.index(self.qclayers.xres)
        except ValueError:
            print "Horizontal resolution error: %.2f"%self.qclayers.xres
        self.inputHorzResBox.setCurrentIndex(n)
        self.inputRepeatsBox.setValue(self.qclayers.repeats)
        self.update_Lp_box()

        self.dirty.emit()
        self.updating = False

    def input_substrate(self, substrateType):
        """ SLOT connected to inputSubstrateBox.currentIndexChanged(QString)
        update substrate chosen """
        # TODO put the following infor mation to a dict
        if substrateType == 'InP':
            self.qclayers.substrate = 'InP'
            self.materialList = ['InGaAs/AlInAs #1', 'InGaAs/AlInAs #2',
                                 'InGaAs/AlInAs #3', 'InGaAs/AlInAs #4']
            self.mtrl_well.setText('<center><b>'
                                   'In<sub>x</sub>Ga<sub>1-x</sub>As'
                                   '</b></center>')
            self.mtrl_barr.setText('<center><b>'
                                   'Al<sub>1-x</sub>In<sub>x</sub>As'
                                   '</b></center')

        elif substrateType == 'GaAs':
            self.qclayers.substrate = 'GaAs'
            self.materialList = ['AlGaAs/AlGaAs #1', 'AlGaAs/AlGaAs #2',
                                 'AlGaAs/AlGaAs #3', 'AlGaAs/AlGaAs #4']
            self.mtrl_well.setText('<center><b>'
                                   'Al<sub>x</sub>Ga<sub>1-x</sub>As'
                                   '</b></center')
            self.mtrl_barr.setText('<center><b>'
                                   'Al<sub>x</sub>Ga<sub>1-x</sub>As'
                                   '</b></center')

        elif substrateType == 'GaSb':
            self.qclayers.substrate = 'GaSb'
            self.materialList = ['InAsSb/AlGaSb #1', 'InAsSb/AlGaSb #2',
                                 'InAsSb/AlGaSb #3', 'InAsSb/AlGaSb #4']
            self.mtrl_well.setText('<center><b>\
                    InAs<sub>y</sub>Sb<sub>1-y</sub>\
                    </b></center')
            self.mtrl_barr.setText('<center><b>\
                    Al<sub>x</sub>Ga<sub>1-x</sub>Sb\
                    </b></center')

        elif substrateType == 'GaN':
            #  self.input_substrate(self.qclayers.substrate)
            QMessageBox.information(
                self, 'ErwinJr Error',
                'III-Nitride substrates have not yet been implemented.')
            self.inputSubstrateBox.setCurrentIndex(
                self.inputSubstrateBox.findText(self.qclayers.substrate))
            return

        else:
            raise TypeError('substrate selection not allowed')
            return

        if not self.updating:
            self.update_Lp_limits()
            self.update_inputBoxes()
            self.layerTable_refresh()
            self.qclayers.populate_x()
            self.update_quantumCanvas()

    @settingslot
    def input_EField(self, *args):
        """ SLOT connected to inputEFieldBox.valueChanged(double)
        update external E field in unit kV/cm """
        self.qclayers.EField = float(self.inputEFieldBox.value())
        self.qclayers.populate_x()
        self.update_quantumCanvas()
        self.dirty.emit()

    @settingslot
    def input_horzRes(self, *args):
        """ SLOT connected to inputHorzResBox.currentIndexChanged(int)
        update position resolution (xres), in angstrom """
        horzRes = self.HorzRes[self.inputHorzResBox.currentIndex()]
        self.qclayers.set_xres(horzRes)
        self.qclayers.populate_x()
        self.update_quantumCanvas()
        self.dirty.emit()

    @settingslot
    def input_vertRes(self, *args):
        """ SLOT connected to inputVertResBox.valueChanged
        Update initial energy resolution for eigensolver. Set this too small
        may result in loss of some eigenvalue. """
        self.qclayers.vertRes = float(self.inputVertResBox.value())
        self.dirty.emit()

    @settingslot
    def input_repeats(self, *args):
        """ SLOT connected to SINGAL self.inputRepeatsBox.valueChanged(int)
        update number of repeats for the whole structure."""
        self.qclayers.repeats = int(self.inputRepeatsBox.value())
        self.qclayers.populate_x()
        self.update_quantumCanvas()
        self.dirty.emit()

    def input_basis(self):
        """ SLOT connected to self.inputARInjectorCheck.stateChanged(int) and
        self.inputInjectorARCheck.stateChanged(int)
        update dividers info
        """
        self.qclayers.basisARInjector = self.inputARInjectorCheck.isChecked()
        self.qclayers.basisInjectorAR = self.inputInjectorARCheck.isChecked()
        self.dirty.emit()

    def update_Lp_limits(self):
        """
        Update Lp select range in the Period Info box (GUI)
        """
        self.LpFirstSpinbox.setRange(1, self.qclayers.layerWidth.size - 1)
        self.LpFirstSpinbox.setValue(1)
        self.LpLastSpinbox.setRange(1, self.qclayers.layerWidth.size - 1)
        self.LpLastSpinbox.setValue(self.qclayers.layerWidth.size - 1)

    def update_Lp_box(self):
        """
        Update Lp box in the Period Info box (GUI):
            Lp:total length
            well: persentage of well material
            nD: average doping (cm-3)
            ns: 2D carrier density in 1E11 cm-2
        """
        LpFirst = self.LpFirstSpinbox.value()
        LpLast = self.LpLastSpinbox.value() + 1
        # +1 because range is not inclusive of last value
        # total length of the layers (1 period)
        Lp = sum(self.qclayers.layerWidth[LpFirst:LpLast]) * self.qclayers.xres
        Lp_string = u"Lp: %g \u212B<br>" % Lp
        # total length of well (1 period)
        Lw = sum((1 - self.qclayers.layerBarriers[LpFirst:LpLast]) *
                 self.qclayers.layerWidth[LpFirst:LpLast]) * self.qclayers.xres
        if Lp == 0:
            Lp_string += u"wells: NA%%<br>"
            # average doping of the layers
            Lp_string += (u"n<sub>D</sub>: NA\u00D710<sup>17</sup>"
                          u"cm<sup>-3</sup><br>")
        else:
            Lp_string += u"wells: %6.1f%%<br>" % (100.0 * Lw / Lp)
            # average doping of the layers
            nD = self.qclayers.xres * sum(
                self.qclayers.layerDopings[LpFirst:LpLast] *
                self.qclayers.layerWidth[LpFirst:LpLast]) / Lp
            Lp_string += (u"n<sub>D</sub>: %6.3f\u00D710<sup>17</sup>"
                          u"cm<sup>-3</sup><br>") % nD
        # 2D carrier density in 1E11cm-2
        ns = self.qclayers.xres * sum(
            self.qclayers.layerDopings[LpFirst:LpLast] *
            self.qclayers.layerWidth[LpFirst:LpLast]) * 1e-2
        Lp_string += (u"n<sub>s</sub>: %6.3f\u00D710<sup>11</sup>"
                      u"cm<sup>-2</sup") % ns
        self.LpStringBox.setText(Lp_string)

    @settingslot
    def input_description(self):
        """ SLOT connected to self.DescriptionBox.textChanged()
        Change description string for the design"""
        self.qclayers.description = self.DescriptionBox.toPlainText()
        self.dirty.emit()

    @settingslot
    def select_mtrl_indx(self, *args):
        """ SLOT connected to mtrl_indxBox.currentIndexChanged(int)
        change material index shown in UI"""
        mtrl_indx = int(unicode(self.mtrl_indxBox.currentText())) - 1
        self.MoleFracWellBox.setValue(
            self.qclayers.moleFrac[2 * mtrl_indx])
        self.MoleFracWellBox.editingFinished.connect(
            partial(self.input_moleFrac, 2 * mtrl_indx))
        self.MoleFracWellBox.setValue(
            self.qclayers.moleFrac[2 * mtrl_indx + 1])
        self.MoleFracBarrBox.editingFinished.connect(
            partial(self.input_moleFrac, 2 * mtrl_indx + 1))

    @settingslot
    def input_moleFrac(self, boxID):
        """ SLOT connected to self.MoleFracWellBox.editingFinished()
            shen mtrl_indx is n
            and MoleFracBarrBox[n].editingFinished()
        Update moleFrac for material."""
        self.qclayers.moleFrac[boxID] = float(
            self.MoleFracWellBox.value() if boxID % 2 == 0
            else self.MoleFracBarrBox.value())
        self.dirty.emit()

        self.update_inputBoxes()
        self.qclayers.update_alloys()
        self.qclayers.update_strain()
        self.qclayers.populate_x()
        self.update_quantumCanvas()

    @settingslot
    def set_temperature(self, T):
        self.qclayers.Temperature = T
        self.qclayers.populate_x()
        self.qclayers.populate_x_band()
        self.update_quantumCanvas()

# ===========================================================================
# Layer Table Control
# ===========================================================================
    def layerNum(self):
        return self.qclayers.layerWidth.size

    def layerTable_refresh(self):
        """Refresh layer table, called every time after data update"""
        # Block itemChanged SIGNAL while refreshing
        #  self.clear_WFs()
        self.layerTable.blockSignals(True)
        self.layerTable.clear()
        self.layerTable.setColumnCount(6)
        self.layerTable.setRowCount(self.qclayers.layerWidth.size + 1)
        self.layerTable.setHorizontalHeaderLabels(
            ['Width', 'ML', 'Brr', 'AR', 'Doping', 'Material'])
        #  vertLabels = []
        #  for n in xrange(self.qclayers.layerWidth.size+1):
        #      vertLabels.append(str(n))
        vertLabels = [str(n) for n in range(self.qclayers.layerWidth.size + 1)]
        self.layerTable.setVerticalHeaderLabels(vertLabels)

        # color for barrier layers
        gray = QColor(230, 230, 240)  # for Barrier layers
        gray2 = QColor(230, 230, 230)  # for unchangable background

        for q, layerWidth in enumerate(self.qclayers.layerWidth):
            # Width Setup
            width = QTableWidgetItem("%5.1f" %
                                     (layerWidth * self.qclayers.xres))
            width.setTextAlignment(Qt.AlignCenter)
            if bool(self.qclayers.layerBarriers[q]):
                width.setBackground(gray)
            self.layerTable.setItem(q, 0, width)
            if q == 0:
                width.setFlags(Qt.NoItemFlags)
                width.setBackground(gray2)

            # ML Setup
            numML = (self.qclayers.xres * layerWidth /
                     self.qclayers.MLThickness[q])
            item = QTableWidgetItem("%5.1f" % numML)
            item.setTextAlignment(Qt.AlignCenter)
            if bool(self.qclayers.layerBarriers[q]):
                item.setBackground(gray)
            self.layerTable.setItem(q, 1, item)
            if q == 0:
                item.setFlags(Qt.NoItemFlags)
                item.setBackground(gray2)

            # Barrier Layer Setup
            item = QTableWidgetItem()
            #  item.setCheckState(int(self.qclayers.layerBarriers[q])*2)
            item.setCheckState(
                Qt.Checked if self.qclayers.layerBarriers[q] == 1
                else Qt.Unchecked)
            if bool(self.qclayers.layerBarriers[q]):
                item.setBackground(gray)
            self.layerTable.setItem(q, 2, item)
            if q == 0:
                item.setFlags(Qt.NoItemFlags)
                item.setBackground(gray2)

            # Active Region Layer Setup
            item = QTableWidgetItem()
            #  item.setCheckState(int(self.qclayers.layerARs[q])*2)
            item.setCheckState(
                Qt.Checked if self.qclayers.layerARs[q] == 1
                else Qt.Unchecked)
            if bool(self.qclayers.layerBarriers[q]):
                item.setBackground(gray)
            self.layerTable.setItem(q, 3, item)
            if q == 0:
                item.setFlags(Qt.NoItemFlags)
                item.setBackground(gray2)

            # Layer Doping Setup
            doping = QTableWidgetItem(unicode(self.qclayers.layerDopings[q]))
            doping.setTextAlignment(Qt.AlignCenter)
            if bool(self.qclayers.layerBarriers[q]):
                doping.setBackground(gray)
            self.layerTable.setItem(q, 4, doping)
            if q == 0:
                doping.setFlags(Qt.NoItemFlags)
                doping.setBackground(gray2)

            # Material Setup
            if q == 0:
                item = QTableWidgetItem(unicode(self.materialList[
                    int(self.qclayers.layerMaterials[q]) - 1]))
                # TODO: reformat layerMaterials to int begin at 0
                item.setBackground(gray2)
                item.setFlags(Qt.NoItemFlags)
                self.layerTable.setItem(q, 5, item)
            else:
                materialWidget = QComboBox()
                materialWidget.addItems(self.materialList)
                materialWidget.setCurrentIndex(
                    self.qclayers.layerMaterials[q] - 1)
                materialWidget.currentIndexChanged.connect(
                    partial(self.layerTable_materialChanged, q))
                self.layerTable.setCellWidget(q, 5, materialWidget)

        self.layerTable.resizeColumnsToContents()

        self.layerTable.blockSignals(False)

    def insert_layerAbove(self):
        """ SLOT connected to self.insertLayerAboveButton.clicked()"""
        row = self.layerTable.currentRow()
        if row == -1:
            return

        self.qclayers.layerWidth = np.insert(
            self.qclayers.layerWidth, row, 0)
        self.qclayers.layerBarriers = np.insert(
            self.qclayers.layerBarriers, row, 0
            if self.qclayers.layerBarriers[row] == 1 else 1)
        self.qclayers.layerARs = np.insert(
            self.qclayers.layerARs, row,
            self.qclayers.layerARs[row])
        self.qclayers.layerMaterials = np.insert(
            self.qclayers.layerMaterials, row,
            self.qclayers.layerMaterials[row])
        self.qclayers.layerDopings = np.insert(
            self.qclayers.layerDopings, row,
            self.qclayers.layerDopings[row])
        self.qclayers.layerDividers = np.insert(
            self.qclayers.layerDividers, row,
            self.qclayers.layerDividers[row])

        self.update_Lp_limits()
        self.update_Lp_box()

        self.layerTable_refresh()
        self.qclayers.populate_x()
        self.layerTable.selectRow(row)
        self.layerTable.setFocus()

        self.dirty.emit()

    def delete_layer(self):
        """ SLOT connected to self.deleteLayerButton.clicked()"""
        # don't delete last layer
        if self.qclayers.layerWidth.size == 1:
            return
        row = self.layerTable.currentRow()
        if row == -1 or row >= self.qclayers.layerWidth.size:
            return

        self.qclayers.layerWidth = np.delete(
            self.qclayers.layerWidth, row)
        self.qclayers.layerBarriers = np.delete(
            self.qclayers.layerBarriers, row)
        self.qclayers.layerARs = np.delete(
            self.qclayers.layerARs, row)
        self.qclayers.layerMaterials = np.delete(
            self.qclayers.layerMaterials, row)
        self.qclayers.layerDopings = np.delete(
            self.qclayers.layerDopings, row)
        self.qclayers.layerDividers = np.delete(
            self.qclayers.layerDividers, row)

        if row == self.qclayers.layerWidth.size:  # if row == last_row
            # make first item the same as last item
            self.qclayers.layerWidth[0] = self.qclayers.layerWidth[-1]
            self.qclayers.layerBarriers[0] = self.qclayers.layerBarriers[-1]
            self.qclayers.layerARs[0] = self.qclayers.layerARs[-1]
            self.qclayers.layerMaterials[0] = self.qclayers.layerMaterials[-1]
            self.qclayers.layerDopings[0] = self.qclayers.layerDopings[-1]
            self.qclayers.layerDividers[0] = self.qclayers.layerDividers[-1]

        self.update_Lp_limits()
        self.update_Lp_box()

        self.qclayers.update_strain()
        self.layerTable_refresh()
        self.qclayers.populate_x()
        self.layerTable.selectRow(row)
        self.layerTable.setFocus()

        self.dirty.emit()

    def layerTable_itemChanged(self, item):
        """SLOT connected to layerTable.itemChanged(QTableWidgetItem*)
        Update layer profile after user input"""
        # TODO: redo illegal input
        #  column = self.layerTable.currentColumn()
        #  row = self.layerTable.currentRow()
        #  print "---debug, itemChanged--- (%d, %d)"%(column, row)
        #  print "--debug, itemChanged (%d, %d)"%(item.column(), item.row())
        #  print item.text()
        #  if column == -1: #column == -1 on GUI initialization
        #      return
        column = item.column()
        row = item.row()
        if column == 0:  # column == 0 for item change in Widths column
            new_width = float(item.text())
            new_width_int = int(np.round(new_width / self.qclayers.xres))
            #  if np.mod(new_width, self.qclayers.xres) != 0 \
            #          and self.qclayers.xres != 0.1:
            if np.abs(new_width_int * self.qclayers.xres - new_width) > 1E-9:
                # TODO: bug to fix, np.mod is not good for xres < 0.5
                # potential solution is to change internal length to int
                # times xres
                QMessageBox.warning(self, "ErwinJr - Warning", (
                    "You entered a width that is not compatible with "
                    "the minimum horizontal resolution. "
                    "%f %% %f = %f" % (new_width, self.qclayers.xres,
                                       np.mod(new_width, self.qclayers.xres))
                ))
                return
            if row == self.qclayers.layerWidth.size:  # add row at end of list
                self.qclayers.layerWidth = np.append(
                    self.qclayers.layerWidth, new_width_int)
                self.qclayers.layerBarriers = np.append(
                    self.qclayers.layerBarriers,
                    0 if self.qclayers.layerBarriers[-1] == 1 else 1)
                self.qclayers.layerARs = np.append(
                    self.qclayers.layerARs,
                    self.qclayers.layerARs[-1])
                self.qclayers.layerMaterials = np.append(
                    self.qclayers.layerMaterials,
                    self.qclayers.layerMaterials[-1])
                self.qclayers.layerDopings = np.append(
                    self.qclayers.layerDopings,
                    self.qclayers.layerDopings[-1])
                self.qclayers.layerDividers = np.append(
                    self.qclayers.layerDividers,
                    self.qclayers.layerDividers[-1])
                row += 1  # used so that last (blank) row is again selected

                # make first item the same as last item
                for LayerD in (self.qclayers.layerWidth,
                               self.qclayers.layerBarriers,
                               self.qclayers.layerARs,
                               self.qclayers.layerMaterials,
                               self.qclayers.layerDopings,
                               self.qclayers.layerDividers):
                    LayerD[0] = LayerD[-1]
                self.update_Lp_limits()

            elif row == self.qclayers.layerWidth.size - 1:
                self.qclayers.layerWidth[row] = new_width_int
                # make first item the same as last item
                self.qclayers.layerWidth[0] = self.qclayers.layerWidth[-1]

            else:  # change Width of selected row in-place
                self.qclayers.layerWidth[row] = new_width_int

        elif column == 1:  # column == 1 for ML
            if self.qclayers.xres != 0.1:
                QMessageBox.warning(self, "ErwinJr - Warning", (
                    u"Horizontal Resolution of 0.1 \u212B required"
                    u"when setting monolayer thicknesses."))
                return
            if row == self.qclayers.layerWidth.size:  # add row at end of list
                pass
            elif row == self.qclayers.layerWidth.size - 1:
                self.qclayers.layerWidth[row] = int(np.round(
                    self.qclayers.MLThickness[row] * float(item.text()) /
                    self.qclayers.xres))

                # make first item the same as last item
                self.qclayers.layerWidth[0] = self.qclayers.layerWidth[-1]

                self.update_Lp_limits()

            else:  # change Width of selected row in-place
                self.qclayers.layerWidth[row] = int(np.round(
                    self.qclayers.MLThickness[row] * float(item.text()) /
                    self.qclayers.xres))
        elif column == 2:  # column == 2 for item change in Barrier column
            if row == self.qclayers.layerWidth.size:
                # don't do anything if row is last row
                return
            #  self.qclayers.layerBarriers[row] = int(item.checkState())//2
            self.qclayers.layerBarriers[row] = (
                item.checkState() == Qt.Checked)
            if row == self.qclayers.layerWidth.size - 1:
                self.qclayers.layerBarriers[0] = \
                    self.qclayers.layerBarriers[-1]

        elif column == 3:  # column == 3 for item change in AR column
            if row == self.qclayers.layerWidth.size:
                # don't do anything if row is last row
                return
            #  self.qclayers.layerARs[row] = int(item.checkState())//2
            self.qclayers.layerARs[row] = (item.checkState() == Qt.Checked)
            if row == self.qclayers.layerWidth.size - 1:
                self.qclayers.layerARs[0] = self.qclayers.layerARs[-1]

        elif column == 4:  # column == 4 for item change in Doping column
            if row == self.qclayers.layerWidth.size:
                # don't do anything if row is last row
                return
            self.qclayers.layerDopings[row] = float(item.text())
            if row == self.qclayers.layerWidth.size - 1:
                self.qclayers.layerDopings[0] = self.qclayers.layerDopings[-1]

        elif column == 5:  # column == 5 for item change in Materials column
            # See layerTable_materialChanged for more information
            # self.qclayers.layerWidth[row] = int(item.text[row])
            if row == self.qclayers.layerWidth.size - 1:
                self.qclayers.layerMaterials[0] =\
                    self.qclayers.layerMaterials[-1]
        else:
            pass

        if not self.updating:
            self.layerTable_refresh()
            self.layerTable.setCurrentCell(row, column)
            self.layerTable.setFocus()
            self.update_Lp_box()
            self.update_quantumCanvas

        self.dirty.emit()

    def layerTable_itemSelectionChanged(self):
        """SLOT connected to layerTable.itemSelectionChanged()"""
        # This is the primary call to update_quantumCanvas
        self.qclayers.layerSelected = self.layerTable.currentRow()
        if not self.updating and self.qclayers.layerSelected >= 0 and \
                self.qclayers.layerSelected < self.qclayers.layerWidth.size:
            self.qclayers.populate_x()
            self.update_quantumCanvas()

    def layerTable_materialChanged(self, row, selection):
        """SLOT as partial(self.layerTable_materialChanged, q)) connected to
        materialWidget.currentIndexChanged(int) """
        self.qclayers.layerMaterials[row] = selection + 1
        # self.layerTable_refresh()
        self.qclayers.populate_x()
        self.layerTable.selectRow(row)

        self.dirty.emit()

    def bump_first_layer(self):
        """Move zeroth layer to first layer"""
        self.qclayers.layerWidth = np.insert(
            self.qclayers.layerWidth, 0, self.qclayers.layerWidth[-1])
        self.qclayers.layerBarriers = np.insert(
            self.qclayers.layerBarriers, 0, self.qclayers.layerBarriers[-1])
        self.qclayers.layerARs = np.insert(
            self.qclayers.layerARs, 0, self.qclayers.layerARs[-1])
        self.qclayers.layerMaterials = np.insert(
            self.qclayers.layerMaterials, 0, self.qclayers.layerMaterials[-1])
        self.qclayers.layerDopings = np.insert(
            self.qclayers.layerDopings, 0, self.qclayers.layerDopings[-1])
        self.qclayers.layerDividers = np.insert(
            self.qclayers.layerDividers, 0, self.qclayers.layerDividers[-1])

        self.update_inputBoxes()
        self.layerTable_refresh()
        self.layerTable.setCurrentCell(1, 0)
        self.layerTable.setFocus()
        self.update_quantumCanvas()
        self.dirty.emit()

    def copy_structure(self):
        clipboard = QApplication.clipboard()
        string = ''
        for layer in self.qclayers.layerWidth[1:]:
            string += '%g\n' % (layer * self.qclayers.xres)
        clipboard.setText(string)

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
        self.quantumCanvas.axes.plot(self.qclayers.xPoints,
                                     self.qclayers.xARs, 'k', linewidth=1.5)
        if self.qclayers.layerSelected >= 0 and \
                self.qclayers.layerSelected < self.qclayers.layerWidth.size:
            self.quantumCanvas.axes.plot(
                self.qclayers.xPoints, self.qclayers.xLayerSelected, 'b',
                linewidth=1.5 if self.qclayers.layerARs[
                    self.qclayers.layerSelected] == 1 else 1)

        if hasattr(self.qclayers, 'EigenE'):
            self.curveWF = []
            #  scale = np.max(np.abs(self.qclayers.xyPsiPlot[:, n]))
            if self.plotType == "mode":
                y = self.qclayers.xyPsiPsi
            elif self.plotType == "wf":
                y = self.qclayers.xyPsiPlot
            for n in range(self.qclayers.EigenE.size):
                curve, = self.quantumCanvas.axes.plot(
                    self.qclayers.xPointsPost,
                    y[:, n] + self.qclayers.EigenE[n],
                    color=self.colors[np.mod(n, len(self.colors))])
                if self.fillplot:
                    self.quantumCanvas.axes.fill_between(
                        self.qclayers.xPointsPost,
                        y[:, n] + self.qclayers.EigenE[n],
                        self.qclayers.EigenE[n],
                        facecolor=self.colors[np.mod(n, len(self.colors))],
                        alpha=self.fillplot)
                self.curveWF.append(curve)
            for n in self.stateHolder:
                curve, = self.quantumCanvas.axes.plot(
                    self.qclayers.xPointsPost,
                    y[:, n] + self.qclayers.EigenE[n],
                    'k', linewidth=2)

        self.quantumCanvas.draw()

    def well_select(self, event):
        """ callback registered in plotControl when it's in wellselect mode.
        It's mpl_connect to button_release_event of quantumCanvas """
        if event.button == 1:  # left button clicked
            x = event.xdata
            xLayerNum = np.argmin((self.qclayers.xPoints - x)**2)
            layerNum = self.qclayers.xLayerNums[xLayerNum]
            self.layerTable.selectRow(layerNum)
            self.layerTable.setFocus()

    def clear_WFs(self):
        if hasattr(self.qclayers, 'EigenE'):
            delattr(self.qclayers, 'EigenE')
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

        if hasattr(self.qclayers, 'xyPsiPsi'):
            # otherwise band structure hasn't been solved yet
            xyPsiPsiEig = np.zeros(self.qclayers.xyPsiPsi.shape)
            for q in xrange(self.qclayers.EigenE.size):
                xyPsiPsiEig[:, q] = (self.qclayers.xyPsiPsi[:, q] +
                                     self.qclayers.EigenE[q])
            np.savetxt(
                fname.split('.')[0] + '_States' + '.csv',
                np.column_stack([self.qclayers.xPointsPost, xyPsiPsiEig]),
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

    def solve_whole(self):  # solves whole structure
        """SLOT connected to solveWholeButton.clicked()
        Whole solver """
        if hasattr(self.qclayers, "EigenE"):
            self.clear_WFs()
        self.pairSelected = False
        self.Calculating(True)

        self.qclayers.populate_x_band()
        try:
            self.stateHolder = []
            self.qclayers.solve_psi()
            self.solveType = 'whole'
            self.update_quantumCanvas()
            self.pairSelectButton.setEnabled(True)
            if DEBUG >= 4:
                with open('qclayer.pkl', 'wb') as f:
                    pickle.dump(self.qclayers, f, pickle.HIGHEST_PROTOCOL)
                print self.qclayers.EigenE
        except (IndexError, TypeError):
            QMessageBox.warning(self, 'ErwinJr - Error',
                                traceback.format_exc())

        self.Calculating(False)

    def solve_basis(self):  # solves structure with basis
        """SLOT connected to solveBasisButton.clicked()
        Basis solver """
        self.Calculating(True)

        try:
            self.stateHolder = []
            self.dCL = self.qclayers.basisSolve()
            self.qclayers.convert_dCL_to_data(self.dCL)
            self.solveType = 'basis'
            self.update_quantumCanvas()
            self.pairSelectButton.setEnabled(True)
        except (ValueError, IndexError):
            QMessageBox.warning(self, "ErwinJr - Error",
                                traceback.format_exc())

        self.Calculating(False)

    def state_pick(self, event):
        """ callback registered in plotControl when it's in pairselect mode.
        It's mpl_connect to button_release_event of quantumCanvas """
        if not hasattr(self.qclayers, 'EigenE'):
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
            xData = np.tile(self.qclayers.xPointsPost,
                            (self.qclayers.xyPsiPsi.shape[1], 1)).T
            if self.plotType in ("mode", "wf"):
                yData = self.qclayers.EigenE + (self.qclayers.xyPsiPsi
                                                if self.plotType == "mode"
                                                else self.qclayers.xyPsiPlot)

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
                self.qclayers.populate_x_band()
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
                self.qclayers.populate_x_band()
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
