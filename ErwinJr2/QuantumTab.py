"""
This file defines the quantum tab of ErwinJr, for simulating electron spectrum
and scattering
"""

# TODO:
# In plot controls, add "show one period"
# save and load pickle for qcLayers
# Reverse layers
# export excel file for growth sheet
# Inverse rotate
# last change history

import sys
import traceback
import numpy as np
from numpy import sqrt
from functools import partial, wraps
# figure is and should be used for gain spectrum plot only!
from matplotlib.pyplot import figure

from .QCLayers import (QCLayers, StateRecognizeError,
                       optimize_layer, optimize_global,
                       QCMaterial, h, hbar, eps0, c0, e0)
from .EJcanvas import EJcanvas, EJplotControl
from .EJcanvas import config as plotconfig

from PyQt5.QtCore import pyqtSignal, pyqtSlot, QObject, QThread, Qt
from PyQt5.QtGui import QPalette, QColor
from PyQt5.QtWidgets import (QWidget, QVBoxLayout, QLabel, QComboBox,
                             QSpinBox, QDoubleSpinBox, QGroupBox,
                             QCheckBox, QTextEdit, QSizePolicy,
                             QGridLayout, QHBoxLayout, QPushButton,
                             QTableWidget, QTableWidgetItem,
                             QApplication, QMessageBox,
                             QDialog, QDialogButtonBox)
from .customQTClass import mtrlComboBox

from .Material import AParam
from .versionAndName import ejError, ejWarning
from .darkDetect import isdark
from .QCPlotter import plotPotential, scaleWF, plotWF


# TODO: this may not be necessary by better designer
def settingslot(fn):
    """ A decorator to ask slots to skip reaction when doing massive
    updating (when self.updating is True)"""
    @wraps(fn)
    def wrapper(self, *args, **kwargs):
        if self.updating:
            return
        else:
            return fn(self, *args, **kwargs)
    return wrapper


class CalculateHolder(QObject):
    finished = pyqtSignal()
    succeed = pyqtSignal()
    failed = pyqtSignal(str)
    warning = pyqtSignal(str)

    def __init__(self, calc, parent=None):
        self.calc = calc
        super().__init__(parent)

    def run(self):
        try:
            self.calc()
            self.succeed.emit()
        except StateRecognizeError as e:
            self.warning.emit('Warning: ' + e.expression)
            self.succeed.emit()
        except (IndexError, TypeError):
            self.failed.emit(traceback.format_exc())
        finally:
            self.finished.emit()


class QuantumTab(QWidget):
    """The Quantum Tab of ErwinJr. This is designed to be a GUI wrapper of
    the class QCLayers
    Member variable (! label public interface):
        !qclayers: A class to describe the physics of QC layers
        mtrlList: for Material column in the layerTable
        colors: colors for different wavefunctions

        customized SIGNAL
            dirty: Show if there's any new changes to qclayers that need to
                   be saved

        --- state selection ---
        stateHolder: the list containing the index of selected states. the
                    index is defined in qclayers
        pairSelected: Boolean flag for if a pair of states is selected
        --- plot flags ---
        !plotVX, !plotVL, !plotLH, !plotSO: show X point of conduction band
                                            (L point, LH band, SO band)

        --- GUI widget ---
        1st column (settingBox):
            descBox
            inputSubstrateBox
            inputEFieldBox
            inputxResBox
            inputEresBox
            inputRepeatsBox
            inputWlBox
            inputARInjectorCheck
            inputInjectorARCheck
            LpFirstSpinbox          LpLastSpinbox

        2nd column (layerBox):
            insertLayerAboveButton  deleteLayerButton
            optimizeLayerButton     globalOptimizeButton
            layerTable

        3rd column (solveBox):
            solveBasisButton
            solveWholeButton
            mtrlTable
            ifrDeltaBox             ifrLambdaBox
            offsetLabel
            netStrainLabel
            LOPhononLabel
            layerSelectButton
            zoomInButton            zoomOutButton
            panButton               clearWFsButton
            pairSelectButton
            FoMButton               fullPopulationButton
            toOpticsButton          gainSpecButton
            stateParamText

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
    toOptics = pyqtSignal(QCLayers)

    def __init__(self, qclayers=None, parent=None):
        super(QuantumTab, self).__init__(parent)
        self.qclayers = qclayers if qclayers else QCLayers()
        self.calcThread = QThread()
        self._worker = None

        # colors for different wavefunctions
        self.colors = ((0.584, 0.450, 0.701), (0.431, 0.486, 0.745),
                       (0.576, 0.694, 0.517), (0.682, 0.780, 0.321),
                       (0.501, 0.501, 0.509), (0.854, 0.741, 0.247),
                       (0.874, 0.607, 0.290), (0.823, 0.341, 0.278),
                       (0.725, 0.321, 0.623), (0.411, 0.741, 0.270),
                       (0.078, 0.078, 0.078), (0.431, 0.803, 0.870),
                       (0.223, 0.321, 0.643))
        # colors for different materials in tables
        mtrlcolors_RGB = ((255, 255, 255), (230, 230, 240), (230, 240, 230),
                          (240, 230, 230), (230, 240, 240), (240, 230, 240),
                          (240, 240, 230), (230, 230, 230))
        self.mtrlcolors = tuple([QColor(*rgb) for rgb in mtrlcolors_RGB])
        if isdark:
            self.mtrlcolors = tuple([QColor(*(2*(255-c) for c in rgb))
                                     for rgb in mtrlcolors_RGB])

        self.updating = False
        self.plotVX = False
        self.plotVL = False
        self.plotLH = False
        self.plotSO = False
        self.layerSelected = None

        self.stateHolder = []
        self.pairSelected = False
        # plotType can be mode, wf or DoS (TODO)
        #  self.plotType = "wf"
        self.plotType = "mode"
        #  self.fillPlot = 0.3  # alpha of fill; False for not fill
        self.fillPlot = False

        # Platform dependent settings, eg. layout size settings
        if sys.platform.startswith('win'):
            settingBoxWidth = 150
            layerBoxWidth = 230
            solveBoxWidth = 170
        elif sys.platform.startswith('darwin'):
            settingBoxWidth = 160
            layerBoxWidth = 250
            solveBoxWidth = 180
        elif sys.platform.startswith('linux'):
            settingBoxWidth = 150
            layerBoxWidth = 250
            solveBoxWidth = 185
        else:
            QMessageBox.warning(self, ejWarning,
                                'Platform %s not tested.' % sys.platform)
            settingBoxWidth = 150
            layerBoxWidth = 400
            solveBoxWidth = 190

        quantumLayout = QHBoxLayout()
        quantumLayout.setSpacing(1)
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
        """ Return a Qt Layout object containing all setting parameters """
        settingBox = QVBoxLayout()

        # set up description box
        settingBox.addWidget(QLabel("<center><b>Description</b></center>"))
        self.descBox = QTextEdit('')
        self.descBox.setReadOnly(False)
        self.descBox.setMinimumHeight(50)
        self.descBox.setMaximumHeight(80)
        self.descBox.setFixedWidth(width)
        self.descBox.setSizePolicy(QSizePolicy(
            QSizePolicy.Fixed, QSizePolicy.Expanding))
        self.descBox.textChanged.connect(self.input_description)
        self.descBox.setToolTip('Description text is not used in modeling.')
        settingBox.addWidget(self.descBox)

        settingBox.addWidget(QLabel(
            "<center><b>Substrate</b></center>"))
        self.inputSubstrateBox = QComboBox()
        self.inputSubstrateBox.setMaximumWidth(width)
        self.inputSubstrateBox.addItems(QCMaterial.keys())
        self.inputSubstrateBox.currentIndexChanged[str].connect(
            self.input_substrate)
        self.inputSubstrateBox.setToolTip(
            'The substrate decides possible choices of materials.')
        settingBox.addWidget(self.inputSubstrateBox)

        settingBox.addWidget(QLabel(
            '<center><b><i>E<sub>field</sub></i></b></center>'))
        self.inputEFieldBox = QDoubleSpinBox()
        self.inputEFieldBox.setMaximumWidth(width)
        self.inputEFieldBox.setDecimals(1)
        self.inputEFieldBox.setSuffix(' kV/cm')
        self.inputEFieldBox.setRange(-250.0, 250.0)
        self.inputEFieldBox.valueChanged[float].connect(self.input_EField)
        self.inputEFieldBox.setToolTip(
            'The Bias setting for the structure. '
            'Negative values may results in artifact results.')
        settingBox.addWidget(self.inputEFieldBox)

        settingBox.addWidget(QLabel(
            '<center><b>Position<br>Resolution</b></center>'))
        self.inputxResBox = QDoubleSpinBox()
        self.inputxResBox.setMaximumWidth(width)
        self.inputxResBox.setDecimals(2)
        self.inputxResBox.setRange(0.01, 1.0)
        self.inputxResBox.setSingleStep(0.1)
        self.inputxResBox.setSuffix(' \u212B')
        self.inputxResBox.valueChanged[float].connect(self.input_xres)
        settingBox.addWidget(self.inputxResBox)

        self.eResLabel = QLabel(
            '<center><b>Energy<br>Resolution</b></center>')
        settingBox.addWidget(self.eResLabel)
        self.inputEresBox = QDoubleSpinBox()
        self.inputEresBox.setMaximumWidth(width)
        self.inputEresBox.setDecimals(2)
        self.inputEresBox.setRange(0.0, 10.0)
        self.inputEresBox.setSingleStep(0.1)
        self.inputEresBox.setSuffix(' meV')
        self.inputEresBox.valueChanged[float].connect(self.input_Eres)
        self.inputEresBox.setToolTip(
            'The energy resolution used for ODE solver root finding. '
            'This number being too large may results in miss of states. \n'
            'The value is not used for matrix solver.')
        settingBox.addWidget(self.inputEresBox)
        self.eCountLabel = QLabel(
            '<center><b>No of States<br>Per Period</b></center>')
        settingBox.addWidget(self.eCountLabel)
        self.inputECountBox = QSpinBox()
        self.inputECountBox.setMaximumWidth(width)
        self.inputECountBox.setRange(1, 100)
        self.inputECountBox.valueChanged.connect(self.input_Ecount)
        self.inputECountBox.setToolTip(
            'Increase this number if you believe the number of states solved'
            'is not enough. This is only used for matrix solver.'
        )
        settingBox.addWidget(self.inputECountBox)
        self.algoParamUpdate()

        settingBox.addWidget(QLabel(
            '<center><b>Repeats</b></center>'))
        self.inputRepeatsBox = QSpinBox()
        self.inputRepeatsBox.setMaximumWidth(width)
        self.inputRepeatsBox.setRange(1, 5)
        self.inputRepeatsBox.valueChanged[int].connect(self.input_repeats)
        settingBox.addWidget(self.inputRepeatsBox)

        settingBox.addWidget(QLabel(
            '<center><b>Wavelength</b></center>'))
        self.inputWlBox = QDoubleSpinBox()
        self.inputWlBox.setMaximumWidth(width)
        self.inputWlBox.setDecimals(1)
        self.inputWlBox.setRange(1.5, 40)
        self.inputWlBox.setSingleStep(1)
        self.inputWlBox.setSuffix(' μm')
        self.inputWlBox.valueChanged[float].connect(self.input_wl)
        settingBox.addWidget(self.inputWlBox)

        # Basis solver divider setting
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
        LpLayoutGroupBox.setFixedWidth(width)
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
            QSizePolicy(QSizePolicy.Fixed, QSizePolicy.Preferred))
        self.LpStringBox.setMinimumHeight(30)
        # self.LpStringBox.setMaximumHeight(120)
        LpLayout = QGridLayout()
        LpLayout.addWidget(QLabel('<b>first</b>'), 0, 0)
        LpLayout.addWidget(QLabel('<b>last</b>'), 0, 1)
        LpLayout.addWidget(self.LpFirstSpinbox, 1, 0)
        LpLayout.addWidget(self.LpLastSpinbox, 1, 1)
        LpLayout.addWidget(self.LpStringBox, 2, 0, 1, 2)
        LpLayout.setSpacing(1)
        LpLayoutGroupBox.setLayout(LpLayout)
        settingBox.addWidget(LpLayoutGroupBox)

        settingBox.addStretch()
        return settingBox
        # _generateSettingBox end

    def _generateLayerBox(self, width):
        """ Return a Qt Layout object containing all layer parameters """
        layerBox = QGridLayout()
        layerBox.setSpacing(5)
        self.insertLayerAboveButton = QPushButton("Insert Layer")
        self.insertLayerAboveButton.clicked.connect(self.insert_layerAbove)
        layerBox.addWidget(self.insertLayerAboveButton, 0, 0)
        self.deleteLayerButton = QPushButton("Delete Layer")
        self.deleteLayerButton.clicked.connect(self.delete_layer)
        layerBox.addWidget(self.deleteLayerButton, 0, 1)
        self.optimizeLayerButton = QPushButton("Optimize Layer")
        self.optimizeLayerButton.clicked.connect(self.optimize_layer)
        self.optimizeLayerButton.setEnabled(False)
        self.globalOptimizeButton = QPushButton("Global Optimize")
        self.globalOptimizeButton.clicked.connect(self.global_optimize)
        self.globalOptimizeButton.setEnabled(False)
        layerBox.addWidget(self.optimizeLayerButton, 1, 0)
        layerBox.addWidget(self.globalOptimizeButton, 1, 1)

        # set up layerTable
        self.layerTable = QTableWidget()
        self.layerTable.setSelectionBehavior(QTableWidget.SelectRows)
        self.layerTable.setSelectionMode(QTableWidget.SingleSelection)
        self.layerTable.setFixedWidth(width)
        self.layerTable.setSizePolicy(
            QSizePolicy(QSizePolicy.Fixed, QSizePolicy.Expanding))
        self.layerTable.itemChanged.connect(self.layerTable_itemChanged)
        self.layerTable.itemSelectionChanged.connect(
            self.layerTable_itemSelectionChanged)
        layerBox.addWidget(self.layerTable, 2, 0, 1, 2)

        return layerBox
        # _generateLayerBox end

    def _generateFigureBox(self):
        """Return:
            A Qt Layout containing a canvas for plotting;
            A QGridLayout containing plot control;"""
        self.quantumCanvas = EJcanvas(xlabel='Position (Å)',
                                      ylabel='Energy (eV)', parent=self)
        self.plotControl = EJplotControl(self.quantumCanvas, self)
        self.quantumCanvas.setSizePolicy(QSizePolicy.Expanding,
                                         QSizePolicy.Expanding)
        figureBox = QVBoxLayout()
        figureBox.addWidget(self.quantumCanvas)

        self.zoomInButton = QPushButton("Zoom")
        self.plotControl.set_action('zoom', self.zoomInButton)
        self.zoomOutButton = QPushButton("Reset")
        self.plotControl.set_action('home', self.zoomOutButton)
        self.panButton = QPushButton("Pan")  # to move
        self.plotControl.set_action('pan', self.panButton)
        self.layerSelectButton = QPushButton("Layer Select")
        self.plotControl.set_custom('layerselect', self.layerSelectButton,
                                    self.layer_select)
        self.layerSelectButton.clicked[bool].connect(self.layerSelectMode)
        self.clearWFsButton = QPushButton("Clear")
        self.clearWFsButton.clicked.connect(self.clear_WFs)
        plotControlGrid = QGridLayout()
        plotControlGrid.addWidget(self.layerSelectButton, 0, 0, 1, 2)
        plotControlGrid.addWidget(self.zoomInButton, 1, 0, 1, 1)
        plotControlGrid.addWidget(self.zoomOutButton, 1, 1, 1, 1)
        plotControlGrid.addWidget(self.panButton, 2, 0, 1, 1)
        plotControlGrid.addWidget(self.clearWFsButton, 2, 1, 1, 1)
        plotControlGrid.setSpacing(5)

        return figureBox, plotControlGrid
        # _generateFigureBox end

    def _generateSolveBox(self, plotControlGrid, width):
        """ Return a Qt Layout containing material information,
        eigensolve control, states properties calculation and plot control"""
        self.solveBox = QVBoxLayout()
        self.solveBasisButton = QPushButton("Solve Basis")
        self.solveBasisButton.clicked.connect(self.solve_basis)
        self.solveBox.addWidget(self.solveBasisButton)
        self.solveWholeButton = QPushButton("Solve Whole")
        self.solveWholeButton.clicked.connect(self.solve_whole)
        self.solveBox.addWidget(self.solveWholeButton)

        # set up material composition inputs
        self.mtrlTable = QTableWidget()
        #  self.mtrlTable.setSelectionBehavior(QTableWidget.SelectItems)
        #  self.mtrlTable.setSelectionMode(QTableWidget.NoSelection)
        self.mtrlTable.setSelectionBehavior(QTableWidget.SelectRows)
        self.mtrlTable.setSelectionMode(QTableWidget.SingleSelection)
        self.mtrlTable.setMaximumWidth(width)
        self.mtrlTable.setMinimumWidth(width)
        self.mtrlTable.setMinimumHeight(100)
        self.mtrlTable.setMaximumHeight(120)
        self.mtrlTable.setSizePolicy(
            QSizePolicy(QSizePolicy.Fixed, QSizePolicy.Minimum))
        self.mtrlTable.itemChanged.connect(self.mtrlTable_itemChanged)
        self.mtrlTable.itemSelectionChanged.connect(
            self.mtrlTable_itemSelectionChanged)
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
        mtrlGrid.setSpacing(5)
        mtrlGroupBox = QGroupBox("Materials")
        mtrlGroupBox.setLayout(mtrlGrid)
        self.solveBox.addWidget(mtrlGroupBox)

        # IFR setting
        self.ifrDefBox = QCheckBox('Constant IFR')
        self.ifrDefBox.stateChanged.connect(self.ifrDefConstant)
        self.ifrDeltaBox = QDoubleSpinBox()
        self.ifrDeltaBox.setMaximumWidth(width)
        self.ifrDeltaBox.setDecimals(1)
        self.ifrDeltaBox.setSuffix(' Å')
        self.ifrDeltaBox.setRange(0.0, 9999.0)
        self.ifrDeltaBox.valueChanged[float].connect(self.input_ifrDelta)
        self.ifrLambdaBox = QDoubleSpinBox()
        self.ifrLambdaBox.setMaximumWidth(width)
        self.ifrLambdaBox.setDecimals(1)
        self.ifrLambdaBox.setSuffix(' Å')
        self.ifrLambdaBox.setRange(0.0, 9999.0)
        ifrGrid = QGridLayout()
        ifrGrid.addWidget(self.ifrDefBox, 0, 0, 1, 2)
        ifrGrid.addWidget(QLabel('<center>IRF Δ</center>'), 1, 0)
        ifrGrid.addWidget(self.ifrDeltaBox, 2, 0)
        ifrGrid.addWidget(QLabel('<center>IRF Λ</center>'), 1, 1)
        ifrGrid.addWidget(self.ifrLambdaBox, 2, 1)
        self.ifrLambdaBox.valueChanged[float].connect(self.input_ifrLambda)
        # This is exposed s.t. it can be hidden when IFR is off.
        self._ifrGroupBox = QGroupBox('Interface Roughness')
        self._ifrGroupBox.setLayout(ifrGrid)
        self.solveBox.addWidget(self._ifrGroupBox)

        # set up plot control inputs
        plotControlGroupBox = QGroupBox("Plot Controls")
        plotControlGroupBox.setLayout(plotControlGrid)
        self.solveBox.addWidget(plotControlGroupBox)

        # set up Calculate controls
        self.pairSelectButton = QPushButton("Pair Select")
        self.pairSelectButton.setEnabled(False)
        self.plotControl.set_custom('pairselect', self.pairSelectButton,
                                    self.state_pick)
        self.pairSelectButton.clicked[bool].connect(self.pairSelectMode)
        self.FoMButton = QPushButton("FoM")
        self.FoMButton.setEnabled(False)
        self.FoMButton.clicked.connect(self.updateFoM)
        self.fullPopulationButton = QPushButton("Population")
        self.fullPopulationButton.setEnabled(False)
        self.fullPopulationButton.clicked.connect(self.fullPopulation)
        self.toOpticsButton = QPushButton("->Optics")
        self.toOpticsButton.setEnabled(False)
        self.gainSpecButton = QPushButton("Gain Spec")
        self.gainSpecButton.setEnabled(False)
        self.gainSpecButton.clicked.connect(self.popGainSpecWindow)
        # signal is processed in ErwinJr main window
        self.stateParamText = QTextEdit('')
        self.stateParamText.setReadOnly(True)
        self.stateParamText.setMaximumWidth(width)
        self.stateParamText.setMinimumWidth(width)
        self.stateParamText.setMinimumHeight(150)
        self.stateParamText.setSizePolicy(QSizePolicy(
            QSizePolicy.Fixed, QSizePolicy.Expanding))
        self.stateParamText.textChanged.connect(
            lambda: self.toOpticsButton.setEnabled(False))
        calculateControlGrid = QGridLayout()
        calculateControlGrid.addWidget(self.pairSelectButton, 0, 0, 1, 2)
        calculateControlGrid.addWidget(self.FoMButton, 1, 0, 1, 1)
        calculateControlGrid.addWidget(self.fullPopulationButton, 1, 1, 1, 1)
        calculateControlGrid.addWidget(self.toOpticsButton, 2, 0, 1, 1)
        calculateControlGrid.addWidget(self.gainSpecButton, 2, 1, 1, 1)
        calculateControlGrid.addWidget(self.stateParamText, 3, 0, 1, 2)
        calculateControlGrid.setSpacing(5)
        # TODO: voltage efficiency, hbar omega / field * Lp
        calculateControlGroupBox = QGroupBox("Calculate")
        calculateControlGroupBox.setLayout(calculateControlGrid)
        self.solveBox.addWidget(calculateControlGroupBox)

        self.solveBox.addStretch()
        return self.solveBox
        # _generateSolveBox end

    def reload(self):
        """ Reload everything from qclayers, typically when it's changed
        externally or loaded"""
        # TODO: reconsider when reload is needed.
        self.qclayers.populate_x()
        self._update_settingBox()
        self._update_mtrlList()
        self._update_mtrlTable()
        self._update_layerTable()
        self._update_Lp_limits()
        self.update_Lp_box()
        self.layerTable.clearSelection()
        self.update_IFR_settings()
        self.IFR_settings_checkConst()
        self.update_quantumCanvas()

    def _update_mtrlList(self):
        """Update self.mtrlList according to self.qclayers material
        information. mtrlList is used for mtrl column in layerTable"""
        self.mtrlList = []
        for n, mtrl in enumerate(self.qclayers.materials):
            name = AParam[mtrl]['name']
            name = name.replace("1-x", str(1-self.qclayers.moleFracs[n]))
            name = name.replace("x", str(self.qclayers.moleFracs[n]))
            name = "#%d " % (n+1)  # + name
            self.mtrlList.append(name)

# ===========================================================================
# SettingBox Controls
# ===========================================================================
    def _update_settingBox(self):
        self.updating = True
        self.descBox.setText(self.qclayers.description)
        self.inputSubstrateBox.setCurrentText(self.qclayers.substrate)
        self.inputEFieldBox.setValue(self.qclayers.EField)
        self.inputxResBox.setValue(self.qclayers.xres)
        self.inputEresBox.setValue(self.qclayers.Eres)
        self.inputECountBox.setValue(self.qclayers.statePerRepeat)
        self.inputWlBox.setValue(self.qclayers.wl)
        self.inputRepeatsBox.setValue(self.qclayers.repeats)
        self.updating = False

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
                '%s substrates have not yet been implemented.' % substrateType)
            self.inputSubstrateBox.setCurrentIndex(
                self.inputSubstrateBox.findText(self.qclayers.substrate))
            return
        self.qclayers.populate_x()
        self.update_mtrl_info()
        self._update_layerTable()
        self._update_mtrlTable()
        self.update_quantumCanvas()
        self.dirty.emit()

    @pyqtSlot(float)
    @settingslot
    def input_EField(self, EField):
        """ SLOT connected to inputEFieldBox.valueChanged(double)
        update external E field in unit kV/cm """
        self.clear_WFs()
        self.qclayers.EField = EField
        self.qclayers.populate_x()
        self.dirty.emit()
        self.update_quantumCanvas()

    @pyqtSlot(float)
    @settingslot
    def input_xres(self, xres):
        """ SLOT connected to inputxResBox.valueChanged
        update position resolution (xres), in angstrom """
        self.clear_WFs()
        self.qclayers.xres = xres
        self.qclayers.populate_x()
        self.dirty.emit()
        self.update_quantumCanvas()

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
    def input_Ecount(self, count: int):
        self.qclayers.statePerRepeat = count
        self.qclayers.matrixEigenCount = count * self.qclayers.repeats
        self.dirty.emit()

    @pyqtSlot(int)
    @settingslot
    def input_repeats(self, repeat):
        """SLOT connected to SINGAL self.inputRepeatsBox.valueChanged(int)
        update number of repeats for the whole structure."""
        self.clear_WFs()
        self.qclayers.repeats = repeat
        self.qclayers.populate_x()
        self.dirty.emit()
        self.update_quantumCanvas()

    @pyqtSlot(float)
    @settingslot
    def input_wl(self, wl):
        self.qclayers.wl = wl

    @pyqtSlot()
    def input_basis(self):
        """SLOT connected to self.inputARInjectorCheck.stateChanged(int) and
        self.inputInjectorARCheck.stateChanged(int)
        update dividers info
        """
        self.qclayers.basisARInjector = self.inputARInjectorCheck.isChecked()
        self.qclayers.basisInjectorAR = self.inputInjectorARCheck.isChecked()

    def _update_Lp_limits(self):
        """Update Lp select range in the Period Info box (GUI)
        This should be called whenever number of layers is changed
        Will trigger LpFirstSpinbox/LpLastSpinBox.valueChanged(int) SIGNAL
        and thus update_Lp_box
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
        LpFirst = self.LpFirstSpinbox.value() - 1
        LpLast = self.LpLastSpinbox.value()
        # +1 because range is not inclusive of last value
        # total length of the layers (1 period)
        Lp = sum(self.qclayers.layerWidths[LpFirst:LpLast])
        Lp_string = "Lp: %.1f \u212B<br>" % Lp
        # average doping of the layers
        ns = sum(self.qclayers.layerDopings[n]
                 * self.qclayers.layerWidths[n]
                 for n in range(LpFirst, LpLast))
        if Lp == 0:
            Lp_string += ("n<sub>D</sub>: NA\u00D710<sup>17</sup>"
                          "cm<sup>-3</sup><br>")
        else:
            nD = ns / Lp
            Lp_string += ("N<sub>D</sub>: %6.3f\u00D710<sup>17</sup>"
                          "cm<sup>-3</sup><br>") % nD
        # 2D carrier density in 1E11cm-2
        ns = ns * 1e-2
        Lp_string += ("N<sub>s</sub>: %6.3f\u00D710<sup>11</sup>"
                      "cm<sup>-2</sup><br>") % ns
        Lp_string += ("n<sub>eff</sub>: %.2f" % self.qclayers.effective_ridx(
            self.inputWlBox.value()))
        self.LpStringBox.setText(Lp_string)

    @pyqtSlot()
    @settingslot
    def input_description(self):
        """ SLOT connected to self.descBox.textChanged()
        Change description string for the design"""
        self.qclayers.description = self.descBox.toPlainText()
        self.dirty.emit()

    def set_temperature(self, T):
        """A public method to wrap temperature settings"""
        self.qclayers.set_temperature(T)
        self.qclayers.populate_x()
        self.dirty.emit()
        self.update_mtrl_info()
        self.update_quantumCanvas()

# ===========================================================================
# Layer Table Control
# ===========================================================================
    def _update_layerTable(self):
        """Refresh layer table, called every time after data update"""
        # Block itemChanged SIGNAL while refreshing
        self.layerTable.blockSignals(True)
        self.layerTable.clear()
        self.layerTable.setColumnCount(5)
        # An extra blank line for adding new layers
        self.layerTable.setRowCount(len(self.qclayers.layerWidths) + 1)
        vertLabels = [str(n+1) for n in range(len(self.qclayers.layerWidths))]
        self.layerTable.setHorizontalHeaderLabels(
            ('Width', 'ML', 'Mtrl', 'AR', 'Doping'))
        self.layerTable.setVerticalHeaderLabels(vertLabels)

        # gray2 = QColor(230, 230, 230)  # for unchangeable background

        for q, layerWidths in enumerate(self.qclayers.layerWidths):
            color = self.mtrlcolors[
                self.qclayers.layerMtrls[q] % len(self.mtrlcolors)]
            # Width Setup
            width = QTableWidgetItem("%5.1f" % layerWidths)
            width.setTextAlignment(Qt.AlignCenter)
            width.setBackground(color)
            self.layerTable.setItem(q, 0, width)

            # "ML" number of monolayer Setup
            # TODO: monolayer <-> thickness should not include
            # temperature/strain correction
            mlThickness = self.qclayers.mtrlAlloys[
                self.qclayers.layerMtrls[q]].a_perp / 2
            # a_perp has two layer of atoms (one III and one V)
            numML = QTableWidgetItem("%5.1f" % (layerWidths / mlThickness))
            numML.setTextAlignment(Qt.AlignCenter)
            numML.setBackground(color)
            self.layerTable.setItem(q, 1, numML)

            # Material Setup
            #  mtrlWidget = QComboBox()
            mtrlWidget = mtrlComboBox()
            mtrlWidget.addItems(self.mtrlList)
            mtrlWidget.setCurrentIndex(self.qclayers.layerMtrls[q])
            mtrlWidget.currentIndexChanged.connect(
                partial(self.layerTable_materialChanged, q))
            self.layerTable.setCellWidget(q, 2, mtrlWidget)

            # Active Region Layer Setup
            ARitem = QTableWidgetItem()
            ARitem.setCheckState(
                Qt.Checked if self.qclayers.layerARs[q] == 1
                else Qt.Unchecked)
            ARitem.setTextAlignment(Qt.AlignCenter)
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
        if row >= len(self.qclayers.layerWidths) or row < 0:
            # Add new lines in the last layer
            row = len(self.qclayers.layerWidths)
            aR = self.qclayers.layerARs[row-1] and self.qclayers.layerARs[0]
            doping = self.qclayers.layerDopings[row-1]
            mtrlIdx = (self.qclayers.layerMtrls[row-1] + 1) % N
        else:
            aR = self.qclayers.layerARs[row] and self.qclayers.layerARs[row-1]
            doping = self.qclayers.layerDopings[row]
            mtrlIdx = (self.qclayers.layerMtrls[row-1] - 1) % N

        self.clear_WFs()
        self.qclayers.add_layer(row, 0.0, mtrlIdx, aR, doping)
        self.qclayers.populate_x()
        self._update_Lp_limits()
        self.update_Lp_box()
        self._update_layerTable()
        self.update_mtrl_info()
        self.layerTable.selectRow(row)
        # Trigger itemSelectionChanged SIGNAL and thus update_quantumCanvas
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

        self.clear_WFs()
        self.qclayers.del_layer(row)
        self.qclayers.populate_x()
        self._update_Lp_limits()
        self.update_Lp_box()
        self._update_layerTable()
        self.update_mtrl_info()
        self.layerTable.selectRow(row)
        # Trigger itemSelectionChanged SIGNAL and thus update_quantumCanvas
        self.dirty.emit()

    def _optimize_layer(self, n, upper, lower):
        optimize_layer(self.qclayers, n, upper, lower)
        self.qclayers.solve_whole()
        self.stateHolder = [upper, lower]
        self.pairSelected = True
        self._calcFoM()

    def showOptmize(self, n=None):
        self._update_layerTable()
        if n is None:
            self._updatePopulation()
        else:
            self.layerTable.setCurrentCell(n, 0)
            # self.updateSelected()
            self._updateFoM
            self.dirty.emit()

    @pyqtSlot()
    def optimize_layer(self):
        """SLOT connected to self.optimizeLayerButton.clicked()"""
        n = self.layerTable.currentRow()
        if self.qclayers.status == 'unsolved':
            QMessageBox.warning(self, ejWarning, "Solve the model first.")
            return
        if n < 0 or n > len(self.qclayers.layerWidths):
            QMessageBox.warning(self, ejError,
                                "Select the layer to optimize.")
            return
        try:
            upper, lower = self.stateHolder
            if self.qclayers.eigenEs[upper] < self.qclayers.eigenEs[lower]:
                upper, lower = lower, upper
        except ValueError:
            QMessageBox.warning(self, ejError,
                                "Select state pair to optimize.")
            return
        self.stateParamText.clear()
        self.clear_WFs()
        self._threadRun(lambda: self._optimize_layer(n, upper, lower),
                        lambda: self.showOptmize(n))

    def _global_optimize(self):
        optimize_global(self.qclayers)
        self._fullPopulation()

    @pyqtSlot()
    def global_optimize(self):
        """SLOT connect to self.globalOptimizeButton.clicked()"""
        message = ("This is an experimental feature "
                   "the result may not be stable and may take a long time."
                   "It's highly recommended to save the structure first.\n"
                   "Do you want to proceed?"
                   )
        if QMessageBox.question(self, ejWarning, message) == QMessageBox.Yes:
            self.stateParamText.clear()
            self.clear_WFs()
            self._threadRun(self._global_optimize, self.showOptmize)

    @pyqtSlot(QTableWidgetItem)
    def layerTable_itemChanged(self, item):
        """SLOT connected to layerTable.itemChanged(QTableWidgetItem*)
        Update layer profile after user input"""
        row = item.row()
        mtrlN = len(self.qclayers.materials)
        if row > len(self.qclayers.layerWidths):
            raise ValueError("Bad layer width input row number")
        column = item.column()
        if column in (0, 1, 4):
            try:
                value = float(item.text())
            except ValueError:
                # invalid input
                QMessageBox.warning(self, ejError,
                                    "This value should be a number")
                self._update_layerTable()
                return
            if column in (0, 1):
                if column == 0:  # column == 0 for Widths column
                    new_width = value
                else:  # column == 1 for "ML" number of monolayers
                    new_width = value * self.qclayers.mtrlAlloys[
                        self.qclayers.layerMtrls[row]].a_perp/2
                    new_width = round(new_width, 1)

                if row == len(self.qclayers.layerWidths):
                    # add row at end of list
                    AR = self.qclayers.layerARs[row-1]
                    doping = self.qclayers.layerDopings[row-1]
                    mtrlIndx = (self.qclayers.layerMtrls[row-1] + 1) % mtrlN
                    self.qclayers.add_layer(row, new_width, mtrlIndx, AR,
                                            doping)
                    row += 1  # used so that last (blank) row is again selected
                    self._update_Lp_limits()
                    self.update_Lp_box()
                else:  # change Width of selected row in-place
                    self.qclayers.layerWidths[row] = new_width
                    self.update_Lp_box()
            else:
                # column == 4 for item change in Doping column
                doping = value
                if row == len(self.qclayers.layerWidths):
                    AR = self.qclayers.layerARs[row-1]
                    mtrlIndx = (self.qclayers.layerMtrls[row-1] + 1) % mtrlN
                    self.qclayers.add_layer(row, 0.0, mtrlIndx, AR, doping)
                    row += 1  # used so that last (blank) row is again selected
                    self._update_Lp_limits()
                    self.update_Lp_box()
                else:
                    self.qclayers.layerDopings[row] = doping

        elif column == 2:
            # column == 2 for item change in mtrl column, should be
            # controlled by layerTable_materialChanged
            raise Exception("Should not be here")

        elif column == 3:
            # column == 3 for item change in AR column
            if row == len(self.qclayers.layerWidths):
                # don't do anything if row is last row
                return
            self.qclayers.layerARs[row] = (item.checkState() == Qt.Checked)

        else:
            raise Exception("Should not be here")

        self._update_layerTable()
        self.clear_WFs()
        self.layerTable.setCurrentCell(row, column)
        self.qclayers.populate_x()
        self.update_mtrl_info()
        self.update_quantumCanvas()
        self.dirty.emit()

    @pyqtSlot()
    def layerTable_itemSelectionChanged(self):
        """SLOT connected to layerTable.itemSelectionChanged()"""
        row = self.layerTable.currentRow()
        if row < len(self.qclayers.layerWidths):
            self.layerSelected = row
        else:
            self.layerSelected = None
            self.layerTable.clearSelection()
        self.update_quantumCanvas()

    def layerTable_materialChanged(self, row, selection):
        """SLOT as partial(self.layerTable_materialChanged, q)) connected to
        materialWidget.currentIndexChanged(int) """
        self.qclayers.layerMtrls[row] = selection
        self.clear_WFs()
        self.qclayers.populate_x()
        self.layerTable.selectRow(row)
        self.update_mtrl_info()
        self._update_layerTable()
        # Trigger itemSelectionChanged SIGNAL and thus update_quantumCanvas
        self.dirty.emit()

    @pyqtSlot()
    def rotate_layer(self):
        """Move last layer to first layer, SLOT open to public"""
        self.clear_WFs()
        self.qclayers.rotate_layer()
        self.qclayers.populate_x()
        self._update_layerTable()
        self.layerTable.setCurrentCell(1, 0)
        self.update_quantumCanvas()
        self.dirty.emit()

    @pyqtSlot()
    def invert_layer(self):
        """Invert the order of the layers"""
        self.clear_WFs()
        self.qclayers.invert_layer()
        self.qclayers.populate_x()
        self._update_layerTable()
        row = self.layerTable.currentRow()
        if row >= 0 and row < len(self.layerWidths):
            self.layerTable.setCurrentCell(
                len(self.qclayers.layerWidths)-1-row, 0)
        self.update_quantumCanvas()
        self.dirty.emit()

    @pyqtSlot()
    def ARonly(self):
        self.qclayers.basisAROnly = not self.qclayers.basisAROnly
        self.clear_WFs()

    @pyqtSlot()
    def copy_structure(self):
        clipboard = QApplication.clipboard()
        string = ''
        for width in self.qclayers.layerWidths:
            string += '%.1f\n' % width
        clipboard.setText(string)

# =========================================================================
# mtrlTable Control
# =========================================================================
    def _update_mtrlTable(self):
        """Set up material table, both format and consistency with qclayers"""
        self.mtrlTable.blockSignals(True)
        self.mtrlTable.clear()
        self.mtrlTable.setColumnCount(3)
        self.mtrlTable.setRowCount(len(self.qclayers.materials))
        self.mtrlTable.setHorizontalHeaderLabels(["#", "material", "x"])
        self.mtrlTable.verticalHeader().hide()
        # TODO: material name support

        possibleMtrl = tuple([AParam[m]['name'] for m in
                              QCMaterial[self.qclayers.substrate]])
        for n, mtrl in enumerate(self.qclayers.materials):
            color = self.mtrlcolors[n % len(self.mtrlcolors)]
            name = QTableWidgetItem(str(n+1))
            name.setTextAlignment(Qt.AlignCenter)
            name.setBackground(color)
            name.setFlags(Qt.NoItemFlags)
            self.mtrlTable.setItem(n, 0, name)

            # Choose from available materials, according to substrate
            #  mtrlItem = QComboBox()
            mtrlItem = mtrlComboBox()
            mtrlItem.addItems(possibleMtrl)
            mtrlItem.setCurrentText(AParam[mtrl]['name'])
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
        self.delMtrlButton.setEnabled(False)
        self.mtrlTable.blockSignals(False)

    @pyqtSlot()
    def mtrlTable_itemSelectionChanged(self):
        """SLOT connected to mtrlTable.itemSelectionChanged()"""
        self.delMtrlButton.setEnabled(len(self.qclayers.materials) > 2)
        self.update_IFR_settings()

    def update_IFR_settings(self):
        if self.qclayers.includeIFR:
            self._ifrGroupBox.setVisible(True)
        else:
            self._ifrGroupBox.setVisible(False)
        if self.qclayers.customIFR:
            self.ifrDefBox.setEnabled(False)
            self.ifrDeltaBox.setEnabled(False)
            self.ifrLambdaBox.setEnabled(False)
            # button to remove custom
        if self.ifrDefBox.isChecked():
            assert(all(x == self.qclayers.mtrlIFRDelta[0]
                       for x in self.qclayers.mtrlIFRDelta[1:]))
            assert(all(x == self.qclayers.mtrlIFRLambda[0]
                       for x in self.qclayers.mtrlIFRLambda[1:]))
            self.ifrDeltaBox.setValue(self.qclayers.mtrlIFRDelta[0])
            self.ifrLambdaBox.setValue(self.qclayers.mtrlIFRLambda[0])
        else:
            row = self.mtrlTable.currentRow()
            if (row < len(self.qclayers.materials)
                    and not self.qclayers.customIFR
                    and self.qclayers.mtrlIFRDelta is not None):
                self.ifrDeltaBox.setValue(self.qclayers.mtrlIFRDelta[row])
                self.ifrLambdaBox.setValue(self.qclayers.mtrlIFRLambda[row])

    def IFR_settings_checkConst(self):
        if (all(x == self.qclayers.mtrlIFRDelta[0]
                for x in self.qclayers.mtrlIFRDelta[1:]) and
            all(x == self.qclayers.mtrlIFRLambda[0]
                for x in self.qclayers.mtrlIFRLambda[1:])):
            self.ifrDefBox.setChecked(True)

    def mtrlTable_mtrlChanged(self, row, selection):
        """SLOT as partial(self.mtrlTable_mrtlChanged, q)) connected to
        mtrlItem.currentIndexChanged(int) """
        # not decorated by pyqt because it's not a real slot but a meta-slot
        self.qclayers.set_mtrl(row, QCMaterial[
            self.qclayers.substrate][selection])
        self.clear_WFs()
        self.qclayers.populate_x()
        self.update_mtrl_info()
        self.update_quantumCanvas()
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
        elif column == 2:
            # change "x" or mole fraction
            try:
                mf = float(item.text())
                assert(mf <= 1)
                self.qclayers.set_mtrl(row, moleFrac=int(100*mf)/100)
            except (ValueError, AssertionError):
                # Input is not a number or is larger than 1
                QMessageBox.warning(self, ejError, "invalid mole Fraction")
            #  item.setText(str(self.qclayers.moleFracs[row]))
            self.clear_WFs()
            self.qclayers.populate_x()
            self.update_mtrl_info()
            self.update_quantumCanvas()
            self.dirty.emit()
        else:
            # Should never be here
            raise ValueError

    @pyqtSlot()
    def add_mtrl(self):
        """SLOT connected to self.addMtrlButton.clicked()"""
        self.qclayers.add_mtrl()
        self._update_mtrlList()
        self._update_mtrlTable()
        self._update_layerTable()

    @pyqtSlot()
    def del_mtrl(self):
        """SLOT connected to self.delMtrlButton.clicked()"""
        self.qclayers.del_mtrl(self.mtrlTable.currentRow())
        self._update_mtrlList()
        self._update_mtrlTable()
        self._update_layerTable()

    @pyqtSlot(int)
    def ifrDefConstant(self, state):
        # state is int from Qt definition but here it's used as a bool
        if state:
            m = len(self.qclayers.materials)
            self.qclayers.mtrlIFRLambda = [self.ifrLambdaBox.value()] * m
            self.qclayers.mtrlIFRDelta = [self.ifrDeltaBox.value()] * m
        self.dirty.emit()

    @pyqtSlot(float)
    def input_ifrLambda(self, ifrLambda):
        if self.ifrDefBox.isChecked():
            m = len(self.qclayers.materials)
            self.qclayers.mtrlIFRLambda = [ifrLambda] * m
        else:
            mn = self.mtrlTable.currentRow()
            self.qclayers.mtrlIFRLambda[mn] = ifrLambda
        self.qclayers.populate_material()
        if self.qclayers.status.startswith('solved'):
            # following also do self.qclayers.status = 'solved'
            self.qclayers.reset_IFR_cache()
            self.gainSpecButton.setEnabled(False)
            self.update_quantumCanvas()
        self.dirty.emit()

    @pyqtSlot(float)
    def input_ifrDelta(self, ifrDelta):
        if self.ifrDefBox.isChecked():
            m = len(self.qclayers.materials)
            self.qclayers.mtrlIFRDelta = [ifrDelta] * m
        else:
            mn = self.mtrlTable.currentRow()
            self.qclayers.mtrlIFRDelta[mn] = ifrDelta
        self.qclayers.populate_material()
        if self.qclayers.status.startswith('solved'):
            # following also do self.qclayers.status = 'solved'
            self.qclayers.reset_IFR_cache()
            self.gainSpecButton.setEnabled(False)
            self.update_quantumCanvas()
        self.dirty.emit()

    def update_mtrl_info(self):
        """Update labels below mtrlTable"""
        self.offsetLabel.setText(
            '<center>ΔE<sub>c</sub>: <b>%6.0f meV </b></center>' %
            (self.qclayers.mtrlOffset * 1000))
        self.netStrainLabel.setText(
            "<center>Net Strain: <b>%6.3f%%</b></center>" %
            self.qclayers.netStrain)
        self.LOPhononLabel.setText(
            "<center>E<sub>LO</sub>: <b>%4.1f meV</b></center>" %
            (1000*self.qclayers.avghwLO))

# =========================================================================
# Quantum Tab Plotting and Plot Control
# =========================================================================
    def update_quantumCanvas(self):
        """Update the canvas to show band diagram, and if states have been
        solved, draw the wavefunctions"""
        if self.plotControl.zoomed:
            xmin, xmax = self.quantumCanvas.axes.get_xlim()
            ymin, ymax = self.quantumCanvas.axes.get_ylim()
        else:
            xmin = xmax = ymin = ymax = None
        self.quantumCanvas.clear()
        axes = self.quantumCanvas.axes
        plotPotential(self.qclayers, axes, self.plotVL, self.plotVX,
                      self.plotLH, self.plotSO)

        if self.layerSelected is not None:
            selectedVc = np.ma.masked_where(
                self.qclayers.xLayerMask(self.layerSelected),
                self.qclayers.xVc)
            axes.plot(self.qclayers.xPoints, selectedVc, 'b',
                      linewidth=1.5 if self.qclayers.layerARs[
                          self.layerSelected] else 1)

        if self.qclayers.status in ('basis', 'solved', 'solved-full'):
            self.wfs = plotWF(self.qclayers, self.plotType,
                              self.fillPlot, self.stateHolder,
                              axes=axes)

        if xmin is not None:
            axes.set_xlim(xmin, xmax)
            axes.set_ylim(ymin, ymax)
        else:
            ymin, ymax = axes.get_ylim()
            if ymax - ymin < 0.4:
                axes.set_ylim(ymin-0.2, ymax+0.2)
        self.quantumCanvas.draw()

    @pyqtSlot(bool)
    def layerSelectMode(self, checked):
        """ SLOT connected to layerSelectButton.clicked[bool]
        to enable layerSelect mode in the canvas.
        """
        self.plotControl.trigger_custom('layerselect')
        if not checked:
            # clear selection
            self.layerSelected = None
            self.layerTable.clearSelection()
            self.update_quantumCanvas()

    @pyqtSlot(bool)
    def pairSelectMode(self, checked):
        """ SLOT connected to pairSelectButton.clicked[bool]
        to enable pairSelect mode in the canvas.
        """
        self.plotControl.trigger_custom('pairselect')
        if not checked:
            self.pairSelected = False
            self.stateHolder = []
            self.update_quantumCanvas()

    def layer_select(self, event):
        """callback registered in plotControl when it's in wellselect mode.
        It's mpl_connect to button_release_event of quantumCanvas """
        if event.button == 1:  # left button clicked
            x = event.xdata
            xLayerNum = np.argmin((self.qclayers.xPoints - x)**2)
            layerNum = self.qclayers.xLayerNums[xLayerNum]
            self.layerTable.selectRow(layerNum)
            # Trigger itemSelectionChanged SIGNAL and thus update_quantumCanvas

    def clear_WFs(self):
        self.qclayers.status = 'unsolved'
        self.stateHolder = []
        self.pairSelectButton.setEnabled(False)
        self.fullPopulationButton.setEnabled(False)
        self.gainSpecButton.setEnabled(False)
        self.optimizeLayerButton.setEnabled(False)
        self.globalOptimizeButton.setEnabled(False)
        self.update_quantumCanvas()

# ===========================================================================
# Export Functions
# ===========================================================================
    def export_quantumCanvas(self, filename=None):
        self.plotControl.save_figure(
            "ErwinJr2 - Export Band Structure Image", filename, 'png')

    def export_band_data(self, fname):
        fnameBase = fname.split('.')
        fnameBase = ''.join(fnameBase[:-1]) if len(fnameBase) > 1 else fname
        np.savetxt(fnameBase + '_CB' + '.csv', np.column_stack([
            self.qclayers.xPoints, self.qclayers.xVc]), delimiter=',')

        if self.qclayers.status.startswith('solved'):
            # otherwise band structure hasn't been solved yet
            np.savetxt(fnameBase + '_WFs' + '.csv', np.column_stack([
                self.qclayers.xPoints, self.qclayers.psis.T]), delimiter=',')
            ys = scaleWF(self.qclayers, self.plotType).T+self.qclayers.eigenEs
            np.savetxt(fnameBase + '.csv', np.column_stack([
                self.qclayers.xPoints, ys]), delimiter=',')
            np.savetxt(fnameBase + '_Es' + '.csv',
                       self.qclayers.eigenEs, delimiter=',')
            if self.qclayers.status == 'solved-full':
                pop = np.array([
                    self.qclayers.state_population(n)
                    if self.qclayers.periodMap[n] is not None else np.NaN
                    for n in range(len(self.qclayers.eigenEs))
                ])
                np.savetxt(fnameBase + '_population' + '.csv',
                           pop, delimiter=',')

# ===========================================================================
# View Band Items
# ===========================================================================
    def view_Band(self, band):
        # TODO: combine following slot to this one
        self.bandFlags[band] = not self.bandFlags[band]
        self.update_quantumCanvas()

    def view_VXBand(self):
        self.plotVX = not self.plotVX
        self.update_quantumCanvas()

    def view_VLBand(self):
        self.plotVL = not self.plotVL
        self.update_quantumCanvas()

    def view_LHBand(self):
        self.plotLH = not self.plotLH
        self.update_quantumCanvas()

    def view_SOBand(self):
        self.plotSO = not self.plotSO
        self.update_quantumCanvas()

    def set_plotwf(self):
        self.plotType = 'wf' if self.plotType != 'wf' else 'mode'
        self.update_quantumCanvas()

    def set_fill(self):
        if self.fillPlot:
            self.fillPlot = False
        else:
            self.fillPlot = 0.3
        self.update_quantumCanvas()

# ===========================================================================
# Model triggering
# ===========================================================================
    def triggerIFR(self):
        self.qclayers.includeIFR = not self.qclayers.includeIFR
        # TODO: dynamically hide IFR block
        if self.qclayers.includeIFR:
            self._ifrGroupBox.setVisible(True)
        else:
            self._ifrGroupBox.setVisible(False)

    def triggerSolver(self, solver):
        if solver != self.qclayers.solver:
            self.qclayers.solver = solver
            self.clear_WFs()
            self.algoParamUpdate()

    def algoParamUpdate(self):
        # This is causing "has active key-value observers (KVO)" when recreate
        # window warning, but we don't see any real issue.
        if self.qclayers.solver == 'matrix':
            self.eResLabel.setVisible(False)
            self.inputEresBox.setVisible(False)
            self.eCountLabel.setVisible(True)
            self.inputECountBox.setVisible(True)
        else:
            self.eResLabel.setVisible(True)
            self.inputEresBox.setVisible(True)
            self.eCountLabel.setVisible(False)
            self.inputECountBox.setVisible(False)

# ===========================================================================
# Calculations
# ===========================================================================
    def _threadRun(self, calc, postCalc=None):
        """ This is a helper function for multiple threading to avoid GUI
        freeze."""
        # TODO: block all change
        self.calcRepaint(True)
        assert(not self.calcThread.isRunning())
        self._worker = CalculateHolder(calc)
        self._worker.moveToThread(self.calcThread)
        self.calcThread.started.connect(self._worker.run)
        self._worker.finished.connect(self.calcThread.quit)
        self._worker.finished.connect(self._worker.deleteLater)
        self._worker.failed.connect(lambda traceInfo: QMessageBox.warning(
            self, ejError, traceInfo))
        self._worker.warning.connect(lambda message: QMessageBox.warning(
            self, ejWarning, message))
        self._worker.failed.connect(self.clear_WFs)
        self.calcThread.finished.connect(lambda: self.calcRepaint(False))
        # self.calcThread.finished.connect(self.thread.deleteLater)
        if postCalc is not None:
            self._worker.succeed.connect(postCalc)
        self.calcThread.start()

    def calcRepaint(self, is_doing):
        """UI repaint for doing calculating """
        for button in (self.solveWholeButton, self.solveBasisButton,
                       self.pairSelectButton):
            button.setEnabled(not is_doing)
            button.repaint()
        if is_doing:
            self.toOpticsButton.setEnabled(False)
        if self.pairSelected:
            self.FoMButton.setEnabled(not is_doing)
            self.FoMButton.repaint()
        if self.qclayers.status.startswith('solved'):
            self.fullPopulationButton.setEnabled(not is_doing)
            self.fullPopulationButton.repaint()

    def _solve_whole(self):
        self.qclayers.solve_whole()
        try:
            self.periodSet = self.qclayers.period_recognize()
        except StateRecognizeError as e:
            self.periodSet = set()
            raise e

    def _solved(self):
        self.optimizeLayerButton.setEnabled(True)
        self.globalOptimizeButton.setEnabled(True)
        self.update_quantumCanvas()

    @pyqtSlot()
    def solve_whole(self):  # solves whole structure
        """SLOT connected to solveWholeButton.clicked(): Whole solver """
        self.clear_WFs()
        self._threadRun(self._solve_whole, self._solved)

    def _solve_basis(self):
        self.qclayers.solve_basis()

    @pyqtSlot()
    def solve_basis(self):  # solves structure with basis
        """SLOT connected to solveBasisButton.clicked(): Basis solver """
        self.clear_WFs()
        self._threadRun(self._solve_basis, self.update_quantumCanvas)

    def state_pick(self, event):
        """Callback registered in plotControl when it's in pairselect mode.
        It's mpl_connect to button_release_event of quantumCanvas """
        if self.qclayers.status == 'unsolved':
            # Not yet solved
            raise ValueError('Structure not solved yet.')

        if event.button == 1:  # left button clicked: select a state
            if len(self.stateHolder) >= 2:
                # start new pair selection
                self.stateHolder = []
                self.pairSelected = False

            # TODO: replace following by matplotlib lines picker
            x = event.xdata
            y = event.ydata
            xData = np.tile(self.qclayers.xPoints,
                            (self.qclayers.psis.shape[0], 1)).T
            if self.plotType in ("mode", "wf"):
                yData = self.qclayers.eigenEs + self.wfs.T
            else:
                yData = self.qclayers.eigenEs

            width, height = self.quantumCanvas.axes.bbox.size
            xmin, xmax = self.quantumCanvas.axes.get_xlim()
            ymin, ymax = self.quantumCanvas.axes.get_ylim()
            xScale = (xmax - xmin) / width
            yScale = (ymax - ymin) / height

            r = np.nanmin(sqrt(((xData - x) / xScale)**2 +
                               ((yData - y) / yScale)**2), axis=0)
            # penalty for x axis away from visible region
            if self.plotType in ("mode", "wf"):
                x0 = np.argmax(
                    np.abs(self.wfs) > plotconfig["wf_almost_zero"], axis=1)
                xL = self.wfs.shape[1] - np.argmax(np.abs(
                    self.wfs[:, ::-1]) > plotconfig["wf_almost_zero"], axis=1)
                outOfSignt = np.arange(self.wfs.shape[0])[
                    (x0*self.qclayers.xres > x+5) |
                    (xL*self.qclayers.xres < x-5)
                    ]
                r[outOfSignt] += np.inf
            ss = np.nanargmin(r)
            if len(self.stateHolder) == 1 and self.stateHolder[0] == ss:
                r[ss] = np.nan
                ss = np.nanargmin(r)
            self.stateHolder.append(ss)
        elif event.button == 3:  # right button clicked: remove last selection
            if len(self.stateHolder) == 2:
                self.pairSelected = False
            self.stateHolder.pop()
        else:
            return
        self.update_quantumCanvas()

        if len(self.stateHolder) == 1:
            self.stateParamText.clear()
            self.pairString = "selected: %d, ..<br>" % self.stateHolder[0]
            if self.qclayers.status == 'solved-full':
                pop = self.qclayers.state_population(self.stateHolder[0])
                if pop is not None:
                    self.pairString += 'population: %.1f%%<br>' % (pop*100)
            self.stateParamText.setText(self.pairString)

        elif len(self.stateHolder) == 2:
            self.pairSelected = True
            self.updateSelected()

    def updateSelected(self):
        """Update easy-to-calculate parameters regarding picked two states"""
        self.FoMButton.setEnabled(True)
        upper, lower = self.stateHolder
        if self.qclayers.eigenEs[upper] < self.qclayers.eigenEs[lower]:
            upper, lower = lower, upper
        self.de = self.qclayers.eigenEs[upper] - self.qclayers.eigenEs[lower]
        self.wavelength = h * c0 / (e0 * self.de) * 1e6  # um

        if self.qclayers.status == 'unsolved':
            self.FoMButton.setEnabled(False)
            self.pairString = ''
        else:
            opticalDipole = self.qclayers.dipole(upper, lower)
            tauLO_ul = 1 / self.qclayers.lo_transition(upper, lower)
            if self.qclayers.includeIFR:
                self.tauIFR_ul = 1 / self.qclayers.ifr_transition(upper, lower)
            self.pairString = (
                "selected: %d, %d<br>"
                "energy diff: <b>%6.1f meV</b> (%6.1f \u00B5m)<br>"
                "dipole: <b>%6.1f \u212B</b><br>") % (
                    self.stateHolder[0], self.stateHolder[1],
                    1000*self.de, self.wavelength, opticalDipole
                )
            self.pairString += "LO scatter: %6.3g ps<br>" % tauLO_ul
            if self.qclayers.includeIFR:
                self.pairString += (
                    "IFR scatter: %6.3g ps<br>" % self.tauIFR_ul)
            if self.qclayers.status == 'solved-full':
                self.pairString += 'population: <br>&nbsp;&nbsp;&nbsp;'
                uPop = self.qclayers.state_population(upper)
                uPop = 'N/A' if uPop is None else '%.1f%%' % (100*uPop)
                lPop = self.qclayers.state_population(lower)
                lPop = 'N/A' if lPop is None else '%.1f%%' % (100*lPop)
                self.pairString += '%d: %s %d: %s<br>' % (
                    upper, uPop, lower, lPop)
        self.stateParamText.clear()
        self.stateParamText.setText(self.pairString)

    def _calcFoM(self):
        upper, lower = self.stateHolder
        if self.qclayers.eigenEs[upper] < self.qclayers.eigenEs[lower]:
            upper, lower = lower, upper
        tauLO_u = self.qclayers.lo_lifetime(upper)
        tauLO_l = self.qclayers.lo_lifetime(lower)
        if self.qclayers.includeIFR:
            tauIFR_u = self.qclayers.ifr_lifetime(upper)
            tauIFR_l = self.qclayers.ifr_lifetime(lower)
        # tau_u = 1/(1/tauLO_u + 1/tauIFR_u)
        # tau_l = 1/(1/tauLO_l + 1/tauIFR_l)
        FoM = self.qclayers.figure_of_merit(upper, lower)
        self.neff = self.qclayers.effective_ridx(self.wavelength)
        # 1E-27 = 1E-32 * 1E5
        # 1E-32 angstrom^2 ps -> m^2 s, 1E5 m/A -> cm/kA
        gamma = self.qclayers.dephasing(upper, lower)
        self.gainCoeff = e0 * FoM * 1E-27 * self.de / (
            gamma * hbar * c0 * eps0 * self.neff * self.qclayers.periodL*1E-10)
        # tauUpperLower is the inverse of transition rate (lifetime)

        self.FoMString = (
            "<i>\u03C4<sup>LO</sup><sub>upper</sub></i> : %6.3f ps<br>"
            "<i>\u03C4<sup>LO</sup><sub>lower</sub></i> : %6.3f ps<br>"
            ) % (tauLO_u, tauLO_l)
        if self.qclayers.includeIFR:
            self.FoMString += (
                "<i>\u03C4<sup>IFR</sup><sub>upper</sub></i> : %6.3f ps<br>"
                "<i>\u03C4<sup>IFR</sup><sub>lower</sub></i> : %6.3f ps<br>"
            ) % (tauIFR_u, tauIFR_l)
        self.FoMString += (
            "FoM: <b>%6.0f ps \u212B<sup>2</sup></b><br>"
            "Gain coefficient:<br>&nbsp;&nbsp; %.2f cm/kA"
        ) % (FoM, self.gainCoeff)

    def _updateFoM(self):
        self.stateParamText.setText(self.pairString + self.FoMString)
        self.stateParamText.verticalScrollBar().setValue(
            self.stateParamText.verticalScrollBar().maximum())
        self.FoMButton.setEnabled(True)
        self.toOpticsButton.setEnabled(True)

    @pyqtSlot()
    def updateFoM(self):
        """SLOT connected to FoMButton.clicked()
        Calculate Figure of merit.  """
        if len(self.stateHolder) < 2:
            print('Warning: FoM button triggered before state selection.')
            return
        self.FoMString = ''
        self._threadRun(self._calcFoM, self._updateFoM)

    def _fullPopulation(self):
        self.stateHolder = []
        self.qclayers.full_population()
        self.wavelength = self.inputWlBox.value()
        self.neff = self.qclayers.effective_ridx(self.wavelength)
        gain = self.qclayers.full_gain_spectrum(self.wavelength)
        self.gainCoeff = gain / self.qclayers.current
        self.specString = (
            "Carrier Leakage: <br>"
            "&nbsp;&nbsp; %.2f%%<br>"
            "Current: %.1f kA/cm<sup>2</sup><br>"
            "Gain coefficient:<br>"
            "&nbsp;&nbsp; %.2f cm/kA"
        ) % (self.qclayers.carrierLeak * 100, self.qclayers.current,
             self.gainCoeff)

    def _updatePopulation(self):
        self.stateParamText.setText(self.specString)
        self.FoMButton.setEnabled(False)
        self.gainSpecButton.setEnabled(True)
        self.toOpticsButton.setEnabled(True)
        self.update_quantumCanvas()

    @pyqtSlot()
    def fullPopulation(self):
        """SLOT connected to fullPopulationButton.clicked()"""
        if not self.qclayers.status.startswith('solved'):
            print('Full_population triggered before solving.')
        self.stateParamText.clear()
        self._threadRun(self._fullPopulation, self._updatePopulation)

    @pyqtSlot()
    def popGainSpecWindow(self):
        """SLOT connected to gainSpecButton.clicked()"""
        wlmin = 1.5
        wlmax = int(4*self.qclayers.wl + 0.5)/2
        wlDialog = WLDialog(self, wlmin, wlmax)
        wlmin, wlmax, buttonRes = wlDialog.exec()
        if buttonRes:
            wlFigure = figure()
            wlFigure.canvas.set_window_title('Gain Spectrum')
            axes = wlFigure.add_subplot(111)
            wls = np.linspace(wlmin, wlmax, 500)
            axes.plot(wls, self.qclayers.full_gain_spectrum(wls))
            axes.axhline(0, ls='--', lw=0.5)
            axes.set_xlabel('Wavelength (μm)')
            axes.set_ylabel('Gain (cm$^{-1}$)')
            if axes.get_ylim()[0] < -150:
                axes.set_ylim(bottom=-150)
            wlFigure.show()


class WLDialog(QDialog):
    def __init__(self, parent, wlMin, wlMax):
        super().__init__(parent)
        self.setModal(True)
        self.setWindowTitle('ErwinJr2: Gain Spectrum Range')
        self.wlMinBox = QDoubleSpinBox()
        self.wlMaxBox = QDoubleSpinBox()
        for box in (self.wlMinBox, self.wlMaxBox):
            box.setDecimals(1)
            box.setRange(0.5, 100)
            box.setSingleStep(1)
            box.setSuffix(' μm')
        self.wlMinBox.setValue(wlMin)
        self.wlMaxBox.setValue(wlMax)
        mainLayout = QGridLayout()
        mainLayout.addWidget(QLabel('min λ: '), 0, 0)
        mainLayout.addWidget(self.wlMinBox, 0, 1)
        mainLayout.addWidget(QLabel('max λ: '), 1, 0)
        mainLayout.addWidget(self.wlMaxBox, 1, 1)
        buttonBox = QDialogButtonBox(
            QDialogButtonBox.Ok | QDialogButtonBox.Cancel, Qt.Horizontal)
        buttonBox.accepted.connect(self.accept)
        buttonBox.rejected.connect(self.reject)
        mainLayout.addWidget(buttonBox, 2, 0, 1, 2)
        self.setLayout(mainLayout)

    def exec(self):
        res = super().exec()
        return self.wlMinBox.value(), self.wlMaxBox.value(), res
