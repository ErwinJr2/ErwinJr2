"""
This file defines the quantum tab of ErwinJr2, for simulating electron spectrum
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
from functools import partial, wraps

import darkdetect
import numpy as np
from matplotlib.pyplot import figure
from numpy import sqrt
from PyQt5.QtCore import QObject, Qt, QThread, pyqtSignal, pyqtSlot
from PyQt5.QtGui import QColor, QPalette
from PyQt5.QtWidgets import (
    QApplication,
    QCheckBox,
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QDoubleSpinBox,
    QGridLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QMessageBox,
    QPushButton,
    QSizePolicy,
    QSpinBox,
    QTableWidget,
    QTableWidgetItem,
    QTextEdit,
    QVBoxLayout,
    QWidget,
)

from ErwinJr2.gui.constants import EJ_ERROR, EJ_WARNING
from ErwinJr2.gui.custom_qt_class import MtrlComboBox
from ErwinJr2.gui.ej_canvas import EJCanvas, EJplotControl, config as plotconfig
from ErwinJr2.material import ALLOY_PARAM
from ErwinJr2.qc_plotter import plot_potential, plot_wf, scale_wf
from ErwinJr2.qclayers import (
    QC_MATERIAL,
    QCLayers,
    StateRecognizeError,
    c0,
    e0,
    eps0,
    h,
    hbar,
    optimize_layer,
)

# colors for different materials in tables
MTRL_COLORS_RGB = (
    (255, 255, 255),
    (230, 230, 240),
    (230, 240, 230),
    (240, 230, 230),
    (230, 240, 240),
    (240, 230, 240),
    (240, 240, 230),
    (230, 230, 230),
)


# TODO: this may not be necessary by better designer
def settingslot(fn):
    """A decorator to ask slots to skip reaction when doing massive
    updating (when self.updating is True)"""

    @wraps(fn)
    def wrapper(self, *args, **kwargs):
        if self.updating:
            return
        else:
            return fn(self, *args, **kwargs)

    return wrapper


class CalculateHolder(QObject):
    """Helper class to hold the calculation function and run it in a separate thread."""

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
            self.warning.emit("Warning: " + e.expression)
            self.succeed.emit()
        except (IndexError, TypeError):
            self.failed.emit(traceback.format_exc())
        finally:
            self.finished.emit()


class QuantumTab(QWidget):
    """The Quantum Tab of ErwinJr2. This is designed to be a GUI wrapper of
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
    to_optics = pyqtSignal(QCLayers)

    def __init__(self, qclayers=None, parent=None):
        super().__init__(parent)
        self.qclayers = qclayers if qclayers else QCLayers()
        self._calc_thread = QThread()
        self._worker = None

        self.mtrlcolors = tuple([QColor(*rgb) for rgb in MTRL_COLORS_RGB])
        if darkdetect.isDark():
            self.mtrlcolors = tuple(
                [QColor(*(2 * (255 - c) for c in rgb)) for rgb in MTRL_COLORS_RGB]
            )

        self.updating = False
        self.plot_vx = False
        self.plot_vl = False
        self.plot_lh = False
        self.plot_so = False
        self.layer_selected = None

        self.state_holder = []
        self.pair_selected = False
        # plotType can be mode, wf or DoS (TODO)
        #  self.plotType = "wf"
        self.plot_type = "mode"
        #  self.fillPlot = 0.3  # alpha of fill; False for not fill
        self.fill_plot = False

        # Platform dependent settings, eg. layout size settings
        platform = sys.platform
        if platform.startswith("win"):
            setting_box_width = 150
            layer_box_width = 230
            solve_box_width = 170
        elif platform.startswith("darwin"):
            setting_box_width = 160
            layer_box_width = 250
            solve_box_width = 180
        elif platform.startswith("linux"):
            setting_box_width = 150
            layer_box_width = 250
            solve_box_width = 185
        else:
            QMessageBox.warning(self, EJ_WARNING, f"Platform {platform} not tested.")
            setting_box_width = 150
            layer_box_width = 400
            solve_box_width = 190

        quantum_layout = QHBoxLayout()
        quantum_layout.setSpacing(1)
        setting_box = self._gen_setting_box(setting_box_width)
        layer_box = self._gen_layer_box(layer_box_width)
        figure_box, plot_control_grid = self._gen_figure_box()
        solve_box = self._gen_solve_box(plot_control_grid, solve_box_width)
        quantum_layout.addLayout(setting_box)
        quantum_layout.addLayout(layer_box)
        quantum_layout.addLayout(solve_box)
        quantum_layout.addLayout(figure_box)
        self.setLayout(quantum_layout)
        self.setAutoFillBackground(True)
        self.setBackgroundRole(QPalette.Window)

        self.reload()

    # __init__ end

    def _gen_setting_box(self, width):
        """Return a Qt Layout object containing all setting parameters"""
        setting_box = QVBoxLayout()

        # set up description box
        setting_box.addWidget(QLabel("<center><b>Description</b></center>"))
        self.desc_box = QTextEdit("")
        self.desc_box.setReadOnly(False)
        self.desc_box.setMinimumHeight(50)
        self.desc_box.setMaximumHeight(80)
        self.desc_box.setFixedWidth(width)
        self.desc_box.setSizePolicy(
            QSizePolicy(QSizePolicy.Fixed, QSizePolicy.Expanding)
        )
        self.desc_box.textChanged.connect(self.input_description)
        self.desc_box.setToolTip("Description text is not used in modeling.")
        setting_box.addWidget(self.desc_box)

        setting_box.addWidget(QLabel("<center><b>Substrate</b></center>"))
        self.input_substrate_box = QComboBox()
        self.input_substrate_box.setMaximumWidth(width)
        self.input_substrate_box.addItems(QC_MATERIAL.keys())
        self.input_substrate_box.currentIndexChanged[str].connect(self.input_substrate)
        self.input_substrate_box.setToolTip(
            "The substrate decides possible choices of materials."
        )
        setting_box.addWidget(self.input_substrate_box)

        setting_box.addWidget(
            QLabel("<center><b><i>E<sub>field</sub></i></b></center>")
        )
        self.input_e_field_box = QDoubleSpinBox()
        self.input_e_field_box.setMaximumWidth(width)
        self.input_e_field_box.setDecimals(1)
        self.input_e_field_box.setSuffix(" kV/cm")
        self.input_e_field_box.setRange(-250.0, 250.0)
        self.input_e_field_box.valueChanged[float].connect(self.input_e_field)
        self.input_e_field_box.setToolTip(
            "The Bias setting for the structure. "
            "Negative values may results in artifact results."
        )
        setting_box.addWidget(self.input_e_field_box)

        setting_box.addWidget(QLabel("<center><b>Position<br>Resolution</b></center>"))
        self.input_x_res_box = QDoubleSpinBox()
        self.input_x_res_box.setMaximumWidth(width)
        self.input_x_res_box.setDecimals(2)
        self.input_x_res_box.setRange(0.01, 1.0)
        self.input_x_res_box.setSingleStep(0.1)
        self.input_x_res_box.setSuffix(" \u212b")
        self.input_x_res_box.valueChanged[float].connect(self.input_xres)
        setting_box.addWidget(self.input_x_res_box)

        self.e_res_label = QLabel("<center><b>Energy<br>Resolution</b></center>")
        setting_box.addWidget(self.e_res_label)
        self.input_e_res_box = QDoubleSpinBox()
        self.input_e_res_box.setMaximumWidth(width)
        self.input_e_res_box.setDecimals(2)
        self.input_e_res_box.setRange(0.0, 10.0)
        self.input_e_res_box.setSingleStep(0.1)
        self.input_e_res_box.setSuffix(" meV")
        self.input_e_res_box.valueChanged[float].connect(self.input_e_res)
        self.input_e_res_box.setToolTip(
            "The energy resolution used for ODE solver root finding. "
            "This number being too large may results in miss of states. \n"
            "The value is not used for matrix solver."
        )
        setting_box.addWidget(self.input_e_res_box)
        self.e_count_label = QLabel(
            "<center><b>No of States<br>Per Period</b></center>"
        )
        setting_box.addWidget(self.e_count_label)
        self.input_e_count_box = QSpinBox()
        self.input_e_count_box.setMaximumWidth(width)
        self.input_e_count_box.setRange(1, 100)
        self.input_e_count_box.valueChanged.connect(self.input_e_count)
        self.input_e_count_box.setToolTip(
            "Increase this number if you believe the number of states solved"
            "is not enough. This is only used for matrix solver."
        )
        setting_box.addWidget(self.input_e_count_box)
        self.algo_param_update()

        setting_box.addWidget(QLabel("<center><b>Repeats</b></center>"))
        self.input_repeats_box = QSpinBox()
        self.input_repeats_box.setMaximumWidth(width)
        self.input_repeats_box.setRange(1, 5)
        self.input_repeats_box.valueChanged[int].connect(self.input_repeats)
        setting_box.addWidget(self.input_repeats_box)

        setting_box.addWidget(QLabel("<center><b>Wavelength</b></center>"))
        self.input_wl_box = QDoubleSpinBox()
        self.input_wl_box.setMaximumWidth(width)
        self.input_wl_box.setDecimals(1)
        self.input_wl_box.setRange(1.5, 40)
        self.input_wl_box.setSingleStep(1)
        self.input_wl_box.setSuffix(" μm")
        self.input_wl_box.valueChanged[float].connect(self.input_wl)
        setting_box.addWidget(self.input_wl_box)

        # Basis solver divider setting
        basis_group_box = QGroupBox("Basis Divisions")
        self.input_ar_injector_check = QCheckBox("AR->Injector")
        self.input_injector_ar_check = QCheckBox("Injector->AR")
        self.input_ar_injector_check.setChecked(True)
        self.input_injector_ar_check.setChecked(True)
        basis_layout = QVBoxLayout()
        basis_layout.addWidget(self.input_ar_injector_check)
        basis_layout.addWidget(self.input_injector_ar_check)
        basis_group_box.setLayout(basis_layout)
        self.input_ar_injector_check.stateChanged.connect(self.input_basis)
        self.input_injector_ar_check.stateChanged.connect(self.input_basis)
        setting_box.addWidget(basis_group_box)

        # Period information groupbox
        lp_layout_group_box = QGroupBox("Period Info")
        lp_layout_group_box.setFixedWidth(width)
        self.lp_first_spin_box = QSpinBox()
        self.lp_first_spin_box.setValue(1)
        self.lp_first_spin_box.setRange(1, 1)
        self.lp_first_spin_box.valueChanged.connect(self.update_lp_box)
        self.lp_last_spin_box = QSpinBox()
        self.lp_last_spin_box.setRange(1, 1)
        self.lp_last_spin_box.valueChanged.connect(self.update_lp_box)
        self.lp_string_box = QTextEdit("")
        self.lp_string_box.setReadOnly(True)
        self.lp_string_box.setSizePolicy(
            QSizePolicy(QSizePolicy.Fixed, QSizePolicy.Preferred)
        )
        self.lp_string_box.setMinimumHeight(30)
        # self.LpStringBox.setMaximumHeight(120)
        lp_layout = QGridLayout()
        lp_layout.addWidget(QLabel("<b>first</b>"), 0, 0)
        lp_layout.addWidget(QLabel("<b>last</b>"), 0, 1)
        lp_layout.addWidget(self.lp_first_spin_box, 1, 0)
        lp_layout.addWidget(self.lp_last_spin_box, 1, 1)
        lp_layout.addWidget(self.lp_string_box, 2, 0, 1, 2)
        lp_layout.setSpacing(1)
        lp_layout_group_box.setLayout(lp_layout)
        setting_box.addWidget(lp_layout_group_box)

        setting_box.addStretch()
        return setting_box
        # _generateSettingBox end

    def _gen_layer_box(self, width):
        """Return a Qt Layout object containing all layer parameters"""
        layer_box = QGridLayout()
        layer_box.setSpacing(5)
        self.insert_layer_above_button = QPushButton("Insert Layer")
        self.insert_layer_above_button.clicked.connect(self.insert_layer_above)
        layer_box.addWidget(self.insert_layer_above_button, 0, 0)
        self.delete_layer_button = QPushButton("Delete Layer")
        self.delete_layer_button.clicked.connect(self.delete_layer)
        layer_box.addWidget(self.delete_layer_button, 0, 1)
        self.optimize_layer_button = QPushButton("Optimize Layer")
        self.optimize_layer_button.clicked.connect(self.optimize_layer)
        self.optimize_layer_button.setEnabled(False)
        # TODO: global optimize is not implemented
        self.global_optimize_button = QPushButton("Global Optimize")
        self.global_optimize_button.setEnabled(False)
        layer_box.addWidget(self.optimize_layer_button, 1, 0)
        layer_box.addWidget(self.global_optimize_button, 1, 1)

        # set up layerTable
        self.layer_table = QTableWidget()
        self.layer_table.setSelectionBehavior(QTableWidget.SelectRows)
        self.layer_table.setSelectionMode(QTableWidget.SingleSelection)
        self.layer_table.setFixedWidth(width)
        self.layer_table.setSizePolicy(
            QSizePolicy(QSizePolicy.Fixed, QSizePolicy.Expanding)
        )
        self.layer_table.itemChanged.connect(self.layer_table_item_changed)
        self.layer_table.itemSelectionChanged.connect(
            self.layer_table_item_selection_changed
        )
        layer_box.addWidget(self.layer_table, 2, 0, 1, 2)

        return layer_box
        # _generateLayerBox end

    def _gen_figure_box(self):
        """Return:
        A Qt Layout containing a canvas for plotting;
        A QGridLayout containing plot control;"""
        self.quantum_canvas = EJCanvas(
            xlabel="Position (Å)", ylabel="Energy (eV)", parent=self
        )
        self.plot_control = EJplotControl(self.quantum_canvas, self)
        self.quantum_canvas.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        figure_box = QVBoxLayout()
        figure_box.addWidget(self.quantum_canvas)

        self.zoom_in_button = QPushButton("Zoom")
        self.plot_control.set_action("zoom", self.zoom_in_button)
        self.zoom_out_button = QPushButton("Reset")
        self.plot_control.set_action("home", self.zoom_out_button)
        self.pan_button = QPushButton("Pan")  # to move
        self.plot_control.set_action("pan", self.pan_button)
        self.layer_select_button = QPushButton("Layer Select")
        self.plot_control.set_custom(
            "layerselect", self.layer_select_button, self.layer_select
        )
        self.layer_select_button.clicked[bool].connect(self.layer_select_mode)
        self.clear_wfs_button = QPushButton("Clear")
        self.clear_wfs_button.clicked.connect(self.clear_wfs)
        plot_control_grid = QGridLayout()
        plot_control_grid.addWidget(self.layer_select_button, 0, 0, 1, 2)
        plot_control_grid.addWidget(self.zoom_in_button, 1, 0, 1, 1)
        plot_control_grid.addWidget(self.zoom_out_button, 1, 1, 1, 1)
        plot_control_grid.addWidget(self.pan_button, 2, 0, 1, 1)
        plot_control_grid.addWidget(self.clear_wfs_button, 2, 1, 1, 1)
        plot_control_grid.setSpacing(5)

        return figure_box, plot_control_grid
        # _generateFigureBox end

    def _gen_solve_box(self, plot_control_grid, width):
        """Return a Qt Layout containing material information,
        eigensolve control, states properties calculation and plot control"""
        self.solve_box = QVBoxLayout()
        self.solve_basis_button = QPushButton("Solve Basis")
        self.solve_basis_button.clicked.connect(self.solve_basis)
        self.solve_box.addWidget(self.solve_basis_button)
        self.solve_whole_button = QPushButton("Solve Whole")
        self.solve_whole_button.clicked.connect(self.solve_whole)
        self.solve_box.addWidget(self.solve_whole_button)

        # set up material composition inputs
        self.mtrl_table = QTableWidget()
        #  self.mtrlTable.setSelectionBehavior(QTableWidget.SelectItems)
        #  self.mtrlTable.setSelectionMode(QTableWidget.NoSelection)
        self.mtrl_table.setSelectionBehavior(QTableWidget.SelectRows)
        self.mtrl_table.setSelectionMode(QTableWidget.SingleSelection)
        self.mtrl_table.setMaximumWidth(width)
        self.mtrl_table.setMinimumWidth(width)
        self.mtrl_table.setMinimumHeight(100)
        self.mtrl_table.setMaximumHeight(120)
        self.mtrl_table.setSizePolicy(
            QSizePolicy(QSizePolicy.Fixed, QSizePolicy.Minimum)
        )
        self.mtrl_table.itemChanged.connect(self.mtrl_table_item_changed)
        self.mtrl_table.itemSelectionChanged.connect(
            self.mtrl_table_item_selection_changed
        )
        self.offset_label = QLabel("")
        self.net_strian_label = QLabel("")
        self.lo_phonon_label = QLabel("")
        self.add_mtrl_button = QPushButton("Add")
        self.add_mtrl_button.clicked.connect(self.add_mtrl)
        self.del_mtrl_button = QPushButton("Del")
        self.del_mtrl_button.clicked.connect(self.del_mtrl)
        mtrl_grid = QGridLayout()
        mtrl_grid.addWidget(self.add_mtrl_button, 0, 0)
        mtrl_grid.addWidget(self.del_mtrl_button, 0, 1)
        mtrl_grid.addWidget(self.mtrl_table, 1, 0, 1, 2)
        mtrl_grid.addWidget(self.offset_label, 2, 0, 1, 2)
        mtrl_grid.addWidget(self.net_strian_label, 3, 0, 1, 2)
        mtrl_grid.addWidget(self.lo_phonon_label, 4, 0, 1, 2)
        mtrl_grid.setSpacing(5)
        mtrl_group_box = QGroupBox("Materials")
        mtrl_group_box.setLayout(mtrl_grid)
        self.solve_box.addWidget(mtrl_group_box)

        # IFR setting
        self.ifr_def_box = QCheckBox("Constant IFR")
        self.ifr_def_box.stateChanged.connect(self.ifr_def_constant)
        self.ifr_delta_box = QDoubleSpinBox()
        self.ifr_delta_box.setMaximumWidth(width)
        self.ifr_delta_box.setDecimals(1)
        self.ifr_delta_box.setSuffix(" Å")
        self.ifr_delta_box.setRange(0.0, 9999.0)
        self.ifr_delta_box.valueChanged[float].connect(self.input_ifr_delta)
        self.ifr_lambda_box = QDoubleSpinBox()
        self.ifr_lambda_box.setMaximumWidth(width)
        self.ifr_lambda_box.setDecimals(1)
        self.ifr_lambda_box.setSuffix(" Å")
        self.ifr_lambda_box.setRange(0.0, 9999.0)
        ifr_grid = QGridLayout()
        ifr_grid.addWidget(self.ifr_def_box, 0, 0, 1, 2)
        ifr_grid.addWidget(QLabel("<center>IRF Δ</center>"), 1, 0)
        ifr_grid.addWidget(self.ifr_delta_box, 2, 0)
        ifr_grid.addWidget(QLabel("<center>IRF Λ</center>"), 1, 1)
        ifr_grid.addWidget(self.ifr_lambda_box, 2, 1)
        self.ifr_lambda_box.valueChanged[float].connect(self.input_ifr_lambda)
        # This is exposed s.t. it can be hidden when IFR is off.
        self._ifr_group_box = QGroupBox("Interface Roughness")
        self._ifr_group_box.setLayout(ifr_grid)
        self.solve_box.addWidget(self._ifr_group_box)

        # set up plot control inputs
        plot_control_group_box = QGroupBox("Plot Controls")
        plot_control_group_box.setLayout(plot_control_grid)
        self.solve_box.addWidget(plot_control_group_box)

        # set up Calculate controls
        self.pair_select_button = QPushButton("Pair Select")
        self.pair_select_button.setEnabled(False)
        self.plot_control.set_custom(
            "pairselect", self.pair_select_button, self.state_pick
        )
        self.pair_select_button.clicked[bool].connect(self.pair_select_mode)
        self.fom_button = QPushButton("FoM")
        self.fom_button.setEnabled(False)
        self.fom_button.clicked.connect(self.update_fom)
        self.full_population_button = QPushButton("Population")
        self.full_population_button.setEnabled(False)
        self.full_population_button.clicked.connect(self.full_population)
        self.to_optics_button = QPushButton("->Optics")
        self.to_optics_button.setEnabled(False)
        self.gain_spec_button = QPushButton("Gain Spec")
        self.gain_spec_button.setEnabled(False)
        self.gain_spec_button.clicked.connect(self.pop_gain_spec_window)
        # signal is processed in ErwinJr main window
        self.state_param_text = QTextEdit("")
        self.state_param_text.setReadOnly(True)
        self.state_param_text.setMaximumWidth(width)
        self.state_param_text.setMinimumWidth(width)
        self.state_param_text.setMinimumHeight(150)
        self.state_param_text.setSizePolicy(
            QSizePolicy(QSizePolicy.Fixed, QSizePolicy.Expanding)
        )
        self.state_param_text.textChanged.connect(
            lambda: self.to_optics_button.setEnabled(False)
        )
        calculate_control_grid = QGridLayout()
        calculate_control_grid.addWidget(self.pair_select_button, 0, 0, 1, 2)
        calculate_control_grid.addWidget(self.fom_button, 1, 0, 1, 1)
        calculate_control_grid.addWidget(self.full_population_button, 1, 1, 1, 1)
        calculate_control_grid.addWidget(self.to_optics_button, 2, 0, 1, 1)
        calculate_control_grid.addWidget(self.gain_spec_button, 2, 1, 1, 1)
        calculate_control_grid.addWidget(self.state_param_text, 3, 0, 1, 2)
        calculate_control_grid.setSpacing(5)
        # TODO: voltage efficiency, hbar omega / field * Lp
        calculate_control_group_box = QGroupBox("Calculate")
        calculate_control_group_box.setLayout(calculate_control_grid)
        self.solve_box.addWidget(calculate_control_group_box)

        self.solve_box.addStretch()
        return self.solve_box
        # _generateSolveBox end

    def reload(self):
        """Reload everything from qclayers, typically when it's changed
        externally or loaded"""
        # TODO: reconsider when reload is needed.
        self.qclayers.populate_x()
        self._update_setting_box()
        self._update_mtrl_list()
        self._update_mtrl_table()
        self._update_layer_table()
        self._update_lp_limits()
        self.update_lp_box()
        self.layer_table.clearSelection()
        self.update_ifr_settings()
        self.ifr_settings_check_const()
        self.update_quantum_canvas()

    def _update_mtrl_list(self):
        """Update self.mtrlList according to self.qclayers material
        information. mtrlList is used for mtrl column in layerTable"""
        self.mtrm_list = []
        for n, mtrl in enumerate(self.qclayers.materials):
            name = ALLOY_PARAM[mtrl]["name"]
            name = name.replace("1-x", str(1 - self.qclayers.mole_fracs[n]))
            name = name.replace("x", str(self.qclayers.mole_fracs[n]))
            name = f"#{n + 1} "  # + name
            self.mtrm_list.append(name)

    # ===========================================================================
    # SettingBox Controls
    # ===========================================================================
    def _update_setting_box(self):
        self.updating = True
        self.desc_box.setText(self.qclayers.description)
        self.input_substrate_box.setCurrentText(self.qclayers.substrate)
        self.input_e_field_box.setValue(self.qclayers.e_field)
        self.input_x_res_box.setValue(self.qclayers.x_step)
        self.input_e_res_box.setValue(self.qclayers.e_step)
        self.input_e_count_box.setValue(self.qclayers.state_per_repeat)
        self.input_wl_box.setValue(self.qclayers.wl)
        self.input_repeats_box.setValue(self.qclayers.repeats)
        self.updating = False

    @pyqtSlot("QString")
    @settingslot
    def input_substrate(self, substrate_type):
        """SLOT connected to inputSubstrateBox.currentIndexChanged(QString)
        update substrate chosen"""
        if substrate_type in ("InP", "GaAs", "GaSb"):
            self.qclayers.set_substrate(substrate_type)
        else:
            QMessageBox.information(
                self,
                EJ_ERROR,
                f"{substrate_type} substrates have not yet been implemented.",
            )
            self.input_substrate_box.setCurrentIndex(
                self.input_substrate_box.findText(self.qclayers.substrate)
            )
            return
        self.qclayers.populate_x()
        self.update_mtrl_info()
        self._update_layer_table()
        self._update_mtrl_table()
        self.update_quantum_canvas()
        self.dirty.emit()

    @pyqtSlot(float)
    @settingslot
    def input_e_field(self, e_field):
        """SLOT connected to inputEFieldBox.valueChanged(double)
        update external E field in unit kV/cm"""
        self.clear_wfs()
        self.qclayers.e_field = e_field
        self.qclayers.populate_x()
        self.dirty.emit()
        self.update_quantum_canvas()

    @pyqtSlot(float)
    @settingslot
    def input_xres(self, xres):
        """SLOT connected to inputxResBox.valueChanged
        update position resolution (xres), in angstrom"""
        self.clear_wfs()
        self.qclayers.x_step = xres
        self.qclayers.populate_x()
        self.dirty.emit()
        self.update_quantum_canvas()

    @pyqtSlot(float)
    @settingslot
    def input_e_res(self, e_res):
        """SLOT connected to inputEresBox.valueChanged
        Update initial energy resolution for eigensolver. Set this too small
        may result in loss of some eigenvalue."""
        self.qclayers.e_step = e_res
        self.dirty.emit()

    @pyqtSlot(int)
    @settingslot
    def input_e_count(self, count: int):
        self.qclayers.state_per_repeat = count
        self.qclayers.matrix_eigen_count = count * self.qclayers.repeats
        self.dirty.emit()

    @pyqtSlot(int)
    @settingslot
    def input_repeats(self, repeat):
        """SLOT connected to SINGAL self.inputRepeatsBox.valueChanged(int)
        update number of repeats for the whole structure."""
        self.clear_wfs()
        self.qclayers.repeats = repeat
        self.qclayers.populate_x()
        self.dirty.emit()
        self.update_quantum_canvas()

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
        self.qclayers.basis_ar_injector = self.input_ar_injector_check.isChecked()
        self.qclayers.basis_injector_ar = self.input_injector_ar_check.isChecked()

    def _update_lp_limits(self):
        """Update Lp select range in the Period Info box (GUI)
        This should be called whenever number of layers is changed
        Will trigger LpFirstSpinbox/LpLastSpinBox.valueChanged(int) SIGNAL
        and thus update_Lp_box
        """
        self.lp_first_spin_box.setRange(1, len(self.qclayers.layer_widths) - 1)
        self.lp_first_spin_box.setValue(1)
        self.lp_last_spin_box.setRange(1, len(self.qclayers.layer_widths))
        self.lp_last_spin_box.setValue(len(self.qclayers.layer_widths))

    @pyqtSlot()
    @settingslot
    def update_lp_box(self):
        """Update Lp box in the Period Info box (GUI):
            Lp:total length
            nD: average doping (cm-3)
            ns: 2D carrier density in 1E11 cm-2
        This needs layerWidths and layerDopings information of self.qclayers
        SLOT connected to LpFirstSpinbox/LpLastSpinBox.valueChanged(int)
        """
        lp_first = self.lp_first_spin_box.value() - 1
        lp_last = self.lp_last_spin_box.value()
        # +1 because range is not inclusive of last value
        # total length of the layers (1 period)
        lp = sum(self.qclayers.layer_widths[lp_first:lp_last])
        lp_string = f"Lp: {lp:.1f} \u212b<br>"
        # average doping of the layers
        ns = sum(
            self.qclayers.layer_dopings[n] * self.qclayers.layer_widths[n]
            for n in range(lp_first, lp_last)
        )
        if lp == 0:
            lp_string += "n<sub>D</sub>: NA\u00d710<sup>17</sup>" "cm<sup>-3</sup><br>"
        else:
            nd = ns / lp
            lp_string += (
                f"N<sub>D</sub>: {nd:6.3f}\u00d710<sup>17</sup>" "cm<sup>-3</sup><br>"
            )
        # 2D carrier density in 1E11cm-2
        ns = ns * 1e-2
        lp_string += (
            f"N<sub>s</sub>: {ns:6.3f}\u00d710<sup>11</sup>" "cm<sup>-2</sup><br>"
        )
        neff = self.qclayers.effective_ridx(self.input_wl_box.value())
        lp_string += f"n<sub>eff</sub>: {neff:.2f}"
        self.lp_string_box.setText(lp_string)

    @pyqtSlot()
    @settingslot
    def input_description(self):
        """SLOT connected to self.descBox.textChanged()
        Change description string for the design"""
        self.qclayers.description = self.desc_box.toPlainText()
        self.dirty.emit()

    def set_temperature(self, temp):
        """A public method to wrap temperature settings"""
        self.qclayers.set_temperature(temp)
        self.qclayers.populate_x()
        self.dirty.emit()
        self.update_mtrl_info()
        self.update_quantum_canvas()

    # ===========================================================================
    # Layer Table Control
    # ===========================================================================
    def _update_layer_table(self):
        """Refresh layer table, called every time after data update"""
        # Block itemChanged SIGNAL while refreshing
        self.layer_table.blockSignals(True)
        self.layer_table.clear()
        self.layer_table.setColumnCount(5)
        # An extra blank line for adding new layers
        self.layer_table.setRowCount(len(self.qclayers.layer_widths) + 1)
        vert_labels = [str(n + 1) for n in range(len(self.qclayers.layer_widths))]
        self.layer_table.setHorizontalHeaderLabels(
            ("Width", "ML", "Mtrl", "AR", "Doping")
        )
        self.layer_table.setVerticalHeaderLabels(vert_labels)

        # gray2 = QColor(230, 230, 230)  # for unchangeable background

        for q, layer_widths in enumerate(self.qclayers.layer_widths):
            color = self.mtrlcolors[self.qclayers.layer_mtrls[q] % len(self.mtrlcolors)]
            # Width Setup
            width = QTableWidgetItem(f"{layer_widths:5.1f}")
            width.setTextAlignment(Qt.AlignCenter)
            width.setBackground(color)
            self.layer_table.setItem(q, 0, width)

            # "ML" number of monolayer Setup
            # TODO: monolayer <-> thickness should not include
            # temperature/strain correction
            ml_thickness = (
                self.qclayers.mtrl_alloys[self.qclayers.layer_mtrls[q]].a_perp / 2
            )
            # a_perp has two layer of atoms (one III and one V)
            num_ml = QTableWidgetItem(f"{layer_widths / ml_thickness:5.1f}")
            num_ml.setTextAlignment(Qt.AlignCenter)
            num_ml.setBackground(color)
            self.layer_table.setItem(q, 1, num_ml)

            # Material Setup
            #  mtrlWidget = QComboBox()
            mtrl_widget = MtrlComboBox()
            mtrl_widget.addItems(self.mtrm_list)
            mtrl_widget.setCurrentIndex(self.qclayers.layer_mtrls[q])
            mtrl_widget.currentIndexChanged.connect(
                partial(self.layer_table_mtrl_changed, q)
            )
            self.layer_table.setCellWidget(q, 2, mtrl_widget)

            # Active Region Layer Setup
            ar_item = QTableWidgetItem()
            ar_item.setCheckState(
                Qt.Checked if self.qclayers.layer_ar[q] == 1 else Qt.Unchecked
            )
            ar_item.setTextAlignment(Qt.AlignCenter)
            ar_item.setBackground(color)
            self.layer_table.setItem(q, 3, ar_item)

            # Layer Doping Setup
            doping = QTableWidgetItem(str(self.qclayers.layer_dopings[q]))
            doping.setTextAlignment(Qt.AlignCenter)
            doping.setBackground(color)
            self.layer_table.setItem(q, 4, doping)

        self.layer_table.resizeColumnsToContents()
        self.layer_table.blockSignals(False)

    @pyqtSlot()
    def insert_layer_above(self):
        """SLOT connected to self.insertLayerAboveButton.clicked()"""
        row = self.layer_table.currentRow()
        mtrl_count = len(self.qclayers.materials)
        if row >= len(self.qclayers.layer_widths) or row < 0:
            # Add new lines in the last layer
            row = len(self.qclayers.layer_widths)
            ar = self.qclayers.layer_ar[row - 1] and self.qclayers.layer_ar[0]
            doping = self.qclayers.layer_dopings[row - 1]
            mtrl_idx = (self.qclayers.layer_mtrls[row - 1] + 1) % mtrl_count
        else:
            ar = self.qclayers.layer_ar[row] and self.qclayers.layer_ar[row - 1]
            doping = self.qclayers.layer_dopings[row]
            mtrl_idx = (self.qclayers.layer_mtrls[row - 1] - 1) % mtrl_count

        self.clear_wfs()
        self.qclayers.add_layer(row, 0.0, mtrl_idx, ar, doping)
        self.qclayers.populate_x()
        self._update_lp_limits()
        self.update_lp_box()
        self._update_layer_table()
        self.update_mtrl_info()
        self.layer_table.selectRow(row)
        # Trigger itemSelectionChanged SIGNAL and thus update_quantumCanvas
        self.dirty.emit()

    @pyqtSlot()
    def delete_layer(self):
        """SLOT connected to self.deleteLayerButton.clicked()"""
        row = self.layer_table.currentRow()
        if row == -1 or row >= len(self.qclayers.layer_widths):
            return
        # don't delete last layer
        if len(self.qclayers.layer_widths) == 1:
            self.qclayers.layer_widths[0] = 0.0
            return

        self.clear_wfs()
        self.qclayers.del_layer(row)
        self.qclayers.populate_x()
        self._update_lp_limits()
        self.update_lp_box()
        self._update_layer_table()
        self.update_mtrl_info()
        self.layer_table.selectRow(row)
        # Trigger itemSelectionChanged SIGNAL and thus update_quantumCanvas
        self.dirty.emit()

    def _optimize_layer(self, n, upper, lower):
        optimize_layer(self.qclayers, n, upper, lower)
        self.qclayers.solve_whole()
        self.state_holder = [upper, lower]
        self.pair_selected = True
        self._calc()

    def show_optimize(self, n=None):
        self._update_layer_table()
        if n is None:
            self._update_population()
        else:
            self.layer_table.setCurrentCell(n, 0)
            # self.updateSelected()
            self._update_fom()
            self.dirty.emit()

    @pyqtSlot()
    def optimize_layer(self):
        """SLOT connected to self.optimizeLayerButton.clicked()"""
        n = self.layer_table.currentRow()
        if self.qclayers.status == "unsolved":
            QMessageBox.warning(self, EJ_WARNING, "Solve the model first.")
            return
        if n < 0 or n > len(self.qclayers.layer_widths):
            QMessageBox.warning(self, EJ_ERROR, "Select the layer to optimize.")
            return
        if len(self.state_holder) != 2:
            QMessageBox.warning(self, EJ_ERROR, "Select state pair to optimize.")
            return
        try:
            upper = self.state_holder[0]
            lower = self.state_holder[1]
            if self.qclayers.eigen_es[upper] < self.qclayers.eigen_es[lower]:
                upper, lower = lower, upper
        except ValueError:
            QMessageBox.warning(self, EJ_ERROR, "Select state pair to optimize.")
            return
        self.state_param_text.clear()
        self.clear_wfs()
        self._thread_run(
            lambda: self._optimize_layer(n, upper, lower), lambda: self.show_optimize(n)
        )

    @pyqtSlot(QTableWidgetItem)
    def layer_table_item_changed(self, item):
        """SLOT connected to layerTable.itemChanged(QTableWidgetItem*)
        Update layer profile after user input"""
        row = item.row()
        mtrl_count = len(self.qclayers.materials)
        if row > len(self.qclayers.layer_widths):
            raise ValueError("Bad layer width input row number")
        column = item.column()
        if column in (0, 1, 4):
            try:
                value = float(item.text())
            except ValueError:
                # invalid input
                QMessageBox.warning(self, EJ_ERROR, "This value should be a number")
                self._update_layer_table()
                return
            if column in (0, 1):
                if column == 0:  # column == 0 for Widths column
                    new_width = value
                else:  # column == 1 for "ML" number of monolayers
                    new_width = (
                        value
                        * self.qclayers.mtrl_alloys[
                            self.qclayers.layer_mtrls[row]
                        ].a_perp
                        / 2
                    )
                    new_width = round(new_width, 1)

                if row == len(self.qclayers.layer_widths):
                    # add row at end of list
                    ar = self.qclayers.layer_ar[row - 1]
                    doping = self.qclayers.layer_dopings[row - 1]
                    mtrl_idx = (self.qclayers.layer_mtrls[row - 1] + 1) % mtrl_count
                    self.qclayers.add_layer(row, new_width, mtrl_idx, ar, doping)
                    row += 1  # used so that last (blank) row is again selected
                    self._update_lp_limits()
                    self.update_lp_box()
                else:  # change Width of selected row in-place
                    self.qclayers.layer_widths[row] = new_width
                    self.update_lp_box()
            else:
                # column == 4 for item change in Doping column
                doping = value
                if row == len(self.qclayers.layer_widths):
                    ar = self.qclayers.layer_ar[row - 1]
                    mtrl_idx = (self.qclayers.layer_mtrls[row - 1] + 1) % mtrl_count
                    self.qclayers.add_layer(row, 0.0, mtrl_idx, ar, doping)
                    row += 1  # used so that last (blank) row is again selected
                    self._update_lp_limits()
                    self.update_lp_box()
                else:
                    self.qclayers.layer_dopings[row] = doping

        elif column == 2:
            # column == 2 for item change in mtrl column, should be
            # controlled by layerTable_materialChanged
            raise RuntimeError("Should not be here")

        elif column == 3:
            # column == 3 for item change in AR column
            if row == len(self.qclayers.layer_widths):
                # don't do anything if row is last row
                return
            self.qclayers.layer_ar[row] = item.checkState() == Qt.Checked

        else:
            raise RuntimeError("Should not be here")

        self._update_layer_table()
        self.clear_wfs()
        self.layer_table.setCurrentCell(row, column)
        self.qclayers.populate_x()
        self.update_mtrl_info()
        self.update_quantum_canvas()
        self.dirty.emit()

    @pyqtSlot()
    def layer_table_item_selection_changed(self):
        """SLOT connected to layerTable.itemSelectionChanged()"""
        row = self.layer_table.currentRow()
        if row < len(self.qclayers.layer_widths):
            self.layer_selected = row
        else:
            self.layer_selected = None
            self.layer_table.clearSelection()
        self.update_quantum_canvas()

    def layer_table_mtrl_changed(self, row, selection):
        """SLOT as partial(self.layerTable_materialChanged, q)) connected to
        materialWidget.currentIndexChanged(int)"""
        self.qclayers.layer_mtrls[row] = selection
        self.clear_wfs()
        self.qclayers.populate_x()
        self.layer_table.selectRow(row)
        self.update_mtrl_info()
        self._update_layer_table()
        # Trigger itemSelectionChanged SIGNAL and thus update_quantumCanvas
        self.dirty.emit()

    @pyqtSlot()
    def rotate_layer(self):
        """Move last layer to first layer, SLOT open to public"""
        self.clear_wfs()
        self.qclayers.rotate_layer()
        self.qclayers.populate_x()
        self._update_layer_table()
        self.layer_table.setCurrentCell(1, 0)
        self.update_quantum_canvas()
        self.dirty.emit()

    @pyqtSlot()
    def invert_layer(self):
        """Invert the order of the layers"""
        self.clear_wfs()
        self.qclayers.invert_layer()
        self.qclayers.populate_x()
        self._update_layer_table()
        row = self.layer_table.currentRow()
        if row >= 0 and row < len(self.layerWidths):
            self.layer_table.setCurrentCell(
                len(self.qclayers.layer_widths) - 1 - row, 0
            )
        self.update_quantum_canvas()
        self.dirty.emit()

    @pyqtSlot()
    def ARonly(self):
        self.qclayers.basis_ar_only = not self.qclayers.basis_ar_only
        self.clear_wfs()

    @pyqtSlot()
    def copy_structure(self):
        clipboard = QApplication.clipboard()
        string = ""
        for width in self.qclayers.layer_widths:
            string += f"{width:.1f}\n"
        clipboard.setText(string)

    # =========================================================================
    # mtrlTable Control
    # =========================================================================
    def _update_mtrl_table(self):
        """Set up material table, both format and consistency with qclayers"""
        self.mtrl_table.blockSignals(True)
        self.mtrl_table.clear()
        self.mtrl_table.setColumnCount(3)
        self.mtrl_table.setRowCount(len(self.qclayers.materials))
        self.mtrl_table.setHorizontalHeaderLabels(["#", "material", "x"])
        self.mtrl_table.verticalHeader().hide()
        # TODO: material name support

        possible_mtrl = tuple(
            [ALLOY_PARAM[m]["name"] for m in QC_MATERIAL[self.qclayers.substrate]]
        )
        for n, mtrl in enumerate(self.qclayers.materials):
            color = self.mtrlcolors[n % len(self.mtrlcolors)]
            name = QTableWidgetItem(str(n + 1))
            name.setTextAlignment(Qt.AlignCenter)
            name.setBackground(color)
            name.setFlags(Qt.NoItemFlags)
            self.mtrl_table.setItem(n, 0, name)

            # Choose from available materials, according to substrate
            #  mtrlItem = QComboBox()
            mtrl_item = MtrlComboBox()
            mtrl_item.addItems(possible_mtrl)
            mtrl_item.setCurrentText(ALLOY_PARAM[mtrl]["name"])
            mtrl_item.currentIndexChanged.connect(
                partial(self.mtrl_table_mtrl_changed, n)
            )
            self.mtrl_table.setCellWidget(n, 1, mtrl_item)

            # Set mole fraction for materials
            # TODO: number only with QItemDelegate
            mtrl_frac = QTableWidgetItem(str(self.qclayers.mole_fracs[n]))
            mtrl_frac.setTextAlignment(Qt.AlignCenter)
            mtrl_frac.setBackground(color)
            self.mtrl_table.setItem(n, 2, mtrl_frac)

        self.mtrl_table.resizeColumnsToContents()
        self.update_mtrl_info()
        self.del_mtrl_button.setEnabled(False)
        self.mtrl_table.blockSignals(False)

    @pyqtSlot()
    def mtrl_table_item_selection_changed(self):
        """SLOT connected to mtrlTable.itemSelectionChanged()"""
        self.del_mtrl_button.setEnabled(len(self.qclayers.materials) > 2)
        self.update_ifr_settings()

    def update_ifr_settings(self):
        if self.qclayers.include_ifr:
            self._ifr_group_box.setVisible(True)
        else:
            self._ifr_group_box.setVisible(False)
        if self.qclayers.custom_ifr:
            self.ifr_def_box.setEnabled(False)
            self.ifr_delta_box.setEnabled(False)
            self.ifr_lambda_box.setEnabled(False)
            # button to remove custom
        if self.ifr_def_box.isChecked():
            assert all(
                x == self.qclayers.mtrl_ifr_delta[0]
                for x in self.qclayers.mtrl_ifr_delta[1:]
            )
            assert all(
                x == self.qclayers.mtrl_ifr_lambda[0]
                for x in self.qclayers.mtrl_ifr_lambda[1:]
            )
            self.ifr_delta_box.setValue(self.qclayers.mtrl_ifr_delta[0])
            self.ifr_lambda_box.setValue(self.qclayers.mtrl_ifr_lambda[0])
        else:
            row = self.mtrl_table.currentRow()
            if (
                row < len(self.qclayers.materials)
                and not self.qclayers.custom_ifr
                and self.qclayers.mtrl_ifr_delta is not None
            ):
                self.ifr_delta_box.setValue(self.qclayers.mtrl_ifr_delta[row])
                self.ifr_lambda_box.setValue(self.qclayers.mtrl_ifr_lambda[row])

    def ifr_settings_check_const(self):
        if all(
            x == self.qclayers.mtrl_ifr_delta[0]
            for x in self.qclayers.mtrl_ifr_delta[1:]
        ) and all(
            x == self.qclayers.mtrl_ifr_lambda[0]
            for x in self.qclayers.mtrl_ifr_lambda[1:]
        ):
            self.ifr_def_box.setChecked(True)

    def mtrl_table_mtrl_changed(self, row, selection):
        """SLOT as partial(self.mtrlTable_mrtlChanged, q)) connected to
        mtrlItem.currentIndexChanged(int)"""
        # not decorated by pyqt because it's not a real slot but a meta-slot
        self.qclayers.set_mtrl(row, QC_MATERIAL[self.qclayers.substrate][selection])
        self.clear_wfs()
        self.qclayers.populate_x()
        self.update_mtrl_info()
        self.update_quantum_canvas()
        self.dirty.emit()

    @pyqtSlot(QTableWidgetItem)
    def mtrl_table_item_changed(self, item):
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
                assert mf <= 1
                self.qclayers.set_mtrl(row, mole_frac=int(100 * mf) / 100)
            except (ValueError, AssertionError):
                # Input is not a number or is larger than 1
                QMessageBox.warning(self, EJ_ERROR, "invalid mole Fraction")
            #  item.setText(str(self.qclayers.moleFracs[row]))
            self.clear_wfs()
            self.qclayers.populate_x()
            self.update_mtrl_info()
            self.update_quantum_canvas()
            self.dirty.emit()
        else:
            # Should never be here
            raise ValueError

    @pyqtSlot()
    def add_mtrl(self):
        """SLOT connected to self.addMtrlButton.clicked()"""
        self.qclayers.add_mtrl()
        self._update_mtrl_list()
        self._update_mtrl_table()
        self._update_layer_table()

    @pyqtSlot()
    def del_mtrl(self):
        """SLOT connected to self.delMtrlButton.clicked()"""
        self.qclayers.del_mtrl(self.mtrl_table.currentRow())
        self._update_mtrl_list()
        self._update_mtrl_table()
        self._update_layer_table()

    @pyqtSlot(int)
    def ifr_def_constant(self, state):
        # state is int from Qt definition but here it's used as a bool
        if state:
            m = len(self.qclayers.materials)
            self.qclayers.mtrl_ifr_lambda = [self.ifr_lambda_box.value()] * m
            self.qclayers.mtrl_ifr_delta = [self.ifr_delta_box.value()] * m
        self.dirty.emit()

    @pyqtSlot(float)
    def input_ifr_lambda(self, ifr_lambda):
        if self.ifr_def_box.isChecked():
            m = len(self.qclayers.materials)
            self.qclayers.mtrl_ifr_lambda = [ifr_lambda] * m
        else:
            mn = self.mtrl_table.currentRow()
            self.qclayers.mtrl_ifr_lambda[mn] = ifr_lambda
        self.qclayers.populate_material()
        if self.qclayers.status.startswith("solved"):
            # following also do self.qclayers.status = 'solved'
            self.qclayers.reset_ifr_cache()
            self.gain_spec_button.setEnabled(False)
            self.update_quantum_canvas()
        self.dirty.emit()

    @pyqtSlot(float)
    def input_ifr_delta(self, ifr_delta):
        if self.ifr_def_box.isChecked():
            m = len(self.qclayers.materials)
            self.qclayers.mtrl_ifr_delta = [ifr_delta] * m
        else:
            mn = self.mtrl_table.currentRow()
            self.qclayers.mtrl_ifr_delta[mn] = ifr_delta
        self.qclayers.populate_material()
        if self.qclayers.status.startswith("solved"):
            # following also do self.qclayers.status = 'solved'
            self.qclayers.reset_ifr_cache()
            self.gain_spec_button.setEnabled(False)
            self.update_quantum_canvas()
        self.dirty.emit()

    def update_mtrl_info(self):
        """Update labels below mtrlTable"""
        mtrl_offset_in_mev = self.qclayers.mtrl_offset * 1000
        self.offset_label.setText(
            f"<center>ΔE<sub>c</sub>: <b>{mtrl_offset_in_mev:6.0f} meV</b></center>"
        )
        self.net_strian_label.setText(
            f"<center>Net Strain: <b>{self.qclayers.net_strain:6.3f}%</b></center>"
        )
        hwlo_in_mev = self.qclayers.avg_hw_lo * 1000
        self.lo_phonon_label.setText(
            f"<center>E<sub>LO</sub>: <b>{hwlo_in_mev:4.1f} meV</b></center>"
        )

    # =========================================================================
    # Quantum Tab Plotting and Plot Control
    # =========================================================================
    def update_quantum_canvas(self):
        """Update the canvas to show band diagram, and if states have been
        solved, draw the wavefunctions"""
        if self.plot_control.zoomed:
            xmin, xmax = self.quantum_canvas.axes.get_xlim()
            ymin, ymax = self.quantum_canvas.axes.get_ylim()
        else:
            xmin = xmax = ymin = ymax = None
        self.quantum_canvas.clear()
        axes = self.quantum_canvas.axes
        plot_potential(
            self.qclayers, axes, self.plot_vl, self.plot_vx, self.plot_lh, self.plot_so
        )

        if self.layer_selected is not None:
            selected_vc = np.ma.masked_where(
                self.qclayers.x_layer_mask(self.layer_selected), self.qclayers.x_vc
            )
            axes.plot(
                self.qclayers.x_points,
                selected_vc,
                "b",
                linewidth=1.5 if self.qclayers.layer_ar[self.layer_selected] else 1,
            )

        if self.qclayers.status in ("basis", "solved", "solved-full"):
            self.wfs = plot_wf(
                self.qclayers,
                self.plot_type,
                self.fill_plot,
                self.state_holder,
                axes=axes,
            )

        if xmin is not None:
            axes.set_xlim(xmin, xmax)
            axes.set_ylim(ymin, ymax)
        else:
            ymin, ymax = axes.get_ylim()
            if ymax - ymin < 0.4:
                axes.set_ylim(ymin - 0.2, ymax + 0.2)
        self.quantum_canvas.draw()

    @pyqtSlot(bool)
    def layer_select_mode(self, checked):
        """SLOT connected to layerSelectButton.clicked[bool]
        to enable layerSelect mode in the canvas.
        """
        self.plot_control.trigger_custom("layerselect")
        if not checked:
            # clear selection
            self.layer_selected = None
            self.layer_table.clearSelection()
            self.update_quantum_canvas()

    @pyqtSlot(bool)
    def pair_select_mode(self, checked):
        """SLOT connected to pairSelectButton.clicked[bool]
        to enable pairSelect mode in the canvas.
        """
        self.plot_control.trigger_custom("pairselect")
        if not checked:
            self.pair_selected = False
            self.state_holder = []
            self.update_quantum_canvas()

    def layer_select(self, event):
        """callback registered in plotControl when it's in wellselect mode.
        It's mpl_connect to button_release_event of quantumCanvas"""
        if event.button == 1:  # left button clicked
            x = event.xdata
            x_layer_num = np.argmin((self.qclayers.x_points - x) ** 2)
            layer_num = self.qclayers.x_layer_nums[x_layer_num]
            self.layer_table.selectRow(layer_num)
            # Trigger itemSelectionChanged SIGNAL and thus update_quantumCanvas

    def clear_wfs(self):
        self.qclayers.status = "unsolved"
        self.state_holder = []
        self.pair_select_button.setEnabled(False)
        self.full_population_button.setEnabled(False)
        self.gain_spec_button.setEnabled(False)
        self.optimize_layer_button.setEnabled(False)
        self.update_quantum_canvas()

    # ===========================================================================
    # Export Functions
    # ===========================================================================
    def export_quantum_canvas(self, filename=None):
        self.plot_control.save_figure(
            "ErwinJr2 - Export Band Structure Image", filename, "png"
        )

    def export_band_data(self, fname):
        fname_base = fname.split(".")
        fname_base = "".join(fname_base[:-1]) if len(fname_base) > 1 else fname
        np.savetxt(
            fname_base + "_CB" + ".csv",
            np.column_stack([self.qclayers.x_points, self.qclayers.x_vc]),
            delimiter=",",
        )

        if self.qclayers.status.startswith("solved"):
            # otherwise band structure hasn't been solved yet
            np.savetxt(
                fname_base + "_WFs" + ".csv",
                np.column_stack([self.qclayers.x_points, self.qclayers.psis.T]),
                delimiter=",",
            )
            ys = scale_wf(self.qclayers, self.plot_type).T + self.qclayers.eigen_es
            np.savetxt(
                fname_base + ".csv",
                np.column_stack([self.qclayers.x_points, ys]),
                delimiter=",",
            )
            np.savetxt(
                fname_base + "_Es" + ".csv", self.qclayers.eigen_es, delimiter=","
            )
            if self.qclayers.status == "solved-full":
                pop = np.array(
                    [
                        (
                            self.qclayers.state_population(n)
                            if self.qclayers.period_map[n] is not None
                            else np.NaN
                        )
                        for n in range(len(self.qclayers.eigen_es))
                    ]
                )
                np.savetxt(fname_base + "_population" + ".csv", pop, delimiter=",")

    # ===========================================================================
    # View Band Items
    # ===========================================================================
    def view_band(self, band):
        # TODO: combine following slot to this one
        self.bandFlags[band] = not self.bandFlags[band]
        self.update_quantum_canvas()

    def view_vx_band(self):
        self.plot_vx = not self.plot_vx
        self.update_quantum_canvas()

    def view_vl_band(self):
        self.plot_vl = not self.plot_vl
        self.update_quantum_canvas()

    def view_lh_band(self):
        self.plot_lh = not self.plot_lh
        self.update_quantum_canvas()

    def view_so_band(self):
        self.plot_so = not self.plot_so
        self.update_quantum_canvas()

    def set_plotwf(self):
        self.plot_type = "wf" if self.plot_type != "wf" else "mode"
        self.update_quantum_canvas()

    def set_fill(self):
        if self.fill_plot:
            self.fill_plot = False
        else:
            self.fill_plot = 0.3
        self.update_quantum_canvas()

    # ===========================================================================
    # Model triggering
    # ===========================================================================
    def trigger_ifr(self):
        self.qclayers.include_ifr = not self.qclayers.include_ifr
        # TODO: dynamically hide IFR block
        if self.qclayers.include_ifr:
            self._ifr_group_box.setVisible(True)
        else:
            self._ifr_group_box.setVisible(False)

    def trigger_solver(self, solver):
        if solver != self.qclayers.solver:
            self.qclayers.solver = solver
            self.clear_wfs()
            self.algo_param_update()

    def algo_param_update(self):
        # This is causing "has active key-value observers (KVO)" when recreate
        # window warning, but we don't see any real issue.
        if self.qclayers.solver == "matrix":
            self.e_res_label.setVisible(False)
            self.input_e_res_box.setVisible(False)
            self.e_count_label.setVisible(True)
            self.input_e_count_box.setVisible(True)
        else:
            self.e_res_label.setVisible(True)
            self.input_e_res_box.setVisible(True)
            self.e_count_label.setVisible(False)
            self.input_e_count_box.setVisible(False)

    # ===========================================================================
    # Calculations
    # ===========================================================================
    def _thread_run(self, calc, post_calc=None):
        """This is a helper function for multiple threading to avoid GUI
        freeze."""
        # TODO: block all change
        self.calc_repaint(True)
        assert not self._calc_thread.isRunning()
        self._worker = CalculateHolder(calc)
        self._worker.moveToThread(self._calc_thread)
        self._calc_thread.started.connect(self._worker.run)
        self._worker.finished.connect(self._calc_thread.quit)
        self._worker.finished.connect(self._worker.deleteLater)
        self._worker.failed.connect(
            lambda traceInfo: QMessageBox.warning(self, EJ_ERROR, traceInfo)
        )
        self._worker.warning.connect(
            lambda message: QMessageBox.warning(self, EJ_WARNING, message)
        )
        self._worker.failed.connect(self.clear_wfs)
        self._calc_thread.finished.connect(lambda: self.calc_repaint(False))
        # self.calcThread.finished.connect(self.thread.deleteLater)
        if post_calc is not None:
            self._worker.succeed.connect(post_calc)
        self._calc_thread.start()

    def calc_repaint(self, is_doing):
        """UI repaint for doing calculating"""
        for button in (
            self.solve_whole_button,
            self.solve_basis_button,
            self.pair_select_button,
        ):
            button.setEnabled(not is_doing)
            button.repaint()
        if is_doing:
            self.to_optics_button.setEnabled(False)
        if self.pair_selected:
            self.fom_button.setEnabled(not is_doing)
            self.fom_button.repaint()
        if self.qclayers.status.startswith("solved"):
            self.full_population_button.setEnabled(not is_doing)
            self.full_population_button.repaint()

    def _solve_whole(self):
        self.qclayers.solve_whole()
        try:
            self.period_set = self.qclayers.period_recognize()
        except StateRecognizeError as e:
            self.period_set = set()
            raise e

    def _solved(self):
        self.optimize_layer_button.setEnabled(True)
        self.update_quantum_canvas()

    @pyqtSlot()
    def solve_whole(self):  # solves whole structure
        """SLOT connected to solveWholeButton.clicked(): Whole solver"""
        self.clear_wfs()
        self._thread_run(self._solve_whole, self._solved)

    def _solve_basis(self):
        self.qclayers.solve_basis()

    @pyqtSlot()
    def solve_basis(self):  # solves structure with basis
        """SLOT connected to solveBasisButton.clicked(): Basis solver"""
        self.clear_wfs()
        self._thread_run(self._solve_basis, self.update_quantum_canvas)

    def state_pick(self, event):
        """Callback registered in plotControl when it's in pairselect mode.
        It's mpl_connect to button_release_event of quantumCanvas"""
        if self.qclayers.status == "unsolved":
            # Not yet solved
            raise ValueError("Structure not solved yet.")

        if event.button == 1:  # left button clicked: select a state
            if len(self.state_holder) >= 2:
                # start new pair selection
                self.state_holder = []
                self.pair_selected = False

            # TODO: replace following by matplotlib lines picker
            x = event.xdata
            y = event.ydata
            x_data = np.tile(self.qclayers.x_points, (self.qclayers.psis.shape[0], 1)).T
            if self.plot_type in ("mode", "wf"):
                y_data = self.qclayers.eigen_es + self.wfs.T
            else:
                y_data = self.qclayers.eigen_es

            width, height = self.quantum_canvas.axes.bbox.size
            xmin, xmax = self.quantum_canvas.axes.get_xlim()
            ymin, ymax = self.quantum_canvas.axes.get_ylim()
            x_scale = (xmax - xmin) / width
            y_scale = (ymax - ymin) / height

            r = np.nanmin(
                sqrt(((x_data - x) / x_scale) ** 2 + ((y_data - y) / y_scale) ** 2),
                axis=0,
            )
            # penalty for x axis away from visible region
            if self.plot_type in ("mode", "wf"):
                x0 = np.argmax(np.abs(self.wfs) > plotconfig["wf_almost_zero"], axis=1)
                xl = self.wfs.shape[1] - np.argmax(
                    np.abs(self.wfs[:, ::-1]) > plotconfig["wf_almost_zero"], axis=1
                )
                out_of_signt = np.arange(self.wfs.shape[0])[
                    (x0 * self.qclayers.x_step > x + 5)
                    | (xl * self.qclayers.x_step < x - 5)
                ]
                r[out_of_signt] += np.inf
            ss = np.nanargmin(r)
            if len(self.state_holder) == 1 and self.state_holder[0] == ss:
                r[ss] = np.nan
                ss = np.nanargmin(r)
            self.state_holder.append(ss)
        elif event.button == 3:  # right button clicked: remove last selection
            if len(self.state_holder) == 2:
                self.pair_selected = False
            self.state_holder.pop()
        else:
            return
        self.update_quantum_canvas()

        if len(self.state_holder) == 1:
            self.state_param_text.clear()
            self.pair_string = f"selected: {self.state_holder[0]}, ..<br>"
            if self.qclayers.status == "solved-full":
                pop = self.qclayers.state_population(self.state_holder[0])
                if pop is not None:
                    self.pair_string += f"population: {pop * 100:.1f}%<br>"
            self.state_param_text.setText(self.pair_string)

        elif len(self.state_holder) == 2:
            self.pair_selected = True
            self.update_selected()

    def update_selected(self):
        """Update easy-to-calculate parameters regarding picked two states"""
        self.fom_button.setEnabled(True)
        upper = self.state_holder[0]
        lower = self.state_holder[1]
        if self.qclayers.eigen_es[upper] < self.qclayers.eigen_es[lower]:
            upper, lower = lower, upper
        self.de = self.qclayers.eigen_es[upper] - self.qclayers.eigen_es[lower]
        self.wavelength = h * c0 / (e0 * self.de) * 1e6  # um

        if self.qclayers.status == "unsolved":
            self.fom_button.setEnabled(False)
            self.pair_string = ""
        else:
            optical_dipole = self.qclayers.dipole(upper, lower)
            tau_lo_ul = 1 / self.qclayers.lo_transition(upper, lower)
            if self.qclayers.include_ifr:
                self.tau_ifr_ul = 1 / self.qclayers.ifr_transition(upper, lower)
            self.pair_string = (
                f"selected: {self.state_holder[0]}, {self.state_holder[1]}<br>"
                f"energy diff: <b>{1000 * self.de:6.1f} meV</b> "
                f"({self.wavelength:6.1f} \u00b5m)<br>"
                f"dipole: <b>{optical_dipole:6.1f} \u212b</b><br>"
            )
            self.pair_string += f"LO scatter: {tau_lo_ul:6.3g} ps<br>"
            if self.qclayers.include_ifr:
                self.pair_string += f"IFR scatter: {self.tau_ifr_ul:6.3g} ps<br>"
            if self.qclayers.status == "solved-full":
                self.pair_string += "population: <br>&nbsp;&nbsp;&nbsp;"
                pip_u = self.qclayers.state_population(upper)
                pip_u = "N/A" if pip_u is None else f"{100 * pip_u:.1f}%"
                pop_l = self.qclayers.state_population(lower)
                pop_l = "N/A" if pop_l is None else f"{100 * pop_l:.1f}%"
                self.pair_string += f"{upper}: {pip_u} {lower}: {pop_l}<br>"
        self.state_param_text.clear()
        self.state_param_text.setText(self.pair_string)

    def _calc_fom(self):
        upper = self.state_holder[0]
        lower = self.state_holder[1]
        if self.qclayers.eigen_es[upper] < self.qclayers.eigen_es[lower]:
            upper, lower = lower, upper
        tau_lo_u = self.qclayers.lo_lifetime(upper)
        tau_lo_l = self.qclayers.lo_lifetime(lower)
        if self.qclayers.include_ifr:
            tau_ifr_u = self.qclayers.ifr_lifetime(upper)
            tau_ifr_l = self.qclayers.ifr_lifetime(lower)
        # tau_u = 1/(1/tauLO_u + 1/tauIFR_u)
        # tau_l = 1/(1/tauLO_l + 1/tauIFR_l)
        fom = self.qclayers.figure_of_merit(upper, lower)
        self.neff = self.qclayers.effective_ridx(self.wavelength)
        # 1E-27 = 1E-32 * 1E5
        # 1E-32 angstrom^2 ps -> m^2 s, 1E5 m/A -> cm/kA
        gamma = self.qclayers.dephasing(upper, lower)
        self.gain_coeff = (
            e0
            * fom
            * 1e-27
            * self.de
            / (gamma * hbar * c0 * eps0 * self.neff * self.qclayers.period_l * 1e-10)
        )
        # tauUpperLower is the inverse of transition rate (lifetime)

        self.fom_string = (
            f"<i>\u03c4<sup>LO</sup><sub>upper</sub></i> : {tau_lo_u:6.3f} ps<br>"
            f"<i>\u03c4<sup>LO</sup><sub>lower</sub></i> : {tau_lo_l:6.3f} ps<br>"
        )
        if self.qclayers.include_ifr:
            self.fom_string += (
                f"<i>\u03c4<sup>IFR</sup><sub>upper</sub></i> : {tau_ifr_u:6.3f} ps<br>"
                f"<i>\u03c4<sup>IFR</sup><sub>lower</sub></i> : {tau_ifr_l:6.3f} ps<br>"
            )
        self.fom_string += (
            f"FoM: <b>{fom:6.0f} ps \u212b<sup>2</sup></b><br>"
            f"Gain coefficient:<br>&nbsp;&nbsp; {self.gain_coeff:.2f} cm/kA"
        )

    def _update_fom(self):
        self.state_param_text.setText(self.pair_string + self.fom_string)
        self.state_param_text.verticalScrollBar().setValue(
            self.state_param_text.verticalScrollBar().maximum()
        )
        self.fom_button.setEnabled(True)
        self.to_optics_button.setEnabled(True)

    @pyqtSlot()
    def update_fom(self):
        """SLOT connected to FoMButton.clicked()
        Calculate Figure of merit."""
        if len(self.state_holder) < 2:
            print("Warning: FoM button triggered before state selection.")
            return
        self.fom_string = ""
        self._thread_run(self._calc_fom, self._update_fom)

    def _full_population(self):
        self.state_holder = []
        self.qclayers.full_population()
        self.wavelength = self.input_wl_box.value()
        self.neff = self.qclayers.effective_ridx(self.wavelength)
        gain = self.qclayers.full_gain_spectrum(self.wavelength)
        self.gain_coeff = gain / self.qclayers.current
        self.spec_string = (
            "Carrier Leakage: <br>"
            f"&nbsp;&nbsp; {self.qclayers.carrier_leak * 100:.2f}%<br>"
            f"Current: {self.qclayers.current:.1f} kA/cm<sup>2</sup><br>"
            f"Gain coefficient:<br>&nbsp;&nbsp; {self.gain_coeff:.2f} cm/kA"
        )

    def _update_population(self):
        self.state_param_text.setText(self.spec_string)
        self.fom_button.setEnabled(False)
        self.gain_spec_button.setEnabled(True)
        self.to_optics_button.setEnabled(True)
        self.update_quantum_canvas()

    @pyqtSlot()
    def full_population(self):
        """SLOT connected to fullPopulationButton.clicked()"""
        if not self.qclayers.status.startswith("solved"):
            print("Full_population triggered before solving.")
        self.state_param_text.clear()
        self._thread_run(self._full_population, self._update_population)

    @pyqtSlot()
    def pop_gain_spec_window(self):
        """SLOT connected to gainSpecButton.clicked()"""
        wlmin = 1.5
        wlmax = int(4 * self.qclayers.wl + 0.5) / 2
        wl_dialog = WLDialog(self, wlmin, wlmax)
        wlmin, wlmax, button_res = wl_dialog.exec()
        if button_res:
            wl_figure = figure()
            wl_figure.canvas.manager.set_window_title("Gain Spectrum")
            axes = wl_figure.add_subplot(111)
            wls = np.linspace(wlmin, wlmax, 500)
            axes.plot(wls, self.qclayers.full_gain_spectrum(wls))
            axes.axhline(0, ls="--", lw=0.5)
            axes.set_xlabel("Wavelength (μm)")
            axes.set_ylabel("Gain (cm$^{-1}$)")
            if axes.get_ylim()[0] < -150:
                axes.set_ylim(bottom=-150)
            wl_figure.show()


class WLDialog(QDialog):
    """Dialog window for setting the wavelength range for gain spectrum"""

    def __init__(self, parent, wl_min, wl_max):
        super().__init__(parent)
        self.setModal(True)
        self.setWindowTitle("ErwinJr2: Gain Spectrum Range")
        self.wl_min_box = QDoubleSpinBox()
        self.wl_max_box = QDoubleSpinBox()
        for box in (self.wl_min_box, self.wl_max_box):
            box.setDecimals(1)
            box.setRange(0.5, 100)
            box.setSingleStep(1)
            box.setSuffix(" μm")
        self.wl_min_box.setValue(wl_min)
        self.wl_max_box.setValue(wl_max)
        main_layout = QGridLayout()
        main_layout.addWidget(QLabel("min λ: "), 0, 0)
        main_layout.addWidget(self.wl_min_box, 0, 1)
        main_layout.addWidget(QLabel("max λ: "), 1, 0)
        main_layout.addWidget(self.wl_max_box, 1, 1)
        button_box = QDialogButtonBox(
            QDialogButtonBox.Ok | QDialogButtonBox.Cancel, Qt.Horizontal
        )
        button_box.accepted.connect(self.accept)
        button_box.rejected.connect(self.reject)
        main_layout.addWidget(button_box, 2, 0, 1, 2)
        self.setLayout(main_layout)

    def exec(self):
        res = super().exec()
        return self.wl_min_box.value(), self.wl_max_box.value(), res
