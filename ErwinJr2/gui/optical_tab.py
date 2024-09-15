"""
This file defines the optical tab of ErwinJr2, for simulation and optimization
of 1D waveguide
"""

import sys
from functools import partial

import numpy as np
from numpy import log, pi
from PyQt5.QtCore import Qt, pyqtSignal, pyqtSlot
from PyQt5.QtGui import QPalette
from PyQt5.QtWidgets import (
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
    QSpinBox,
    QTableWidget,
    QTableWidgetItem,
    QTextEdit,
    QVBoxLayout,
    QWidget,
)

from ErwinJr2.gui.constants import EJ_ERROR
from ErwinJr2.gui.custom_qt_class import MtrlComboBox
from ErwinJr2.gui.ej_canvas import EJCanvas
from ErwinJr2.opt_strata import (
    ALLOY_MAP,
    DOPABLE_MATERIALS,
    OptStrata,
    optimize_opt_strata,
    rIdx,
)

mtrlList = list(rIdx.keys()) + list(ALLOY_MAP.keys())
facetList = ("cleaved", "perfect AR", "perfect HR", "custom")


class OpticalTab(QWidget):
    """The Optical Tab of ErwinJr2. This is designed to be a GUI wrapper of
    the class Stratum
    Member variable
        stratum: A class to describe the physics of 1D optical slab waveguide

        customized SIGNAL
            dirty: Show if there's any new changes to qclayers that need to
                   be saved

        plotImag: bool, whether to plot the imaginary part of the refractive
            index

        redActive: bool, whether to plot the active part red

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
        super().__init__(parent)
        self.stratum = stratum if stratum else OptStrata(3.0)
        self.periods = {}
        self.ridge_length = 3.0
        self.ridge_loss = 0.0
        self.facet = [facetList[0], facetList[0]]
        self.beta = None
        self.select = None
        self.red_active = False
        self.plot_imag = True

        self.setAutoFillBackground(True)
        self.setBackgroundRole(QPalette.Window)

        if sys.platform.startswith("win"):
            setting_box_width = 200
            strate_box_width = 350
        elif sys.platform.startswith("darwin"):
            setting_box_width = 350
            strate_box_width = 370
        elif sys.platform.startswith("linux"):
            setting_box_width = 350
            strate_box_width = 350
        else:
            setting_box_width = 350
            strate_box_width = 375

        optical_layout = QHBoxLayout()
        optical_layout.setSpacing(0)
        optical_layout.addLayout(self._setting_box(setting_box_width))
        optical_layout.addLayout(self._strata_box(strate_box_width))
        figure_box = QVBoxLayout()
        self.opt_canvas = EJCanvas(
            xlabel="Position $x$ (μm)", ylabel="Refractive index $n$"
        )
        self.ridx_axes = self.opt_canvas.axes
        self.mode_axes = self.ridx_axes.twinx()
        self.mode_axes.set_frame_on(False)
        self.mode_axes.get_yaxis().set_visible(False)
        # self.modeAxis = self.ridxAxis
        # self.modeAxis.set_ylabel("TM Mode $E_y$ (a.u.)", color='C0',
        #                          fontsize=canvasConfig["fontsize"])
        # self.modeAxis.tick_params(axis='y', labelcolor='C0')
        figure_box.addWidget(self.opt_canvas)
        optical_layout.addLayout(figure_box)

        self.setLayout(optical_layout)
        self.dirty.connect(self.update_xs)
        self.update_xs()
        self.update_loss()

    def _setting_box(self, width):
        """Return a Qt Layout object containing all setting widgets"""
        setting_box = QVBoxLayout()

        setting_box.addWidget(QLabel("<b><center>Wavelength</center></b>"))
        self.wl_box = QDoubleSpinBox()
        self.wl_box.setDecimals(1)
        self.wl_box.setValue(self.stratum.wl)
        self.wl_box.setSuffix(" μm")
        self.wl_box.setRange(0.0, 100.0)
        self.wl_box.setMaximumWidth(width)
        setting_box.addWidget(self.wl_box)

        mtrl_group_box = QGroupBox("Custom Material")
        mtrl_group_box.setMaximumWidth(width)
        mtrl_layout = QGridLayout()
        mtrl_layout.addWidget(QLabel("<center>Material Name</center>"), 0, 0, 1, 2)
        self.mtrl_box = QComboBox()
        self.mtrl_box.addItems(self.stratum.cstm_idx.keys())
        mtrl_layout.addWidget(self.mtrl_box, 1, 0, 1, 2)
        mtrl_layout.addWidget(QLabel("<center>index n<sub>eff</sub></center>"), 2, 0)
        self.ridx_real_box = QDoubleSpinBox()
        self.ridx_real_box.setSingleStep(0.1)
        mtrl_layout.addWidget(self.ridx_real_box, 3, 0)
        mtrl_layout.addWidget(QLabel("<center>passive loss k</center>"), 2, 1)
        self.ridx_imag_box = QDoubleSpinBox()
        self.ridx_imag_box.setSingleStep(0.1)
        mtrl_layout.addWidget(self.ridx_imag_box, 3, 1)
        mtrl_layout.addWidget(QLabel("<center>Period Length</center>"), 4, 0)
        self.period_box = QDoubleSpinBox()
        self.period_box.setSuffix(" Å")
        self.period_box.setMaximum(19999.99)
        self.repeat_box = QSpinBox()
        self.repeat_box.setMaximum(999)
        mtrl_layout.addWidget(self.period_box, 5, 0)
        mtrl_layout.addWidget(QLabel("<center>Periods</center>"), 4, 1)
        mtrl_layout.addWidget(self.repeat_box, 5, 1)
        mtrl_layout.addWidget(QLabel("<center>Active property</center>"), 6, 0, 1, 2)
        mtrl_layout.addWidget(QLabel("<center>Bias Field</center>"), 7, 0)
        self.field_box = QDoubleSpinBox()
        self.field_box.setSuffix(" kV/cm")
        self.field_box.setEnabled(False)
        self.field_box.setDecimals(1)
        self.field_box.setMaximum(500)
        mtrl_layout.addWidget(self.field_box, 8, 0)
        mtrl_layout.addWidget(QLabel("<center>Gain Coeff.</center>"), 7, 1)
        self.gain_coeff_box = QDoubleSpinBox()
        self.gain_coeff_box.setSuffix(" cm/kA")
        self.gain_coeff_box.setEnabled(False)
        self.gain_coeff_box.setDecimals(1)
        mtrl_layout.addWidget(self.gain_coeff_box, 8, 1)
        # TODO
        self.update_mtrl()
        mtrl_layout.setHorizontalSpacing(1)
        mtrl_layout.setVerticalSpacing(3)
        mtrl_group_box.setLayout(mtrl_layout)
        setting_box.addWidget(mtrl_group_box)

        self.wl_box.valueChanged[float].connect(self.input_wl)
        self.mtrl_box.currentIndexChanged[int].connect(self.update_mtrl)
        self.ridx_real_box.valueChanged[float].connect(self.input_ridx)
        self.ridx_imag_box.valueChanged[float].connect(self.input_alpha)
        self.period_box.valueChanged[float].connect(self.input_period)
        self.repeat_box.valueChanged[int].connect(self.input_repeats)
        self.gain_coeff_box.valueChanged[float].connect(self.input_gain_coeff)

        ridge_box = QGroupBox("Ridge Geometry")
        ridge_box.setMaximumWidth(width)
        ridge_layout = QGridLayout()
        self.facet_box = [None] * 2
        self.ref_box = [None] * 2
        for n in (0, 1):
            ridge_layout.addWidget(QLabel(f"<center>Facet{n+1}</center>"), 0, n)
            self.facet_box[n] = QComboBox()
            self.facet_box[n].addItems(facetList)
            self.facet_box[n].setCurrentText(self.facet[n])
            ridge_layout.addWidget(self.facet_box[n], 1, n)
            self.ref_box[n] = QDoubleSpinBox()
            self.ref_box[n].setDecimals(1)
            self.ref_box[n].setRange(0.0, 100.0)
            self.ref_box[n].setSuffix(" %")
            self.ref_box[n].setValue(self.facet_ref_ct(n))
            ridge_layout.addWidget(self.ref_box[n], 3, n)
            self.facet_box[n].currentIndexChanged[int].connect(
                partial(self.input_facet, n)
            )
            self.ref_box[n].valueChanged[float].connect(partial(self.input_ref, n))
            if self.facet[n] == "custom":
                self.ref_box[n].setEnabled(True)
                self.ref_box[n].blockSignals(False)
            else:
                self.ref_box[n].setEnabled(False)
                self.ref_box[n].blockSignals(True)
        ridge_layout.addWidget(QLabel("<center>Reflectivity</center>"), 2, 0, 1, 2)

        ridge_layout.addWidget(QLabel("<center>Ridge Length</center>"), 4, 0)
        self.ridge_length_box = QDoubleSpinBox()
        self.ridge_length_box.setDecimals(1)
        self.ridge_length_box.setRange(0.0, 20.0)
        self.ridge_length_box.setSingleStep(1)
        self.ridge_length_box.setValue(self.ridge_length)
        self.ridge_length_box.setSuffix(" mm")
        self.ridge_length_box.valueChanged[float].connect(self.input_ridge_l)
        ridge_layout.addWidget(self.ridge_length_box, 5, 0)
        ridge_layout.addWidget(QLabel("<center>Mirror Loss: </center>"), 4, 1)
        self.mirror_loss = QLabel("<center>0.0 cm<sup>-1</sup></center>")
        ridge_layout.addWidget(self.mirror_loss, 5, 1)
        ridge_layout.setHorizontalSpacing(1)
        ridge_layout.setVerticalSpacing(3)
        ridge_box.setLayout(ridge_layout)
        setting_box.addWidget(ridge_box)

        setting_box.addStretch()
        return setting_box
        # _settingBox end

    def _strata_box(self, width):
        """Return a Qt Layout object containing all strata table widgets"""
        strata_box = QGridLayout()
        strata_box.setSpacing(5)
        self.insert_button = QPushButton("Insert Strata")
        self.insert_button.clicked.connect(self.insert_strata)
        self.delete_button = QPushButton("Delete Strata")
        self.delete_button.clicked.connect(self.delete_strata)
        strata_box.addWidget(self.insert_button, 0, 0)
        strata_box.addWidget(self.delete_button, 0, 1)
        self.strata_table = QTableWidget()
        self.strata_table.setMaximumWidth(width)
        self.strata_table.setMinimumWidth(width)
        self.strata_table.setMinimumHeight(375)
        self.strata_table.setMaximumHeight(500)
        self.strata_table.setSelectionBehavior(QTableWidget.SelectRows)
        self.strata_table.setSelectionMode(QTableWidget.SingleSelection)
        self.strata_table.itemChanged.connect(self.strata_table_item)
        self.strata_table.itemSelectionChanged.connect(self.strata_table_select)
        self.strata_table_refresh()
        strata_box.addWidget(self.strata_table, 1, 0, 1, 2)
        self.solve_button = QPushButton("Solve Mode")
        self.solve_button.clicked.connect(self.solve)
        self.optimize_button = QPushButton("Optimize Strata")
        self.optimize_button.clicked.connect(self.optimize_strata)
        self.optimize_button.setEnabled(False)
        strata_box.addWidget(self.solve_button, 2, 0)
        strata_box.addWidget(self.optimize_button, 2, 1)
        strata_box.addWidget(QLabel("Performance"), 3, 0, 1, 2)
        self.info_box = QTextEdit()
        self.info_box.setReadOnly(True)
        self.info_box.setMinimumWidth(width)
        self.info_box.setMaximumWidth(width)
        self.info_box.setMinimumHeight(120)
        strata_box.addWidget(self.info_box, 4, 0, 1, 2)

        strata_box.setRowStretch(5, 0)
        return strata_box

    @pyqtSlot(float)
    def input_wl(self, wl):
        """SLOT connected to self.wlBox.valueChanged[float]"""
        self.stratum.set_wl(wl)
        self.dirty.emit()

    @pyqtSlot(int)
    def update_mtrl(self, idx=0):
        """SLOT connect to mtrlsBox.currentIndexChanged[int]"""
        if self.stratum.cstm_idx:
            self.mtrl_box.setCurrentIndex(idx)
            mtrl = self.mtrl_box.currentText()
            mtrl_idx = self.stratum.cstm_idx[mtrl]
            self.ridx_real_box.setValue(mtrl_idx.real)
            self.ridx_imag_box.setValue(mtrl_idx.imag)
            if mtrl in self.stratum.cstm_prd and self.stratum.cstm_prd[mtrl][0] > 0:
                self.period_box.setValue(self.stratum.cstm_prd[mtrl][0])
                self.repeat_box.setValue(self.stratum.cstm_prd[mtrl][1])
            else:
                self.period_box.setValue(0.0)
                self.repeat_box.setValue(0)
                self.repeat_box.setEnabled(False)
            if mtrl in self.stratum.cstm_gain:
                self.gain_coeff_box.setValue(self.stratum.cstm_gain[mtrl])
        else:
            self.ridx_real_box.setValue(1.0)
            self.ridx_real_box.setEnabled(False)
            self.ridx_imag_box.setValue(0.0)
            self.ridx_imag_box.setEnabled(False)
            self.period_box.setValue(0.0)
            self.period_box.setEnabled(False)
            return

    @pyqtSlot(float)
    def input_ridx(self, value):
        """SLOT connected to self.rIdxRealBox.valueChanged[float]"""
        mtrl = self.mtrl_box.currentText()
        alpha = self.stratum.cstm_idx[mtrl].imag
        self.stratum.cstm_idx[mtrl] = value + 1j * alpha
        self.dirty.emit()

    @pyqtSlot(float)
    def input_alpha(self, value):
        """SLOT connected to self.rIdxImagBox.valueChanged[float]"""
        mtrl = self.mtrl_box.currentText()
        neff = self.stratum.cstm_idx[mtrl].real
        self.stratum.cstm_idx[mtrl] = neff + 1j * value
        self.dirty.emit()

    @pyqtSlot(float)
    def input_period(self, value):
        """SLOT connected to self.periodBox.valueChanged[float]"""
        mtrl = self.mtrl_box.currentText()
        if mtrl not in self.stratum.cstm_prd:
            self.stratum.cstm_prd[mtrl] = [value, 1]
        else:
            self.stratum.cstm_prd[mtrl][0] = value
        if value == 0:
            self.repeat_box.setEnabled(False)
        else:
            self.repeat_box.setEnabled(True)
            self.update_custom_length()

    @pyqtSlot(int)
    def input_repeats(self, value):
        """SLOT connected to self.repeatBox.valueChanged[int]"""
        mtrl = self.mtrl_box.currentText()
        if mtrl not in self.stratum.cstm_prd:
            self.stratum.cstm_prd[mtrl] = [1, value]
        else:
            self.stratum.cstm_prd[mtrl][1] = value
        self.update_custom_length()

    @pyqtSlot(float)
    def input_gain_coeff(self, value):
        """SLOT connected to self.gainCoeffBox.valueChanged[float]"""
        mtrl = self.mtrl_box.currentText()
        self.stratum.cstm_gain[mtrl] = value

    def setup_active(self, wl, e_field, gain_coeff, ridx, lp):
        """Interface to get parameters from quantum tab"""
        self.wl_box.setValue(wl)
        mtrl = "Active Core"
        # self.stratum.wl = wl
        self.stratum.cstm_idx[mtrl] = ridx
        if "Active Core" in self.stratum.cstm_prd:
            self.stratum.cstm_prd[mtrl][0] = lp
        else:
            self.stratum.cstm_prd[mtrl] = [lp, 1]
        self.stratum.cstm_gain[mtrl] = gain_coeff

        # set up active core
        if self.mtrl_box.findText(mtrl) == -1:
            self.mtrl_box.addItem(mtrl)
        self.mtrl_box.setCurrentText(mtrl)
        self.field_box.setValue(e_field)
        self.period_box.setValue(lp)
        self.gain_coeff_box.setValue(gain_coeff)
        self.update_custom_length()
        self.update_mtrl()

    def update_custom_length(self):
        length = self.period_box.value() * self.repeat_box.value() / 1e4  # to um
        for n, mtrl in enumerate(self.stratum.materials):
            if mtrl == self.mtrl_box.currentText():
                self.stratum.layer_widths[n] = length
        self.dirty.emit()

    def input_facet(self, n, idx):
        """SLOT as partial(self.input_facet, n) connected to
        facetBox[n].currentIndexChanged(int)"""
        self.facet[n] = facetList[idx]
        self.ref_box[n].setValue(self.facet_ref_ct(n) * 100)
        if self.facet[n] == "custom":
            self.ref_box[n].setEnabled(True)
            self.ref_box[n].blockSignals(False)
        else:
            self.ref_box[n].setEnabled(False)
            self.ref_box[n].blockSignals(True)
        self.update_loss()
        self.dirty.emit()

    def input_ref(self, n, ref):
        """SLOT as partial(self.input_facet, n) connected to
        facetBox[n].currentIndexChanged(int)"""
        _ = n, ref
        self.update_loss()
        self.dirty.emit()

    def facet_ref_ct(self, n):
        """Return the reflectivity of facet n"""
        if self.facet[n] == "cleaved":
            if self.beta is None:
                return -1
            reflectivity = abs((self.beta.real - 1) / (self.beta.real + 1)) ** 2
            self.ref_box[n].setValue(100 * reflectivity)
            return reflectivity
        if self.facet[n] == "perfect AR":
            return 1e-9
        if self.facet[n] == "perfect HR":
            return 1
        if self.facet[n] == "custom":
            return self.ref_box[n].value() / 100
        else:
            raise ValueError(f"Wrong facet {self.facet[n]}")

    @pyqtSlot(float)
    def input_ridge_l(self, value):
        """SLOT connected to ridgeLengthBox.valueChanged[float]"""
        self.ridge_length = value
        self.update_loss()
        self.dirty.emit()

    def update_loss(self):
        """Update the mirror loss, should be called whenever the facet
        settings are changed"""
        per_fun_loss = self.facet_ref_ct(1) * self.facet_ref_ct(0)
        self.alpha_m = -log(per_fun_loss) / (2 * self.ridge_length / 10)  # to cm-1
        self.mirror_loss.setText(f"<center>{self.alpha_m:.1f} cm<sup>-1</sup></center>")

    @pyqtSlot()
    def strata_table_select(self):
        """SLOT connected to SIGNAL strataTable.itemSelectionChanged"""
        row = self.strata_table.currentRow()
        if row < len(self.stratum.materials):
            self.select = row
        else:
            self.select = None
            self.stratTable.clearSelection()
        self.update_canvas()

    @pyqtSlot(QTableWidgetItem)
    def strata_table_item(self, item: QTableWidgetItem):
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
                QMessageBox.warning(
                    self, EJ_ERROR, "This value should be a non-negative number"
                )
                self.strata_table_refresh()
                return
            if column == 1:
                # mole fraction
                if value > 1:
                    QMessageBox.warning(
                        self, EJ_ERROR, "Mole Fraction must be between 0 and 1"
                    )
                    self.strata_table.setItem(
                        row,
                        column,
                        QTableWidgetItem(f"{self.stratum.mole_fracs[row]:.2f}"),
                    )
                    return
                self.stratum.mole_fracs[row] = value
            if column == 2:
                # width
                self.stratum.layer_widths[row] = value
            if column == 3:
                # doping
                self.stratum.dopings[row] = value

        self.stratum.update_indices()
        self.dirty.emit()

    @pyqtSlot()
    def insert_strata(self):
        """SLOT connected to self.insertButton.clicked()"""
        row = self.strata_table.currentRow()
        if row >= len(self.stratum.materials) or row < 0:
            # Add new lines in the last strata
            row = len(self.stratum.materials) - 1
        elif row == 0:
            row = 1
        # len(self.stratum.materials) is at least 2 for top and substrate
        # s.t. row >= 1 and <= len(self.stratum.materials) - 1
        self.stratum.insert(row)
        self.dirty.emit()
        self.strata_table.selectRow(row)

    @pyqtSlot()
    def delete_strata(self):
        """SLOT connected to self.deleteButton.clicked()"""
        row = self.strata_table.currentRow()
        if row >= len(self.stratum.materials) - 1 or row < 1:
            return
        # s.t. len(self.stratum.materials) > 2 and Ls is not empty
        self.stratum.delete(row)
        self.dirty.emit()
        self.strata_table.selectRow(row)

    def strata_table_refresh(self):
        """Update strataTable content, should be called whenever the strata
        structure is changed"""
        self.strata_table.blockSignals(True)
        self.stratum.update_indices()
        self.strata_table.clear()
        self.strata_table.setColumnCount(5)
        self.strata_table.setHorizontalHeaderLabels(
            ["Material", "x", "Width", "Doping", "n"]
        )
        self.strata_table.setRowCount(len(self.stratum.materials))
        self.strata_table.setVerticalHeaderLabels(
            str(n + 1) for n in range(len(self.stratum.materials))
        )

        for q, mtrl in enumerate(self.stratum.materials):
            # material name
            mtrl_name = MtrlComboBox()
            mtrl_name.addItems(mtrlList)
            mtrl_name.addItems(self.stratum.cstm_idx.keys())
            mtrl_name.setCurrentText(mtrl)
            mtrl_name.currentTextChanged.connect(
                partial(self.strata_table_mtrl_changed, q)
            )
            self.strata_table.setCellWidget(q, 0, mtrl_name)

            # Mole Frac
            if mtrl not in ALLOY_MAP:
                mole_frac = QTableWidgetItem("N/A")
                mole_frac.setFlags(Qt.ItemIsSelectable)
            else:
                mole_frac = QTableWidgetItem(f"{self.stratum.mole_fracs[q]:.2f}")
            # moleFrac.setTextAlignment(Qt.AlignCenter)
            self.strata_table.setItem(q, 1, mole_frac)

            # Thickness
            if q == 0 or q == len(self.stratum.materials) - 1:
                # TODO: To check the best default width for substrate and air
                thickness = QTableWidgetItem("1.0" if q == 0 else "3.0")
                thickness.setFlags(Qt.ItemIsSelectable)
            else:
                thickness = QTableWidgetItem(f"{self.stratum.layer_widths[q]:.2f}")
                if mtrl in self.stratum.cstm_idx:
                    thickness.setFlags(Qt.ItemIsSelectable)
            self.strata_table.setItem(q, 2, thickness)

            # doping
            if mtrl in DOPABLE_MATERIALS:
                doping = QTableWidgetItem(f"{self.stratum.dopings[q]:.1f}")
            else:
                doping = QTableWidgetItem("N/A")
                doping.setFlags(Qt.ItemIsSelectable)
            self.strata_table.setItem(q, 3, doping)

            # refractive index
            if q == 0:
                ridx = self.stratum.index_0
            elif q == len(self.stratum.materials) - 1:
                ridx = self.stratum.index_s
            else:
                ridx = self.stratum.indices[q - 1]
            ridx = QTableWidgetItem(f"{ridx.real:.3f} + {ridx.imag:.3f}i")
            ridx.setFlags(Qt.ItemIsSelectable)
            self.strata_table.setItem(q, 4, ridx)

        self.strata_table.resizeColumnsToContents()
        self.strata_table.blockSignals(False)

    def strata_table_mtrl_changed(self, row, selection):
        """SLOT as partial(self.strataTable_mtrlChanged, q) connected to
        mtrlName.currentTextChanged(str)"""
        self.stratum.materials[row] = selection
        if selection == "Active Core":
            self.update_custom_length()
        self.strata_table.selectRow(row)
        self.dirty.emit()

    @pyqtSlot()
    def update_xs(self):
        """Update position vector for plotting and confinement factor integral,
        also as SLOT to self.dirty SGINAL"""
        self.beta = None
        self.optimize_button.setEnabled(False)
        # self.stratum.updateIndices()
        self.strata_table_refresh()
        self.xs = np.linspace(-1, sum(self.stratum.layer_widths[1:]), 5000)
        self.update_canvas()

    @pyqtSlot()
    def solve(self):
        """SLOT connected to self.solveButton.clicked

        Yield
        -----
        Ey : np.ndarray(complex)
            The field to plot, normalized to max(abs(Ey)) = 1

        confinement : float
            The confinement factor in unit 1 (percentage 100*confinement)

        alphaW : float
            The waveguide loss in unit cm-1
        """
        try:
            self.beta = self.stratum.bound_mode_tm()
        except (TimeoutError, ValueError):
            QMessageBox.warning(self, EJ_ERROR, "Failed to solve for modes")
            return
        self.e_y, _, _ = self.stratum.populate_mode(self.beta, self.xs)
        # TODO
        self.confinement = self.stratum.confinementy(self.beta, self.xs, self.e_y)
        self.alpha_w = 4 * pi / (self.stratum.wl / 1e4) * self.beta.imag  # cm^-1
        self.update_canvas()
        self.update_loss()
        self.update_info()
        self.optimize_button.setEnabled(True)
        return self.e_y

    def update_info(self):
        """Update information in the info box"""
        info = ""
        if self.beta is not None:
            info += "Effective refractive index:\n"
            info += f"  β = {self.beta.real:.3f} + ({self.beta.imag:.3g})i\n"
            info += "Waveguide loss:\n"
            info += f"  α<sub>w</sub> = {self.alpha_w:.3f} cm<sup>-1</sup>\n"
            info += "Confinement factor: "
            info += f"  Γ = {self.confinement * 100:.1f}%\n"
            info += "Threshold gain:\n"
            gth = (self.alpha_m + self.alpha_w) / self.confinement
            info += f" g<sub>th</sub> = {gth:.1f} cm<sup>-1</sup>\n"
            try:
                self.jth = gth / self.stratum.cstm_gain["Active Core"]
                info += "Threshold current:\n"
                info += f"  J<sub>th</sub> = {self.jth:.1f} kA/cm<sup>-2</sup>"
            except (AttributeError, ZeroDivisionError, KeyError):
                info += "\nUse the quantum tab to define Active region"
        self.info_box.setHtml(info.replace("\n", "<br>"))

    def update_canvas(self):
        """Update figure in optCanvas"""
        self.ridx_axes.clear()
        self.mode_axes.clear()
        nx = self.stratum.populate_indices(self.xs)
        self.ridx_axes.plot(self.xs, nx.real, "k", lw=1)
        self.ridx_axes.set_xlabel("Position (μm)")
        self.ridx_axes.set_ylabel("Mode Intensity (a.u.) or Refractive Index")
        if self.red_active:
            for ar in self.stratum.populate_mtrl(self.xs):
                self.ridx_axes.plot(self.xs[ar], nx.real[ar], "r", lw=2)
        if np.max(nx.real) > 5:
            self.ridx_axes.set_ylim(top=5)
        if self.plot_imag:
            self.ridx_axes.plot(self.xs, nx.imag, "orange", lw=1)
            if np.max(nx.imag) > 5:
                self.ridx_axes.set_ylim(top=5)
        # plot select strata
        if self.select is not None:
            if self.select == 0:
                idx = self.xs < 0
            elif self.select == len(self.stratum.layer_widths) - 1:
                idx = self.xs > sum(self.stratum.layer_widths[1:-1])
            else:
                lsum = sum(self.stratum.layer_widths[1 : self.select])
                idx = (self.xs >= lsum) & (
                    self.xs < lsum + self.stratum.layer_widths[self.select]
                )
            self.ridx_axes.plot(self.xs[idx], nx.real[idx], "b", lw=1.5)
        if self.beta is not None:
            self.mode_axes.plot(self.xs, np.abs(self.e_y) ** 2, color="C0")
        # self.optCanvas.figure.tight_layout()
        self.opt_canvas.draw()

    @pyqtSlot()
    def view_red_active(self):
        self.red_active = not self.red_active
        self.update_canvas()

    @pyqtSlot()
    def optimize_strata(self):
        # TODO use threadRun in QuantumTab
        # TODO: block the button before solve
        nmax = len(self.stratum.layer_widths) - 2
        to_optimize, lmax, button_res = OptimizeInfoDialog(
            self,
            [
                self.stratum.materials[q + 1] not in self.stratum.cstm_idx
                for q in range(nmax)
            ],
            nmax,
            sum(self.stratum.layer_widths[1:-1]),
        ).exec()
        if button_res:
            optimize_opt_strata(self.stratum, self.alpha_m, to_optimize, lmax)
            self.strata_table_refresh()
            self.solve()


class OptimizeInfoDialog(QDialog):
    """Dialog to get information for optimization"""

    def __init__(self, parent, optimizable, nmax, lmax):
        super().__init__(parent)
        self.setModal(True)
        self.setWindowTitle("ErwinJr2: Optimize Optical Stratum")
        main_layout = QGridLayout()
        main_layout.addWidget(
            QLabel("Select the layer to optimize: "), 0, 0, 1, nmax + 1
        )
        main_layout.addWidget(QLabel("Layer:"), 1, 0, 1, 1)
        self.check_boxes = [QCheckBox() for _ in range(nmax)]
        for n in range(nmax):
            main_layout.addWidget(QLabel(str(n + 2)), 1, n + 1)
            main_layout.addWidget(self.check_boxes[n], 2, n + 1)
            if optimizable[n]:
                self.check_boxes[n].setChecked(True)
            else:
                self.check_boxes[n].setEnabled(False)
        length_layout = QHBoxLayout()
        length_layout.addWidget(QLabel("Maximum Length: "))
        self.length_box = QDoubleSpinBox()
        self.lengthBox.setSuffix(" μm")
        self.lengthBox.setRange(0.0, 20)
        self.lengthBox.setSingleStep(0.1)
        self.lengthBox.setDecimals(0.01)
        self.lengthBox.setValue(lmax)
        length_layout.addWidget(self.lengthBox)
        main_layout.addLayout(length_layout, 3, 0, 1, nmax + 1)
        button_box = QDialogButtonBox(
            QDialogButtonBox.Ok | QDialogButtonBox.Cancel, Qt.Horizontal
        )
        button_box.accepted.connect(self.accept)
        button_box.rejected.connect(self.reject)
        main_layout.addWidget(button_box, 4, 0, 1, nmax + 1)
        self.setLayout(main_layout)

    def exec(self):
        res = super().exec()
        checked = [n for n, box in enumerate(self.check_boxes) if box.isChecked()]
        return checked, self.lengthBox.value(), res
