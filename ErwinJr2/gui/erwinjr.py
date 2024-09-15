#!/usr/bin/env python3
# -*- coding:utf-8 -*-

import os
import sys
import traceback
from functools import partial

from PyQt5.QtCore import (
    QDir,
    QFile,
    QFileInfo,
    QSettings,
    QStandardPaths,
    QUrl,
    QVariant,
)
from PyQt5.QtGui import QDesktopServices, QIcon, QKeySequence, QPixmap
from PyQt5.QtWidgets import (
    QAction,
    QApplication,
    QFileDialog,
    QInputDialog,
    QMainWindow,
    QMessageBox,
    QSplashScreen,
    QTabWidget,
)

from ErwinJr2 import save_load
from ErwinJr2.gui.constants import ABOUT_MSG, COPYRIGHT, EJ_FULLNAME, EJ_NAME
from ErwinJr2.gui.optical_tab import OpticalTab
from ErwinJr2.gui.quantum_tab import QuantumTab
from ErwinJr2.opt_strata import OptStrata
from ErwinJr2.qc_layers import QCLayers, onedq
from ErwinJr2.version import VERSION

PROJ_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_FILE_DIR = os.path.join(
    QStandardPaths.writableLocation(QStandardPaths.DocumentsLocation), "QCLDesign"
)


class MainWindow(QMainWindow):
    """Main window for ErwinJr2"""

    def __init__(self, fname=None, parent=None):
        super().__init__(parent)
        self.qsettings = QSettings(
            QSettings.IniFormat, QSettings.UserScope, "ErwinJr", EJ_NAME, self
        )

        if self.qsettings.value("firstRun", True, type=bool):
            self.qsettings.setValue("firstRun", False)
            if not fname:
                first_run_box = QMessageBox(
                    QMessageBox.Question,
                    EJ_FULLNAME,
                    (
                        "Welcome to ErwinJr2!\n"
                        "Since this is your first time running the program, "
                        "would you like to open an example or a blank file?"
                    ),
                    parent=self,
                )
                first_run_box.addButton("Blank File", QMessageBox.NoRole)
                first_run_box.addButton("Example File", QMessageBox.YesRole)
                ansr = first_run_box.exec_()
                if ansr:
                    if not QDir(DEFAULT_FILE_DIR).exists():
                        os.makedirs(DEFAULT_FILE_DIR)
                    fname = os.path.join(DEFAULT_FILE_DIR, "PQLiu.json")
                    if not QFile(fname).exists():
                        example = os.path.join(PROJ_ROOT, "example", "PQLiu.json")
                        with open(example, "r") as fin:
                            with open(fname, "w") as fout:
                                fout.write(fin.read())
        elif not fname:
            # Load last file
            fname = self.qsettings.value("LastFile", None)
        self.dirty = True
        qclayers = None
        stratum = None
        rf = self.qsettings.value("RecentFiles")
        self.recent_files = rf if rf else []
        if fname and QFile.exists(fname):
            try:
                with open(fname, "r") as f:
                    qclayers, stratum = save_load.load_both(f)
                self.filename = fname
                self.add_recent_file(fname)
                self.dirty = False
            except Exception:  # pylint: disable=broad-except
                QMessageBox.warning(
                    self,
                    "ErwinJr2 - Warning",
                    "Could not load the *.json file.\n" + traceback.format_exc(),
                )
                qclayers = None
                stratum = None
                self.filename = None
        else:
            self.filename = None

        self.main_tab_widget = QTabWidget()
        self._set_tabs(qclayers, stratum)
        # self.mainTabWidget.setCurrentIndex(1)

        self.setCentralWidget(self.main_tab_widget)

        if self.qsettings.value("MainWindow/Geometry"):
            self.restoreGeometry(self.qsettings.value("MainWindow/Geometry"))
        if self.qsettings.value("MainWindow/State"):
            self.restoreState(self.qsettings.value("MainWindow/State"))
        self.update_file_menu()
        self.update_window_title()

    def _set_tabs(self, qclayers, stratum):
        self.main_tab_widget.clear()
        # ==========================
        # Quantum Tab
        # ==========================
        self.qtab = QuantumTab(qclayers, self)
        self.qtab.dirty.connect(self.things_changed)
        self.main_tab_widget.addTab(self.qtab, "Quantum")
        # ==========================
        # Optical Tab
        # ==========================
        self.otab = OpticalTab(stratum)
        self.otab.dirty.connect(self.things_changed)
        self.main_tab_widget.addTab(self.otab, "Optics")
        self.qtab.to_optics_button.clicked.connect(self.q2o)
        self.otab.fieldBox.setValue(self.qtab.qclayers.e_field)
        self.create_menu()
        self.main_tab_widget.currentChanged.connect(self.create_menu)

    def q2o(self):
        self.otab.setupActive(
            self.qtab.wavelength,
            self.qtab.qclayers.e_field,
            self.qtab.gain_coeff,
            self.qtab.neff,
            sum(self.qtab.qclayers.layer_widths),
        )
        self.main_tab_widget.setCurrentIndex(1)

    def things_changed(self):
        """SLOT connected to self.qtab.dirty"""
        self.dirty = True
        self.update_window_title()

    # ===========================================================================
    # General Menu Functions
    # ===========================================================================
    def add_actions(self, target, actions):
        for action in actions:
            if action is None:
                target.addSeparator()
            else:
                target.addAction(action)

    def create_action(
        self,
        text,
        slot=None,
        shortcut=None,
        icon=None,
        tip=None,
        checkable=False,
        ischecked=False,
        setting_sync=False,
    ):
        action = QAction(text, self)
        if icon:
            action.setIcon(QIcon(os.path.join(PROJ_ROOT, "images", f"{icon}.png")))
        if shortcut:
            action.setShortcut(shortcut)
        if tip:
            action.setToolTip(tip)
            action.setStatusTip(tip)
        if slot:
            action.triggered.connect(slot)
        if checkable:
            action.setCheckable(True)
        if ischecked:
            action.setChecked(True)
        if setting_sync:
            assert checkable
            action.triggered.connect(lambda: self.action_setting(text))
            # initialize
            if self.qsettings.value(text, ischecked, type=bool) != ischecked:
                action.setChecked(not ischecked)
                if slot:
                    slot()
        return action

    def action_setting(self, name):
        self.qsettings.setValue(name, not self.qsettings.value(name, False, type=bool))

    def create_menu(self, tab_index=0):
        self.menuBar().clear()
        # file menu
        self.file_menu = self.menuBar().addMenu("&File")
        new_file_action = self.create_action(
            "&New...",
            slot=self.file_new,
            shortcut=QKeySequence.New,
            icon="filenew",
            tip="New ErwinJr2 file",
        )
        open_file_action = self.create_action(
            "&Open",
            shortcut="Ctrl+O",
            slot=self.file_open,
            tip="Open ErwinJr2 file",
            icon="fileopen",
        )
        save_file_action = self.create_action(
            "&Save",
            shortcut="Ctrl+S",
            slot=self.file_save,
            tip="Save ErwinJr2 file",
            icon="filesave",
        )
        save_as_file_action = self.create_action(
            "S&ave As",
            shortcut="Ctrl+Alt+S",
            slot=self.file_save_as,
            tip="Save ErwinJr2 file as",
            icon="filesaveas",
        )
        export_quantum_canvas_action = self.create_action(
            "Export Band Diagram Image",
            slot=self.export_band_diagram,
            tip="Export Band Diagram Image",
        )
        export_band_csv_action = self.create_action(
            "Export Band Diagram Data",
            slot=self.export_band_diagram_data,
            tip="Export Band Diagram Data",
        )
        quit_action = self.create_action(
            "&Quit",
            slot=self.close,
            shortcut="Ctrl+W",
            tip="Close the application",
            icon="filequit",
        )
        self.file_menu_actions = (
            new_file_action,
            open_file_action,
            save_file_action,
            save_as_file_action,
            None,
            export_band_csv_action,
            export_quantum_canvas_action,
            None,
            quit_action,
        )
        self.add_actions(self.file_menu, self.file_menu_actions)
        self.file_menu.aboutToShow.connect(self.update_file_menu)

        # edit menu
        if tab_index == 0:
            self.edit_menu = self.menuBar().addMenu("&Edit")
            temperature_action = self.create_action(
                "&Temperature", slot=self.set_temperature, tip="Set temperature"
            )
            rotate_layer_action = self.create_action(
                "&Rotate Layer Table",
                slot=self.qtab.rotate_layer,
                tip="Move zeroth layer to first layer",
            )
            rotate_layer_action.setShortcut("Ctrl+T")
            invert_layer_action = self.create_action(
                "&Invert Layer Table",
                slot=self.qtab.invert_layer,
                tip="Invert the layer order",
            )
            solve_ar_only = self.create_action(
                "&Solve Active Only",
                checkable=True,
                ischecked=self.qtab.qclayers.basis_ar_only,
                slot=self.qtab.ARonly,
            )
            copy_structure_action = self.create_action(
                "&Copy Structure",
                slot=self.qtab.copy_structure,
                tip="Copy Layer Structure to Clipboard",
            )
            self.add_actions(
                self.edit_menu,
                (
                    temperature_action,
                    rotate_layer_action,
                    invert_layer_action,
                    None,
                    solve_ar_only,
                    None,
                    copy_structure_action,
                ),
            )

        # view menu
        self.view_menu = self.menuBar().addMenu("&View")
        if tab_index == 0:
            vx_band_action = self.create_action(
                "X Valley Conduction Band",
                checkable=True,
                ischecked=self.qtab.plot_vx,
                slot=self.qtab.view_vx_band,
            )
            vl_band_action = self.create_action(
                "L Valley Conduction Band",
                checkable=True,
                ischecked=self.qtab.plot_vl,
                slot=self.qtab.view_vl_band,
            )
            lh_band_action = self.create_action(
                "Light Hole Valence Band",
                checkable=True,
                ischecked=self.qtab.plot_lh,
                slot=self.qtab.view_lh_band,
            )
            so_band_action = self.create_action(
                "Split Off Valence Band",
                checkable=True,
                ischecked=self.qtab.plot_so,
                slot=self.qtab.view_so_band,
            )
            plotwf = self.create_action(
                "Plot Wave function",
                checkable=True,
                ischecked=self.qtab.plot_type == "wf",
                slot=self.qtab.set_plotwf,
                setting_sync=True,
            )
            plot_fill = self.create_action(
                "Fill wave function curve",
                checkable=True,
                ischecked=bool(self.qtab.fill_plot),
                slot=self.qtab.set_fill,
                setting_sync=True,
            )
            self.add_actions(
                self.view_menu,
                (
                    vx_band_action,
                    vl_band_action,
                    lh_band_action,
                    so_band_action,
                    None,
                    plotwf,
                    plot_fill,
                ),
            )
        else:
            optical_active_action = self.create_action(
                "Optical: Show Active Core",
                checkable=True,
                ischecked=self.otab.redActive,
                slot=self.otab.view_redActive,
                setting_sync=True,
            )
            self.add_actions(self.view_menu, (optical_active_action,))

        if tab_index == 0:
            # model menu
            self.model_menu = self.menuBar().addMenu("&Model")
            self.solver_actions = {}
            for solver in ("ODE", "matrix"):
                self.solver_actions[solver] = self.create_action(
                    solver,
                    checkable=True,
                    ischecked=self.qtab.qclayers.solver == solver,
                    slot=partial(self.choose_solver, solver),
                )
            if onedq is None:
                # C library does not exist
                self.solver_actions["ODE"].setEnabled(False)
            solver_menu = self.model_menu.addMenu("Eigen Solver")
            solver_menu.addActions(self.solver_actions.values())
            ifr_action = self.create_action(
                "IFR scattering",
                checkable=True,
                ischecked=self.qtab.qclayers.include_ifr,
                slot=self.qtab.trigger_ifr,
            )
            self.add_actions(self.model_menu, (None, ifr_action))

        # help menu
        # MacOS automatically remove the about item...
        about_text = (
            "&Who write this?" if sys.platform.startswith("darwin") else "&About"
        )
        self.help_menu = self.menuBar().addMenu("&Help")
        about_action = self.create_action(about_text, slot=self.on_about)
        licenses_action = self.create_action("&License", slot=self.on_licenses)
        tutorial_action = self.create_action(
            "&Documents", shortcut="F1", slot=self.on_tutorial
        )
        self.add_actions(
            self.help_menu, (tutorial_action, about_action, licenses_action)
        )

    # ===========================================================================
    # File Menu Items
    # ===========================================================================
    def update_file_menu(self):
        """SLOT connected to self.file_menu.aboutToShow()
        Update for recent files"""
        self.file_menu.clear()
        self.add_actions(self.file_menu, self.file_menu_actions[:-1])
        recent_files = []
        for fname in self.recent_files:
            if fname != self.filename and QFile.exists(fname):
                recent_files.append(fname)
        if recent_files:
            self.file_menu.addSeparator()
            for i, fname in enumerate(recent_files):
                action = QAction(f"&{i + 1}  {QFileInfo(fname).fileName()}", self)
                action.setData(QVariant(fname))
                action.triggered.connect(partial(self.file_open, fname))
                self.file_menu.addAction(action)
        self.file_menu.addSeparator()
        self.file_menu.addAction(self.file_menu_actions[-1])

    # ===========================================================================
    # Save, Load, RecentFiles
    # ===========================================================================
    def update_window_title(self):
        title = "ErwinJr2 " + VERSION
        if self.filename:
            title += " - " + os.path.basename(self.filename)
        if self.dirty:
            title += "[*]"
        self.setWindowTitle(title)
        self.setWindowModified(self.dirty)

    def add_recent_file(self, fname):
        if fname in self.recent_files:
            self.recent_files.pop(self.recent_files.index(fname))
        self.recent_files.insert(0, fname)
        while len(self.recent_files) > 9:
            self.recent_files.pop()

    def unsave_confirm(self):
        """Confirm if unsaved data should be saved. This returns False if the
        user cancels the operation"""
        if self.dirty:
            reply = QMessageBox.question(
                self,
                EJ_FULLNAME + " - Unsaved Changes",
                "Save unsaved changes?",
                QMessageBox.Yes | QMessageBox.No | QMessageBox.Cancel,
            )
            if reply == QMessageBox.Cancel:
                return False
            elif reply == QMessageBox.Yes:
                self.file_save()
        return True

    def file_new(self):
        """Start a new file, confirm if there's unsaved data"""
        if not self.unsave_confirm():
            return False
        self._set_tabs(QCLayers(), OptStrata())
        self.filename = None
        self.dirty = True
        self.update_window_title()
        return True

    def file_open(self, fname=None):
        """Clear all old data and load a new file. This will check
        self.unsaveConfirm and return False if user cancels it.
        This is used as user action and SLOT to fileOpen related signals."""
        if not self.unsave_confirm():
            return False
        if not fname:
            fname, _ = QFileDialog.getOpenFileName(
                self,
                "ErwinJr2 - Choose file",
                os.path.dirname(self.filename) if self.filename else ".",
                "ErwinJr2 files (*.json)\nAll files (*.*)",
            )
        if fname:
            self.qcl_load(fname)
            return True
        else:
            # User cancels the OpenFile dialog
            return False

    def qcl_load(self, fname):
        """Load from file "fname", and update everything for consistency."""
        try:
            with open(fname, "r") as f:
                qclayers, stratum = save_load.load_both(f)
                if stratum is None:
                    stratum = OptStrata(3.0)
        except Exception:  # pylint: disable=broad-except
            QMessageBox.warning(
                self,
                "ErwinJr2 - Warning",
                "Could not load the *.json file.\n" + traceback.format_exc(),
            )
            return
        self._set_tabs(qclayers, stratum)
        self.filename = fname
        self.add_recent_file(fname)
        self.dirty = False
        self.update_window_title()

    def file_save(self):
        if self.filename is None:
            return self.file_save_as()
        if self.filename.split(".")[-1] == "json":
            try:
                with open(self.filename, "w") as f:
                    save_load.save_json(f, self.qtab.qclayers, self.otab.stratum)
            except OSError:
                QMessageBox.warning(
                    self,
                    "ErwinJr2 - Warning",
                    "Could not save *.json file.\n" + traceback.format_exc(),
                )
                return False
            self.dirty = False
            self.update_window_title()
            return True
        else:
            raise IOError(
                "The *." + self.filename.split(".")[-1] + " extension is not supported."
            )

    def file_save_as(self):
        fname = self.filename if self.filename is not None else "."
        type_string = "ErwinJr2 file (*.json)\nAll files (*.*)"
        fname, _ = QFileDialog.getSaveFileName(
            self, "ErwinJr2 - Save File", fname, type_string
        )
        if fname:
            if "." not in fname:
                fname += ".json"
            self.add_recent_file(fname)
            self.filename = fname
            return self.file_save()
        return False

    def close_event(self, event):
        if self.unsave_confirm():
            filename = QVariant(self.filename) if self.filename else QVariant()
            self.qsettings.setValue("LastFile", filename)
            recent_files = (
                QVariant(self.recent_files) if self.recent_files else QVariant()
            )
            self.qsettings.setValue("RecentFiles", recent_files)
            self.qsettings.setValue(
                "MainWindow/Geometry", QVariant(self.saveGeometry())
            )
            self.qsettings.setValue("MainWindow/State", QVariant(self.saveState()))
        else:
            event.ignore()

    # ===========================================================================
    # Export Functions
    # ===========================================================================
    def export_band_diagram(self):
        if self.filename is None:
            filename = "new_qcl"
        else:
            filename = self.filename.split(".")[0]
        self.qtab.export_quantum_canvas(filename)

    def export_band_diagram_data(self):
        if self.filename is None:
            filename = "new_qcl"
        else:
            filename = self.filename.split(".")[0]
        fname, _ = QFileDialog.getSaveFileName(
            self,
            f"{EJ_NAME} - Export Band Structure Data",
            filename,
            "Comma-Separated Value file (*.csv)",
        )
        if fname:
            # if user doesn't click cancel
            self.qtab.export_band_data(fname)

    # ===========================================================================
    # Edit Menu Items
    # ===========================================================================
    def set_temperature(self):
        temp_now = self.qtab.qclayers.temperature
        temp_new, button_response = QInputDialog.getDouble(
            self, "ErwinJr2 Input Temperature", "Set Temperature", value=temp_now, min=0
        )
        if button_response:
            self.qtab.set_temperature(temp_new)

    # ===========================================================================
    # Model Menu Items
    # ===========================================================================
    def choose_solver(self, solver):
        if solver not in ("matrix", "ODE"):
            QMessageBox.warning(self, f"Unknown solver: {solver}")
            return
        if solver != self.qtab.qclayers.solver:
            for action in self.solver_actions.values():
                action.setChecked(False)
            self.solver_actions[solver].setChecked(True)
            self.qtab.trigger_solver(solver)
            self.dirty = True
            self.things_changed()

    # ===========================================================================
    # Help Menu Items
    # ===========================================================================
    def on_about(self):
        QMessageBox.about(self, EJ_FULLNAME, ABOUT_MSG.strip())

    def on_licenses(self):
        QMessageBox.about(self, EJ_FULLNAME, COPYRIGHT.strip())

    def on_tutorial(self):
        path = os.path.join(
            os.path.dirname(PROJ_ROOT), "docs", "_build", "html", "index.html"
        )
        if os.path.exists(path):
            QDesktopServices.openUrl(QUrl("file://" + os.path.abspath(path)))
        else:
            # TODO change URL according to version?
            QDesktopServices.openUrl(QUrl("https://erwinjr2.readthedocs.io/en/stable/"))


def main(filename=None):
    app = QApplication(sys.argv)
    app.setOrganizationName("ErwinJr")
    app.setOrganizationDomain("github.com/ErwinJr2")
    app.setApplicationName("ErwinJr2")
    app.setWindowIcon(QIcon(os.path.join(PROJ_ROOT, "images", "EJpng256.png")))

    # Create and display the splash screen
    splash_pix = QPixmap(os.path.join(PROJ_ROOT, "images", "splash.png"))
    splash = QSplashScreen(splash_pix)
    splash.setMask(splash_pix.mask())
    splash.show()

    form = MainWindow(filename)
    form.show()
    splash.finish(form)
    if sys.platform.startswith("win"):
        # The default font for win is not consistent with the OS
        app.setFont(QApplication.font("QMenu"))
    return app.exec_()
