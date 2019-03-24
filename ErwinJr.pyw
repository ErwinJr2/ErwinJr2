#!/usr/bin/env python
# -*- coding:utf-8 -*-

# TODO:
# Ctrl+z support
# add status bar

import os, sys
import traceback
from functools import partial

from QCLayers import QCLayers
import SaveLoad

from PyQt5.QtCore import (QSettings, QFile, QUrl,
                          QFileInfo, QVariant, Qt)
from PyQt5.QtGui import QIcon, QKeySequence, QPalette, QDesktopServices
from PyQt5.QtWidgets import (QApplication, QMainWindow, QTabWidget,
                             QAction, QMessageBox, QFileDialog,
                             QInputDialog, QSplashScreen)

from QuantumTab import QuantumTab

Version = '0.1'

class MainWindow(QMainWindow):
    def __init__(self, fname=None, parent=None):
        super(MainWindow, self).__init__(parent)
        self.qsettings = QSettings(QSettings.IniFormat, QSettings.UserScope, 
                                    "ErwinJr", "ErwinJr2", self)

        if self.qsettings.value('firstRun', True, type=bool):
            self.qsettings.setValue("firstRun", False)
            if not fname:
                firstRunBox = QMessageBox(
                    QMessageBox.Question, 'EwrinJr2 ' + Version,
                    ("Welcome to ErwinJr2!\n"
                     "Since this is your first time running the program, "
                     "would you like to open an example file or a blank file?"),
                    parent=self)
                firstRunBox.addButton("Blank File", QMessageBox.NoRole)
                firstRunBox.addButton("Example File", QMessageBox.YesRole)
                ansr = firstRunBox.exec_()
                if ansr:
                    fname = './example/PQLiu.json'
        elif not fname:
            # Load last file
            fname = self.qsettings.value("LastFile", None)
        self.dirty = True
        qclayers = None
        rf = self.qsettings.value("RecentFiles")
        self.recentFiles = rf if rf else []
        if fname and QFile.exists(fname):
            try:
                with open(fname, 'r') as f:
                    qclayers = SaveLoad.qclLoad(f)
                self.filename = fname
                self.addRecentFile(fname)
                self.dirty = False
            except Exception:
                QMessageBox.warning(self, "ErwinJr2 - Warning",
                                    "Could not load the *.json file.\n" +
                                    traceback.format_exc())
                qclayers = None
                self.filename = None
        else:
            self.filename = None

        self.mainTabWidget = QTabWidget()
        # ==========================
        # Quantum Tab
        # ==========================
        self.qtab = QuantumTab(qclayers, self)
        self.qtab.dirty.connect(self.thingsChanged)
        self.mainTabWidget.addTab(self.qtab, 'Quantum')
        # ==========================
        # Optical Tab
        # ==========================
        # TODO

        self.setCentralWidget(self.mainTabWidget)
        self.create_menu()

        if self.qsettings.value("MainWindow/Geometry"):
            self.restoreGeometry(self.qsettings.value("MainWindow/Geometry"))
        if self.qsettings.value("MainWindow/State"): 
            self.restoreState(self.qsettings.value("MainWindow/State"))
        self.updateFileMenu()
        self.updateWindowTitle()

    def thingsChanged(self):
        """ SLOT connected to self.qtab.dirty"""
        self.dirty = True
        self.updateWindowTitle()

# ===========================================================================
# General Menu Functions
# ===========================================================================
    def add_actions(self, target, actions):
        for action in actions:
            if action is None:
                target.addSeparator()
            else:
                target.addAction(action)

    def create_action(self, text, slot=None, shortcut=None, icon=None,
                      tip=None, checkable=False, ischecked=False):
        action = QAction(text, self)
        if icon:
            action.setIcon(QIcon("images/%s.png" % icon))
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
        return action

    def create_menu(self):
        # file menu
        self.file_menu = self.menuBar().addMenu("&File")
        newFileAction = self.create_action(
            "&New...", slot=self.fileNew, shortcut=QKeySequence.New,
            icon="filenew", tip="New ErwinJr2 file")
        openFileAction = self.create_action(
            "&Open", shortcut="Ctrl+O", slot=self.fileOpen,
            tip="Open ErwinJr2 file", icon="fileopen")
        saveFileAction = self.create_action(
            "&Save", shortcut="Ctrl+S", slot=self.fileSave,
            tip="Save ErwinJr2 file", icon="filesave")
        saveAsFileAction = self.create_action(
            "S&ave As", shortcut="Ctrl+W", slot=self.fileSaveAs,
            tip="Save ErwinJr2 file as", icon="filesaveas")
        exportQuantumCanvasAction = self.create_action(
            "Export Band Diagram Image", slot=self.exportBandDiagram,
            tip="Export Band Diagram Image")
        exportBandCSVAction = self.create_action(
            "Export Band Diagram Data", slot=self.export_band_diagram_data,
            tip="Export Band Diagram Data")
        quit_action = self.create_action(
            "&Quit", slot=self.close, shortcut="Ctrl+Q",
            tip="Close the application", icon="filequit")
        self.fileMenuActions = (
            newFileAction, openFileAction, saveFileAction, saveAsFileAction,
            None, exportBandCSVAction, exportQuantumCanvasAction, None,
            quit_action)
        self.file_menu.aboutToShow.connect(self.updateFileMenu)

        # edit menu
        self.edit_menu = self.menuBar().addMenu("&Edit")
        temperatureAction = self.create_action(
            "&Temperature", slot=self.set_temperature, tip="Set temperature")
        rotateLayerAction = self.create_action(
            "&Rotate Layer Table", slot=self.qtab.rotate_layer,
            tip="Move zeroth layer to first layer")
        solveARonly = self.create_action(
            "&Solve Active Only", checkable=True, 
            ischecked=self.qtab.qclayers.basisARonly, 
            slot=self.qtab.ARonly)
        copyStructureAction = self.create_action(
            "&Copy Structure", slot=self.qtab.copy_structure,
            tip="Copy Layer Structure to Clipboard")
        self.add_actions(self.edit_menu, (temperatureAction,
                                          rotateLayerAction, None,
                                          solveARonly, None,
                                          copyStructureAction))

        # view menu
        self.view_menu = self.menuBar().addMenu("&View")
        VXBandAction = self.create_action(
            "X Valley Conduction Band", checkable=True,
            ischecked=self.qtab.plotVX,
            slot=self.qtab.view_VXBand)
        VLBandAction = self.create_action(
            "L Valley Conduction Band",
            checkable=True, ischecked=self.qtab.plotVL,
            slot=self.qtab.view_VLBand)
        LHBandAction = self.create_action(
            "Light Hole Valence Band",
            checkable=True, ischecked=self.qtab.plotLH,
            slot=self.qtab.view_LHBand)
        SOBandAction = self.create_action(
            "Split Off Valence Band",
            checkable=True, ischecked=self.qtab.plotSO,
            slot=self.qtab.view_SOBand)
        plotwf = self.create_action(
            "Plot Wave function", 
            checkable=True, ischecked=self.qtab.plotType == 'wf', 
            slot=self.qtab.set_plotwf)
        plotFill = self.create_action(
            "Fill wave function curve", 
            checkable=True, ischecked=self.qtab.fillplot, 
            slot=self.qtab.set_fill)
        self.add_actions(self.view_menu, (VXBandAction,
                                          VLBandAction,
                                          LHBandAction,
                                          SOBandAction, 
                                          None, 
                                          plotwf, plotFill))

        # help menu
        self.help_menu = self.menuBar().addMenu("&Help")
        about_action = self.create_action("&About", 
                                          slot=self.on_about)
        licenses_action = self.create_action("&License",
                                             slot=self.on_licenses)
        tutorialAction = self.create_action("&Documents", shortcut='F1',
                                            slot=self.on_tutorial)
        self.add_actions(self.help_menu, (tutorialAction,
                                          about_action,
                                          licenses_action))

# ===========================================================================
# File Menu Items
# ===========================================================================
    def updateFileMenu(self):
        """SLOT connected to self.file_menu.aboutToShow()
        Update for recent files"""
        self.file_menu.clear()
        self.add_actions(self.file_menu, self.fileMenuActions[:-1])
        recentFiles = []
        for fname in self.recentFiles:
            if fname != self.filename and QFile.exists(fname):
                recentFiles.append(fname)
        if recentFiles:
            self.file_menu.addSeparator()
            for i, fname in enumerate(recentFiles):
                action = QAction(
                    "&{1}  {1}".format(i + 1,
                                       QFileInfo(fname).fileName()), self)
                action.setData(QVariant(fname))
                action.triggered.connect(partial(self.fileOpen, fname))
                self.file_menu.addAction(action)
        self.file_menu.addSeparator()
        self.file_menu.addAction(self.fileMenuActions[-1])


# ===========================================================================
# Save, Load, RecentFiles 
# ===========================================================================
    def updateWindowTitle(self):
        title = "ErwinJr2 " + Version
        if self.filename:
            title += " - %s"%os.path.basename(self.filename)
        if self.dirty:
            title += "[*]"
        self.setWindowTitle(title)
        self.setWindowModified(self.dirty)

    def addRecentFile(self, fname):
        if fname in self.recentFiles: 
            self.recentFiles.pop(self.recentFiles.index(fname))
        self.recentFiles.insert(0, fname)
        while len(self.recentFiles) > 9:
            self.recentFiles.pop()

    def unsaveConfirm(self):
        """Confirm if unsaved data should be saved. This returns False if the
        user canels the operation"""
        if self.dirty:
            reply = QMessageBox.question(
                self, "ErwinJr2 " + Version + " - Unsaved Changes",
                "Save unsaved changes?",
                QMessageBox.Yes | QMessageBox.No | QMessageBox.Cancel)
            if reply == QMessageBox.Cancel:
                return False
            elif reply == QMessageBox.Yes:
                self.fileSave()
        return True

    def fileNew(self):
        """Start a new file, confirm if there's unsaved data"""
        if not self.unsaveConfirm():
            return False
        self.qtab.qclayers = QCLayers()
        self.qtab.reload()
        self.filename = None
        self.dirty = True
        self.updateWindowTitle()
        return True

    def fileOpen(self, fname=None):
        """ Clear all old data and load a new file. This will check
        self.unsaveConfirm and return False if user cancels it. 
        This is used as user action and SLOT to fileOpen related signals."""
        if not self.unsaveConfirm():
            return False
        if not fname:
            fname, flt = QFileDialog.getOpenFileName(
                self, "ErwinJr2 - Choose file",
                os.path.dirname(self.filename) if self.filename else ".",
                "ErwinJr2 files (*.json)\nAll files (*.*)")
        if fname: 
            self.qclLoad(fname)
            return True
        else:
            # User cancels the OpenFile dialog
            return False

    def qclLoad(self, fname):
        """Load from file "fname", and update everything for consistency."""
        try:
            with open(fname, 'r') as f:
                self.qtab.qclayers = SaveLoad.qclLoad(f)
        except Exception:
            QMessageBox.warning(self, "ErwinJr2 - Warning",
                                "Could not load the *.json file.\n" +
                                traceback.format_exc())
            return
        self.qtab.reload()
        self.filename = fname
        self.addRecentFile(fname)
        self.dirty = False
        self.updateWindowTitle()

    def fileSave(self):
        if self.filename is None:
            return self.fileSaveAs()
        if self.filename.split('.')[-1] == 'json':
            try:
                with open(self.filename, 'w') as f:
                    SaveLoad.qclSaveJSON(f, self.qtab.qclayers)
            except:
                QMessageBox.warning(self, "ErwinJr2 - Warning",
                                    "Could not save *.json file.\n" +
                                    traceback.format_exc())
                return False
            self.dirty = False
            self.updateWindowTitle()
            return True
        else:
            raise IOError('The *.' + self.filename.split('.')[-1] +
                          ' extension is not supported.')

    def fileSaveAs(self):
        fname = self.filename if self.filename is not None else "."
        typeString = "ErwinJr2 file (*.json)\nAll files (*.*)"
        fname, flt = QFileDialog.getSaveFileName(
            self, "ErwinJr2 - Save File", fname, typeString)
        if fname:
            if "." not in fname:
                fname += ".json"
            self.addRecentFile(fname)
            self.filename = fname
            return self.fileSave()
        return False

    def closeEvent(self, event):
        if self.unsaveConfirm():
            filename = (QVariant(self.filename) if self.filename
                        else QVariant())
            self.qsettings.setValue("LastFile", filename)
            recentFiles = (QVariant(self.recentFiles) if self.recentFiles
                           else QVariant())
            self.qsettings.setValue("RecentFiles", recentFiles)
            self.qsettings.setValue(
                "MainWindow/Geometry", QVariant(self.saveGeometry()))
            self.qsettings.setValue(
                "MainWindow/State", QVariant(self.saveState()))
        else:
            event.ignore()


# ===========================================================================
# Export Functions
# ===========================================================================
    def exportBandDiagram(self):
        self.qtab.export_quantumCanvas(self.filename.split('.')[0])

    def export_band_diagram_data(self):
        fname, flt = QFileDialog.getSaveFileName(
            self, "ErwinJr2 - Export Band Structure Data",
            self.filename.split('.')[0],
            "Comma-Separated Value file (*.csv)")
        if fname:
            # if user doesn't click cancel
            self.qtab.export_band_data(fname)


# ===========================================================================
# Edit Menu Items
# ===========================================================================
    def set_temperature(self):
        nowTemp = self.qtab.qclayers.Temperature
        newTemp, buttonResponse = QInputDialog.getDouble(
            self, 'ErwinJr2 Input Dialog', 'Set Temperature',
            value=nowTemp, min=0)
        if buttonResponse:
            self.qtab.set_temperature(nowTemp)

# ===========================================================================
# Help Menu Items
# ===========================================================================
    def on_about(self):
        msg = """ ErwinJr2 0.x Authors and Contributors

         * Ming Lyu
            minglyu@princeton.edu
         * Yaofeng (Desmond) Zhong
         * Xiaowen Chen

        ErwinJr 2.x Authors and Contributors

         * Kale J. Franz, PhD (Jet Propulsion Laboratory)
            kfranz@alumni.princeton.edu
            www.kalefranz.com

With contributions from:
         * Yamac Dikmelik (Johns Hopkins University)
         * Yu Song (Princeton University)
        """
        QMessageBox.about(self, "ErwinJr2 " + Version, msg.strip())

    def on_licenses(self):
        copyright = """
ErwinJr2 is a simulation program for quantum semiconductor lasers.
Copyright (C) 2017 Ming Lyu
Copyright (C) 2012 Kale J. Franz, PhD

A portion of this code is Copyright (c) 2011, California Institute of
Technology ("Caltech"). U.S. Government sponsorship acknowledged.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""
        QMessageBox.about(self, "ErwinJr2 " + Version,
                          copyright.strip())

    def on_tutorial(self):
        QDesktopServices.openUrl(QUrl("./docs/_build/html/index.html"))


def main(filename=None):
    app = QApplication(sys.argv)
    app.setOrganizationName("ErwinJr")
    app.setOrganizationDomain("princetonuniversity.github.io/ErwinJr2")
    app.setApplicationName("ErwinJr2") 
    app.setWindowIcon(QIcon('images/EJpng48x48.png'))

    form = MainWindow(fileName)
    form.show()
    app.exec_()

if __name__ == "__main__":
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    # this block handles a filename passed in by command line
    try:
        fileName = os.path.abspath(sys.argv[1])
    except IndexError:
        fileName = None
    main(fileName)

# vim: ts=4 sw=4 sts=4 expandtab
