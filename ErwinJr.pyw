#!/usr/bin/env python
# -*- coding:utf-8 -*-

# ===========================================================================
# ErwinJr2 is a simulation program for quantum semiconductor lasers.
# Copyright (C) 2012 Kale J. Franz, PhD
# Copyright (C) 2017 Ming Lyu (CareF)
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
# find replacement for psyco
# check unnecessary function call
# Ctrl+z support
# add status bar

import os
import sys
import traceback
from functools import partial
import time

from QCLayers import QCLayers
import SaveLoad

from PyQt5.QtCore import (QSettings, SIGNAL, QString, QFile,
                          QFileInfo, QVariant, Qt)
from PyQt5.QtGui import QIcon, QKeySequence, QPalette, QPixmap
from PyQt5.QtWidgets import (QApplication, QMainWindow, QTabWidget,
                             QAction, QMessageBox, QFileDialog,
                             QInputDialog, QSplashScreen)

from QuantumTab import QuantumTab

# ===========================================================================
# Version
# ===========================================================================
ejVersion = 171109
majorVersion = '3.4.0'

class MainWindow(QMainWindow):
    def __init__(self, fileName=None, parent=None):
        super(MainWindow, self).__init__(parent)

        self.filename = fileName

        self.mainTabWidget = QTabWidget()
        # ==========================
        # Quantum Tab
        # ==========================
        self.qtab = QuantumTab(self)
        self.qtab.dirty.connect(self.winUpdate)
        self.mainTabWidget.addTab(self.qtab, 'Quantum')
        # ==========================
        # Optical Tab
        # ==========================
        # TODO
        self.setCentralWidget(self.mainTabWidget)

        self.create_menu()

        qsettings = QSettings(parent=self)
        self.recentFiles = qsettings.value("RecentFiles").toStringList()
        self.restoreGeometry(
            qsettings.value("MainWindow/Geometry").toByteArray())
        self.restoreState(qsettings.value("MainWindow/State").toByteArray())
        self.updateFileMenu()

        self.dirty = False

        if self.filename:
            self.fileOpen(self.filename)
        else:
            self.loadInitialFile()

        self.dirty = False
        self.update_windowTitle()

    def winUpdate(self):
        """ SLOT connected to self.qtab.dirty"""
        self.dirty = True
        self.update_windowTitle()

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
        if icon is not None:
            action.setIcon(QIcon("images/%s.png" % icon))
        if shortcut is not None:
            action.setShortcut(shortcut)
        if tip is not None:
            action.setToolTip(tip)
            action.setStatusTip(tip)
        if slot is not None:
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
        bumpLayerAction = self.create_action(
            "&Bump First Layer", slot=self.qtab.bump_first_layer,
            tip="Move zeroth layer to first layer")
        copyStructureAction = self.create_action(
            "&Copy Structure", slot=self.qtab.copy_structure,
            tip="Copy Layer Structure to Clipboard")
        self.add_actions(self.edit_menu, (temperatureAction,
                                          bumpLayerAction, None,
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
        self.add_actions(self.view_menu, (VXBandAction,
                                          VLBandAction,
                                          LHBandAction,
                                          SOBandAction))

        # help menu
        self.help_menu = self.menuBar().addMenu("&Help")
        about_action = self.create_action("&About", shortcut='F1',
                                          slot=self.on_about)
        licenses_action = self.create_action("&License",
                                             slot=self.on_licenses)
        tutorialAction = self.create_action("&Tutorial",
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
        current = (QString(self.filename)
                   if self.filename is not None else None)
        recentFiles = []
        for fname in self.recentFiles:
            if fname != current and QFile.exists(fname):
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

    def addRecentFile(self, fname):
        if fname is None:
            return
        if fname not in self.recentFiles:
            self.recentFiles.insert(0, fname)
            while self.recentFiles.count() > 9:
                self.recentFiles.pop()

    def loadInitialFile(self):
        qsettings = QSettings()
        fname = qsettings.value("LastFile").toString()
        if fname and QFile.exists(fname):
            self.qclLoad(fname)
            self.qtab.reload()

        self.filename = fname
        self.addRecentFile(fname)
        self.dirty = False
        self.update_windowTitle()

    def fileNew(self):
        """Start a new file, confirm if there's unsaved data"""
        if not self.unsaveConfirm():
            return False

        self.filename = None
        self.qtab.qclayers = QCLayers()
        self.qtab.reload()

        self.dirty = False
        self.update_windowTitle()

        return True

    def unsaveConfirm(self):
        """Confirm if unsaved data should be saved"""
        if self.dirty:
            reply = QMessageBox.question(
                self, "ErwinJr2 " + str(majorVersion) + " - Unsaved Changes",
                "Save unsaved changes?",
                QMessageBox.Yes | QMessageBox.No | QMessageBox.Cancel)
            if reply == QMessageBox.Cancel:
                return False
            elif reply == QMessageBox.Yes:
                self.fileSave()
        return True

    def update_windowTitle(self):
        if self.filename is not None:
            self.setWindowTitle("ErwinJr2 " + str(majorVersion) + " - %s[*]" %
                                os.path.basename(str(self.filename)))
        else:
            self.setWindowTitle("ErwinJr2 " + str(majorVersion) + "[*]")
        self.setWindowModified(self.dirty)

    def fileOpen(self, fname=None):
        # clear all old data, also calls self.unsaveConfirm()
        if not self.fileNew():
            return False
        if not fname:
            dir = os.path.dirname(str(self.filename)) if self.filename else "."
            fname = QFileDialog.getOpenFileName(
                self, "ErwinJr2 - Choose file", dir,
                "ErwinJr2 files (*.qcl)\nAll files (*.*)")
        # open file and determine if it is from the Matlab version of ErwinJr
        filehandle = open(fname, 'r')
        firstLine = filehandle.readline()
        filehandle.close()
        if fname:
            if firstLine.split(':')[0] == 'Description':
                QMessageBox.warning(
                    self, 'ErwinJr2 Error',
                    'Older .qcl format is no longer supported for Ver>3.0.')
                #  self.qclPtonLoad(fname)
            elif firstLine == 'ErwinJr2 Data File\n':
                self.qclLoad(fname)
            else:
                QMessageBox.warning(self, 'ErwinJr2 Error',
                                    'Could not recognize input file.')
                return
            self.qtab.reload()

        self.filename = fname
        self.addRecentFile(fname)
        self.dirty = False
        self.update_windowTitle()

        return True

    def qclLoad(self, fname):
        try:
            with open(fname, 'r') as f:
                self.qtab.qclayers = SaveLoad.qclLoad(f)
        except Exception:
            QMessageBox.warning(self, "ErwinJr2 - Warning",
                                "Could not load *.qcl file.\n" +
                                traceback.format_exc())

    def fileSave(self):
        if self.filename is None:
            return self.fileSaveAs()
        else:
            # os.path.extsep
            if self.filename.split('.')[-1] == 'qcl':
                if self.qclSave(self.filename):
                    self.dirty = False
                    self.update_windowTitle()
                    return True
                else:
                    return False
            else:
                raise IOError('The *.' + self.filename.split('.')[-1] +
                              ' extension is not supported.')
                return False

    def fileSaveAs(self):
        fname = self.filename if self.filename is not None else "."
        typeString = "ErwinJr 2.x file (*.qcl)\nAll files (*.*)"
        fname = QFileDialog.getSaveFileName(
            self, "ErwinJr2 - Save File", QString(fname), typeString)
        if fname:
            if "." not in fname:
                fname += ".qcl"
            self.addRecentFile(fname)
            self.filename = fname
            return self.fileSave()
        return False

    def qclSave(self, fname):
        try:
            with open(fname, 'w') as f:
                SaveLoad.qclSaveJSON(f, self.qtab.qclayers)
        except Exception:
            QMessageBox.warning(self, "ErwinJr2 - Warning",
                                "Could not save *.qcl file.\n" +
                                traceback.format_exc())
        return True

    def closeEvent(self, event):
        if self.unsaveConfirm():
            qsettings = QSettings()
            filename = (QVariant(QString(self.filename)) if self.filename
                        else QVariant())
            qsettings.setValue("LastFile", filename)
            recentFiles = (QVariant(self.recentFiles) if self.recentFiles
                           else QVariant())
            qsettings.setValue("RecentFiles", recentFiles)
            qsettings.setValue(
                "MainWindow/Geometry", QVariant(self.saveGeometry()))
            qsettings.setValue(
                "MainWindow/State", QVariant(self.saveState()))
        else:
            event.ignore()


# ===========================================================================
# Export Functions
# ===========================================================================
    def exportBandDiagram(self):
        if __USE_MATPLOTLIB__:
            self.qtab.export_quantumCanvas(
                self.filename.split('.')[0])
        else:
            fname = unicode(QFileDialog.getSaveFileName(
                self, "ErwinJr2 - Export Band Structure Image",
                self.filename.split('.')[0],
                "Portable Network Graphics file (*.png)"))
            if not fname:
                return

            # set background color to white and save presets
            bgRole = self.mainTabWidget.backgroundRole()
            self.mainTabWidget.setBackgroundRole(QPalette.Base)
            self.qtab.export_quantumCanvas(fname)
            self.mainTabWidget.setBackgroundRole(bgRole)

    def export_band_diagram_data(self):
        fname = unicode(QFileDialog.getSaveFileName(
            self, "ErwinJr2 - Export Band Structure Data",
            self.filename.split('.')[0],
            "Comma-Separated Value file (*.csv)"))
        if fname != '':
            # if user doesn't click cancel
            self.qtab.export_band_data(fname)


# ===========================================================================
# Edit Menu Items
# ===========================================================================
    def set_temperature(self):
        # TODO
        nowTemp = 300
        newTemp, buttonResponse = QInputDialog.getDouble(
            self, 'ErwinJr2 Input Dialog', 'Set Temperature',
            value=nowTemp, min=0)
        if buttonResponse:
            pass

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
        QMessageBox.about(self, "ErwinJr2 " + str(ejVersion), msg.strip())

    def on_licenses(self):
        copyright1 = """
#=======================================
# ErwinJr2 is a simulation program for quantum semiconductor lasers.
# Copyright (C) 2017 Ming Lyu
# Copyright (C) 2012 Kale J. Franz, PhD
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
#=======================================
"""
        QMessageBox.about(self, "ErwinJr2 " + str(ejVersion),
                          copyright1.strip())

    def on_tutorial(self):
        if os.name == "nt":
            os.startfile("tutorial.pdf")
        elif os.name == "posix":
            os.system("/usr/bin/xdg-open tutorial.pdf")


def main():
    app = QApplication(sys.argv)
    app.setOrganizationName("ErwinJr")
    app.setOrganizationDomain("princetonuniversity.github.io/ErwinJr2")
    app.setApplicationName("ErwinJr2")
    qsettingsSystem = QSettings(QSettings.SystemScope, 
                                "Princeton", "ErwinJr2")
    installDirectory = qsettingsSystem.value('installDirectory').toString()
    if installDirectory:
        os.chdir(installDirectory)

    # Create and display the splash screen
    splash_pix = QPixmap('images/erwinjr_splash.png')
    splash = QSplashScreen(splash_pix, Qt.WindowStaysOnTopHint)
    splash.setMask(splash_pix.mask())
    splash.show()
    app.processEvents()

    time.sleep(1)

    app.setWindowIcon(QIcon('images/EJpng48x48.png'))

    # this block handles a filename passed in by command line
    try:
        fileName = sys.argv[1]
        name, ext = os.path.splitext(fileName)
        assert ext == ".qcl"
        assert os.path.exists(fileName)
        fileName = os.path.abspath(fileName)
    except (IndexError, AssertionError):
        fileName = None

    form = MainWindow(fileName)
    form.show()
    splash.finish(form)

    qsettings = QSettings()
    if not qsettings.value('firstRun').toInt()[1]:
        if not installDirectory:
            qsettingsSystem.setValue("installDirectory", QVariant(os.getcwd()))
        firstRunBox = QMessageBox(
            QMessageBox.Question, 'EwrinJr2 ' + str(majorVersion),
            ("Welcome to ErwinJr2!\n"
             "Since this is your first time running the program, "
             "would you like to open an example file or a blank file?"),
            parent=form)
        firstRunBox.addButton("Blank File", QMessageBox.NoRole)
        firstRunBox.addButton("Example File", QMessageBox.YesRole)
        ansr = firstRunBox.exec_()
        if ansr:
            form.fileOpen('examples/NPhoton PQLiu.qcl')
        else:
            form.fileNew()
        qsettings.setValue("firstRun", 1)

    app.exec_()


if __name__ == "__main__":
    main()

# vim: ts=4 sw=4 sts=4 expandtab
