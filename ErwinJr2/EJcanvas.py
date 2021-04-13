"""
This file defines class for plotting figures in ErwinJr
"""

import numpy as np
import os

import matplotlib
from matplotlib.backends.backend_qt5agg import (FigureCanvasQTAgg as
                                                FigureCanvas)
from matplotlib.backends.backend_qt5 import cursord  # , NavigationToolbar2QT
from matplotlib.backend_bases import (NavigationToolbar2, cursors)
# import matplotlib.backends.qt_editor.figureoptions as figureoptions
from matplotlib.figure import Figure

from PyQt5.QtCore import QObject, Signal
from PyQt5.QtWidgets import QSizePolicy, QMessageBox, QFileDialog

import sys
import json
from warnings import warn
matplotlib.use('Qt5Agg')
config = {
    "PlotMargin": {'l': 0.9, 'r': 0.12, 'b': 0.6, 't': 0.09},
    "fontsize": 12,
    "wfscale": 0.3,
    "modescale": 3,
    # This number not necessarily but is chosen to be consistent with the
    # default value of tol for QCLayers.periodRecognize as
    # wf_almost_zero / modescale = tol
    "wf_almost_zero": 1.5e-4,
}
if sys.platform == 'win32':
    config["PlotMargin"] = {'l': 0.9, 'r': 0.12, 'b': 0.55, 't': 0.09}
elif sys.platform == 'darwin':
    config["PlotMargin"] = {'l': 0.90, 'r': 0.12, 'b': 0.55, 't': 0.09}


def LoadConfig(fname="plotconfig.json"):
    try:
        with open(fname, 'r') as f:
            userConfig = json.load(f)
    except FileNotFoundError:
        userConfig = {}
        warn("Cannot load plot config file %s. Use default config." % fname,
             UserWarning)
    for k in userConfig:
        config[k] = userConfig[k]


class EJcanvas(FigureCanvas):
    """EJcanvas is designed for ErwinJr as a canvas for energy band and
    wavefunctions"""
    def __init__(self, xlabel='x', ylabel='y', parent=None):
        self.figure = Figure()
        super(EJcanvas, self).__init__(self.figure)
        #  NavigationToolbar2.__init__(self, self)
        self.setParent(parent)
        self.setMinimumWidth(200)
        self.setMinimumHeight(200)
        self.setSizePolicy(QSizePolicy.MinimumExpanding,
                           QSizePolicy.MinimumExpanding)
        self.updateGeometry()
        self.axes = self.figure.add_subplot(111)
        self.xlabel = xlabel
        self.ylabel = ylabel

    def set_axes(self, fsize=config["fontsize"]):
        self.axes.tick_params(axis='both', which='major', labelsize=fsize)
        self.axes.spines['top'].set_color('none')
        self.axes.spines['right'].set_color('none')
        self.axes.spines['bottom'].set_position(('outward', 5))
        self.axes.set_xlabel(self.xlabel, fontsize=fsize)
        self.axes.spines['left'].set_position(('outward', 5))
        self.axes.set_ylabel(self.ylabel, fontsize=fsize)
        self.axes.autoscale(enable=True, axis='x', tight=True)
        self.axes.autoscale(enable=True, axis='y', tight=True)

    def test(self):
        """A test function, plotting sin(x)"""
        x = np.linspace(0, 10, 100)
        self.axes.plot(x, np.sin(x))

    def resizeEvent(self, event):
        super(EJcanvas, self).resizeEvent(event)
        height = self.figure.get_figheight()
        width = self.figure.get_figwidth()
        margin = config["PlotMargin"]
        self.figure.subplots_adjust(
            left=margin['l'] / width, bottom=margin['b'] / height,
            right=1 - margin['r'] / width, top=1 - margin['t'] / height,
            wspace=0, hspace=0)

    def clear(self):
        self.axes.clear()
        self.set_axes()


class EJplotControl(NavigationToolbar2, QObject):
    """This class is an implementation of NavigationToolbar2 for controlling
    the canvas of ErwinJr plotting. The class is mainly inspired by
    `NavigationToolbar2QT` in backend_qt5.py
    Critical APIs are:
    - set_action to set necessary actions that are part of `toolitems` of the
      controller.
    - set_custom to set customized actions, basically 'layerselect' and
      `pairselect` for quantum wells/barriers and quantum states selection.
    - trigger_custom to trigger custom actions
    """
    message = Signal(str)

    def __init__(self, canvas, parent):
        # Not QToolBar to not draw any toolbar icons.
        QObject.__init__(self, parent)
        NavigationToolbar2.__init__(self, canvas)
        self._actions = {}
        self._custom_active = {}
        self._custom_cursor = {}
        self.zoomed = False

    def _init_toolbar(self):
        # override for backward compatibility.
        pass

    def set_action(self, callback, button):
        """According to matplotlib.backend_bases, supported actions are:
        toolitems = (
            ('Home', 'Reset original view', 'home', 'home'),
            ('Back', 'Back to  previous view', 'back', 'back'),
            ('Forward', 'Forward to next view', 'forward', 'forward'),
            (None, None, None, None),
            ('Pan',
            'Left button pans, Right button zooms\n'
            'x/y fixes axis, CTRL fixes aspect',
            'move', 'pan'),
            ('Zoom', 'Zoom to rectangle \nx/y fixes axis, CTRL fixes aspect',
            'zoom_to_rect', 'zoom'),
            ('Subplots', 'Configure subplots', 'subplots',
                'configure_subplots'),
            (None, None, None, None),
            ('Save', 'Save the figure', 'filesave', 'save_figure'),
          )
        format (text, tooltip_text, image_file, callback)
        """
        if callback in (self.toolitems[n][-1] for n in
                        range(len(self.toolitems))):
            button.clicked.connect(getattr(self, callback))
            self._actions[callback] = button
            if callback in ['zoom', 'pan']:
                button.setCheckable(True)
        else:
            raise TypeError("%s not supported." % callback)

    def set_custom(self, name, button, callback, cursor=cursors.HAND):
        """customized callback action,
        s.t. this class manages the button's check status
        and its interaction with canvas.
        Note that callback is called when clicked on canvas"""
        self._custom_active[name] = callback
        self._actions[name] = button
        self._custom_cursor[name] = cursor
        self._custom_mode = None
        self._custom_callBack = None
        button.setCheckable(True)

    def _get_mode_name(self):
        # This is a work around for backward compatibility due to
        # https://github.com/matplotlib/matplotlib/pull/17135
        if self._custom_mode:
            return self._custom_mode
        try:
            name = self.mode.name
        except AttributeError:
            name = self._active
        return name

    def _update_buttons_checked(self):
        # sync button check-states to match active mode
        if 'pan' in self._actions:
            self._actions['pan'].setChecked(self._get_mode_name() == 'PAN')
        if 'zoom' in self._actions:
            self._actions['zoom'].setChecked(self._get_mode_name() == 'ZOOM')
        for mode in self._custom_active:
            self._actions[mode].setChecked(self._get_mode_name() == mode)

    def _reset_mode(self):
        # This is a work around for self.mode value type being private.
        name = self._get_mode_name()
        if name == 'PAN':
            self.pan()
            return
        if name == 'ZOOM':
            self.zoom()
            return

    def pan(self, *args):
        # override
        super(EJplotControl, self).pan(*args)
        self._reset_custom()
        self._update_buttons_checked()

    def zoom(self, *args):
        # override
        super(EJplotControl, self).zoom(*args)
        self._reset_custom()
        self.zoomed = True
        self._update_buttons_checked()

    def home(self, *args):
        # override
        super(EJplotControl, self).home(*args)
        self.zoomed = False

    def _reset_custom(self):
        self._custom_mode = None
        if self._custom_callBack is not None:
            self.canvas.mpl_disconnect(self._custom_callBack)
            self._custom_callBack = None

    def trigger_custom(self, mode):
        self._reset_mode()
        if self._custom_mode != mode:
            self._reset_custom()
            self._custom_mode = mode
            self._custom_callBack = self.canvas.mpl_connect(
                'button_release_event', self._custom_active[mode])
            self.canvas.widgetlock(self)
        else:
            self._reset_custom()
            self.canvas.widgetlock.release(self)

        for a in self.canvas.figure.get_axes():
            a.set_navigate_mode(None)

        self._update_buttons_checked()

    def _update_cursor(self, event):
        # override backend_bases.NavigationToolbar2._update_cursor
        if (event.inaxes and self._custom_mode and
                self._lastCursor != self._custom_cursor[self._custom_mode]):
            self.set_cursor(self._custom_cursor[self._custom_mode])
            return
        super(EJplotControl, self)._update_cursor(event)

    def set_cursor(self, cursor):
        # implement
        self.canvas.setCursor(cursord[cursor])

    def draw_rubberband(self, event, x0, y0, x1, y1):
        # implement
        height = self.canvas.figure.bbox.height
        y1 = height - y1
        y0 = height - y0
        rect = [int(val) for val in (x0, y0, x1 - x0, y1 - y0)]
        self.canvas.drawRectangle(rect)

    def remove_rubberband(self):
        # implement
        self.canvas.drawRectangle(None)

    def save_figure(self, caption="Choose a filename to save to",
                    filename=None, default_filetype=None):
        # This method is override to support customized filename.
        filetypes = self.canvas.get_supported_filetypes_grouped()
        sorted_filetypes = sorted(filetypes.items())
        if default_filetype is None:
            default_filetype = self.canvas.get_default_filetype()

        if not filename:
            startpath = os.path.expanduser(
                matplotlib.rcParams['savefig.directory'])
            filename = os.path.join(startpath,
                                    self.canvas.get_default_filename())
        filters = []
        selectedFilter = None
        for name, exts in sorted_filetypes:
            exts_list = " ".join(['*.%s' % ext for ext in exts])
            filter = '%s (%s)' % (name, exts_list)
            #  filter = '(*.%s) %s' % (exts_list, name)
            if default_filetype in exts:
                selectedFilter = filter
            filters.append(filter)
        filters = ';;'.join(filters)
        print(filters, selectedFilter)

        fname, filter = QFileDialog.getSaveFileName(
            self.parent(), caption, filename, filters, selectedFilter)
        if fname:
            try:
                self.canvas.figure.savefig(fname, dpi=300,
                                           bbox_inches='tight')
            except Exception as e:
                QMessageBox.critical(
                    self.parent(), "Error saving band diagram",
                    e, QMessageBox.Ok, QMessageBox.NoButton)
