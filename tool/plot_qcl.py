#!/usr/bin/env python3
# -*- coding:utf-8 -*-
import os

import matplotlib.pyplot as plt

from ErwinJr2 import save_load
from ErwinJr2.qc_layers import QCLayers
from ErwinJr2.qc_plotter import plotPotential, plotWF

_REPO_ROOT = os.path.dirname(os.path.dirname(__file__))
_TEST_SAMPLE_FILE = os.path.join(_REPO_ROOT, "ErwinJr2/example/PQLiu.json")

with open(_TEST_SAMPLE_FILE) as f:
    qcl = save_load.qclLoad(f)

qcl = QCLayers(
    substrate="InP",
    EField=51.0,
    wl=8.0,
    layerWidths=[
        34.0,
        14.0,
        33.0,
        13.0,
        32.0,
        15.0,
        31.0,
        19.0,
        29.0,
        23.0,
        27.0,
        25.0,
        27.0,
        44.0,
        18.0,
        9.0,
        57.0,
        11.0,
        54.0,
        12.0,
        45.0,
        25.0,
    ],
    layerMtrls=[0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1],
    layerDopings=[
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        2.0,
        2.0,
        2.0,
        2.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
    ],
    layerARs=[
        False,
        False,
        False,
        False,
        False,
        False,
        False,
        False,
        False,
        False,
        False,
        False,
        False,
        True,
        True,
        True,
        True,
        True,
        True,
        True,
        True,
        False,
    ],
    mtrlIFRDelta=[1.2, 1.2],
    mtrlIFRLambda=[90.0, 90.0],
)
qcl.includeIFR = True

qcl.populate_x()
qcl.solve_whole()
plotPotential(qcl)
plotWF(qcl)
plt.xlabel("Position (Å)")
plt.ylabel("Energy (eV)")
plt.show()

with open("testSave.json", "w") as f:
    save_load.EJSaveJSON(f, qcl)
