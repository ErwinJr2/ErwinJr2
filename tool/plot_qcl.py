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
    qcl = save_load.qcl_load(f)

qcl = QCLayers(
    substrate="InP",
    e_field=51.0,
    wl=8.0,
    layer_widths=[
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
    layer_matrls=[0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1],
    layer_dopings=[
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
    layer_ar=[
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
    mtrl_ifr_delta=[1.2, 1.2],
    mtrl_ifr_lambda=[90.0, 90.0],
)
qcl.include_ifr = True

qcl.populate_x()
qcl.solve_whole()
plotPotential(qcl)
plotWF(qcl)
plt.xlabel("Position (Å)")
plt.ylabel("Energy (eV)")
plt.show()

with open("testSave.json", "w") as f:
    save_load.save_json(f, qcl)
