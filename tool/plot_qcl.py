#!/usr/bin/env python3
# -*- coding:utf-8 -*-
from context import *  # type: ignore # noqa: F401, F403
from ErwinJr2 import SaveLoad
from ErwinJr2.QCPlotter import plotPotential, plotWF
from ErwinJr2.QCLayers import QCLayers
import matplotlib.pyplot as plt

with open("../ErwinJr2/example/PQLiu.json") as f:
    qcl = SaveLoad.qclLoad(f)

qcl = QCLayers(
    substrate='InP',
    EField=51.0,
    wl=8.0,
    layerWidths=[
        34.0, 14.0, 33.0, 13.0, 32.0, 15.0, 31.0, 19.0, 29.0, 23.0, 27.0, 25.0,
        27.0, 44.0, 18.0, 9.0, 57.0, 11.0, 54.0, 12.0, 45.0, 25.0],
    layerMtrls=[
        0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1],
    layerDopings=[
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 2.0, 2.0, 2.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    layerARs=[
        False, False, False, False, False, False, False, False, False, False,
        False, False, False, True, True, True, True, True, True, True, True,
        False],
    mtrlIFRDelta=[1.2, 1.2],
    mtrlIFRLambda=[90.0, 90.0]
)
qcl.includeIFR = True

qcl.populate_x()
qcl.solve_whole()
plotPotential(qcl)
plotWF(qcl)
plt.xlabel('Position (Å)')
plt.ylabel('Energy (eV)')
plt.show()

with open('testSave.json', 'w') as f:
    SaveLoad.EJSaveJSON(f, qcl)
