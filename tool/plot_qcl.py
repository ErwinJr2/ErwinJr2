#!/usr/bin/env python3
# -*- coding:utf-8 -*-
from context import *  # type: ignore # noqa: F401, F403
from ErwinJr2 import SaveLoad
from ErwinJr2.QCPlotter import plotPotential, plotWF
import matplotlib.pyplot as plt

with open("../ErwinJr2/example/PQLiu.json") as f:
    qcl = SaveLoad.qclLoad(f)

qcl.populate_x()
qcl.solve_whole()
plotPotential(qcl)
plotWF(qcl)
plt.xlabel('Position (Å)')
plt.ylabel('Energy (eV)')
plt.show()
