#!/usr/bin/env python3
# -*- coding:utf-8 -*-
from plotter import plot_band
from context import *
import SaveLoad
import numpy as np
import matplotlib.pyplot as plt




if __name__ == "__main__":
    with open("../example/PQLiu.json") as f:
        qcl = SaveLoad.qclLoad(f)

    qcl.layerSelected = 3
    qcl.NonParabolic = False
    qcl.populate_x()
    qcl.solve_whole()
    qcl.loTransition(5, 19)
    qcl.loTransition(4, 18)
    axes = plt.axes()
    plot_band(axes, qcl)
    plt.show()
