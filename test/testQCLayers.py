#!/usr/bin/env python
# -*- coding:utf-8 -*-
from context import *
import SaveLoad
import numpy as np
import matplotlib.pyplot as plt

def plot_band(axes, qcLayers):
    """ Plot potential (quantum barriers and wells) and other band parameters 
    of the layer scturecture on axes, assuming already populated"""
    # for xv, conf in ((qcLayers.xVL, 'g--'),
    #                  (qcLayers.xVX, 'm-.'), 
    #                  (qcLayers.xVLH, 'k'), 
    #                  (qcLayers.xVSO, 'r--')):
    #     axes.plot(qcLayers.xPoints, xv, conf, linewidth=1)

    fsize = 12
    axes.set_xlabel('Position (Å)', fontsize=fsize)
    axes.set_ylabel('Energy (eV)', fontsize=fsize)
    axes.plot(qcLayers.xPoints, qcLayers.xVc, 'k', linewidth=1)

    if qcLayers.layerSelected is not None:
        for r in range(qcLayers.repeats):
            axes.plot(qcLayers.xPoints[qcLayers.indicesSelected[r, :]],
                      qcLayers.xVc[qcLayers.indicesSelected[r, :]], 'b',
                      linewidth=2)

    # if hasattr(qcLayers, 'eigenEs'): 
    #     for n in range(qcLayers.eigenEs.size): 
    #         axes.plot(qcLayers.xPoints, 
    #                   10*qcLayers.psis[n, :]**2 + qcLayers.eigenEs[n])

if __name__ == "__main__":
    with open("../example/PQLiu.json") as f:
         qcl = SaveLoad.qclLoad(f)

    qcl.layerSelected = 3
    qcl.populate_x()
    qcl.solve_whole()
    axes = plt.axes()
    plot_band(axes, qcl)
    plt.show()
    
