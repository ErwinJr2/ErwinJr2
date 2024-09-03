#!/usr/bin/env python3
# -*- coding:utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np

from ErwinJr2 import QCLayers, SaveLoad


def plot_band(axes, qcLayers: QCLayers):
    """ Plot potential (quantum barriers and wells) and other band parameters
    of the layer scturecture on axes, assuming already populated"""

    fsize = 12
    axes.set_xlabel('Position (Å)', fontsize=fsize)
    axes.set_ylabel('Energy (eV)', fontsize=fsize)
    axes.plot(qcLayers.xPoints, qcLayers.xVc, 'k', linewidth=1)
    # for xv, conf in ((qcLayers.xVL, 'g--'),
    #                  (qcLayers.xVX, 'm-.'),
    #                  (qcLayers.xVLH, 'k'),
    #                  (qcLayers.xVSO, 'r--')):
    #     axes.plot(qcLayers.xPoints, xv, conf, linewidth=1)

    if hasattr(qcLayers, 'xlayerSelected'):
        axes.plot(qcLayers.xPoints, qcLayers.xlayerSelected, 'b', linewidth=1)

    if hasattr(qcLayers, 'eigenEs'):
        for n in range(qcLayers.eigenEs.size):
            psi = qcLayers.psis[n, :]**2
            tol = 1e-6
            start = np.argmax(psi > tol)
            end = len(psi) - np.argmax(psi[::-1] > tol)
            axes.plot(qcLayers.xPoints[start:end],
                      5*psi[start:end] + qcLayers.eigenEs[n])


if __name__ == '__main__':
    with open("../example/new16umgaas-cladding.json") as f:
        qcl = SaveLoad.qclLoad(f)
    qcl.populate_x()
    qcl.solve_whole()

    axes = plt.axes()
    plot_band(axes, qcl)
    plt.show()
