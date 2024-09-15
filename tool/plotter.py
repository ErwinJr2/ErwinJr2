#!/usr/bin/env python3
# -*- coding:utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np

from ErwinJr2 import qc_layers, save_load


def plot_band(axes, qclayer: qc_layers.QCLayers):
    """Plot potential (quantum barriers and wells) and other band parameters
    of the layer scturecture on axes, assuming already populated"""

    fsize = 12
    axes.set_xlabel("Position (Å)", fontsize=fsize)
    axes.set_ylabel("Energy (eV)", fontsize=fsize)
    axes.plot(qclayer.x_points, qclayer.x_vc, "k", linewidth=1)
    # for xv, conf in ((qcLayers.xVL, 'g--'),
    #                  (qcLayers.xVX, 'm-.'),
    #                  (qcLayers.xVLH, 'k'),
    #                  (qcLayers.xVSO, 'r--')):
    #     axes.plot(qcLayers.xPoints, xv, conf, linewidth=1)

    if hasattr(qclayer, "x_layer_selected"):
        axes.plot(qclayer.x_points, qclayer.x_layer_selected, "b", linewidth=1)

    if hasattr(qclayer, "eigen_es"):
        for n in range(qclayer.eigen_es.size):
            psi = qclayer.psis[n, :] ** 2
            tol = 1e-6
            start = np.argmax(psi > tol)
            end = len(psi) - np.argmax(psi[::-1] > tol)
            axes.plot(
                qclayer.x_points[start:end], 5 * psi[start:end] + qclayer.eigen_es[n]
            )


if __name__ == "__main__":
    with open("../example/PQLiue.json") as f:
        qcl = save_load.qcl_load(f)
    qcl.populate_x()
    qcl.solve_whole()

    axes = plt.axes()
    plot_band(axes, qcl)
    plt.show()
