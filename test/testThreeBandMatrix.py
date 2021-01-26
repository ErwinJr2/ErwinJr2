#!/usr/bin/env python3
# -*- coding:utf-8 -*-
from context import *  # type: ignore # noqa: F401, F403
from QCLayers import SchrodingerLayer, QCLayers
import numpy as np
from numpy import sqrt
from scipy.constants import e as e0, hbar as hbar, electron_mass as m0
import unittest

"""This unit test is to compare simulation with PhysRevB.50.8663"""

class TestThreeBandMatrix(unittest.TestCase):
    def test_three_well(self):
        layers = [300, 46, 10, 20, 10, 19, 300]
        mtrls = [1, 0, 1, 0, 1, 0, 1]
        qcLayers = QCLayers(xres=0.01, layerWidths=layers,
                            layerMtrls=mtrls, repeats=1, T=10.0)
        qcLayers.populate_x()
        qcLayers.solve_whole()
        psis_ode, eigenEs_ode = qcLayers.psis, qcLayers.eigenEs
        dipole_ode = qcLayers.dipole(0, 1)
        qcLayers.solver = 'matrix'
        qcLayers.matrixEigenCount = 4
        qcLayers.solve_whole()
        psis_mtx, eigenEs_mtx = qcLayers.psis, qcLayers.eigenEs
        dipole_mtx = qcLayers.dipole(0, 1)
        np.testing.assert_almost_equal(eigenEs_ode, eigenEs_mtx, decimal=8)
        np.testing.assert_almost_equal(psis_ode, psis_mtx, decimal=6)
        self.assertAlmostEqual(dipole_ode, dipole_mtx, places=3)


def plot_debugger():
    import matplotlib.pyplot as plt
    layers = [300, 46, 10, 20, 10, 19, 300]
    mtrls = [1, 0, 1, 0, 1, 0, 1]
    qcLayers = QCLayers(xres=0.01, layerWidths=layers,
                        layerMtrls=mtrls, repeats=1, T=10.0)
    qcLayers.populate_x()
    qcLayers.solve_whole()
    psis_ode, eigenEs_ode = qcLayers.psis, qcLayers.eigenEs
    qcLayers.solver = 'matrix'
    qcLayers.matrixEigenCount = 4
    qcLayers.solve_whole()
    # validate_ODE(qcLayers)
    print(qcLayers.eigenEs)
    print(np.diff(qcLayers.eigenEs))
    print(qcLayers.dipole(0, 1))
    # [0.12559815 0.24226628 0.38040656 0.48960184]
    # [0.11666813 0.13814028 0.10919528]
    # 17.735693339174034

    # qcLayers = singleWell
    plt.xlabel('Position (Å)')
    plt.ylabel('Energy (eV)')
    xEg, xF, xEp, xESO = qcLayers.bandParams
    plt.plot(qcLayers.xPoints, qcLayers.xVc, 'k', linewidth=1)
    # plt.plot(qcLayers.xPoints, qcLayers.xVc-xEg, 'k', linewidth=1)
    # plt.plot(qcLayers.xPoints, qcLayers.xVc-xEg-xESO, 'k', linewidth=1)
    for n in range(qcLayers.eigenEs.size):
        plt.plot(qcLayers.xPoints,
                 10*qcLayers.psis[n, :] + qcLayers.eigenEs[n])
        plt.plot(qcLayers.xPoints,
                 10*psis_ode[n, :] + qcLayers.eigenEs[n], '--')
    plt.show()


if __name__ == "__main__":
    unittest.main()
    # plot_debugger()
