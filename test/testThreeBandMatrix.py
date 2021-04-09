#!/usr/bin/env python3
# -*- coding:utf-8 -*-
from context import *  # type: ignore # noqa: F401, F403
from ErwinJr2.QCLayers import QCLayers
import numpy as np
import unittest

"""This unit test is to compare simulation with PhysRevB.50.8663"""


class TestThreeBandMatrix(unittest.TestCase):
    def test_three_well(self):
        layers = [300, 46, 10, 20, 10, 19, 300]
        mtrls = [1, 0, 1, 0, 1, 0, 1]
        qcLayers = QCLayers(xres=0.01, layerWidths=layers,
                            layerMtrls=mtrls, repeats=1, T=10.0)
        qcLayers.solver = 'ODE'
        qcLayers.populate_x()
        qcLayers.solve_whole()
        psis_ode, eigenEs_ode = qcLayers.psis, qcLayers.eigenEs
        philh_ode, phiso_ode = qcLayers.philh, qcLayers.phiso
        dipole_ode = qcLayers.dipole(0, 1)
        qcLayers.solver = 'matrix'
        qcLayers.matrixEigenCount = 4
        qcLayers.solve_whole()
        psis_mtx, eigenEs_mtx = qcLayers.psis, qcLayers.eigenEs
        philh_mtx, phiso_mtx = qcLayers.philh, qcLayers.phiso
        dipole_mtx = qcLayers.dipole(0, 1)
        np.testing.assert_almost_equal(eigenEs_ode, eigenEs_mtx, decimal=8)
        np.testing.assert_almost_equal(psis_ode, psis_mtx, decimal=6)
        np.testing.assert_almost_equal(philh_ode, philh_mtx, decimal=3)
        np.testing.assert_almost_equal(phiso_ode, phiso_mtx, decimal=3)
        self.assertAlmostEqual(dipole_ode, dipole_mtx, places=3)

    def test_three_well_biased(self):
        layers = [300, 46, 10, 20, 10, 19, 300]
        mtrls = [1, 0, 1, 0, 1, 0, 1]
        qcLayers = QCLayers(xres=0.01, layerWidths=layers, EField=30,
                            layerMtrls=mtrls, repeats=1, T=10.0)
        qcLayers.solver = 'ODE'
        qcLayers.populate_x()
        qcLayers.solve_whole()
        psis_ode, eigenEs_ode = qcLayers.psis[:4], qcLayers.eigenEs[:4]
        dipole_ode = qcLayers.dipole(0, 1)
        qcLayers.solver = 'matrix'
        qcLayers.matrixEigenCount = 4
        qcLayers.solve_whole()
        psis_mtx, eigenEs_mtx = qcLayers.psis, qcLayers.eigenEs
        dipole_mtx = qcLayers.dipole(0, 1)
        np.testing.assert_almost_equal(eigenEs_ode, eigenEs_mtx, decimal=5)
        np.testing.assert_almost_equal(psis_ode, psis_mtx, decimal=5)
        self.assertAlmostEqual(dipole_ode, dipole_mtx, places=3)


if __name__ == "__main__":
    unittest.main()
