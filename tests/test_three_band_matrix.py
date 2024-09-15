#!/usr/bin/env python3
# -*- coding:utf-8 -*-
import unittest

import numpy as np

from ErwinJr2.qc_layers import QCLayers

"""This unit test is to compare simulation with PhysRevB.50.8663"""


class TestThreeBandMatrix(unittest.TestCase):
    def test_three_well(self):
        layers = [300, 46, 10, 20, 10, 19, 300]
        mtrls = [1, 0, 1, 0, 1, 0, 1]
        qcLayers = QCLayers(
            x_res=0.01, layer_widths=layers, layer_matrls=mtrls, repeats=1, temp=10.0
        )
        qcLayers.solver = "ODE"
        qcLayers.populate_x()
        qcLayers.solve_whole()
        psis_ode, eigenEs_ode = qcLayers.psis, qcLayers.eigen_es
        philh_ode, phiso_ode = qcLayers.philh, qcLayers.phiso
        dipole_ode = qcLayers.dipole(0, 1)
        qcLayers.solver = "matrix"
        qcLayers.matrix_eigen_count = 4
        qcLayers.solve_whole()
        psis_mtx, eigenEs_mtx = qcLayers.psis, qcLayers.eigen_es
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
        qcLayers = QCLayers(
            x_res=0.01,
            layer_widths=layers,
            e_field=30,
            layer_matrls=mtrls,
            repeats=1,
            temp=10.0,
        )
        qcLayers.solver = "ODE"
        qcLayers.populate_x()
        qcLayers.solve_whole()
        psis_ode, eigenEs_ode = qcLayers.psis[:4], qcLayers.eigen_es[:4]
        dipole_ode = qcLayers.dipole(0, 1)
        qcLayers.solver = "matrix"
        qcLayers.matrix_eigen_count = 4
        qcLayers.solve_whole()
        psis_mtx, eigenEs_mtx = qcLayers.psis, qcLayers.eigen_es
        dipole_mtx = qcLayers.dipole(0, 1)
        np.testing.assert_almost_equal(eigenEs_ode, eigenEs_mtx, decimal=5)
        np.testing.assert_almost_equal(psis_ode, psis_mtx, decimal=5)
        self.assertAlmostEqual(dipole_ode, dipole_mtx, places=3)


if __name__ == "__main__":
    unittest.main()
