#!/usr/bin/env python3
# -*- coding:utf-8 -*-
import unittest

import numpy as np

from ErwinJr2.qclayers import QCLayers


class TestThreeBandMatrix(unittest.TestCase):
    """This unit test is to compare simulation with PhysRevB.50.8663"""

    def test_three_well(self):
        layers = [300, 46, 10, 20, 10, 19, 300]
        mtrls = [1, 0, 1, 0, 1, 0, 1]
        qclayers = QCLayers(
            x_res=0.01, layer_widths=layers, layer_matrls=mtrls, repeats=1, temp=10.0
        )
        qclayers.solver = "ODE"
        qclayers.populate_x()
        qclayers.solve_whole()
        psis_ode, eigen_es_ode = qclayers.psis, qclayers.eigen_es
        philh_ode, phiso_ode = qclayers.philh, qclayers.phiso
        dipole_ode = qclayers.dipole(0, 1)
        qclayers.solver = "matrix"
        qclayers.matrix_eigen_count = 4
        qclayers.solve_whole()
        psis_mtx, eigen_es_mtx = qclayers.psis, qclayers.eigen_es
        philh_mtx, phiso_mtx = qclayers.philh, qclayers.phiso
        dipole_mtx = qclayers.dipole(0, 1)
        np.testing.assert_almost_equal(eigen_es_ode, eigen_es_mtx, decimal=8)
        np.testing.assert_almost_equal(psis_ode, psis_mtx, decimal=6)
        np.testing.assert_almost_equal(philh_ode, philh_mtx, decimal=3)
        np.testing.assert_almost_equal(phiso_ode, phiso_mtx, decimal=3)
        self.assertAlmostEqual(dipole_ode, dipole_mtx, places=3)

    def test_three_well_biased(self):
        layers = [300, 46, 10, 20, 10, 19, 300]
        mtrls = [1, 0, 1, 0, 1, 0, 1]
        qclayers = QCLayers(
            x_res=0.01,
            layer_widths=layers,
            e_field=30,
            layer_matrls=mtrls,
            repeats=1,
            temp=10.0,
        )
        qclayers.solver = "ODE"
        qclayers.populate_x()
        qclayers.solve_whole()
        psis_ode, eigen_es_ode = qclayers.psis[:4], qclayers.eigen_es[:4]
        dipole_ode = qclayers.dipole(0, 1)
        qclayers.solver = "matrix"
        qclayers.matrix_eigen_count = 4
        qclayers.solve_whole()
        psis_mtx, eigen_es_mtx = qclayers.psis, qclayers.eigen_es
        dipole_mtx = qclayers.dipole(0, 1)
        np.testing.assert_almost_equal(eigen_es_ode, eigen_es_mtx, decimal=5)
        np.testing.assert_almost_equal(psis_ode, psis_mtx, decimal=5)
        self.assertAlmostEqual(dipole_ode, dipole_mtx, places=3)


if __name__ == "__main__":
    unittest.main()
