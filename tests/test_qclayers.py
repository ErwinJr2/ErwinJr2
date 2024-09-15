#!/usr/bin/env python3
# -*- coding:utf-8 -*-
import os
import unittest

import numpy as np

from ErwinJr2 import save_load

_REPO_ROOT = os.path.dirname(os.path.dirname(__file__))
_TEST_SAMPLE_FILE = os.path.join(_REPO_ROOT, "ErwinJr2/example/PQLiu.json")


class TestQCLayers(unittest.TestCase):
    """Test QCLayers class"""

    # pylint: disable=protected-access
    def test_solve_whole_matrix(self):
        with open(_TEST_SAMPLE_FILE) as f:
            qcl = save_load.qcl_load(f)
        qcl.solver = "matrix"
        qcl.populate_x()
        qcl.solve_whole()
        # 0.264 is approximately 4.7um
        self.assertAlmostEqual(qcl.eigen_es[41] - qcl.eigen_es[31], 0.264, 2)
        qcl.period_recognize()
        qcl.period_map_build()
        self.assertAlmostEqual(qcl.eigen_es[41] - qcl.eigen_es[31], 0.26, 2)
        for i, j in ((25, 41), (24, 40), (26, 42), (31, 49)):
            self.assertEqual(qcl.period_map[i][1], 1)
            self.assertEqual(qcl.period_idx[qcl.period_map[i][0]], j)
            np.testing.assert_almost_equal(
                qcl.psi_overlap(44, j, 1), qcl.psi_overlap(44, i), decimal=6
            )
            self.assertAlmostEqual(qcl._dipole(44, j, 1), qcl._dipole(44, i), 3)
            self.assertAlmostEqual(
                qcl._lo_transition(44, j, 1), qcl._lo_transition(44, i), 4
            )
            self.assertAlmostEqual(
                qcl._ifr_transition(44, j, 1), qcl._ifr_transition(44, i), 5
            )
        # Test overlapping
        np.testing.assert_almost_equal(
            qcl.psi_overlap(41, 49, 1), qcl.psi_overlap(41, 31), decimal=6
        )

    def test_solve_whole_ode(self):
        with open(_TEST_SAMPLE_FILE) as f:
            qcl = save_load.qcl_load(f)
        qcl.solver = "ODE"
        qcl.populate_x()
        qcl.solve_whole()
        # 0.264 is approximately 4.7um
        self.assertAlmostEqual(qcl.eigen_es[31] - qcl.eigen_es[21], 0.264, 2)
        qcl.period_recognize()
        qcl.period_map_build()
        for i, j in ((21, 39), (14, 30), (15, 31)):
            self.assertEqual(qcl.period_map[i][1], 1)
            self.assertEqual(qcl.period_idx[qcl.period_map[i][0]], j)
            np.testing.assert_almost_equal(
                qcl.psi_overlap(33, j, 1), qcl.psi_overlap(33, i), decimal=6
            )
            self.assertAlmostEqual(qcl._dipole(33, j, 1), qcl._dipole(33, i), 2)

    def test_solve_basis(self):
        with open(_TEST_SAMPLE_FILE) as f:
            qcl = save_load.qcl_load(f)
        qcl.populate_x()
        qcl.solve_basis()
        self.assertAlmostEqual(qcl.eigen_es[32] - qcl.eigen_es[31], 0.27, 2)
        self.assertEqual(qcl.eigen_es.shape, (96,))
        self.assertEqual(qcl.psis.shape, (96, 1384))

    def test_cache_consistency(self):
        with open(_TEST_SAMPLE_FILE) as f:
            qcl = save_load.qcl_load(f)
        qcl.include_ifr = True
        qcl.mtrl_ifr_lambda = [5.0] * 2
        qcl.mtrl_ifr_delta = [5.0] * 2
        qcl.repeats = 4
        qcl.populate_x()
        qcl.solve_whole()
        self.assertFalse(hasattr(qcl, "period_map"))
        self.assertEqual(qcl.status, "solved")
        taulo = qcl.lo_lifetime(33)
        tauifr = qcl.ifr_lifetime(33)
        tau = qcl.lifetime(33)
        gamma = qcl.ifr_broadening(31, 21)
        qcl.period_recognize()
        qcl.period_map_build()
        qcl.full_population()
        self.assertEqual(qcl.status, "solved-full")
        self.assertTrue(hasattr(qcl, "period_map"))
        # A different cache is used
        self.assertNotEqual(taulo, qcl.lo_lifetime(33), 4)
        self.assertNotEqual(tauifr, qcl.ifr_lifetime(33), 4)
        self.assertNotEqual(tau, qcl.lifetime(33), 4)
        self.assertNotEqual(gamma, qcl.ifr_broadening(39, 31), 9)
        # But result should be approximately same
        self.assertAlmostEqual(taulo, qcl.lo_lifetime(33), 4)
        self.assertAlmostEqual(tauifr, qcl.ifr_lifetime(33), 4)
        self.assertAlmostEqual(tau, qcl.lifetime(33), 4)
        self.assertAlmostEqual(gamma, qcl.ifr_broadening(31, 21), 0)

    def test_revert_layer(self):
        with open(_TEST_SAMPLE_FILE) as f:
            qcl = save_load.qcl_load(f)
        qcl.solver = "matrix"
        qcl.repeats = 2
        qcl.x_step = 0.02
        qcl.populate_x()
        e1 = qcl.solve_whole() - np.min(qcl.x_vc)
        qcl.invert_layer()
        qcl.e_field = -qcl.e_field
        qcl.populate_x()
        e2 = qcl.solve_whole() - np.min(qcl.x_vc)
        np.testing.assert_almost_equal(e1, e2, decimal=3)


if __name__ == "__main__":
    unittest.main()
