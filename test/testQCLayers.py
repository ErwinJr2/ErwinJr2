#!/usr/bin/env python3
# -*- coding:utf-8 -*-
from context import *  # type: ignore # noqa: F401, F403
import numpy as np
from ErwinJr2 import SaveLoad
import unittest


class TestQCLayers(unittest.TestCase):
    def test_solve_whole_matrix(self):
        with open("../ErwinJr2/example/PQLiu.json") as f:
            qcl = SaveLoad.qclLoad(f)
        qcl.solver = 'matrix'
        qcl.populate_x()
        qcl.solve_whole()
        # 0.264 is approximately 4.7um
        self.assertAlmostEqual(qcl.eigenEs[41] - qcl.eigenEs[31], 0.264, 2)
        qcl.period_recognize()
        qcl.period_map_build()
        self.assertAlmostEqual(qcl.eigenEs[41] - qcl.eigenEs[31], 0.26, 2)
        for i, j in ((25, 41), (24, 40), (26, 42), (31, 49)):
            self.assertEqual(qcl.periodMap[i][1], 1)
            self.assertEqual(qcl.periodIdx[qcl.periodMap[i][0]], j)
            np.testing.assert_almost_equal(
                qcl.psi_overlap(44, j, 1), qcl.psi_overlap(44, i), decimal=6)
            self.assertAlmostEqual(
                qcl._dipole(44, j, 1), qcl._dipole(44, i), 3)
            self.assertAlmostEqual(
                qcl._lo_transition(44, j, 1), qcl._lo_transition(44, i), 4)
            self.assertAlmostEqual(
                qcl._ifr_transition(44, j, 1), qcl._ifr_transition(44, i), 5)
        # Test overlapping
        np.testing.assert_almost_equal(
            qcl.psi_overlap(41, 49, 1), qcl.psi_overlap(41, 31), decimal=6)

    def test_solve_whole_ode(self):
        with open("../ErwinJr2/example/PQLiu.json") as f:
            qcl = SaveLoad.qclLoad(f)
        qcl.solver = 'ODE'
        qcl.populate_x()
        qcl.solve_whole()
        # 0.264 is approximately 4.7um
        self.assertAlmostEqual(qcl.eigenEs[31] - qcl.eigenEs[21], 0.264, 2)
        qcl.period_recognize()
        qcl.period_map_build()
        for i, j in ((21, 39), (14, 30), (15, 31)):
            self.assertEqual(qcl.periodMap[i][1], 1)
            self.assertEqual(qcl.periodIdx[qcl.periodMap[i][0]], j)
            np.testing.assert_almost_equal(
                qcl.psi_overlap(33, j, 1), qcl.psi_overlap(33, i), decimal=6)
            self.assertAlmostEqual(
                qcl._dipole(33, j, 1), qcl._dipole(33, i), 2)

    def test_solve_basis(self):
        with open("../ErwinJr2/example/PQLiu.json") as f:
            qcl = SaveLoad.qclLoad(f)
        qcl.populate_x()
        qcl.solve_basis()
        self.assertAlmostEqual(qcl.eigenEs[32] - qcl.eigenEs[31], 0.27, 2)
        self.assertEqual(qcl.eigenEs.shape, (96,))
        self.assertEqual(qcl.psis.shape, (96, 1384))

    def test_cache_consistency(self):
        with open("../ErwinJr2/example/PQLiu.json") as f:
            qcl = SaveLoad.qclLoad(f)
        qcl.includeIFR = True
        qcl.mtrlIFRLambda = [5.0] * 2
        qcl.mtrlIFRDelta = [5.0] * 2
        qcl.repeats = 4
        qcl.populate_x()
        qcl.solve_whole()
        self.assertFalse(hasattr(qcl, 'periodMap'))
        self.assertEqual(qcl.status, 'solved')
        taulo = qcl.lo_lifetime(33)
        tauifr = qcl.ifr_lifetime(33)
        tau = qcl.lifetime(33)
        gamma = qcl.ifr_broadening(31, 21)
        qcl.period_recognize()
        qcl.period_map_build()
        qcl.full_population()
        self.assertEqual(qcl.status, 'solved-full')
        self.assertTrue(hasattr(qcl, 'periodMap'))
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
        with open("../ErwinJr2/example/PQLiu.json") as f:
            qcl = SaveLoad.qclLoad(f)
        qcl.solver = 'matrix'
        qcl.repeats = 2
        qcl.xres = 0.02
        qcl.populate_x()
        e1 = qcl.solve_whole() - np.min(qcl.xVc)
        qcl.invert_layer()
        qcl.EField = -qcl.EField
        qcl.populate_x()
        e2 = qcl.solve_whole() - np.min(qcl.xVc)
        np.testing.assert_almost_equal(e1, e2, decimal=3)


if __name__ == "__main__":
    unittest.main()
