#!/usr/bin/env python3
# -*- coding:utf-8 -*-
from context import *  # type: ignore # noqa: F401, F403
import SaveLoad
import unittest


class TestQCLayers(unittest.TestCase):
    def test_solve_whole(self):
        with open("../example/PQLiu.json") as f:
            qcl = SaveLoad.qclLoad(f)
        qcl.solver = 'matrix'
        qcl.populate_x()
        qcl.solve_whole()
        qcl.period_recognize()
        qcl.period_map()
        self.assertAlmostEqual(qcl.eigenEs[41] - qcl.eigenEs[31], 0.26, 2)
        for i, j in ((25, 41), (24, 40), (26, 42), (31, 49)):
            self.assertEqual(qcl.periodMap[i][1], 1)
            self.assertEqual(qcl.singlePeriodIdx[qcl.periodMap[i][0]], j)
        # This is just to make sure the program runs.
        # The consistency is checked in test/testThreeBandMatrix.py
        qcl.solver = 'ODE'
        qcl.populate_x()
        qcl.solve_whole()

    def test_solve_basis(self):
        with open("../example/PQLiu.json") as f:
            qcl = SaveLoad.qclLoad(f)
        qcl.populate_x()
        qcl.solve_basis()
        self.assertAlmostEqual(qcl.eigenEs[32] - qcl.eigenEs[31], 0.27, 2)
        self.assertEqual(qcl.eigenEs.shape, (103,))
        self.assertEqual(qcl.psis.shape, (103, 1384))


if __name__ == "__main__":
    unittest.main()
