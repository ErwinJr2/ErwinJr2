#!/usr/bin/env python3
# -*- coding:utf-8 -*-
from context import *  # type: ignore # noqa: F401, F403
from ErwinJr2 import SchrodingerLayer
import numpy as np
from scipy.constants import e as e0, hbar as hbar, electron_mass as m0
import scipy.sparse as sparse
import scipy.sparse.linalg as splg
import unittest

"""This unit test is to compare simulation with PhysRevB.50.8663"""


class GaAs_Layer(SchrodingerLayer):
    """Sadly the parameters given in PhysRevB.50.8663 is not complete"""
    def __init__(self, xres, layerWidths):
        # fixed offset 0.51 eV
        super().__init__(xres, layerWidths=layerWidths,
                         layerVc=np.zeros(len(layerWidths)))
        self.crystalType = 'ZincBlende'

    def populate_material(self):
        N = self.xPoints.size
        ones = np.ones(N)
        self.xF = -1.94*ones
        self.xEp = 28.8*ones
        self.xESO = 0.341*ones
        self.xEg = 1.519*ones
        self.bandParams = (self.xEg, self.xF, self.xEp, self.xESO)
        self.gamma1 = 6.98*ones
        self.gamma2 = 2.06*ones
        self.gamma3 = 2.93*ones
        self.luttinger = (self.gamma1, self.gamma2, self.gamma3)


class TestEffectiveMass(unittest.TestCase):
    def test_single_well(self):
        L = 1E4
        layer = GaAs_Layer(5, [L])
        ks = np.pi/L * np.arange(1, 4)
        layer.solver = 'matrix'
        layer.populate_x()
        layer.populate_Kane_matrix()
        eigen_c = splg.eigsh(layer.Hsparse, 3, sigma=0,
                             return_eigenvectors=False)
        diag = sparse.diags(np.ones(3*layer.xPoints.size))
        eigen_lh = splg.eigsh(layer.Hsparse+1.519*diag, 3,
                              sigma=0, return_eigenvectors=False)
        eigen_so = splg.eigsh(layer.Hsparse+(1.519+0.341)*diag, 3,
                              sigma=0, return_eigenvectors=False)
        kunit = hbar**2/(2*m0*e0*1E-20)*ks**2
        mc = 1/(1 - 2*1.94 + 2*28.8/1.519/3 + 28.8/(1.519+0.341)/3)
        self.assertAlmostEqual(mc, 0.067, 4)
        np.testing.assert_allclose(eigen_c, kunit/mc, rtol=1.5E-2)

        mso = 1/(6.98 - 28.8*0.341/(3*1.519*(1.519+0.341)))
        mlh = 1/(6.98 + 2*2.06)
        self.assertAlmostEqual(mso, 0.172, 3)
        np.testing.assert_allclose(eigen_so[::-1], -kunit/mso, rtol=1E-2)
        np.testing.assert_allclose(eigen_lh[::-1], -kunit/mlh, rtol=1E-2)


if __name__ == "__main__":
    unittest.main()
