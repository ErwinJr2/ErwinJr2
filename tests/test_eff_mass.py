#!/usr/bin/env python3
# -*- coding:utf-8 -*-
import unittest

import numpy as np
import scipy.sparse as sparse
import scipy.sparse.linalg as splg
from scipy.constants import e as e0, electron_mass as m0, hbar

from ErwinJr2.qc_layers import SchrodingerLayer


class GaAsLayer(SchrodingerLayer):
    """Sadly the parameters given in PhysRevB.50.8663 is not complete"""

    def __init__(self, xres, layer_width):
        # fixed offset 0.51 eV
        super().__init__(
            xres, layer_widths=layer_width, layer_vc=np.zeros(len(layer_width))
        )
        self.crystal_type = "ZincBlende"

    def populate_material(self):
        length = self.x_points.size
        ones = np.ones(length)
        self.x_f = -1.94 * ones
        self.x_ep = 28.8 * ones
        self.x_eso = 0.341 * ones
        self.x_eg = 1.519 * ones
        self.band_params = (self.x_eg, self.x_f, self.x_ep, self.x_eso)
        self.gamma1 = 6.98 * ones
        self.gamma2 = 2.06 * ones
        self.gamma3 = 2.93 * ones
        self.luttinger = (self.gamma1, self.gamma2, self.gamma3)


class TestEffectiveMass(unittest.TestCase):
    """This unit test is to compare simulation with PhysRevB.50.8663"""

    def test_single_well(self):
        layer_width = 1e4
        layer = GaAsLayer(5, [layer_width])
        ks = np.pi / layer_width * np.arange(1, 4)
        layer.solver = "matrix"
        layer.populate_x()
        layer.populate_kane_matrix()
        eigen_c = splg.eigsh(layer.h_sparse, 3, sigma=0, return_eigenvectors=False)
        diag = sparse.diags(np.ones(3 * layer.x_points.size))
        eigen_lh = splg.eigsh(
            layer.h_sparse + 1.519 * diag, 3, sigma=0, return_eigenvectors=False
        )
        eigen_so = splg.eigsh(
            layer.h_sparse + (1.519 + 0.341) * diag,
            3,
            sigma=0,
            return_eigenvectors=False,
        )
        kunit = hbar**2 / (2 * m0 * e0 * 1e-20) * ks**2
        mc = 1 / (1 - 2 * 1.94 + 2 * 28.8 / 1.519 / 3 + 28.8 / (1.519 + 0.341) / 3)
        self.assertAlmostEqual(mc, 0.067, 4)
        np.testing.assert_allclose(eigen_c, kunit / mc, rtol=1.5e-2)

        mso = 1 / (6.98 - 28.8 * 0.341 / (3 * 1.519 * (1.519 + 0.341)))
        mlh = 1 / (6.98 + 2 * 2.06)
        self.assertAlmostEqual(mso, 0.172, 3)
        np.testing.assert_allclose(eigen_so[::-1], -kunit / mso, rtol=1e-2)
        np.testing.assert_allclose(eigen_lh[::-1], -kunit / mlh, rtol=1e-2)


if __name__ == "__main__":
    unittest.main()
