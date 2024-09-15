#!/usr/bin/env python3
# -*- coding:utf-8 -*-
import unittest

import numpy as np

from ErwinJr2.opt_strata import MaxwellLayer, MaxwellLayerAnisotropic


class TestTransferMatrix(unittest.TestCase):
    """Test 1D Maxwell solver"""

    def test_real_bound_mode(self):
        """TM0 mode in [John Chilwell and Ian Hodgkinson,
        J. Opt. Soc. Am. A 1, 742-753 (1984)"""
        hs = np.array([500] * 4)
        indices = np.array([1.0, 1.66, 1.53, 1.60, 1.66, 1.5])
        wl = 632.8
        stratum = MaxwellLayer(wl, hs, indices)
        beta = stratum.bound_mode_tm(max(indices))
        self.assertAlmostEqual(beta, 1.620031, 6)

        xs = np.linspace(-2000, 4000, 20000)
        n = stratum.populate_indices(xs)
        _, h_x, e_z = stratum.populate_mode(beta, xs)
        # dHx/dy = i k n^2 Ez
        lhs = (h_x[2:] - h_x[:-2]) / (xs[2] - xs[0])
        rhs = 1j * 2 * np.pi / wl * n**2 * e_z
        np.testing.assert_almost_equal(lhs, rhs[1:-1], decimal=3)
        # n^2 d/dy n^-2 d/dy Hx = -(n^2k^2 - beta^2) Hx
        lhs = (
            1j * 2 * np.pi / wl * n[1:-1] ** 2 * (e_z[2:] - e_z[:-2]) / (xs[2] - xs[0])
        )
        rhs = -(n**2 - beta**2) * (2 * np.pi / wl) ** 2 * h_x
        np.testing.assert_almost_equal(lhs, rhs[1:-1], decimal=3)
        self.assertAlmostEqual(e_z[0], 0, 5)
        self.assertAlmostEqual(e_z[-1], 0, 5)

    def test_anisotropy(self):
        hs = np.array([500] * 4)
        indices = np.array([1.0, 1.66, 1.53, 1.60, 1.66, 1.5])
        wl = 632.8
        stratum = MaxwellLayerAnisotropic(wl, hs, indices)
        stratum.indexy[2] = 1.8
        beta = stratum.bound_mode_tm(max(indices))

        xs = np.linspace(-500, 2500, 10000)
        nz, ny = stratum.populate_indices(xs)
        np.testing.assert_raises(AssertionError, np.testing.assert_almost_equal, ny, nz)
        _, h_x, e_z = stratum.populate_mode(beta, xs)
        # dHx/dy = i k nz^2 Ez
        lhs = (h_x[2:] - h_x[:-2]) / (xs[2] - xs[0])
        rhs = 1j * 2 * np.pi / wl * nz**2 * e_z
        np.testing.assert_almost_equal(lhs, rhs[1:-1], decimal=3)
        # ny^2 d/dy nz^-2 d/dy Hx = -(ny^2k^2 - beta^2) Hx
        lhs = (
            1j * 2 * np.pi / wl * ny[1:-1] ** 2 * (e_z[2:] - e_z[:-2]) / (xs[2] - xs[0])
        )
        rhs = -(ny**2 - beta**2) * (2 * np.pi / wl) ** 2 * h_x
        np.testing.assert_almost_equal(lhs, rhs[1:-1], decimal=3)


if __name__ == "__main__":
    unittest.main()
