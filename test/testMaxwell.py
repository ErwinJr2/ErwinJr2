#!/usr/bin/env python3
# -*- coding:utf-8 -*-
import context  # type: ignore # noqa: F401
from ErwinJr2.OptStrata import MaxwellLayer, MaxwellLayer_anisotropic
import numpy as np
import unittest


class TestTransferMatrix(unittest.TestCase):
    def test_real_bound_mode(self):
        """TM0 mode in [John Chilwell and Ian Hodgkinson,
        J. Opt. Soc. Am. A 1, 742-753 (1984)"""
        hs = np.array([500]*4)
        indices = np.array([1.0, 1.66, 1.53, 1.60, 1.66, 1.5])
        wl = 632.8
        stratum = MaxwellLayer(wl, hs, indices)
        beta = stratum.boundModeTM(max(indices))
        self.assertAlmostEqual(beta, 1.620031, 6)

        xs = np.linspace(-2000, 4000, 20000)
        n = stratum.populateIndices(xs)
        Ey, Hx, Ez = stratum.populateMode(beta, xs)
        # dHx/dy = i k n^2 Ez
        lhs = (Hx[2:] - Hx[:-2])/(xs[2]-xs[0])
        rhs = 1j*2*np.pi/wl*n**2*Ez
        np.testing.assert_almost_equal(lhs, rhs[1:-1], decimal=3)
        # n^2 d/dy n^-2 d/dy Hx = -(n^2k^2 - beta^2) Hx
        lhs = 1j*2*np.pi/wl*n[1:-1]**2*(Ez[2:] - Ez[:-2])/(xs[2]-xs[0])
        rhs = -(n**2 - beta**2)*(2*np.pi/wl)**2*Hx
        np.testing.assert_almost_equal(lhs, rhs[1:-1], decimal=3)
        self.assertAlmostEqual(Ez[0], 0, 5)
        self.assertAlmostEqual(Ez[-1], 0, 5)

    def test_anisotropy(self):
        hs = np.array([500]*4)
        indices = np.array([1.0, 1.66, 1.53, 1.60, 1.66, 1.5])
        wl = 632.8
        stratum = MaxwellLayer_anisotropic(wl, hs, indices)
        stratum.indexy[2] = 1.8
        beta = stratum.boundModeTM(max(indices))

        xs = np.linspace(-500, 2500, 10000)
        nz, ny = stratum.populateIndices(xs)
        np.testing.assert_raises(AssertionError,
                                 np.testing.assert_almost_equal, ny, nz)
        Ey, Hx, Ez = stratum.populateMode(beta, xs)
        # dHx/dy = i k nz^2 Ez
        lhs = (Hx[2:] - Hx[:-2])/(xs[2]-xs[0])
        rhs = 1j*2*np.pi/wl*nz**2*Ez
        np.testing.assert_almost_equal(lhs, rhs[1:-1], decimal=3)
        # ny^2 d/dy nz^-2 d/dy Hx = -(ny^2k^2 - beta^2) Hx
        lhs = 1j*2*np.pi/wl*ny[1:-1]**2*(Ez[2:] - Ez[:-2])/(xs[2]-xs[0])
        rhs = -(ny**2 - beta**2)*(2*np.pi/wl)**2*Hx
        np.testing.assert_almost_equal(lhs, rhs[1:-1], decimal=3)


if __name__ == "__main__":
    unittest.main()
