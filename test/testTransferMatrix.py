#!/usr/bin/env python3
# -*- coding:utf-8 -*-
import context
from OptStrata import OptStrata
import numpy as np
import unittest


class TestTransferMatrix(unittest.TestCase):
    def test_real_bound_mode(self):
        """TM0 mode in [John Chilwell and Ian Hodgkinson,
        J. Opt. Soc. Am. A 1, 742-753 (1984)"""
        hs = np.array([500]*4)
        indices = np.array([1.66, 1.53, 1.60, 1.66])
        n0 = 1.0
        ns = 1.5
        wl = 632.8
        NLayer = len(hs) + 2
        stratum = OptStrata(wl, ["test"]*NLayer, [0.0]*NLayer,
                            [0.0]*NLayer, hs)
        stratum.index0 = n0
        stratum.indexs = ns
        stratum.indices = indices
        beta = stratum.boundModeTM(max(indices))
        self.assertTrue(abs(beta - 1.620031) < 1e-6)

        xs = np.linspace(-500, 2500, 10000)
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


if __name__ == "__main__":
    unittest.main()
# vim: ts=4 sw=4 sts=4 expandtab
