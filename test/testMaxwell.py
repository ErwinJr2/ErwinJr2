#!/usr/bin/env python
# -*- coding:utf-8 -*-
from context import *
from OneDQuantum import *
import numpy as np
import unittest

class TestMaxwell(unittest.TestCase):
    def test_real_bound_mode(self):
        """TM0 mode in [John Chilwell and Ian Hodgkinson, 
        J. Opt. Soc. Am. A 1, 742-753 (1984)"""
        hs = np.array([500]*4)
        indices = np.array([1.66, 1.53, 1.60, 1.66])
        n0 = 1.0
        ns = 1.5
        wl = 632.8
        beta = boundModeTM(max(indices), wl, hs, indices, n0, ns)
        self.assertTrue(abs(beta - 1.620031) < 1e-6)

        xs = np.linspace(-500, 2500, 10000)
        lsum = np.zeros(len(hs)+1)
        lsum[1:] = np.cumsum(hs)
        n = np.piecewise(xs, [(xs>=lsum[i]) & (xs<lsum[i+1]) 
            for i in range(len(hs))], indices)
        n[xs<0] = n0
        n[xs>=lsum[-1]] = ns

        Ey, Hx, Ez = genModeTM(beta, xs, wl, hs, indices, n0, ns)
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