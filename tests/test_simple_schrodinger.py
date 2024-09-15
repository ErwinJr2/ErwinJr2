#!/usr/bin/env python3
# -*- coding:utf-8 -*-
import unittest

import numpy as np
from scipy.constants import e, hbar, m_e, pi
from scipy.special import ai_zeros, airy

from ErwinJr2.oned_quantum.c_schrodinger import cSimpleFillPsi, cSimpleSolve1D

ANG = 1e-10
# test precision
WF_DECIMAL = 3
E_DECIMAL = 6


class TestSimpleSchrodinger(unittest.TestCase):
    """Test simple Schrodinger solver"""

    def test_inf_sq_well(self):
        xmax = 200
        x = np.linspace(0, xmax, 200)
        mass = 0.067
        step = x[1] - x[0]
        es = np.linspace(0, 0.4, 300)
        v = np.zeros(x.shape)
        eigen_es = cSimpleSolve1D(step, es, v, mass)
        psis = cSimpleFillPsi(step, eigen_es, v, mass)

        # Theoretical result
        an = np.arange(1, 5)
        eigen_es_th = hbar**2 / (2 * m_e * mass) * an**2 * pi**2 / (xmax * ANG) ** 2 / e
        psis_th = np.sqrt(2 / xmax) * np.sin(np.outer(an, x) * pi / xmax)

        np.testing.assert_array_almost_equal(
            eigen_es[: len(an)], eigen_es_th, decimal=E_DECIMAL, verbose=True
        )
        np.testing.assert_array_almost_equal(
            psis[: len(an), :], psis_th, decimal=WF_DECIMAL, verbose=True
        )

    def test_trangle_well(self):
        f = 2e-4
        xmax = 1e3
        mass = 0.067  # Typical parameters for GaAs
        x = np.linspace(0, xmax, 200)
        step = x[1] - x[0]
        v = f * (xmax - x)
        es = np.linspace(0, 0.2, 300)
        eigen_es = cSimpleSolve1D(step, es, v, mass)
        psis = cSimpleFillPsi(step, eigen_es, v, mass)

        # Theoretical result
        nmax = 8
        an = ai_zeros(nmax)[0]
        eigen_es_th = (
            -((hbar**2 * f**2 / (2 * m_e * mass * e * ANG**2)) ** (1 / 3)) * an
        )
        psis_th = np.array(
            [
                airy(
                    (2 * mass * m_e * e * ANG**2 * f / hbar**2) ** (1 / 3)
                    * (x - eigen_e / f)
                )[0]
                for eigen_e in eigen_es_th
            ]
        )
        psis_th /= (np.linalg.norm(psis_th, axis=1) * np.sqrt(step))[:, None]

        np.testing.assert_array_almost_equal(
            eigen_es[:nmax], eigen_es_th, decimal=E_DECIMAL, verbose=True
        )
        np.testing.assert_array_almost_equal(
            psis[:nmax, ::-1], psis_th, decimal=WF_DECIMAL, verbose=True
        )


if __name__ == "__main__":
    unittest.main()
