#!/usr/bin/env python3
# -*- coding:utf-8 -*-
import unittest

import numpy as np
from scipy.constants import e, hbar, m_e, pi
from scipy.special import ai_zeros, airy

from ErwinJr2.OneDQuantum.c_schrodinger import cSimpleFillPsi, cSimpleSolve1D

ANG = 1E-10
# test precision
wf_dec = 3
E_dec = 6


class TestSimpleSchrodinger(unittest.TestCase):
    def test_inf_sq_well(self):
        xmax = 200
        x = np.linspace(0, xmax, 200)
        mass = 0.067
        step = x[1] - x[0]
        Es = np.linspace(0, 0.4, 300)
        V = np.zeros(x.shape)
        EigenEs = cSimpleSolve1D(step, Es, V, mass)
        psis = cSimpleFillPsi(step, EigenEs, V, mass)

        # Theoretical result
        an = np.arange(1, 5)
        EigenEs_th = hbar**2/(2*m_e*mass)*an**2*pi**2/(xmax*ANG)**2/e
        psis_th = np.sqrt(2/xmax) * np.sin(np.outer(an, x) * pi / xmax)

        np.testing.assert_array_almost_equal(
            EigenEs[:len(an)], EigenEs_th, decimal=E_dec, verbose=True)
        np.testing.assert_array_almost_equal(
            psis[:len(an), :], psis_th, decimal=wf_dec, verbose=True)

    def test_trangle_well(self):
        F = 2e-4
        xmax = 1E3
        mass = 0.067  # Typical parameters for GaAs
        x = np.linspace(0, xmax, 200)
        step = x[1] - x[0]
        V = F * (xmax - x)
        Es = np.linspace(0, 0.2, 300)
        EigenEs = cSimpleSolve1D(step, Es, V, mass)
        psis = cSimpleFillPsi(step, EigenEs, V, mass)

        # Theoretical result
        nmax = 8
        an = ai_zeros(nmax)[0]
        EigenEs_th = -(hbar**2*F**2/(2*m_e*mass*e*ANG**2))**(1/3)*an
        psis_th = np.array([airy((2*mass*m_e*e*ANG**2*F/hbar**2)**(1/3)
                                 * (x-E/F))[0] for E in EigenEs_th])
        psis_th /= (np.linalg.norm(psis_th, axis=1) * np.sqrt(step))[:, None]

        np.testing.assert_array_almost_equal(EigenEs[:nmax], EigenEs_th,
                                             decimal=E_dec, verbose=True)
        np.testing.assert_array_almost_equal(psis[:nmax, ::-1], psis_th,
                                             decimal=wf_dec, verbose=True)


if __name__ == "__main__":
    unittest.main()
