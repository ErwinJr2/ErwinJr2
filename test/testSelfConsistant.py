#!/usr/bin/env python
# -*- coding:utf-8 -*-
from context import *
from OneDQuantum import *
from pylab import *
from scipy.constants import hbar, e, m_e, pi
import numpy as np
import unittest
import numpy.testing as npt
ANG=1E-8 # angstrom in cm  

#test precision
wf_dec = 3
E_dec = 6
class TestSelfConstistant(unittest.Testcase):
    def test_triangle_well(self):
        F = ??
        Vmax = 0.4
        Emax = 0.15
        xmax = Vmax / F
        x = np.linspace(0, xmax, 5000)
        step = x[1]-x[0]
        V = F * (xmax - x)
        mass = 0.067
        EigenEs = []
        Es = np.linspace(0, Emax, 100)
        EigenEs = cSimpleSolve1D(step, Es, V, mass)
        print(EigenEs)
        psis = cSimpleFillPsi(step, EigenEs, V, mass)
        eDensity, EF = cFermiDirac0N(sheet, EigenEs,
                                     0.067, psis, x[1]-x[0])
        Vc = -cCoulombField(x[1]-x[0], -eDensity, 13)
        yield x, np.ascontiguousarray(V[::-1]) , EigenEs,\
            np.ascontiguousarray(psis[:,::-1]),\
            np.ascontiguousarray(Vc[::-1]),\
            np.ascontiguousarray(eDensity[::-1]), EF
        xmax = (Vmax+Vc[-1]) / F
        x = np.linspace(0, xmax, 5000)
        V = F * (xmax - x) + Vc
        Vc -= Vc[-1]
        step = x[1]-x[0]
        Es = np.linspace(np.min(V), np.min(V)+Emax, 100)
        
        # Want result to be self consistent with what?
        npt.assert_array_almost_equal(???, ???, decimal=E_dec, verbose=True)
        npt.assert_array_almost_equal(???, ???, decimal=wf_dec, verbose=True)


if __name__ == "__main__":
    import sys
    # if the argument has MP in it, test for MP
    if len(sys.argv) > 1 and sys.argv.pop() == "MP":
        OneDSchrodinger.bindOpenMP(True)
    else:
        OneDSchrodinger.bindOpenMP(False)
    unittest.main()
# vim: ts=4 sw=4 sts=4 expandtab                                               

