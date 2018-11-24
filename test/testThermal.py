#!/usr/bin/env python
# -*- coding:utf-8 -*-
from context import *
from OneDQuantum import *
from pylab import *
from scipy.constants import hbar, e, m_e, pi
import unittest
ANG=1E-10

#test precision

class TestTriangleWellThermal(unittest.TestCase):
    sheet = 1e17 * 50e-8 * (1e-8)**2
    F = 2.02e-4
    xmax = 1E3
    x = np.linspace(0, xmax, 5000)
    step = x[1]-x[0]
    V = F * (xmax - x)
    mass = 0.067
    Es = np.linspace(0, 0.15, 100)
    EigenEs = cSimpleSolve1D(step, Es, V, mass)
    psis = cSimpleFillPsi(step, EigenEs, V, mass)
    V = np.ascontiguousarray(V[::-1])
    psis = np.ascontiguousarray(psis[:,::-1])
    
    def testTriangleWellBoltzmann(self):
        T = 3000
        eDensity, EF = cBoltzmannN(T, self.sheet, self.EigenEs, 
                                   self.mass, self.psis, self.step)
        # Need to think about tests for high T Boltzmann distribution
        self.assertAlmostEqual(EF, np.mean(self.EigenEs))
        # Why does the average energy < 0? all eigen energy > 0
        self.assertAlmostEqual(self.step * np.sum(eDensity), self.sheet)


    # Zero-temperature Fermi-Dirac
    def testTriangleWellFermiDirac0T(self):
        eDensity, EF = cFermiDirac0N(self.sheet, self.EigenEs, 
                                     self.mass, self.psis, self.step)
        print(eDensity)
        self.assertAlmostEqual(EF, self.EigenEs[0], places=3)
        self.assertAlmostEqual(self.step * np.sum(eDensity), self.sheet)

    # Test for finite temperature Fermi-Dirac
    def testTriangleWellFermiDiracFT(self):
        pass

if __name__ == "__main__":
    unittest.main()

