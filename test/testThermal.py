#!/usr/bin/env python
# -*- coding:utf-8 -*-
from context import *
from OneDQuantum import *
from pylab import *
from scipy.constants import hbar, e, m_e, pi
import unittest
ANG=1E-10

#test precision
eDensity_dec = 5


class TestUniformPsis(unittest.TestCase):
    EigenEs = np.linspace(0, 1, 2)
    xmax = 1
    l = 3
    x = np.linspace(0, xmax, l)
    step = x[1] - x[0]
    psis = np.ones((2,l))/np.sqrt(l)/np.sqrt(step)
    mass = 1
    sheet = 0.04179

    # Zero-temperature Fermi-Dirac = low temperature Fermi-Dirac
    def testUniformFermiDirac0T(self):
        EF = 0.5
        eDensity0T = cFermiDirac0(EF, self.EigenEs, self.mass, self.psis,
                                self.step)
        eDensityLowT = cFermiDirac(0.01, EF, self.EigenEs, self.mass, self.psis, 
                                   self.step)
        self.assertAlmostEqual(np.sum(eDensity0T), np.sum(eDensityLowT))
        np.testing.assert_array_almost_equal(eDensity0T, eDensityLowT,
                                             decimal=eDensity_dec,
                                             verbose=True)

    def testUniformFTFD(self):
        T = 300
        eDensityN, EF = cFermiDiracN(T, self.sheet, self.EigenEs,
                                     self.mass, self.psis, self.step)
        eDensity = cFermiDirac(T, EF, self.EigenEs, self.mass, self.psis,
                                self.step)
        self.assertAlmostEqual(self.step * np.sum(eDensityN), self.sheet)
        self.assertAlmostEqual(self.step * np.sum(eDensity), self.sheet)
        np.testing.assert_array_almost_equal(eDensityN, eDensity,
                                             decimal=eDensity_dec,
                                             verbose=True)

    def testUniformHighT(self):
        T = 300
        eDensityFD, EF_FD = cFermiDiracN(T, self.sheet, self.EigenEs, 
                                     self.mass, self.psis, self.step)
        eDensityB, EF_B = cBoltzmannN(T, self.sheet, self.EigenEs,
                                      self.mass, self.psis, self.step)
        self.assertAlmostEqual(self.step * np.sum(eDensityFD), 
                               self.step * np.sum(eDensityB))


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
    
    # Boltzmann distribution self-consistency test (rho -> EF -> rho)
    def testTriangleWellBoltzmann(self):
        T = 300
        eDensityN, EF = cBoltzmannN(T, self.sheet, self.EigenEs, 
                                   self.mass, self.psis, self.step)
        eDensity = cBoltzmann(T, EF, self.EigenEs, self.mass, self.psis, 
                              self.step)
        self.assertAlmostEqual(self.step * np.sum(eDensity), self.sheet)
        self.assertAlmostEqual(self.step * np.sum(eDensityN), self.sheet)
        np.testing.assert_array_almost_equal(eDensityN, eDensity, 
                                             decimal=eDensity_dec,
                                             verbose=True)

    # Zero-temperature Fermi-Dirac (rho -> EF -> rho)
    def testTriangleWellFermiDirac0T(self):
        eDensityN, EF = cFermiDirac0N(self.sheet, self.EigenEs, 
                                     self.mass, self.psis, self.step)
        eDensity = cFermiDirac0(EF, self.EigenEs, self.mass, self.psis,
                                self.step)
        self.assertAlmostEqual(self.step * np.sum(eDensityN), self.sheet)
        self.assertAlmostEqual(self.step * np.sum(eDensity), self.sheet)
        np.testing.assert_array_almost_equal(eDensityN, eDensity, 
                                             decimal=eDensity_dec, 
                                             verbose=True)

    # test finite temp F-D distribution at T = 0
    def testTriangleWellLowT(self):
        T = 0.01
        eDensityN, EF = cFermiDirac0N(self.sheet, self.EigenEs, 
                                      self.mass, self.psis, self.step)
        eDensityFD = cFermiDirac(T, EF, self.EigenEs, self.mass, self.psis, 
                                 self.step)
        eDensityFD0 = cFermiDirac0(EF, self.EigenEs, self.mass, self.psis, 
                                    self.step)
        np.testing.assert_array_almost_equal(eDensityFD, eDensityFD0,
                                             decimal=eDensity_dec,
                                             verbose=True)

    # Test for finite temperature Fermi-Dirac 
    # (at high T, density of states = Boltzmann)
    def testTriangleWellHighT(self):
        T = 3000
        eDensityFD, EF_FD = cFermiDiracN(T, self.sheet, self.EigenEs, 
                                      self.mass, self.psis, self.step)
        eDensityB, EF_B = cBoltzmannN(T, self.sheet, self.EigenEs,
                                      self.mass, self.psis, self.step)
        self.assertAlmostEqual(np.sum(eDensityFD), np.sum(eDensityB))
        np.testing.assert_array_almost_equal(eDensityFD, eDensityB, 
                                             decimal=eDensity_dec, 
                                             verbose=True)

if __name__ == "__main__":
    unittest.main()

