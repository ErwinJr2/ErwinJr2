#!/usr/bin/env python
# -*- coding:utf-8 -*-
from pylab import *
from scipy.constants import hbar, e, m_e, pi
from OneDSchrodinger import * 
from OneDThermal import *
ANG=1E-8 # angstrom in cm

def triangle_well(F, sheet, xmax = 1E3): 
    x = np.linspace(0, xmax, 5000)
    step = x[1]-x[0]
    V = F * (xmax - x)
    mass = 0.067
    EigenEs = []
    for i in range(10):
        print(i, "-th interation")
        Es = np.linspace(0, 0.15, 100)
        EigenEs_old = EigenEs
        EigenEs = cSimpleSolve1D(step, Es, V, mass)
        print(EigenEs)
        psis = cFillPsi(step, EigenEs, V, mass)
        eDensity, EF = cFermiDirac0N(sheet, EigenEs, 
                                     0.067, psis, x[1]-x[0])
        Vc = cCoulombField(x[1]-x[0], -eDensity, 3)
        plot(x, V[::-1])
        for n in range(0, EigenEs.size):
            plot(x, EigenEs[n] + psis[n,::-1]/np.max(psis[n,:])*0.01)
        fill_between(x, EF + eDensity[::-1]/np.max(eDensity)*0.1, 
                     EF, alpha=0.3)
        plot(x, Vc[::-1])
        if i > 0:
            print("Diff", EigenEs - EigenEs_old)
        show()
        #  print(EF-EigenEs[0], (x[1]-x[0])*np.sum(eDensity)-sheet)
        V = F * (xmax - x) + Vc
    V = np.ascontiguousarray(V[::-1])
    psis = np.ascontiguousarray(psis[:,::-1])
    Vc = np.ascontiguousarray(Vc[::-1])
    eDensity = np.ascontiguousarray(eDensity[::-1])
    return x, V, EigenEs, psis, Vc, eDensity, EF

if __name__ == "__main__":
    sheet = 1E11 * (ANG)**2
    x, V, EigenEs, psis, Vc, eDensity, EF = triangle_well(2.02e-4, sheet=sheet)
    
# vim: ts=4 sw=4 sts=4 expandtab
