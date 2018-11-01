#!/usr/bin/env python
# -*- coding:utf-8 -*-
from context import *
from OneDQuantum import *
from pylab import *
from scipy.constants import hbar, e, m_e, pi
ANG=1E-10

def triangle_well(F, xmax = 1E3): 
    x = np.linspace(0, xmax, 5000)
    step = x[1]-x[0]
    V = F * (xmax - x)
    mass = 0.067
    Es = np.linspace(0, 0.15, 100)
    EigenEs = cSimpleSolve1D(step, Es, V, mass)
    psis = cSimpleFillPsi(step, EigenEs, V, mass)
    V = np.ascontiguousarray(V[::-1])
    psis = np.ascontiguousarray(psis[:,::-1])
    plot(x, V)
    for n in range(0, EigenEs.size):
        plot(x, EigenEs[n] + psis[n,:]/np.max(psis[n,:])*0.01)
    return x, V, EigenEs, psis

if __name__ == "__main__":
    x, V, EigenEs, psis = triangle_well(2.02e-4)
    print(EigenEs)
    sheet = 1e17 * 50e-8 * (1e-8)**2
    #  eDensity, EF = cBoltzmannN(300, sheet, EigenEs, 
    #                               0.067, psis, x[1]-x[0])
    eDensity, EF = cFermiDirac0N(sheet, EigenEs, 
                                 0.067, psis, x[1]-x[0])
    #  plot(x, EF + eDensity/np.max(eDensity)*0.1, '--')
    fill_between(x, EF + eDensity/np.max(eDensity)*0.1, EF, alpha=0.3)
    print(EF-EigenEs[0], (x[1]-x[0])*np.sum(eDensity)-sheet)
    show()
    
# vim: ts=4 sw=4 sts=4 expandtab
