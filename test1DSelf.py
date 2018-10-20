#!/usr/bin/env python
# -*- coding:utf-8 -*-
from pylab import *
from scipy.constants import hbar, e, m_e, pi
from OneDSchrodinger import * 
from OneDThermal import *
from OneDMaxwell import *
ANG=1E-8 # angstrom in cm

def triangle_well(F, sheet, Vmax=0.4, Emax=0.15): 
    xmax = Vmax / F
    x = np.linspace(0, xmax, 5000)
    step = x[1]-x[0]
    V = F * (xmax - x)
    mass = 0.067
    EigenEs = []
    Es = np.linspace(0, Emax, 100)
    for i in range(40):
        print(i, "-th interation")
        EigenEs_old = EigenEs
        EigenEs = cSimpleSolve1D(step, Es, V, mass)
        print(EigenEs)
        psis = cSimpleFillPsi(step, EigenEs, V, mass)
        eDensity, EF = cFermiDirac0N(sheet, EigenEs, 
                                     0.067, psis, x[1]-x[0])
        plot(x, V[::-1])
        for n in range(0, EigenEs.size):
            plot(x, EigenEs[n] + psis[n,::-1]/np.max(psis[n,:])*0.01)
        fill_between(x, EF + eDensity[::-1]/np.max(eDensity)*0.1, 
                     EF, alpha=0.3)
        if i > 0:
            nmax = min(EigenEs.size, EigenEs_old.size)
            Diff = EigenEs[0:nmax] - EigenEs_old[0:nmax]
            print("Diff", Diff)
            if np.max(Diff**2) < 1E-12: 
                print("converged")
                break
        Vc = -cCoulombField(x[1]-x[0], -eDensity, 13)
        xmax = (Vmax+Vc[-1]) / F
        Vc -= Vc[-1]
        plot(x, Vc[::-1])
        print("xmax", xmax)
        show()
        #  print(EF-EigenEs[0], (x[1]-x[0])*np.sum(eDensity)-sheet)
        x = np.linspace(0, xmax, 5000)
        step = x[1]-x[0]
        V = F * (xmax - x) + Vc
        Es = np.linspace(np.min(V), np.min(V)+Emax, 100)
    V = np.ascontiguousarray(V[::-1])
    psis = np.ascontiguousarray(psis[:,::-1])
    Vc = np.ascontiguousarray(Vc[::-1])
    eDensity = np.ascontiguousarray(eDensity[::-1])
    return x, V, EigenEs, psis, Vc, eDensity, EF

if __name__ == "__main__":
    sheet = 5E11 * (ANG)**2
    x, V, EigenEs, psis, Vc, eDensity, EF = triangle_well(2.02e-4, 
                                                          sheet=sheet, 
                                                          Emax=0.2)
    plot(x, V)
    for n in range(0, EigenEs.size):
        plot(x, EigenEs[n] + psis[n,:]/np.max(psis[n,:])*0.01)
    plot(x, Vc)
    show()
    
# vim: ts=4 sw=4 sts=4 expandtab
