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
    while True: 
        EigenEs = cSimpleSolve1D(step, Es, V, mass)
        if test: 
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
        Vc -= Vc[-1]
        x = np.linspace(0, xmax, 5000)
        V = F * (xmax - x) + Vc
        step = x[1]-x[0]
        Es = np.linspace(np.min(V), np.min(V)+Emax, 100)

if __name__ == "__main__":
    sheet = 5E11 * (ANG)**2
    EigenEs_old = np.zeros(1)
    for x, V, EigenEs, psis, Vc, eDensity, EF in triangle_well(2.02e-4, 
                                                               sheet=sheet, 
                                                               Emax=0.2): 
        nmax = min(EigenEs.size, EigenEs_old.size)
        Diff = EigenEs[0:nmax] - EigenEs_old[0:nmax]
        print("Diff", Diff)
        plot(x, V)
        for n in range(0, EigenEs.size):
            plot(x, EigenEs[n] + psis[n,:]/np.max(psis[n,:])*0.01)
        plot(x, Vc)
        fill_between(x, EF + eDensity/np.max(eDensity)*0.1, 
                     EF, alpha=0.3)
        show()
        if np.max(Diff**2) < 1E-12: 
            print("converged")
            break
        EigenEs_old = EigenEs
    plot(x, V)
    for n in range(0, EigenEs.size):
        plot(x, EigenEs[n] + psis[n,:]/np.max(psis[n,:])*0.01)
    plot(x, Vc)
    fill_between(x, EF + eDensity/np.max(eDensity)*0.1, 
                 EF, alpha=0.3)
    show()
    
# vim: ts=4 sw=4 sts=4 expandtab
