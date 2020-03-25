#!/usr/bin/env python3
# -*- coding:utf-8 -*-
from context import *
from OneDQuantum import *
from pylab import *
from scipy.constants import hbar, e, m_e, pi 
EgG = {'AlAs': 3.099, 'GaAs': 1.519}
VBO = {'AlAs': -1.33, 'GaAs': -0.80}
F   = {'AlAs': -0.46, 'GaAs': -1.94}
Ep  = {'AlAs': 21.1,  'GaAs': 28.8}
DSO = {'AlAs': 0.28,  'GaAs': 0.341}

def qwells(layers):
    xs = [0.0] + [x for m, x in layers]
    xx = np.linspace(0, xs[-1], 5000, endpoint=False)
    step = xx[1]-xx[0]
    xEg = np.empty(xx.shape)
    xVc = np.empty(xx.shape)
    xm = np.empty(xx.shape)
    xF = np.empty(xx.shape)
    xEp = np.empty(xx.shape)
    xESO = np.empty(xx.shape)
    for n in range(len(layers)):
        xEg[(xx>=xs[n])&(xx<xs[n+1])] = EgG[layers[n][0]]
        xVc[(xx>=xs[n])&(xx<xs[n+1])] = VBO[layers[n][0]]
        xF[(xx>=xs[n])&(xx<xs[n+1])] = F[layers[n][0]]
        xEp[(xx>=xs[n])&(xx<xs[n+1])] = Ep[layers[n][0]]
        xESO[(xx>=xs[n])&(xx<xs[n+1])] = DSO[layers[n][0]]
    xVc += xEg
    band = Band('ZincBlende', xEg, xF, xEp, xESO)
    Es = np.linspace(np.min(xVc), np.min(xVc)+0.2, 100)
    EigenEs = cBandSolve1D(step, Es, xVc, band)
    print(EigenEs)
    plot(xx, xVc)
    psis = cBandFillPsi(step, EigenEs, xVc, band)
    for n in range(EigenEs.size): 
        plot(xx, EigenEs[n] + psis[n, :]/np.max(psis[n,:])*0.02)
    #  plot(xx, xm)
    #  plot(xx, np.ctypeslib.as_array(band.contents.m, (band.contents.N,)))


if __name__ == "__main__":
    qwells([('AlAs', 200), ('GaAs', 700), ('AlAs', 750)])
    show()

# vim: ts=4 sw=4 sts=4 expandtab
