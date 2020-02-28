#!/usr/bin/env python
# -*- coding:utf-8 -*-
from context import *
from OneDQuantum.OneDMaxwell import chiMTM, boundModeTM, genModeTM
from pylab import *

if __name__ == '__main__':
    # Plot the TM0 mode in [John Chilwell and Ian Hodgkinson, 
    # J. Opt. Soc. Am. A 1, 742-753 (1984)
    hs = np.array([500]*4)
    indices = np.array([1.66, 1.53, 1.60, 1.66])
    n0 = 1.0
    ns = 1.5
    wl = 632.8
    beta = boundModeTM(max(indices), wl, hs, indices, n0, ns)

    xs = np.linspace(-500, 2500, 10000)
    lsum = np.zeros(len(hs)+1)
    lsum[1:] = np.cumsum(hs)
    n = np.piecewise(xs, [(xs>=lsum[i]) & (xs<lsum[i+1]) 
    	for i in range(len(hs))], indices)
    n[xs<0] = n0
    n[xs>=lsum[-1]] = ns

    Ey, Hx, Ez = genModeTM(beta, xs, wl, hs, indices, n0, ns)
    plot(xs, n)
    plot(xs, Ey.real)
    show()