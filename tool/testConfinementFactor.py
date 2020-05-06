#!/usr/bin/env python3
# -*- coding:utf-8 -*-
import context
import numpy as np
import matplotlib.pyplot as plt
from OptStrata import MaxwellLayer


def firstGaAs():
    wl = 9.4
    hs = [0.7,  1.0,  1.55, 1.4,  1.4,  0.6]
    ns = np.array([1.0, 3.23, 2.87, 3.23, 3.21, 3.23, 2.87, 2.68])
    stratum = MaxwellLayer(wl, ns, hs)
    beta = stratum.boundModeTM(max(ns.real))
    xs = np.linspace(-1, sum(hs)+1, 10000)
    n = stratum.populateIndices(xs)
    Ey, Hx, Ez = stratum.populateMode(beta, xs)
    plt.plot(xs, n.real, label="index")
    plt.plot(xs, abs(Ey)**2, label="field")
    plt.plot(xs, np.angle(Ey)/np.pi, label="phase")
    plt.legend()
    plt.show()


if __name__ == "__main__":
    firstGaAs()
