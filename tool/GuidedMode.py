#!/usr/bin/env python3
# -*- coding:utf-8 -*-
import context
import numpy as np
import matplotlib.pyplot as plt
from OptStrata import Stratum


def firstGaAs():
    stratum = Stratum(
        9.4,
        ['Air', 'GaAs', 'AlGaAs', 'GaAs', 'Active', 'GaAs', 'AlGaAs', 'GaAs'],
        [0, 0, 0.9, 0, 0, 0, 0.9, 0], [0, 0, 0, 0, 0, 0, 0, 0],
        [0.7, 1.0, 1.55, 1.4, 1.4, 0.6])
    stratum.addCustomMtrl('Active', 3.21)
    stratum.updateIndices()
    print("No doping", stratum.index0, stratum.indices, stratum.indexs)
    stratum.dopings = [0, 90, 6, 0.4, 0, 0.4, 6, 30]
    stratum.updateIndices()
    print(stratum.index0, stratum.indices, stratum.indexs)
    beta = stratum.boundModeTM()
    print(beta)

    xs = np.linspace(-2, 10, 10000)
    n = stratum.populateIndices(xs)
    Ey, Hx, Ez = stratum.populateMode(beta, xs)
    plt.plot(xs, n.real)
    plt.plot(xs, abs(Ey)**2)
    plt.show()


def firstPlasmon():
    hs = np.array([25*(31+24)/1000, 0.7])
    stratum = Stratum(
        11.5, ['Pt', 'Active', 'InGaAs', 'InP'],
        [0, 0, 0.53, 0], [0, 0, 0.6, 2], hs)
    stratum.addCustomMtrl('Pt', 3.85+49.2j)
    stratum.addCustomMtrl('Active', 3.38)
    stratum.updateIndices()
    print(stratum.index0, stratum.indices, stratum.indexs)
    beta = stratum.boundModeTM()
    print(beta)

    xs = np.linspace(-2, 6, 1000)
    n = stratum.populateIndices(xs)
    Ey, Hx, Ez = stratum.populateMode(beta, xs)
    plt.plot(xs, n.real)
    plt.plot(xs, abs(Ey)**2)
    plt.show()


if __name__ == '__main__':
    firstPlasmon()
    # firstGaAs()

# vim: ts=4 sw=4 sts=4 expandtab
