#!/usr/bin/env python3
# -*- coding:utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np
from numpy import pi

from ErwinJr2.opt_strata import OptStrata


def firstGaAs():
    wl = 9.4
    stratum = OptStrata(
        wl,
        ["Air", "GaAs", "AlGaAs", "GaAs", "Active", "GaAs", "AlGaAs", "GaAs"],
        [0, 0, 0.9, 0, 0, 0, 0.9, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0.7, 1.0, 1.55, 1.4, 1.4, 0.6],
        cstm_idx={"Active": 3.21},
    )
    # print("No doping", stratum.index0, stratum.indices, stratum.indexS)
    stratum.dopings = [0, 90, 6, 0.4, 0, 0.4, 6, 30]
    stratum.update_indices()
    print(stratum.index_0, stratum.indices, stratum.index_s)
    beta = stratum.bound_mode_tm()
    print(beta)
    print(4 * pi / (wl * 1e-4) * beta.imag)

    xs = np.linspace(-2, 10, 10000)
    n = stratum.populate_indices(xs)
    Ey, Hx, Ez = stratum.populate_mode(beta, xs)
    plt.plot(xs, n.real)
    plt.plot(xs, abs(Ey) ** 2)
    for xx in stratum.populate_mtrl(xs):
        plt.plot(xs[xx], n.real[xx], "r")
    plt.show()


def firstPlasmon():
    wl = 11.5
    hs = np.array([0.1, 25 * (31 + 24) / 1000, 0.7])
    stratum = OptStrata(
        wl,
        ["Air", "Pt", "Active", "InGaAs", "InP"],
        [0, 0, 0, 0.53, 0],
        [0, 0, 0, 0.6, 2],
        hs,
        cstm_idx={"Pt": 3.85 + 49.2j, "Active": 3.38},
    )
    # print(stratum.indices[1])
    stratum.indices[2] += 0.1
    # stratum.indexS += 0.2
    print(stratum.index_0, stratum.indices, stratum.index_s)
    beta = stratum.bound_mode_tm()
    print(beta)
    print(4 * pi / (wl * 1e-4) * beta.imag)

    xs = np.linspace(-0.5, 6, 1000)
    n = stratum.populate_indices(xs)
    Ey, Hx, Ez = stratum.populate_mode(beta, xs)
    print(stratum.confinementy(beta, xs, Ey) * 100)
    plt.plot(xs, n.real, label="n")
    plt.plot(xs, n.imag, label="k")
    plt.plot(xs, abs(Ey) ** 2, label="Field")
    # plt.plot(xs, np.angle(Ey), label="Phase")
    plt.ylim(-1, 4)
    plt.legend()
    plt.show()


if __name__ == "__main__":
    firstPlasmon()
    firstGaAs()
