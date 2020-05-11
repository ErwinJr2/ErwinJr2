#!/usr/bin/env python3
# -*- coding:utf-8 -*-
import context
import numpy as np
from numpy import sqrt, pi, abs
import matplotlib.pyplot as plt
from OptStrata import MaxwellLayer_anisotropic as MaxwellLayer


def firstGaAs():
    wl = 9.4
    hs = [0.7,  1.0,  1.55, 1.4,  1.4,  0.6]
    ns = np.array([1.0, 3.28, 2.87, 3.28, 3.21, 3.28, 2.87, 3.28])

    # plasmon loss
    plasmonUnit = 8.9698E-5 * wl**2
    nue = 1.0
    gammaUnit = 5.3088E-3
    me = np.array([1, 0.067, 0.15, 0.067, 0.067, 0.067, 0.15, 0.067])
    dopings = np.array([0, 90, 6, 0.4, 0, 0.4, 6, 30])
    ns = sqrt(ns**2 - plasmonUnit * dopings /
              me / (1 + 1j * gammaUnit * wl * nue))
    # print(ns)
    stratum = MaxwellLayer(wl, hs, ns)
    beta = stratum.boundModeTM(max(ns.real))
    print(4*pi/(wl*1E-4)*beta.imag)
    xs = np.linspace(-3, sum(hs)+3, 500000)
    nn = stratum.populateIndices(xs)
    Ey, Hx, Ez = stratum.populateMode(beta, xs)
    # plt.plot(xs, nx.real, label="index")
    # plt.plot(xs, abs(Ey)**2, label="field")
    # plt.plot(xs, np.angle(Ey)/np.pi, label="phase")

    ac = ((xs >= sum(hs[:3])) & (xs <= sum(hs[:4])))
    Gamma = beta * np.trapz(nn[ac] * (Ey[ac]**2 - Ez[ac]**2), xs[ac])/(
                np.trapz((nn*Ey)**2, xs))
    Gamma = Gamma.real
    Gamma0 = beta * np.trapz(nn[ac] * Ey[ac]**2, xs[ac])/(
                 np.trapz((nn*Ey)**2, xs))
    Gamma0 = Gamma0.real
    Gamma1 = beta.real * np.trapz(nn[ac].real * abs(Ey[ac])**2, xs[ac])/(
        np.trapz(abs(nn.real*Ey)**2, xs))
    Gamma2 = np.trapz(nn[ac].real * abs(Ey[ac])**2, xs[ac])/(
        np.trapz(nn.real * np.abs(Ey)**2, xs))
    Gamma3 = np.trapz(abs(Ey[ac])**2, xs[ac])/(
        np.trapz(np.abs(Ey)**2+np.abs(Ez)**2, xs))

    chis = np.linspace(0, 0.03, 100)
    g0 = 2*pi/(wl*1E-4*3.21)*chis
    betas = np.empty(chis.shape, dtype=np.complex128)
    dbeta = np.empty(chis.shape)
    stratum.indexy = np.copy(stratum.indices)
    for n, chi in enumerate(chis):
        stratum.indexy[3] = sqrt(3.21**2 - 1j*chi)
        betas[n] = stratum.boundModeTM()
    dbeta = (betas.imag - beta.imag)
    gs = -4*pi/(wl*1E-4) * dbeta
    # plt.plot(chis, gs, label="strict")
    # gs = 0
    gs[0] = 0
    plt.plot(chis, g0*Gamma / gs - 1, label="Gamma")
    plt.plot(chis, g0*Gamma0 / gs - 1, label="Gamma0")
    plt.plot(chis, g0*Gamma1 / gs - 1, label="Gamma1")
    plt.axhline(0, ls='--', color='k', lw=0.5)
    # plt.plot(chis, g0*Gamma2 - gs, label="Gamma2")
    # plt.plot(chis, g0*Gamma3 - gs, label="Gamma3")
    plt.legend()
    plt.show()


if __name__ == "__main__":
    firstGaAs()
