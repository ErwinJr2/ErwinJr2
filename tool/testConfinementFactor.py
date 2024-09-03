#!/usr/bin/env python3
# -*- coding:utf-8 -*-
import context  # type: ignore # noqa: F401, F403
import numpy as np
from numpy import sqrt, pi, abs
import matplotlib.pyplot as plt
from ErwinJr2.OptStrata import MaxwellLayer, MaxwellLayer_anisotropic


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
    stratum = MaxwellLayer_anisotropic(wl, hs, ns)
    beta = stratum.boundModeTM(max(ns.real))
    print(4*pi/(wl*1E-4)*beta.imag)
    xs = np.linspace(-3, sum(hs)+3, 500000)
    nn, _ = stratum.populateIndices(xs)
    Ey, Hx, Ez = stratum.populateMode(beta, xs)
    # plt.plot(xs, nx.real, label="index")
    # plt.plot(xs, abs(Ey)**2, label="field")
    # plt.plot(xs, np.angle(Ey)/np.pi, label="phase")

    ac = ((xs >= sum(hs[:3])) & (xs <= sum(hs[:4])))
    Gamma = beta * np.trapezoid(nn[ac] * (Ey[ac]**2 - Ez[ac]**2), xs[ac])/(
                np.trapezoid((nn*Ey)**2, xs))
    Gamma = Gamma.real
    Gamma0 = beta * np.trapezoid(nn[ac] * Ey[ac]**2, xs[ac])/(
                 np.trapezoid((nn*Ey)**2, xs))
    Gamma0 = Gamma0.real
    Gamma1 = beta.real * np.trapezoid(nn[ac].real * abs(Ey[ac])**2, xs[ac])/(
        np.trapezoid(abs(nn.real*Ey)**2, xs))
    Gamma2 = np.trapezoid(nn[ac].real * abs(Ey[ac])**2, xs[ac])/(
        np.trapezoid(nn.real * np.abs(Ey)**2, xs))
    Gamma3 = np.trapezoid(abs(Ey[ac])**2, xs[ac])/(
        np.trapezoid(np.abs(Ey)**2+np.abs(Ez)**2, xs))

    chis = np.linspace(0, 2, 100)
    g0 = 2*pi/(wl*1E-4*3.21)*chis
    betas = np.empty(chis.shape, dtype=np.complex128)
    dbeta = np.empty(chis.shape)
    stratum.indexy = np.copy(stratum.indices)
    for n, chi in enumerate(chis):
        stratum.indexy[3] = sqrt(3.21**2 - 1j*chi)
        betas[n] = stratum.boundModeTM()
    dbeta = (betas.imag - beta.imag)
    gs = -4*pi/(wl*1E-4) * dbeta

    plt.plot(chis, g0*Gamma, label="Isotropic")
    plt.plot(chis, g0*Gamma0, label="Anisotropic")
    plt.plot(chis, g0*Gamma1, label="Anisotropic-real")
    # plt.axhline(0, ls='--', color='k', lw=0.5)
    plt.plot(chis, g0*Gamma2, label="Thesis")
    plt.plot(chis, g0*Gamma3, label="APR 1, 011307 (2014)")
    plt.plot(chis, gs, label="strict")
    plt.xlabel(r"imag $\chi$")
    plt.ylabel("$g$ (cm-1)")
    plt.legend()
    plt.show()

    gs[0] = 0
    plt.plot(chis, g0*Gamma / gs - 1, label="Isotropic")
    plt.plot(chis, g0*Gamma0 / gs - 1, label="Anisotropic")
    plt.plot(chis, g0*Gamma1 / gs - 1, label="Anisotropic-real")
    plt.axhline(0, ls='--', color='k', lw=0.5)
    # plt.plot(chis, g0*Gamma2 / gs - 1, label="Thesis")
    # plt.plot(chis, g0*Gamma3 / gs - 1, label="APR 1, 011307 (2014)")
    plt.xlabel(r"imag $\chi$")
    plt.ylabel(r"$\Delta g/g_0$")
    plt.ticklabel_format(axis="x", style="sci", scilimits=(0, 0))
    plt.legend()
    plt.show()


def findOutOfPhase():
    wl = 9.4
    N = 5
    hs = [0.7,  1.0,  1.55] + [0.2, 0.1]*N + [1.4,  0.6]
    ns = np.array(
        [1.0, 3.28, 2.87, 3.28] + [3.21, 3.28+1j]*N + [3.28, 2.87, 3.28])

    # plasmon loss
    plasmonUnit = 8.9698E-5 * wl**2
    nue = 1.0
    gammaUnit = 5.3088E-3
    me = np.array([1, 0.067, 0.15, 0.067] + [0.67, 0.67]*N
                  + [0.067, 0.15, 0.067])
    dopings = np.array([0, 90, 6, 0.4] + [0, 0]*N + [0.4, 6, 30])
    ns = sqrt(ns**2 - plasmonUnit * dopings /
              me / (1 + 1j * gammaUnit * wl * nue))
    # print(ns)
    stratum = MaxwellLayer(wl, hs, ns)
    beta = stratum.boundModeTM(max(ns.real))
    print(beta, 4*pi/(wl*1E-4)*beta.imag)

    xs = np.linspace(-3, sum(hs)+3, 5000)
    nn = stratum.populateIndices(xs)
    Ey, Hx, Ez = stratum.populateMode(beta, xs)

    acs = [((xs >= sum(hs[:3+2*n])) & (xs <= sum(hs[:4+2*n])))
           for n in range(N)]
    # ac1 = ((xs >= sum(hs[:3])) & (xs <= sum(hs[:4])))
    # ac2 = ((xs >= sum(hs[:5])) & (xs <= sum(hs[:6])))

    Gamma0 = beta * sum(np.trapezoid(nn[ac] * Ey[ac]**2, xs[ac])
                        for ac in acs)/(
                 np.trapezoid((nn*Ey)**2, xs))
    Gamma1 = beta.real * sum(np.trapezoid(nn[ac].real * abs(Ey[ac])**2, xs[ac])
                             for ac in acs)/(
        np.trapezoid(abs(nn.real*Ey)**2, xs))
    print(Gamma0, Gamma1)

    plt.plot(xs, nn.real, label="n")
    plt.plot(xs, nn.imag, label="k")
    plt.plot(xs, abs(Ey)**2, label="field")
    plt.plot(xs, (Ey**2).real, label="field real")
    # plt.plot(xs, np.angle(Ey), label="phase")
    plt.ylim(-1, 4)
    plt.legend()
    plt.show()


if __name__ == "__main__":
    # firstGaAs()
    findOutOfPhase()
