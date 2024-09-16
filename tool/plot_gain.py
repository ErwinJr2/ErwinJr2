#!/usr/bin/env python3
# -*- coding:utf-8 -*-

import logging
import os

import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import c as c0, e as e0, epsilon_0 as eps0, h as h, hbar as hbar

from ErwinJr2 import save_load
from ErwinJr2.qclayers import QCLayers

_REPO_ROOT = os.path.dirname(os.path.dirname(__file__))
_TEST_SAMPLE_FILE = os.path.join(_REPO_ROOT, "ErwinJr2/example/test.json")


def main():
    with open(_TEST_SAMPLE_FILE) as f:
        qcl = save_load.qcl_load(f)
    logging.info("Start.")
    qcl.populate_x()
    qcl.solve_whole()
    qcl.period_recognize()
    # plotPotential(qcl)
    # plotWF(qcl)
    # plt.xlabel('Position (Å)')
    # plt.ylabel('Energy (eV)')
    # plt.show()
    logging.info("WF solved.")
    qcl.full_population()
    logging.info("Population solved.")
    wls = np.linspace(3.0, 10.0, 1001)
    gain = qcl.full_gain_spectrum(wls)
    fig = plt.figure()
    axes = fig.add_subplot(111)
    axes.plot(wls, gain)
    axes.axhline(0, ls="--", lw=0.5)
    axes.set_xlabel("Wavelength (μm)")
    axes.set_ylabel("Gain (cm$^{-1}$)")
    fig.show()


def full_gain_spectrum(qcl: QCLayers, wl=None):
    """Perform fully automatic calculation for the gain on wavelength(s)."""
    if wl is None:
        wl = qcl.wl
    neff = qcl.effective_ridx(wl)
    de0 = h * c0 / (wl * 1e-6) / e0
    gain = 0
    for i in range(len(qcl.period_idx)):
        upper = qcl.period_idx[i]
        for j in range(i + 1, len(qcl.period_idx)):
            lower = qcl.period_idx[j]
            for shift in (-1, 0, 1):
                dipole = qcl._dipole(upper, lower, shift) * 1e-8  # to cm
                Eu = qcl.eigen_es[upper]
                El = qcl.eigen_es[lower] - shift * qcl.e_shift
                de = Eu - El
                dpop = qcl.population[i] - qcl.population[j]
                if de < 0:
                    de, dpop = -de, -dpop
                if qcl.include_ifr:
                    gamma_para = qcl._p_gamma[shift][j][i]
                else:
                    gamma_para = qcl.dephasing(upper, lower)
                    gamma_para /= 1e12 * hbar / e0  # unit to Hz
                gamma = (
                    (gamma_para + (qcl.decay_rates[i] + qcl.decay_rates[j]) / 2)
                    * 1e12
                    * hbar
                    / e0
                )
                gainul = dpop * dipole**2 * gamma / (gamma**2 + (de - de0) ** 2)
                # if np.max(np.abs(gainul)) > 4E-13:
                #     qcl.gainul[(upper, lower)] = gainul
                #     print(dpop * dipole**2, h*c0/(de*e0)*1E6, gamma/de,
                #           upper, lower, shift)
                gain = gain + gainul
    # e0^2 / (hbar * c0 * eps0) is dimension 1.. 1E-8 makes unit cm^-1
    gain *= (
        e0**2
        * de0
        * qcl.sheet_density
        / (hbar * neff * c0 * eps0 * qcl.period_l * 1e-8)
    )
    return gain


if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG)
    main()
