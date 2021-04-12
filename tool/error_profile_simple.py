#!/usr/bin/env python3
# -*- coding:utf-8 -*-
# from plotter import plot_band
from context import *  # type: ignore # noqa: F401, F403
from plotter import plot_band
from ErwinJr2.QCLayers import SchrodingerLayer, QCLayers
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams.update({'font.size': 12})
import pickle
from scipy.special import ai_zeros
from scipy.constants import hbar, e, m_e
ANG = 1E-10

mass = 0.067  # Typical parameters for GaAs
F = 2e-4
xmax = 1E3
nmax = 8
an = ai_zeros(nmax)[0]
EigenEs_th = -(hbar**2*F**2/(2*m_e*mass*e*ANG**2))**(1/3)*an
de_th = EigenEs_th[1:] - EigenEs_th[0]
qcl = SchrodingerLayer(EField=F*1E5, layerWidths=[xmax], layerVc=[xmax*F],
                       layerMc=[mass])
qcl = QCLayers('GaAs', ['AlGaAs'], [0.0], layerWidths=[xmax],
               EField=F*1E5, repeats=1)


def error_test_simple(reses):
    odeError = {}
    matrixError = {}

    for res in reses:
        qcl.xres = res
        qcl.solver = 'ODE'
        qcl.populate_x()
        qcl.Es = np.arange(0.15, 0.4, qcl.Eres/1E3)
        # odeError[res] = qcl.solve_whole()[:nmax] - EigenEs_th
        eigen = qcl.solve_whole()[:nmax]
        odeError[res] = eigen[1:] - eigen[0]
        qcl.solver = 'matrix'
        qcl.populate_x()
        qcl.matrixSigma = 0.2
        qcl.matrixEigenCount = 10
        # matrixError[res] = qcl.solve_whole() - EigenEs_th
        eigen = qcl.solve_whole()
        eigen = eigen[qcl.psis[:, 0] < 1E-2][:nmax]
        matrixError[res] = eigen[1:] - eigen[0]
        print(res, 'complete.')
    return odeError, matrixError


if __name__ == "__main__":

    testSet = (100, 50, 20, 10, 5, 1, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01)
    testSet = np.logspace(-1.7, 1.7)[::-1]
    # print(testSet)
    # with open('errorPorf_simple.pkl', 'wb') as f:
    #     pickle.dump(error_test_simple(testSet), f)

    # plot Error
    with open('errorPorf_simple.pkl', 'rb') as f:
        odeError, matrixError = pickle.load(f)
        # print(odeEigen, matrixEigen)
    o = np.array([odeError[r] for r in testSet])*1E3
    m = np.array([matrixError[r] for r in testSet])*1E3
    testSet = np.array(testSet)
    plt.loglog(testSet, np.abs(o[:, 0]-o[-1, 0]), 'C0.', ms=2, label='shooting')
    plt.loglog(testSet, np.abs(m[:, 0]-m[-1, 0]), 'C1.', ms=2, label='matrix')
    plt.loglog(testSet, np.abs(o[:, 1:]-o[-1, 1:]), 'C0.', ms=2)
    plt.loglog(testSet, np.abs(m[:, 1:]-m[-1, 1:]), 'C1.', ms=2)
    plt.loglog(testSet, testSet**2*7E-3, 'C2--', label='$O(h^2)$')
    plt.loglog(testSet, testSet**4/5E5, 'C3--', label='$O(h^4)$')
    plt.legend()
    plt.xlabel("Step size (Å)")
    plt.ylabel("Eigen Energy Error (meV)")

    qcl.xres = 0.5
    qcl.populate_x()
    qcl.Es = np.arange(0.15, 0.35, qcl.Eres/1E3)
    qcl.solve_whole()
    eigen = qcl.eigenEs[:nmax]
    inset_axes = plt.gca().inset_axes([0.62, 0.15, 0.36, 0.3])
    plot_band(inset_axes, qcl)
    inset_axes.set_xlabel('Position (Å)', fontsize=10)
    inset_axes.set_ylabel('Energy (eV)', fontsize=10)
    inset_axes.tick_params(axis='both', labelsize=10)
    plt.gcf().set_size_inches(8, 6)
    plt.tight_layout()
    # plt.show()
    plt.savefig('trw_error_prof.pdf')
