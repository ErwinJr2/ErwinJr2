#!/usr/bin/env python3
# -*- coding:utf-8 -*-
from context import *  # type: ignore # noqa: F401, F403
from ErwinJr2 import SaveLoad
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams.update({'font.size': 12})
import pickle
import timeit

with open("../example/new16umgaas-cladding.json") as f:
    qcl = SaveLoad.qclLoad(f)
testRes = (5, 2, 1, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001)
testRes = (5, 2, 1, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005)


def error_test_qcl(qcl, reses, repeats=1):
    odeEigen = {}
    matrixEigen = {}
    odeTimer = {}
    matrixTimer = {}
    for res in reses:
        qcl.xres = res
        qcl.solver = 'ODE'
        qcl.populate_x()
        qcl.Es = np.arange(0.1, 0.3, qcl.Eres/1E3)
        odeTimer[res] = timeit.timeit(
            'qcl.solve_whole()',
            globals=globals(), number=repeats)
        odeEigen[res] = qcl.eigenEs
        qcl.solver = 'matrix'
        qcl.matrixSigma = 0.2
        qcl.matrixEigenCount = 20
        qcl.populate_x()
        matrixTimer[res] = timeit.timeit(
            'qcl.solve_whole()',
            globals=globals(), number=repeats)
        matrixEigen[res] = qcl.eigenEs
        print(res, 'calculated')
        print('timer: ', odeTimer[res], matrixTimer[res])
    return odeEigen, matrixEigen, odeTimer, matrixTimer


def plot_error():
    with open('errorPorf.pkl', 'rb') as f:
        odeEigen, matrixEigen, _, _ = pickle.load(f)
        # print(odeEigen, matrixEigen)
    testSet = np.array(testRes)
    # o = np.array([odeEigen[r] for r in testSet])
    o = np.zeros((len(testSet), 20))
    for n, r in enumerate(testSet):
        eigen = odeEigen[r]
        o[n, :len(eigen)] = eigen
    m = np.array([matrixEigen[r] for r in testSet])
    plt.loglog(testSet, np.abs(o[:, 0]-o[-1, 0]), 'C0.', ms=2,
               label='shooting')
    plt.loglog(testSet*1.05, np.abs(m[:, 0]-m[-1, 0]), 'C1.', ms=2,
               label='matrix')
    plt.loglog(testSet, np.abs(o[:, 1:]-o[-1, 1:]), 'C0.', ms=2)
    plt.loglog(testSet*1.05, np.abs(m[:, 1:]-m[-1, 1:]), 'C1.', ms=2)
    plt.loglog(testSet, testSet**2*3E-3, 'C2--', label='$O(h^2)$')
    plt.loglog(testSet, testSet*5E-4, 'C3--', label='$O(h)$')
    plt.legend()
    plt.xlabel("Step size (Å)")
    plt.ylabel("Eigen Energy Error (meV)")
    plt.gcf().set_size_inches(10, 6)
    plt.tight_layout()
    # plt.show()
    plt.savefig('qcl_error_prof.pdf')


def plot_time():
    with open('errorPorf.pkl', 'rb') as f:
        _, _, odeTimer, matrixTimer = pickle.load(f)
        # print(odeEigen, matrixEigen)
    testSet = np.array(testRes)
    length = np.sum(qcl.layerWidths)*qcl.repeats // testSet
    o = np.array([odeTimer[r] for r in testSet]) / 50
    m = np.array([matrixTimer[r] for r in testSet]) / 50
    plt.loglog(length, o, 'C0.', label='shooting')
    plt.loglog(length, m, 'C1.', label='matrix')
    from scipy.stats import linregress
    oslope, ointercept, _, _, _ = linregress(length, o)
    mslope, mintercept, _, _, _ = linregress(length, m)
    print(oslope, ointercept)
    print(mslope, mintercept)
    print(mslope/oslope)
    plt.loglog(length, oslope*length, 'C0--')
    plt.loglog(length, mslope*length, 'C1--')
    plt.legend()
    plt.xlabel("Number of points")
    plt.ylabel("Computing time (s)")
    plt.gcf().set_size_inches(8, 6)
    plt.tight_layout()
    # plt.show()
    plt.savefig('qcl_time_prof.pdf')


if __name__ == "__main__":
    # with open('errorPorf.pkl', 'wb') as f:
    #     pickle.dump(error_test_qcl(qcl, testRes, 1), f)

    plot_error()
    # plot_time()
