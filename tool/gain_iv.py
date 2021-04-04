#!/usr/bin/env python3
# -*- coding:utf-8 -*-
from context import *
import SaveLoad
import numpy as np
import matplotlib.pyplot as plt


def gainSpec(qcl):
    # qcl.solver = 'matrix'
    # qcl.EField = 50.5
    wls = np.linspace(1, 15, 1000)
    qcl.populate_x()
    qcl.solve_whole()
    qcl.period_recognize()
    qcl.full_population()
    print(1/qcl.flow, qcl.current, qcl.full_gain_spectrum(8.0))
    gain = qcl.full_gain_spectrum(wls)
    plt.plot(wls, gain)
    # gg = 0
    # for key in qcl.gainul:
    #     plt.plot(wls, qcl.gainul[key], label=str(key))
    #     gg += qcl.gainul[key]
    # plt.legend()
    # plt.show()


def iv(qcl):
    qcl.solver = 'matrix'
    qcl.repeats = 4
    bias = np.linspace(1, 1.5*qcl.EField)
    current = np.empty(bias.shape)
    for n, v in enumerate(bias):
        qcl.EField = v
        qcl.populate_x()
        qcl.solve_whole()
        qcl.period_recognize()
        qcl.full_population()
        current[n] = qcl.current
    plt.plot(current, bias)
    plt.show()


if __name__ == "__main__":
    # with open("../example/PQLiu.json") as f:
    with open("../example/std8um.json") as f:
        qcl = SaveLoad.qclLoad(f)

    gainSpec(qcl)
    # iv(qcl)
