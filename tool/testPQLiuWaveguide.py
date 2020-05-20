#!/usr/bin/env python3
# -*- coding:utf-8 -*-
from context import *
import SaveLoad
import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":
    with open("../example/PQLiu.json") as f:
        q, s = SaveLoad.loadBoth(f)

    s.updateIndices()
    # betas = np.linspace(2, 3.5, 1000)
    # chi = np.fromiter((s.chiMTM(b) for b in betas), dtype=np.complex)
    # plt.plot(betas, chi.imag)
    # plt.show()
    beta = s.boundModeTM()
    print(beta)
