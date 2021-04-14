#!/usr/bin/env python3
# -*- coding:utf-8 -*-
from context import *  # type: ignore # noqa: F401, F403
from ErwinJr2 import SaveLoad
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

matplotlib.rcParams.update({'font.size': 12})
with open("../ErwinJr2/example/PQLiu.json", 'r') as f:
    strata = SaveLoad.optLoad(f)

xs = np.linspace(-1, sum(strata.Ls[1:]), 5000)
nx = strata.populateIndices(xs)
beta = strata.boundModeTM()
Ey, _, _ = strata.populateMode(beta, xs)
print("beta = ", beta)
print("confinement: ", strata.confinementy(beta, xs, Ey))

rIdxAxes = plt.gca()
modeAxes = rIdxAxes.twinx()
modeAxes.set_frame_on(False)
modeAxes.get_yaxis().set_visible(False)
rIdxAxes.set_xlabel('Position (μm)')
rIdxAxes.set_ylabel('Refractive index $n$')

lNReal = rIdxAxes.plot(xs, nx.real, 'k', label=r'$\mathrm{Re}[n]$')
lNImag = rIdxAxes.plot(xs, nx.imag, 'orange', label=r'$\mathrm{Im}[n]$')
lMode = modeAxes.plot(xs, np.abs(Ey)**2, color='C0', label=r'Mode')

lns = lNReal + lNImag + lMode
labs = [line.get_label() for line in lns]
plt.legend(lns, labs)
rIdxAxes.set_ylim(0, 5)
modeAxes.set_ylim(bottom=0)
# plt.xlim(right=20)

plt.gcf().set_size_inches(8, 6)
plt.tight_layout()
# plt.savefig('qcl_waveguide.pdf')
plt.show()
