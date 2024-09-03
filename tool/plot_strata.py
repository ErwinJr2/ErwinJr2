#!/usr/bin/env python3
# -*- coding:utf-8 -*-
from ErwinJr2 import SaveLoad
from ErwinJr2.OptStrata import OptStrata
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

matplotlib.rcParams.update({'font.size': 12})
with open("../ErwinJr2/example/PQLiu.json", 'r') as f:
    strata = SaveLoad.optLoad(f)

strata = OptStrata(
    wl=8.0,
    materials=[
        "Air", "InP", "InP", "InP", "InP", "InxGa1-xAs", "Active Core",
        "InxGa1-xAs", "InP", "InP"],
    moleFracs=[0.0, 0.0, 0.0, 0.0, 0.0, 0.53, 0.53, 0.53, 0.53, 0.0],
    dopings=[0.0, 1000.0, 80.0, 2.0, 1.0, 0.5, 0, 0.5, 0, 0.0],
    Ls=[1.0, 0.01, 0.35, 0.5, 2.5, 0.5, 2.0895, 0.5, 5.0, 1.0],
    # the properties for the active core
    cstmIndx={"Active Core": 3.28+0j},
    cstmPrd={"Active Core": [597.0, 35]},
    cstmGain={"Active Core": 39.6}
)

beta = strata.boundModeTM()
xs = np.linspace(-1, sum(strata.Ls[1:]), 5000)
nx = strata.populateIndices(xs)
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
