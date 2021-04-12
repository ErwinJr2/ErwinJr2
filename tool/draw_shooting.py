#!/usr/bin/env python3
# -*- coding:utf-8 -*-
from context import *  # type: ignore # noqa: F401, F403
from ErwinJr2 import SchrodingerLayer
from ErwinJr2.OneDQuantum import OneDSchrodinger as oneds
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from ctypes import c_int, c_double


class GaInAs_AlInAs_Layer(SchrodingerLayer):
    """Sadly the parameters given in PhysRevB.50.8663 is not complete"""
    def __init__(self, xres, layerWidths, layerMtrl, EField=0):
        # fixed offset 0.51 eV
        Vc = {0: 0, 1: 0.6}
        layerVc = [Vc[m] for m in layerMtrl]
        self.layerMtrl = layerMtrl
        super().__init__(xres, layerWidths=layerWidths, EField=EField,
                         layerVc=layerVc)
        self.crystalType = 'ZincBlende'

    def layerMc(self, n):
        Mc = {0: 0.043, 1: 0.072}
        return Mc[self.layerMtrl[n]]

    def populate_material(self):
        N = self.xPoints.size
        self.xF = -0.5*np.ones(N)  # s.t. 1+2F = 0
        self.xEp = 18.3*np.ones(N)
        self.xESO = np.zeros(N)
        self.xEg = 0.79*np.ones(N)
        self.bandParams = (self.xEg, self.xF, self.xEp, self.xESO)


def plot_shooting():
    layers = [100, 45, 8, 30, 8, 10, 100]
    mtrls = [1, 0, 1, 0, 1, 0, 1]
    Mc = {0: 0.043, 1: 0.12}
    Vc = {0: 0, 1: 0.8}
    doubleWell = SchrodingerLayer(0.1, layerWidths=layers,
                                  layerVc=[Vc[m] for m in mtrls],
                                  layerMc=[Mc[m] for m in mtrls],
                                  EField=50)
    doubleWell.populate_x()
    doubleWell.solve_whole()
    print(doubleWell.eigenEs)

    matplotlib.rcParams.update({'font.size': 12})
    ax1 = plt.subplot(121)
    ax1.set_xlabel('Position $z$')
    ax1.set_ylabel('Energy $E$')
    ax1.tick_params(axis='x', which='both', bottom=False, labelbottom=False)
    # ax1.tick_params(axis='y', which='both', left=False, labelleft=False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax1.plot(1, -0.15, ">k", transform=ax1.get_yaxis_transform(), clip_on=False)
    ax1.plot(0, 1, "^k", transform=ax1.get_xaxis_transform(), clip_on=False)

    ax2 = plt.subplot(122)
    ax2.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.spines["left"].set_position(("data", 0))
    plt.setp(ax2.get_xticklabels(), visible=False)
    ax2.tick_params(axis='x', which='both', bottom=False, labelbottom=False)
    ax2.tick_params(axis='y', which='both', left=False, labelleft=False)
    ax2.set_xlabel(r'Normalized $\phi^{\mathrm{end}}$')
    ax2.plot(1, -0.15, ">k", transform=ax2.get_yaxis_transform(), clip_on=False)
    ax2.plot(0, 1, "^k", transform=ax2.get_xaxis_transform(), clip_on=False)

    area = (doubleWell.xPoints > 50) & (doubleWell.xPoints < 300)
    xx = doubleWell.xPoints[area] - 50
    ax1.plot(xx, doubleWell.xVc[area], 'k', linewidth=1)
    color = ['C0', 'C1', 'C2', 'C3']

    Es = np.linspace(-0.1, 0.65, 500)
    N = len(doubleWell.xPoints)
    phi = np.empty((len(Es), N))
    doubleArray = np.ctypeslib.ndpointer(
        dtype=np.float64, ndim=1, flags="C_CONTIGUOUS")
    oneds._clib.rk4.argtypes = [c_double, c_int, c_double, c_double, c_double,
                                doubleArray, doubleArray, doubleArray]
    oneds._clib.rk4.restype = c_double
    phiend = np.array([oneds._clib.rk4(
            c_double(doubleWell.xres), N,
            c_double(0), c_double(1), c_double(Es[n]),
            doubleWell.xVc, doubleWell.xMc, phi[n]
        ) for n in range(len(Es))])
    ax2.plot(phiend/1E14, Es, 'C0',
             label=r'$\phi^{\mathrm{end}}$')
    ax2.plot(phiend*np.exp(-42*np.sqrt(0.8-Es))*2, Es, 'C0',
             ls=(0, (2, 1)),
             label=r'$\mathrm{e}^{-\sqrt{V - E}}\phi^{\mathrm{end}}$')
    for n in range(4):
        E = doubleWell.eigenEs[n]
        y = 0.5*doubleWell.psis[n, :] + E
        c = color[n]
        ax1.plot(xx, y[area], color=c)
        ax1.fill_between(xx, y[area], E,
                         facecolor=c, alpha=0.3)
        ax2.annotate('', xy=(0, E), xycoords=ax1.transData,
                     xytext=(0.01, E), textcoords=ax2.transData,
                     ha="left", va="center",
                     arrowprops=dict(
                        arrowstyle="-", linestyle='--', linewidth=0.5)
                     )
    ax1.set_yticks(doubleWell.eigenEs[:3])
    ax1.set_yticklabels(('$E_1$', '$E_2$', '$E_3$'))
    ax1.set_ylim(-0.15, 0.8)
    ax1.set_xlim(0, 210)
    ax2.set_ylim(-0.15, 0.8)
    ax2.set_xlim(-0.6, 1)
    ax2.legend()

    plt.gcf().set_size_inches(8, 4)
    plt.tight_layout()
    plt.show()
    # plt.savefig('shooting_demo.pdf')


if __name__ == "__main__":
    plot_shooting()
