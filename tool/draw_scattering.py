#!/usr/bin/env python3
# -*- coding:utf-8 -*-
from context import *  # type: ignore # noqa: F401, F403
from ErwinJr2.QCLayers import SchrodingerLayer
import numpy as np
import matplotlib
import matplotlib.pyplot as plt


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


def plot_scattering():
    layers = [100, 45, 10, 35, 100]
    mtrls = [1, 0, 1, 0, 1]
    doubleWell = GaInAs_AlInAs_Layer(0.01, layers, mtrls, 50)
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
    ax1.plot(1, -0.1, ">k", transform=ax1.get_yaxis_transform(), clip_on=False)
    ax1.plot(0, 1, "^k", transform=ax1.get_xaxis_transform(), clip_on=False)

    ax2 = plt.subplot(122)
    ax2.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.spines["left"].set_position(("data", 0))
    plt.setp(ax2.get_xticklabels(), visible=False)
    ax2.tick_params(axis='x', which='both', bottom=False, labelbottom=False)
    ax2.tick_params(axis='y', which='both', left=False, labelleft=False)
    ax2.set_xlabel(r'In-plane momentum $k_{\parallel}$')
    ax2.plot(1, -0.1, ">k", transform=ax2.get_yaxis_transform(), clip_on=False)
    ax2.plot(0, 1, "^k", transform=ax2.get_xaxis_transform(), clip_on=False)

    area = (doubleWell.xPoints > 50) & (doubleWell.xPoints < 230)
    xx = doubleWell.xPoints[area] - 50
    ax1.plot(xx, doubleWell.xVc[area], 'k', linewidth=1)
    color = ['C0', 'C1', 'C2']

    def eofk(k):
        return 0.25 * k**2
    ks = np.linspace(-1, 1, 500)
    for n in range(3):
        E = doubleWell.eigenEs[n]
        y = 0.5*doubleWell.psis[n, :] + E
        c = color[n]
        ax1.plot(xx, y[area], color=c)
        ax1.fill_between(xx, y[area], E,
                         facecolor=c, alpha=0.3)
        ax2.plot(ks, eofk(ks) + E, color=c)
        ax2.annotate('', xy=(0, E), xycoords=ax1.transData,
                     xytext=(0.05, E), textcoords=ax2.transData,
                     ha="left", va="center",
                     arrowprops=dict(
                        arrowstyle="-", linestyle='--', linewidth=0.5)
                     )
    ax1.set_yticks(doubleWell.eigenEs[:3])
    ax1.set_yticklabels(('$E_1$', '$E_2$', '$E_3$'))
    ax1.set_ylim(-0.1, 0.6)
    ax1.set_xlim(0, 200)
    ax2.set_ylim(-0.1, 0.6)

    # scatter
    ki = 0.2
    kf = 0.4
    ei = eofk(0.2)+doubleWell.eigenEs[2]
    ef = eofk(0.4)+doubleWell.eigenEs[1]
    ax2.annotate('', xy=(ki, ei),
                 xycoords=ax2.transData,
                 xytext=(kf, ef),
                 textcoords=ax2.transData,
                 ha="center", va="center",
                 arrowprops=dict(
                    arrowstyle="<-", linewidth=1)
                 )
    ax2.plot((kf, kf), (ei, ef), 'k--', lw=0.5)
    ax2.plot((ki, kf), (ei, ei), 'k--', lw=0.5)
    ax2.text(kf, (ei+ef)/2, r'$\Delta E$',
             va='center', ha='left')
    ax2.text((ki+kf)/2, ei, r'$\Delta k$',
             va='bottom', ha='center')
    plt.gcf().set_size_inches(8, 4)
    plt.tight_layout()
    # plt.show()
    plt.savefig('subband_demo.pdf')


if __name__ == "__main__":
    plot_scattering()
