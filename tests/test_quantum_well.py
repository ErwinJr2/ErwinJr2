#!/usr/bin/env python3
# -*- coding:utf-8 -*-
import typing
import unittest

import numpy as np
from scipy.constants import e as e0, electron_mass as m0, hbar as hbar

from ErwinJr2.qc_layers import QCLayers, SchrodingerLayer

"""This unit test is to compare simulation with PhysRevB.50.8663"""


class GaInAs_AlInAs_Layer(SchrodingerLayer):
    """Sadly the parameters given in PhysRevB.50.8663 is not complete"""

    def __init__(self, xres, layerWidths, layerMtrl, EField=0):
        # fixed offset 0.51 eV
        Vc = {0: 0, 1: 0.51}
        layerVc = [Vc[m] for m in layerMtrl]
        self.layerMtrl = layerMtrl
        super().__init__(
            xres, layer_widths=layerWidths, e_field=EField, layer_vc=layerVc
        )
        self.crystal_type = "ZincBlende"

    def layer_mc(self, n):
        Mc = {0: 0.043, 1: 0.072}
        return Mc[self.layerMtrl[n]]

    def populate_material(self):
        N = self.x_points.size
        self.x_f = -0.5 * np.ones(N)  # s.t. 1+2F = 0
        self.x_ep = 18.3 * np.ones(N)
        self.x_eso = np.zeros(N)
        self.x_eg = 0.79 * np.ones(N)
        self.band_params = (self.x_eg, self.x_f, self.x_ep, self.x_eso)


def validate_ode(q: typing.Union[GaInAs_AlInAs_Layer, QCLayers]):
    interface = np.abs(np.diff(q.x_vc)) > 0.1
    interface = interface[1:] | interface[:-1]
    residules = []
    for t in range(len(q.eigen_es)):
        E = q.eigen_es[t]
        psi = q.psis[t]
        x_eg, x_f, x_ep, x_eso = q.band_params
        meff = 1 / (
            1
            + 2 * x_f
            + x_ep / 3 * (1 / (E - q.x_vc + x_eg + x_eso) + 2 / (E - q.x_vc + x_eg))
        )
        # interpolate half grid
        meff = (meff[1:] + meff[:-1]) / 2
        residule = (
            -np.diff(hbar**2 / (2 * m0 * meff * e0) * np.diff(psi))
            / (q.x_step * 1e-10) ** 2
            + q.x_vc[1:-1] * psi[1:-1]
            - E * psi[1:-1]
        )
        residules.append(residule)
        inter_res = residule[interface]
        # The interface has larger error
        np.testing.assert_almost_equal(inter_res, np.zeros(inter_res.shape), decimal=2)
        other = residule[~interface]
        np.testing.assert_almost_equal(other, np.zeros(other.shape))
    return np.array(residules)


class TestQuantumWell(unittest.TestCase):
    def test_single_well(self):
        # Test PhysRevB.50.8663 original model
        layers = [100, 52, 100]
        mtrls = [1, 0, 1]
        singleWell = GaInAs_AlInAs_Layer(0.01, layers, mtrls)
        singleWell.populate_x()
        singleWell.solve_whole()
        validate_ode(singleWell)
        Ei = singleWell.eigen_es[1]
        Ej = singleWell.eigen_es[0]
        self.assertAlmostEqual((Ei - Ej) / (0.381 - 0.123), 1, 1)
        self.assertAlmostEqual(abs(singleWell.dipole(1, 0)) / 15.3, 1, 0)

        # Test QCLayers
        qcLayers = QCLayers(
            x_res=0.01, layer_widths=layers, layer_matrls=mtrls, repeats=1, temp=10.0
        )
        qcLayers.populate_x()
        qcLayers.solve_whole()
        validate_ode(qcLayers)
        Ei = qcLayers.eigen_es[1]
        Ej = qcLayers.eigen_es[0]
        self.assertAlmostEqual((Ei - Ej) / (0.381 - 0.123), 1, 1)
        self.assertAlmostEqual(abs(qcLayers.dipole(1, 0)) / 15.3, 1, 0)

    def test_double_well(self):
        layers = [200, 59, 13, 24, 200]
        mtrls = [1, 0, 1, 0, 1]
        qcLayers = QCLayers(
            x_res=0.01, layer_widths=layers, layer_matrls=mtrls, repeats=1, temp=10.0
        )
        qcLayers.populate_x()
        qcLayers.solve_whole()
        validate_ode(qcLayers)
        self.assertAlmostEqual(
            qcLayers.eigen_es[1] - qcLayers.eigen_es[0], 0.252 - 0.102, 2
        )
        self.assertAlmostEqual(
            qcLayers.eigen_es[2] - qcLayers.eigen_es[0], 0.373 - 0.102, 1
        )
        self.assertAlmostEqual(abs(qcLayers.dipole(1, 0)), 16.4, 0)

    def test_three_well(self):
        layers = [300, 46, 10, 20, 10, 19, 300]
        mtrls = [1, 0, 1, 0, 1, 0, 1]
        qcLayers = QCLayers(
            x_res=0.01, layer_widths=layers, layer_matrls=mtrls, repeats=1, temp=10.0
        )
        qcLayers.populate_x()
        qcLayers.solve_whole()
        validate_ode(qcLayers)
        self.assertAlmostEqual(
            qcLayers.eigen_es[1] - qcLayers.eigen_es[0], 0.242 - 0.126, 2
        )
        self.assertAlmostEqual(
            qcLayers.eigen_es[2] - qcLayers.eigen_es[0], 0.383 - 0.126, 2
        )
        self.assertAlmostEqual(
            qcLayers.eigen_es[3] - qcLayers.eigen_es[0], 0.494 - 0.126, 2
        )
        self.assertAlmostEqual(abs(qcLayers.dipole(1, 0)), 18.6, -1)


if __name__ == "__main__":
    unittest.main()
