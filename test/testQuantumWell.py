#!/usr/bin/env python3
# -*- coding:utf-8 -*-
from context import *  # type: ignore # noqa: F401, F403
from ErwinJr2.QCLayers import SchrodingerLayer, QCLayers
import numpy as np
from scipy.constants import e as e0, hbar as hbar, electron_mass as m0
import unittest

"""This unit test is to compare simulation with PhysRevB.50.8663"""


class GaInAs_AlInAs_Layer(SchrodingerLayer):
    """Sadly the parameters given in PhysRevB.50.8663 is not complete"""
    def __init__(self, xres, layerWidths, layerMtrl, EField=0):
        # fixed offset 0.51 eV
        Vc = {0: 0, 1: 0.51}
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


def validate_ODE(q):
    interface = (np.abs(np.diff(q.xVc)) > 0.1)
    interface = (interface[1:] | interface[:-1])
    residules = []
    for t in range(len(q.eigenEs)):
        E = q.eigenEs[t]
        psi = q.psis[t]
        xEg, xF, xEp, xESO = q.bandParams
        meff = 1/(1 + 2*xF + xEp/3*(1/(E-q.xVc+xEg+xESO) + 2/(E-q.xVc+xEg)))
        # interpolate half grid
        meff = (meff[1:] + meff[:-1])/2
        residule = (
            -np.diff(hbar**2/(2*m0*meff*e0)*np.diff(psi))/(q.xres*1E-10)**2
            + q.xVc[1:-1] * psi[1:-1] - E*psi[1:-1])
        residules.append(residule)
        interRes = residule[interface]
        # The interface has larger error
        np.testing.assert_almost_equal(interRes, np.zeros(interRes.shape),
                                       decimal=2)
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
        validate_ODE(singleWell)
        Ei = singleWell.eigenEs[1]
        Ej = singleWell.eigenEs[0]
        self.assertAlmostEqual((Ei-Ej) / (0.381-0.123), 1, 1)
        self.assertAlmostEqual(abs(singleWell.dipole(1, 0)) / 15.3, 1, 0)

        # Test QCLayers
        qcLayers = QCLayers(xres=0.01, layerWidths=layers,
                            layerMtrls=mtrls, repeats=1, T=10.0)
        qcLayers.populate_x()
        qcLayers.solve_whole()
        validate_ODE(qcLayers)
        Ei = qcLayers.eigenEs[1]
        Ej = qcLayers.eigenEs[0]
        self.assertAlmostEqual((Ei-Ej) / (0.381-0.123), 1, 1)
        self.assertAlmostEqual(abs(qcLayers.dipole(1, 0)) / 15.3, 1, 0)

    def test_double_well(self):
        layers = [200, 59, 13, 24, 200]
        mtrls = [1, 0, 1, 0, 1]
        qcLayers = QCLayers(xres=0.01, layerWidths=layers,
                            layerMtrls=mtrls, repeats=1, T=10.0)
        qcLayers.populate_x()
        qcLayers.solve_whole()
        validate_ODE(qcLayers)
        self.assertAlmostEqual(qcLayers.eigenEs[1] - qcLayers.eigenEs[0],
                               0.252 - 0.102, 2)
        self.assertAlmostEqual(qcLayers.eigenEs[2] - qcLayers.eigenEs[0],
                               0.373 - 0.102, 1)
        self.assertAlmostEqual(abs(qcLayers.dipole(1, 0)), 16.4, 0)

    def test_three_well(self):
        layers = [300, 46, 10, 20, 10, 19, 300]
        mtrls = [1, 0, 1, 0, 1, 0, 1]
        qcLayers = QCLayers(xres=0.01, layerWidths=layers,
                            layerMtrls=mtrls, repeats=1, T=10.0)
        qcLayers.populate_x()
        qcLayers.solve_whole()
        validate_ODE(qcLayers)
        self.assertAlmostEqual(qcLayers.eigenEs[1] - qcLayers.eigenEs[0],
                               0.242 - 0.126, 2)
        self.assertAlmostEqual(qcLayers.eigenEs[2] - qcLayers.eigenEs[0],
                               0.383 - 0.126, 2)
        self.assertAlmostEqual(qcLayers.eigenEs[3] - qcLayers.eigenEs[0],
                               0.494 - 0.126, 2)
        self.assertAlmostEqual(abs(qcLayers.dipole(1, 0)), 18.6, -1)


if __name__ == "__main__":
    unittest.main()
