#!/usr/bin/env python3
# -*- coding:utf-8 -*-
from context import *  # type: ignore # noqa: F401, F403
import unittest
from ErwinJr2 import Material


class TestMaterial(unittest.TestCase):
    def test_GaAs(self):
        Material.main("GaAs")

    def test_AlGaAs(self):
        Material.main("AlGaAs")

    def test_Alloy_AlGaAs(self):
        algaas0 = Material.Alloy("AlGaAs", 0.0)  # so that is pure GaAs
        algaas1 = Material.Alloy("AlGaAs", 1.0)  # so that is pure GaAs
        algaas = Material.Alloy("AlGaAs", 0.33)  # so that is pure GaAs
        gaas = Material.Material("GaAs")
        alas = Material.Material("AlAs")
        for key in gaas.param:
            self.assertAlmostEqual(algaas0.param[key], gaas.param[key])
            self.assertAlmostEqual(algaas1.param[key], alas.param[key])
            self.assertTrue(
                (algaas0.param[key] <= algaas.param[key] <= algaas1.param[key])
                or
                (algaas1.param[key] <= algaas.param[key] <= algaas0.param[key])
                or
                # Equal with machine precision
                (abs(algaas1.param[key]-algaas.param[key]) < 1E-18
                 and abs(algaas0.param[key]-algaas.param[key]) < 1E-18)
            )

    def test_zero_strain(self):
        gaas = Material.Material("GaAs")
        gaas.set_strain(gaas.param["alc"])
        self.assertAlmostEqual(gaas.a_perp, gaas.param["alc"])


if __name__ == "__main__":
    unittest.main()
