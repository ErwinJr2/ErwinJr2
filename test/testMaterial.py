#!/usr/bin/env python3
# -*- coding:utf-8 -*-
from context import *  # type: ignore # noqa: F401, F403
import unittest
import Material


class TestMaterial(unittest.TestCase):
    def test_GaAs(self):
        Material.main("GaAs")

    def test_AlGaAs(self):
        Material.main("AlGaAs")

    def test_Alloy_AlGaAs(self):
        algaas = Material.Alloy("AlGaAs", 0.0)  # so that is pure GaAs
        gaas = Material.Material("GaAs")
        #  self.assertEqual(gaas.param.pop("Crystal"), "ZincBlende")
        for key in gaas.param:
            self.assertAlmostEqual(algaas.param[key], gaas.param[key])

    def test_zero_strain(self):
        gaas = Material.Material("GaAs")
        gaas.set_strain(gaas.param["alc"])
        self.assertAlmostEqual(gaas.a_perp, gaas.param["alc"])


if __name__ == "__main__":
    unittest.main()
# vim: ts=4 sw=4 sts=4 expandtab
