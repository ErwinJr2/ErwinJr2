#!/usr/bin/env python
# -*- coding:utf-8 -*-
from context import *
import unittest
import Material
class TestMaterial(unittest.TestCase):
    def test_GaAs(self):
        Material.main("GaAs")

    def test_AlGaAs(self):
        Material.main("AlGaAs")

    def test_Alloy_AlGaAs(self):
        algaas = Material.Alloy("AlGaAs", 0.0) # so that is pure GaAs
        gaas = Material.Material("GaAs")
        #  self.assertEqual(gaas.parm.pop("Crystal"), "ZincBlende")
        for key in gaas.parm:
            self.assertAlmostEqual(algaas.parm[key], gaas.parm[key])

    def test_zero_strain(self):
        gaas = Material.Material("GaAs")
        gaas.set_strain(gaas.parm["alc"])
        self.assertAlmostEqual(gaas.a_perp, gaas.parm["alc"])

if __name__ == "__main__":
    unittest.main()
# vim: ts=4 sw=4 sts=4 expandtab
