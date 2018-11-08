#!/usr/bin/env python
# -*- coding:utf-8 -*-
import numpy as np
from numpy import sqrt, exp, pi
import scipy
from scipy.constants import (e as e0, epsilon_0 as eps0, h as h, 
                             hbar as hbar, electron_mass as m0, c as c0)
import OneDQuantum as onedq

qcMaterial = {
    "InP":  ["InGaAs", "AlInAs"], 
    "GaAs": ["AlGaAs", "AlGaAs"], 
    "GaSb": ["InAsSb", "AlGaSb"]
}

class QCLayers(object):
    """Class for QCLayers"""
    def __init__(self, substrate="InP", materials=["InGaAs", "AlInAs"], 
                 moleFracs=[0.53, 0.52], xres=0.5, Eres=0.5, 
                 layerWidths=[0.0], layerMaterials=[0], layerDopings=[0.0], 
                 layerARs=[True], EField=0, repeats=3, T=300.0, Solver="ODE", 
                 description=""):
        self.substrate = substrate
        self.materials = materials
        self.moleFracs = moleFracs
        self.xres = xres
        self.Eres = Eres
        self.layerWidths = layerWidths
        self.layerMaterials = layerMaterials
        self.layerDopings = layerDopings
        self.layerARs = layerARs
        self.EField = EField
        self.repeats = repeats
        self.Temperature = T
        self.Solver = Solver
        self.description = description

# vim: ts=4 sw=4 sts=4 expandtab
