#!/usr/bin/env python
# -*- coding:utf-8 -*-
import numpy as np
from numpy import sqrt, exp, pi
import scipy
from scipy.constants import (e as e0, epsilon_0 as eps0, h as h, 
                             hbar as hbar, electron_mass as m0, c as c0)
import OneDQuantum as onedq
import Material

qcMaterial = {
    "InP":  ["InGaAs", "AlInAs"], 
    "GaAs": ["AlGaAs", "AlGaAs"], 
    "GaSb": ["InAsSb", "AlGaSb"]
}

class QCLayers(object):
    """Class for QCLayers
    Member variables:
        parameters for each layer, np.array type, with len = No. of layers:
            layerWidths - width of each layer, float
            layerMaterialIdxs - label of materials, binary int
            layerDopings - 
            layerARs - if the layer is activ eor not, binary int
    """
    def __init__(self, substrate="InP", materials=["InGaAs", "AlInAs"], 
                 moleFracs=[0.53, 0.52], xres=0.5, Eres=0.5, 
                 layerWidths=[0.0], layerMaterialIdxs=[0], layerDopings=[0.0], 
                 layerARs=[True], EField=0, repeats=3, T=300.0, Solver="ODE", 
                 description=""):
        self.substrate = substrate
        self.materials = materials
        self.moleFracs = moleFracs
        self.xres = xres
        self.Eres = Eres
        self.layerWidths = layerWidths
        self.layerMaterialIdxs = layerMaterialIdxs
        self.layerDopings = layerDopings
        self.layerARs = layerARs
        self.EField = EField
        self.repeats = repeats
        self.Temperature = T
        self.Solver = Solver
        self.description = description

        self.subM = Material.Material(self.substrate, self.Temperature)

    def update_strain(self):
        self.layerMaterials = [Material.Alloy(self.materials[idx],
                                              self.moleFracs[idx],
                                              self.Temperature)
                               for idx in self.layerMaterialIdxs]
        self.a_parallel = self.subM.parm['alc']
        for alloy in self.layerMaterials:
            alloy.set_strain(self.a_parallel)

    def populate_x(self):
        self.update_strain()
        self.xPoints = np.arange(0, sum(self.layerWidths), self.xres)
        # largest index of grid
        N = self.xPoints.size
        self.xMaterialsIdxs = np.empty(N, dtype=int)
        self.xDopings = np.empty(N)
        self.xARs = np.empty(N, dtype=bool)
        self.xLayerNums = np.empty(N, dtype=int)

        self.xVc = np.empty(N)
        self.xEg = np.empty(N)

        layerNumCumSum = [0] + np.cumsum(self.layerWidths).tolist()
        for n in range(len(self.layerWidths)):
            Indices = np.logical_and(self.xPoints >= layerNumCumSum[n],
                                     self.xPoints < layerNumCumSum[n+1])
            
            self.xLayerNums[Indices] = n
            self.xMaterialsIdxs[Indices] = self.layerMaterialIdxs[n]
            self.xDopings[Indices] = self.layerDopings[n]
            self.xARs[Indices] = self.layerARs[n]

            for (p, key) in ((self.xVc, 'EcG'), (self.xEg, 'EgLH')):
                p[Indices] = self.layerMaterials[n].parm[key]

        self.xVc -= self.xPoints * self.EField

# vim: ts=4 sw=4 sts=4 expandtab
