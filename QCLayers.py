#!/usr/bin/env python
# -*- coding:utf-8 -*-
import numpy as np
from numpy import sqrt, exp, pi
import scipy
from scipy.constants import (e as e0, epsilon_0 as eps0, h as h, 
                             hbar as hbar, electron_mass as m0, c as c0)
import OneDQuantum as onedq
import Material
import copy

EUnit = 1e-5    # E field unit from kV/cm to V/Angtrom

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
        self.NonParabolic = True
        self.layerSelected = None

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
        layerNumCumSum = [0] + np.cumsum(self.layerWidths).tolist()
        periodL = layerNumCumSum[-1]
        self.xPoints = np.arange(0, periodL* self.repeats, self.xres)
        # largest index of grid
        N = self.xPoints.size
        self.xMaterialsIdxs = np.empty(N, dtype=int)
        self.xDopings = np.empty(N)
        self.xARs = np.empty(N, dtype=bool)
        self.xLayerNums = np.empty(N, dtype=int)

        self.xVc = np.empty(N)
        self.xVX = np.empty(N)
        self.xVL = np.empty(N)
        self.xVLH = np.empty(N)
        self.xVSO = np.empty(N)
        
        self.xEg = np.empty(N)
        self.xMc = np.empty(N)
        self.xESO = np.empty(N)
        self.xEp = np.empty(N)
        self.xF = np.empty(N)

        for n in range(len(self.layerWidths)):
            #  Indices = np.logical_and(self.xPoints >= layerNumCumSum[n],
            #                           self.xPoints < layerNumCumSum[n+1])
            Indices = np.logical_or.reduce([
                (self.xPoints >= layerNumCumSum[n] + k * periodL) 
                & (self.xPoints < layerNumCumSum[n+1] + k * periodL)
                for k in range(self.repeats)])
            
            self.xLayerNums[Indices] = n
            self.xMaterialsIdxs[Indices] = self.layerMaterialIdxs[n]
            self.xDopings[Indices] = self.layerDopings[n]
            self.xARs[Indices] = self.layerARs[n]

            for (p, key) in ((self.xVc, 'EcG'), (self.xEg, 'EgLH'),
                             (self.xVX, 'EcX'), (self.xMc, 'me0'),
                             (self.xVL, 'EcL'), (self.xESO, 'ESO'),
                             (self.xVLH, 'EvLH'), (self.xEp, 'Ep'),
                             (self.xVSO, 'EvSO'), (self.xF, 'F')):
                p[Indices] = self.layerMaterials[n].parm[key]

        ExtField = self.xPoints * self.EField * EUnit
        for p in (self.xVc, self.xVX, self.xVL, self.xVLH, self.xVSO):
            p -= ExtField

        #  if self.layerSelected:
        xSlt = (self.xLayerNums == self.layerSelected)
        xSlt = (xSlt | np.roll(xSlt, 1) | np.roll(xSlt, -1))
        self.xlayerSelected = np.ma.masked_where(~xSlt , self.xVc)

    def solve_whole(self):
        mass = self.xMc[np.argmin(self.xVc)]
        Emin = 2.33810741 * (hbar**2*(self.EField*EUnit)**2/(
            2*m0*mass*e0**2))**(1/3)
        #  print(Emin)
        Es = np.linspace(np.min(self.xVc), np.max(self.xVc)+Emin, 1000)
        #  Es = np.linspace(-1.35, -0.95, 1000)
        if self.NonParabolic:
            band = onedq.Band("ZincBlende", self.xEg, self.xF, self.xEp,
                              self.xESO)
            self.eigenEs = onedq.cBandSolve1D(self.xres, Es, self.xVc, band)
            self.psis = onedq.cBandFillPsi(self.xres, self.eigenEs, 
                                           self.xVc, band)
        else:
            self.eigenEs = onedq.cSimpleSolve1D(self.xres, Es,
                                                self.xVc, self.xMc)
            self.psis = onedq.cSimpleFillPsi(self.xres, self.eigenEs,
                                             self.xVc, self.xMc)

    def solve_basis(self):
        IndSep = np.nonzero(np.array(self.layerARs[1:]) !=
                            np.array(self.layerARs[:-1]))[0] + 1
        IndSep = np.insert(IndSep, len(IndSep), len(self.layerARs))
        if self.layerARs[0]:
            IndSep = np.insert(IndSep, 0, 0)

        dCL = []
        for n in range(0, int(len(IndSep)/2)):
            dCL.append(copy.deepcopy(self))
            dCL[n].repeats = 1
            
            dCL[n].layerWidths = self.layerWidths[IndSep[2*n]:IndSep[2*n+1]]
            print("n = {0}, dCL[n].layerWidths = {1}".format(n, dCL[n].layerWidths))
            
            dCL[n].layerMaterialIdxs = self.layerMaterialIdxs[IndSep[2*n]:
                                                              IndSep[2*n+1]]
            dCL[n].layerDopings = self.layerDopings[IndSep[2*n]:IndSep[2*n+1]]
            dCL[n].layerARs = self.layerARs[IndSep[2*n]:IndSep[2*n+1]]

            dCL[n].populate_x()
            dCL[n].solve_whole()
            
            print("dCL[n].eigenEs.shape = {0}".format(dCL[n].eigenEs.shape))
            print("dCL[n].psis.shape = {0}".format(dCL[n].psis.shape))

        
        
# vim: ts=4 sw=4 sts=4 expandtab
