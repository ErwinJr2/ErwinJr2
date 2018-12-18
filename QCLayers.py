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
    "GaAs": ["AlGaAs"], 
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
        self.update_strain()

    def update_strain(self): 
        self.a_parallel = self.subM.parm['alc']
        self.mtrlAlloys = [Material.Alloy(self.materials[idx], 
                                          self.moleFracs[idx], 
                                          self.Temperature)
                           for idx in range(len(self.materials))]

        for al in self.mtrlAlloys:
            al.set_strain(self.a_parallel)

    def set_mtrl(self, n, mtrl=None, moleFrac=None):
        """Set material[n] to new material (mtrl) and/or moleFrac"""
        if not mtrl:
            self.mtrlAlloys[n].set_molefrac(moleFrac)
        if not moleFrac:
            moleFrac = self.moleFracs[n]
            self.mtrlAlloys[n] = Material.Alloy(
                mtrl, moleFrac, self.Temperature)
        self.mtrlAlloys[n].set_strain(self.a_parallel)

    def add_layer(self, n, width, mtrlIdx, AR, doping):
        self.layerWidths.insert(n, width)
        self.layerMaterialIdxs.insert(n, mtrlIdx)
        self.layerARs.insert(n, AR)
        self.layerDopings.insert(n, doping)
    
    def del_layer(self, n):
        for layerList in (self.layerWidths, self.layerARs,
                          self.layerMaterialIdxs, self.layerDopings):
            layerList.pop(n)

    def rotate_layer(self):
        for layerList in (self.layerWidths, self.layerARs,
                          self.layerMaterialIdxs, self.layerDopings):
            layerList.insert(0, layerList.pop())

    def set_substrate(self, subs):
        if subs in qcMaterial:
            self.substrate = subs
            self.update_strain()
        else:
            raise TypeError("Substrate %s not supported"%subs)

    def set_temperature(self, T):
        self.qclayers.Temperature = T
        self.subM.set_temperature(T)
        self.update_strain()

    def offset(self):
        """Return the conduction band offset (difference between highest
        conduction band and lowest conduction band energy) of materials, 
        in unit eV"""
        ecgs = [alloy.parm['EcG'] for alloy in self.mtrlAlloys]
        return max(ecgs) - min(ecgs)

    def netStrain(self):
        """Return average strain perpendicular to the layer plane, in
        percentage."""
        totalStrain = sum(self.mtrlAlloys[self.layerMaterialIdxs[n]].eps_perp
                          * self.layerWidths[n] 
                          for n in range(len(self.layerWidths)))
        return 100 * totalStrain / sum(self.layerWidths)

    def avghwLO(self):
        """Return average LO phonont energy in unit eV"""
        sumhwlo = sum(self.mtrlAlloys[self.layerMaterialIdxs[n]].parm['hwLO']
                          * self.layerWidths[n] 
                          for n in range(len(self.layerWidths)))
        return sumhwlo / sum(self.layerWidths)

    def populate_x(self):
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
                p[Indices] = self.mtrlAlloys[self.layerMaterialIdxs[n]].parm[key]

        ExtField = self.xPoints * self.EField * EUnit
        for p in (self.xVc, self.xVX, self.xVL, self.xVLH, self.xVSO):
            p -= ExtField

        self.xRepeats = np.empty((0), dtype=int)
        for k in range(self.repeats):
            self.xRepeats = np.concatenate(
                (self.xRepeats, k*np.ones(int(N/self.repeats), dtype=int)))

    def xLayerSelected(self):
        xSlt = (self.xLayerNums == self.layerSelected)
        xSlt = (xSlt | np.roll(xSlt, 1) | np.roll(xSlt, -1))
        return np.ma.masked_where(~xSlt, self.xVc)

    def xLayerAR(self):
        return np.ma.masked_where(~self.xARs, self.xVc)

    def solve_whole(self):
        mass = self.xMc[np.argmin(self.xVc)]
        Emin = 2.33810741 * (hbar**2*(self.EField*EUnit)**2/(
            2*m0*mass*e0**2))**(1/3)
        #  print(Emin)
        Es = np.linspace(np.min(self.xVc)+Emin, np.max(self.xVc), 1000)
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
        if self.layerARs[0]:
            IndSep = np.concatenate(([0], IndSep))
        if len(IndSep) % 2:
            IndSep = np.concatenate((IndSep, len(self.layerARs)))

        StartInd = IndSep[0::2]
        EndInd = IndSep[1::2]
        # create two lists to store start ind and end ind

        self.eigenEs = np.empty((0))
        self.psis = np.empty((0, self.xPoints.size))
        for n in range(0, len(StartInd)):
            dCL = copy.deepcopy(self)
            dCL.repeats = 1
            dCL.layerWidths = self.layerWidths[StartInd[n]:EndInd[n]]
            dCL.layerMaterialIdxs = self.layerMaterialIdxs[StartInd[n]:
                                                           EndInd[n]]
            dCL.layerDopings = self.layerDopings[StartInd[n]:EndInd[n]]
            dCL.layerARs = self.layerARs[StartInd[n]:EndInd[n]]

            dCL.populate_x()
            dCL.solve_whole()

            xInd_all = (self.xLayerNums >= StartInd[n]) & \
                (self.xLayerNums < EndInd[n])

            psis_re = np.zeros((0, self.xPoints.size))
            eigenEs_re = np.zeros(0)
            for j in range(0, self.repeats):
                xInd = xInd_all & (self.xRepeats == j)
                psis_fill = np.zeros((dCL.eigenEs.shape[0], self.xPoints.size))
                psis_fill[:, xInd] = dCL.psis
                psis_re = np.concatenate((psis_re, psis_fill), axis=0)

                eigenEs_re = \
                    np.concatenate((eigenEs_re,
                                    dCL.eigenEs -
                                    j*sum(self.layerWidths)*self.EField*EUnit))

            self.eigenEs = np.concatenate((self.eigenEs, eigenEs_re), axis=0)
            self.psis = np.concatenate((self.psis, psis_re), axis=0)

# vim: ts=4 sw=4 sts=4 expandtab
