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
    """
    Class for QCLayers

    Parameters
    ----------
    substrate : str
        The substrate material for the device, which determines the well and
        barrier material

        =========    ==================================   ==================================
        substrate             well                        barrier
        =========    ==================================   ==================================
        InP          In :math:`_x` Ga :math:`_{1-x}` As   Al :math:`_{1-x}` In :math:`_x` As
        GaAs         Al :math:`_x` Ga :math:`_{1-x}` As   Al :math:`_x` Ga :math:`_{1-x}` As
        GaSb         InAs :math:`_y` Sb :math:`_{1-y}`    Al :math:`_x` Ga :math:`_{1-x}` Sb
        =========    ==================================   ==================================

    materials : list of str, len = 2
        Name of alloys
    moleFrac : list of float
        mole fraction for each possible layer material, 
        in format [well, barrier]*4
    xres : float
        Position resolution, in armstrong
    Eres : float
        External (static) electrical field, in kV/cm = 1e5 V/m

    layerWidths : np.array of float, len = No. of layers
        Width of each layer, in mm
    layerMtrls : np.array of binary int, len = No. of layers
        Label of materials, depending on substrate
    layerDopings : np.array of float, len = No. of layers
        Doping per volumn in unit 1e17 cm-3
    layerARs : np.array of binary, len = No. of layers
        Binaries indicating if the layer is active(True) or not(False), 
        only affects basis solver

    EField : float
        External (static) electrical field, in kV/cm = 1e5 V/m
    repeats : int
        Number of repeat times for the given structure
    T : float
        Temperature of the device, affecting material property
    Solver : str
        Name of solver
    description : str
        Description of the data

    
    Other Parameters
    --------------------
    NonParabolic : bool
        TBD
    layerSelected : int
        a label indicating which layer is selected in GUI, with default None 
        indicating no layer is selected
        
    """
    def __init__(self, substrate="InP", materials=["InGaAs", "AlInAs"], 
                 moleFracs=[0.53, 0.52], xres=0.5, Eres=0.5, 
                 layerWidths=[0.0], layerMtrls=[0], layerDopings=[0.0], 
                 layerARs=[True], EField=0, repeats=3, T=300.0, Solver="ODE", 
                 description=""):
        self.substrate = substrate
        self.materials = materials
        self.moleFracs = moleFracs
        self.xres = xres
        self.Eres = Eres
        self.layerWidths = layerWidths
        self.layerMtrls = layerMtrls
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
        self.layerMtrls.insert(n, mtrlIdx)
        self.layerARs.insert(n, AR)
        self.layerDopings.insert(n, doping)
    
    def del_layer(self, n):
        for layerList in (self.layerWidths, self.layerARs,
                          self.layerMtrls, self.layerDopings):
            layerList.pop(n)

    def rotate_layer(self):
        for layerList in (self.layerWidths, self.layerARs,
                          self.layerMtrls, self.layerDopings):
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
        if sum(self.layerWidths) <= 1e-5: 
            return -1
        totalStrain = sum(self.mtrlAlloys[self.layerMtrls[n]].eps_perp
                          * self.layerWidths[n] 
                          for n in range(len(self.layerWidths)))
        return 100 * totalStrain / sum(self.layerWidths)

    def avghwLO(self):
        """Return average LO phonont energy in unit eV"""
        if sum(self.layerWidths) <= 1e-5: 
            return -1
        sumhwlo = sum(self.mtrlAlloys[self.layerMtrls[n]].parm['hwLO']
                          * self.layerWidths[n] 
                          for n in range(len(self.layerWidths)))
        return sumhwlo / sum(self.layerWidths)

    def populate_x(self):
        """Calculate the properties in terms of position

        Yield
        -----
        xPoints : np.array of float
            position grid
        xMaterialsIdxs : np.array of binary int
            label of materials at each position
        xDopings : np.array of float
            Doping per volumn at each position
        xARs : np.array of bool
            Binaries indicating if the layer is active(True) or not(False) at
            each position
        xLayerNums : np.array of int
            at xPoints[i] it's xLayerNums[i]-th layer


        Also update following band structure parameters (with type *np.array of float*):
        **xVc, xVX, xVL, xVLH, xVSO, xEg, xMc, xESP, xF**

        """
        layerNumCumSum = [0] + np.cumsum(self.layerWidths).tolist()
        periodL = layerNumCumSum[-1]
        self.xPoints = np.arange(0, periodL* self.repeats, self.xres)
        """ np.array """
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
            self.xMaterialsIdxs[Indices] = self.layerMtrls[n]
            self.xDopings[Indices] = self.layerDopings[n]
            self.xARs[Indices] = self.layerARs[n]

            for (p, key) in ((self.xVc, 'EcG'), (self.xEg, 'EgLH'),
                             (self.xVX, 'EcX'), (self.xMc, 'me0'),
                             (self.xVL, 'EcL'), (self.xESO, 'ESO'),
                             (self.xVLH, 'EvLH'), (self.xEp, 'Ep'),
                             (self.xVSO, 'EvSO'), (self.xF, 'F')):
                p[Indices] = self.mtrlAlloys[self.layerMtrls[n]].parm[key]

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
        """
        solve eigen modes for the whole structure

        Yield
        -----
        eigenEs : np.array of float
            the eigenenergy of the layer structure
        psis : np.array of float
            the wave function

        """
        mass = self.xMc[np.argmin(self.xVc)]
        # ground state for triangular well
        Emin = 2.33810741 * (hbar**2*(self.EField*EUnit)**2/(
            2*m0*mass*e0**2))**(1/3)
        Es = np.linspace(np.min(self.xVc)+Emin, np.max(self.xVc), 1000)
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
        #  self.stateFilter()

    def solve_basis(self):
        """
        solve basis for the QC device, with each basis being the eigen mode of 
        a seperate part of the layer structure

        Yield
        -----
        eigenEs : np.array of float
            the eigenenergy of the layer structure
        psis : np.array of float
            the wave function
        
        """
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
            dCL.layerMtrls = self.layerMtrls[StartInd[n]:
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

    def stateFilter(self):
        """Filter unbounded states: 
        States with energy higher than potential at the end and wf[-2] higher
        than a threshold (1e-4) is considered unbounded. """
        bounded = []
        semiBounded = []
        for n in range(self.eigenEs.size):
            #  E = self.eigenEs[n]
            #  k = sqrt((E-self.xVc[-1])*2*self.xMc[-1])/hbar
            wf = self.psis[n,:]
            if self.eigenEs[n] < self.xVc[-1] or abs(wf[-2]) < 1e-4:
                bounded.append(n)
            elif np.all(self.xVc[wf**2 > 1e-4 ] > self.eigenEs[n]):
                semiBounded.append(n)
        ss = sorted(bounded + semiBounded)
        self.eigenEs = self.eigenEs[ss]
        self.psis = self.psis[ss,:]

    def dipole(self, upper, lower):
        """Return Electrical dipole between upper and lower states, 
        update self.dipole. Should be called for any other related physics
        quantities."""
        if upper < lower:
            upper, lower = lower, upper
        psi_i = self.psis[upper, :]
        psi_j = self.psis[lower, :]
        Ei = self.eigenEs[upper]
        Ej = self.eigenEs[lower]
        #TODO: eff mass for non-parabolic
        avgpsi_i = (psi_i[:-1] + psi_i[1:])/2
        avgxMc = (self.xMc[:-1] + self.xMc[1:])/2
        self.z = np.sum(avgpsi_i * np.diff(psi_j/self.xMc) 
                   + 1/avgxMc * (avgpsi_i * np.diff(psi_j)))
        self.z *= hbar**2 / (2 * (Ei - Ej) * e0 * m0) / (1E-10)**2 #Angtrom
        return self.z

    def loTransition(self, upper, lower):
        INV_INF = 1e-20  # for infinit small decay rate (ns-1)
        if upper < lower:
            upper, lower = lower, upper

        psi_i = self.psis[upper, :]
        psi_j = self.psis[lower, :]
        Ei = self.eigenEs[upper]
        Ej = self.eigenEs[lower]
        hwLO = self.avghwLO()

        if Ei - Ej - hwLO < 0:
            # energy difference is smaller than a LO phonon
            # LO phonon scatering doesn't happen
            return INV_INF

        # TODO: imporve expression for eff mass
        mass = m0 * sqrt(np.sum(self.xMc * psi_i**2 * self.xres)
                         * np.sum(self.xMc * psi_j**2 * self.xres))
        kl = sqrt(2 * mass / hbar**2 * (Ei-Ej-hwLO) * e0)
        dIij = np.empty(self.xPoints.size)
        for n in range(self.xPoints.size):
            x1 = self.xPoints[n]
            x2 = self.xPoints
            dIij[n] = np.sum(psi_i * psi_j * exp(-kl * abs(x1 - x2)*1e-10)
                             * psi_i[n] * psi_j[n] * self.xres**2)
        Iij = np.sum(dIij)
        epsInf = np.array([a.parm["epsInf"] for a in self.mtrlAlloys])
        epss = np.array([a.parm["epss"] for a in self.mtrlAlloys])
        epsrho = 1 / (1/epsInf - 1/epss)
        epsrho = (np.sum(epsrho[self.layerMtrls] * self.layerWidths) 
                  / sum(self.layerWidths))
        inv_tau = (mass * e0**2 * self.avghwLO() * e0 / hbar * Iij 
                   / (4 * hbar**2 * epsrho * eps0 * kl))
        return inv_tau / 1e12 #unit ps

    def loLifeTime(self, state):
        """ return the life time due to LO phonon scattering of the
        given state(label)"""
        rate = [self.loTransition(state, q) for q in range(state)]
        return 1/sum(rate)

    def FoM(self, upper, lower):
        tauLower = self.loLifeTime(lower)
        tauUpper = self.loLifeTime(upper)
        tauUpperLower = 1/self.loTransition(upper, lower)
        self.FoM = self.z**2 * tauUpper * (1 - tauLower / tauUpperLower)
        return self.FoM

    def coupleBroadening(self, upper, lower):
        """Boradening of nergy difference between upper state and lower state
        due to coupling to other states"""
        #TODO
        return 0

    def ifrBroadening(self, upper, lower):
        """Interface roughness induced broadening"""
        #TODO
        return 0

    def alphaISB(self, upper, lower):
        #TODO
        return 0

# vim: ts=4 sw=4 sts=4 expandtab
