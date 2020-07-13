"""
This file defines the QCLayer class for simulating QC structure
"""
import numpy as np
from numpy import sqrt
from scipy.constants import (e as e0, epsilon_0 as eps0, h as h,
                             hbar as hbar, electron_mass as m0, c as c0)
import scipy.linalg as slg
import OneDQuantum as onedq
import Material
from OptStrata import rIdx
import copy

# onedq.OneDSchrodinger.bindOpenMP(rk4=True)

EUnit = 1e-5    # E field unit from kV/cm to V/Angtrom
# upper and lower bond extension for periodic solver
PeriodU = 1
PeriodL = 0.2

qcMaterial = {
    "InP":  ["InGaAs", "AlInAs"],
    "GaAs": ["AlGaAs"],
    "GaSb": ["InAsSb", "AlGaSb"]
}


class QCLayers(object):
    r"""Class for Quantum Cascade Layers

Parameters
----------
substrate : str
    The substrate material for the device, which determines the well and
    barrier material

    ========= ================================ ================================
    substrate              well                         barrier
    ========= ================================ ================================
    InP       In\ :sub:`x`\ Ga\ :sub:`1-x`\ As Al\ :sub:`1-x`\ In\ :sub:`x`\ As
    GaAs      Al\ :sub:`x`\ Ga\ :sub:`1-x`\ As Al\ :sub:`x`\ Ga\ :sub:`1-x`\ As
    GaSb      InAs\ :sub:`y`\ Sb\ :sub:`1-y`   Al\ :sub:`x`\ Ga\ :sub:`1-x`\ Sb
    ========= ================================ ================================

materials : list of str, len = 2
    Name of alloys
moleFrac : list of float
    mole fraction for each possible layer material,
    in format [well, barrier]*4
xres : float
    Position resolution, in armstrong
Eres : float
    External (static) electrical field, in kV/cm = 1e5 V/m
wl : float
    The wavelength for the design, in unit um, for book keeping and
    optimization, but doesn't go into quantum solver

layerWidths : list of float, len = No. of layers
    Width of each layer, in mm
layerMtrls : list of int, len = No. of layers
    Label of materials, depending on substrate
layerDopings : list of float, len = No. of layers
    Doping per volumn in unit 1e17 cm-3
layerARs : list of bool, len = No. of layers
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
    A flag to decide whether to use a constant effective mass, assuming
    perfect parabolic dispersion or an energy dependent effective mass
    to correct non-parabolic dispersion relation of electrons.
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
        self.crystalType = Material.MParm[substrate]["Crystal"]
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
        self.periodic = False
        self.solveMatrix = False
        self.layerSelected = None

        self.subM = Material.Material(self.substrate, self.Temperature)
        self.update_strain()

        self.basisARInjector = True
        self.basisInjectorAR = True
        self.basisARonly = False

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
        if mtrl is None and moleFrac is None:
            raise Exception("Nothing changed")
        if mtrl is None:
            mtrl = self.materials[n]
        if moleFrac is None:
            moleFrac = self.moleFracs[n]
        self.moleFracs[n] = moleFrac
        self.materials[n] = mtrl
        self.mtrlAlloys[n] = Material.Alloy(mtrl, moleFrac, self.Temperature)
        self.mtrlAlloys[n].set_strain(self.a_parallel)

    def add_mtrl(self, mtrl=None, moleFrac=None):
        """Add a new material possiblities"""
        self.materials.append(mtrl if mtrl else
                              qcMaterial[self.substrate][0])
        self.moleFracs.append(moleFrac if moleFrac else 0.0)

    def del_mtrl(self, n):
        """Delete materials labeled n. All layers of this material will
        become previous materials[n-1 if n >0 else 1].  There should be
        at least two materials otherwise there will be error"""
        if len(self.materials) <= 2:
            raise ValueError("There should be at least 2 materials")
        self.materials.pop(n)
        self.moleFracs.pop(n)
        for i in range(len(self.layerMtrls)):
            if self.layerMtrls[i] >= n:
                self.layerMtrls[i] = (self.layerMtrls[i] - 1
                                      if self.layerMtrls[i] > 0
                                      else 0)

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
            self.crystalType = Material.MParm[subs]["Crystal"]
            matlN = len(self.materials)
            self.materials = (qcMaterial[subs]*matlN)[0:matlN]
            self.update_strain()
        else:
            raise TypeError("Substrate %s not supported" % subs)

    def set_temperature(self, T):
        self.Temperature = T
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


        Also update following band structure parameters (with type *np.array
        of float*): **xVc, xVX, xVL, xVLH, xVSO, xEg, xMc, xESP, xF**

        """
        layerCumSum = [0] + np.cumsum(self.layerWidths).tolist()
        periodL = layerCumSum[-1]
        self.xPoints = np.arange(0, periodL * self.repeats, self.xres)
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
            Indices = np.logical_or.reduce([
                (self.xPoints >= layerCumSum[n] + k * periodL)
                & (self.xPoints < layerCumSum[n+1] + k * periodL)
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

        # populate mass half grid self.xMc[i] is m at i+-0.5
        # for matrix solver only
        self.xMcplus = np.empty(self.xPoints.size)
        self.xMcminus = np.empty(self.xPoints.size)
        for n in range(len(self.layerWidths)):
            Indices = np.logical_or.reduce([
                (self.xPoints >= layerCumSum[n] + k * periodL - self.xres/2)
                & (self.xPoints < layerCumSum[n+1] + k * periodL - self.xres/2)
                for k in range(self.repeats)])
            self.xMcplus[Indices] = self.mtrlAlloys[
                                        self.layerMtrls[n]].parm['me0']
        if N > 1:
            self.xMcplus[-1] = self.xMcplus[-2]
            self.xMcminus[1:] = self.xMcplus[:-1]
            self.xMcminus[0] = self.xMcminus[1]

        ExtField = self.xPoints * self.EField * EUnit
        for p in (self.xVc, self.xVX, self.xVL, self.xVLH, self.xVSO):
            p -= ExtField

        self.xRepeats = np.empty(N, dtype=int)
        for k in range(self.repeats):
            self.xRepeats[(self.xPoints < (k+1)*periodL)
                          & (self.xPoints >= k*periodL)] = k
        #  self.xRepeats = np.empty((0), dtype=int)
        #  for k in range(self.repeats):
        #      self.xRepeats = np.concatenate(
        #          (self.xRepeats, k*np.ones(int(N/self.repeats), dtype=int)))

    def xLayerSelected(self):
        xSlt = (self.xLayerNums == self.layerSelected)
        xSlt = (xSlt | np.roll(xSlt, 1) | np.roll(xSlt, -1))
        return np.ma.masked_where(~xSlt, self.xVc)

    def xLayerAR(self):
        return np.ma.masked_where(~self.xARs, self.xVc)

    def solve_whole(self):
        if self.solveMatrix:
            return self.solve_whole_matrix()
        else:
            return self.solve_whole_ode()

    def solve_whole_ode(self):
        """
        solve eigen modes for the whole structure

        Yield
        -----
        eigenEs : np.array of float
            the eigenenergy of the layer structure
        psis : np.array of float
            the wave function

        """
        # mass = self.xMc[np.argmin(self.xVc)]
        # # ground state for triangular well
        # Emin = 2.33810741 * (hbar**2*(self.EField*EUnit)**2/(
        #     2*m0*mass*e0**2))**(1/3)
        # Es = np.linspace(np.min(self.xVc)+Emin, np.max(self.xVc), 1024)
        Es = np.arange(np.min(self.xVc), np.max(self.xVc), self.Eres/1E3)
        if self.periodic:
            band = onedq.Band(self.crystalType, self.xEg, self.xF, self.xEp,
                              self.xESO)
            offset = self.offset()
            self.Emin = (min(m.parm['EcG'] for m in self.mtrlAlloys)
                         - PeriodL*offset)
            self.Emax = (max(m.parm['EcG'] for m in self.mtrlAlloys)
                         + PeriodU*offset)
            Eshift = self.EField*EUnit*sum(self.layerWidths)
            Es = np.arange(self.Emin - Eshift, self.Emin, self.Eres/1E3)
            eigenEs = onedq.cBandSolve1DBonded(
                self.xres, Es, self.Emin, self.Emax,
                self.EField, self.xVc, band)
            psis = onedq.cBandFillPsi(
                self.xres, eigenEs, self.xVc, band,
                Elower=self.Emin, Eupper=self.Emax, field=self.EField)
            self.psis, self.eigenEs = self.shiftPeriod(
                                          (-1, 0, 1, 2), psis, eigenEs)
        elif self.NonParabolic:
            band = onedq.Band(self.crystalType, self.xEg, self.xF, self.xEp,
                              self.xESO)
            self.eigenEs = onedq.cBandSolve1D(self.xres, Es, self.xVc, band)
            self.psis = onedq.cBandFillPsi(self.xres, self.eigenEs,
                                           self.xVc, band)
        else:
            self.eigenEs = onedq.cSimpleSolve1D(self.xres, Es,
                                                self.xVc, self.xMc)
            self.psis = onedq.cSimpleFillPsi(self.xres, self.eigenEs,
                                             self.xVc, self.xMc)
        self.loMatrix = [[None]*len(self.eigenEs)
                         for _ in range(len(self.eigenEs))]
        return self.eigenEs
        #  self.stateFilter()

    def solve_whole_matrix(self):
        """
        solve eigen modes for the whole structure with matrix eigen-solver

        Yield
        -----
        eigenEs : np.array of float
            the eigenenergy of the layer structure
        psis : np.array of float
            the wave function

        """
        # populate mass half grid self.xMc[i] is m at i+-0.5
        # for matrix solver only
        layerCumSum = [0] + np.cumsum(self.layerWidths).tolist()
        periodL = layerCumSum[-1]
        self.xMcplus = np.empty(self.xPoints.size)
        self.xMcminus = np.empty(self.xPoints.size)
        for n in range(len(self.layerWidths)):
            Indices = np.logical_or.reduce([
                (self.xPoints >= layerCumSum[n] + k * periodL - self.xres/2)
                & (self.xPoints < layerCumSum[n+1] + k * periodL - self.xres/2)
                for k in range(self.repeats)])
            self.xMcplus[Indices] = self.mtrlAlloys[
                                        self.layerMtrls[n]].parm['me0']
        self.xMcplus[-1] = self.xMcplus[-2]
        self.xMcminus[1:] = self.xMcplus[:-1]
        self.xMcminus[0] = self.xMcminus[1]

        # unit eV/step^2
        unit = hbar**2/(2*e0*m0*(1E-10*self.xres)**2)
        # diagonal and subdiagonal of Hamiltonian
        self.Hdiag = unit*(1/self.xMcplus + 1/self.xMcminus) + self.xVc
        self.Hsubd = -unit / self.xMcplus[:-1]
        self.eigenEs, self.psis = slg.eigh_tridiagonal(
            self.Hdiag, self.Hsubd, select='v',
            select_range=(np.min(self.xVc), np.max(self.xVc)))
        self.psis /= sqrt(self.xres)
        self.psis = self.psis.T
        self.loMatrix = [[None]*len(self.eigenEs)
                         for _ in range(len(self.eigenEs))]
        return self.eigenEs

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
        StartInd = []
        EndInd = []
        if self.basisARonly:
            if self.layerARs[0]:
                StartInd.append(0)
            for n in range(1, len(self.layerARs)):
                if not self.layerARs[n-1] and self.layerARs[n]:
                    StartInd.append(n)
                if self.layerARs[n-1] and not self.layerARs[n]:
                    EndInd.append(n)
            if self.layerARs[-1]:
                EndInd.append(len(self.layerARs))
        else:
            StartInd.append(0)
            for n in range(1, len(self.layerARs)):
                if ((self.basisInjectorAR and
                     not self.layerARs[n-1] and self.layerARs[n]) or
                    (self.basisARInjector and
                     self.layerARs[n-1] and not self.layerARs[n])):
                    StartInd.append(n)
                    EndInd.append(n)
            EndInd.append(len(self.layerARs))

        self.eigenEs = np.empty((0))
        self.psis = np.empty((0, self.xPoints.size))
        for n in range(0, len(StartInd)):
            dCL = copy.deepcopy(self)
            dCL.repeats = 1
            dCL.periodic = False
            # solve for Active region wavefunctions
            dCL.layerWidths = self.layerWidths[StartInd[n]:EndInd[n]]
            dCL.layerMtrls = self.layerMtrls[StartInd[n]: EndInd[n]]
            dCL.layerDopings = self.layerDopings[StartInd[n]:EndInd[n]]
            dCL.layerARs = self.layerARs[StartInd[n]:EndInd[n]]

            dCL.populate_x()
            dCL.solve_whole()

            # map dCL result back to self
            shift = sum(self.layerWidths[:StartInd[n]])
            psis, eigenEs = self.shiftPeriod(
                range(self.repeats), dCL.psis,
                dCL.eigenEs - shift*self.EField*EUnit,
                dCL.xPoints + shift)
            self.eigenEs = np.concatenate((self.eigenEs, eigenEs))
            self.psis = np.concatenate((self.psis, psis))
        self.loMatrix = [[None]*len(self.eigenEs)
                         for _ in range(len(self.eigenEs))]
        return self.eigenEs

    def shiftPeriod(self, ns, psis0, eigenEs0, xPoints=None):
        """Shift wavefunctions n (in ns) period(s) and return coorelated
        wavefunctions and EigenEs. psis0 corresponds to xPoints[0:]"""
        if xPoints is None:
            xPoints = self.xPoints
        period = sum(self.layerWidths)
        Eshift = period * self.EField * EUnit
        psis = np.empty((0, self.xPoints.size))
        eigenEs = np.empty(0)
        for n in sorted(ns, reverse=True):
            psisn = np.array([np.interp(self.xPoints, xPoints+period*n,
                                        psi) for psi in psis0])
            # filter almost zero sols)
            if len(psisn) > 0:
                idx = np.sum(psisn**2, axis=1)*self.xres > 0.1
                psis = np.concatenate((psis, psisn[idx]))
                eigenEs = np.concatenate((eigenEs, eigenEs0[idx] - Eshift*n))
        return psis, eigenEs

    def stateFilter(self):
        """Filter unbounded states:
        States with energy higher than potential at the end and wf[-2] higher
        than a threshold (1e-4) is considered unbounded. """
        bounded = []
        semiBounded = []
        for n in range(self.eigenEs.size):
            #  E = self.eigenEs[n]
            #  k = sqrt((E-self.xVc[-1])*2*self.xMc[-1])/hbar
            wf = self.psis[n, :]
            if self.eigenEs[n] < self.xVc[-1] or abs(wf[-2]) < 1e-4:
                bounded.append(n)
            elif np.all(self.xVc[wf**2 > 1e-4] > self.eigenEs[n]):
                semiBounded.append(n)
        ss = sorted(bounded + semiBounded)
        self.eigenEs = self.eigenEs[ss]
        self.psis = self.psis[ss, :]

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
        self.deltaE = Ei - Ej
        # TODO: eff mass for non-parabolic
        avgpsi_i = (psi_i[:-1] + psi_i[1:])/2
        avgxMc = (self.xMc[:-1] + self.xMc[1:])/2
        self.z = np.sum(avgpsi_i * np.diff(psi_j/self.xMc)
                        + 1/avgxMc * (avgpsi_i * np.diff(psi_j)))
        self.z *= hbar**2 / (2 * (Ei - Ej) * e0 * m0) / (1E-10)**2  # Angtrom
        return self.z

    def loTransition(self, upper, lower):
        INV_INF = 1e-20  # for infinit small decay rate (ns-1)
        if upper < lower:
            upper, lower = lower, upper

        if self.loMatrix[upper][lower] is None:

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
            # to imporve this by adding the knowledge of zero's of psi
            # convpsi = fft.irfft(
            #     np.abs(fft.rfft(psi_i*psi_j, 2*len(psi_i)))**2)[:len(psi_i)]
            # Iij = 2*self.xres**2*np.trapz(
            #     exp(-kl*self.xPoints*1E-10)*convpsi)
            # Python implementation
            # dIij = np.empty(self.xPoints.size)
            # for n in range(self.xPoints.size):
            #     x1 = self.xPoints[n]
            #     x2 = self.xPoints
            #     dIij[n] = np.sum(psi_i * psi_j * exp(-kl*abs(x1 - x2)*1e-10)
            #                      * psi_i[n] * psi_j[n] * self.xres**2)
            # Iij = np.sum(dIij)
            # C implementation
            Iij = onedq.OneDSchrodinger.cLOphononScatter(self.xres, kl,
                                                         psi_i, psi_j)
            # print(upper, lower, Iij, Iijc)
            epsInf = np.array([a.parm["epsInf"] for a in self.mtrlAlloys])
            epss = np.array([a.parm["epss"] for a in self.mtrlAlloys])
            epsrho = 1 / (1/epsInf - 1/epss)
            epsrho = (np.sum(epsrho[self.layerMtrls] * self.layerWidths)
                      / sum(self.layerWidths))
            self.loMatrix[upper][lower] = (
                mass * e0**2 * hwLO * e0 / hbar * Iij
                / (4 * hbar**2 * epsrho * eps0 * kl))
        return self.loMatrix[upper][lower] / 1e12  # unit ps^-1

    def loLifeTime(self, state):
        """ Return the life time due to LO phonon scattering of the
        given state(label)"""
        return 1/sum(self.loTransition(state, q) for q in range(state))
        #  Ei = self.eigenEs[state]
        #  psi_i = self.psis[state]
        #  hwLO = self.avghwLO()
        #  idxs = self.eigenEs <= Ei - hwLO
        #  psi_js = self.psis[idxs]
        #  Ejs = self.eigenEs[idxs]
        #  masses = m0 * sqrt(np.sum(self.xMc*psi_i**2*self.xres) *
        #                     np.sum(self.xMc*psi_js**2*self.xres, axis=1))
        #  kls = sqrt(2 * masses / hbar**2 * (Ei - Ejs - hwLO) * e0)
        #  epsInf = np.array([a.parm["epsInf"] for a in self.mtrlAlloys])
        #  epss = np.array([a.parm["epss"] for a in self.mtrlAlloys])
        #  epsrho = 1 / (1/epsInf - 1/epss)
        #  epsrho = (np.sum(epsrho[self.layerMtrls] * self.layerWidths)
        #            / sum(self.layerWidths))
        #  fjs = (masses * e0**2 * hwLO * e0 / hbar
        #             / (4 * hbar**2 * epsrho * eps0 * kls))
        #  Iijtotal = onedq.OneDSchrodinger.cLOtotal(
        #      self.xres, kls, psi_i, psi_js, fjs)
        #  return 1e12 / Iijtotal if Iijtotal > 0 else 1E20

    def dephasing(self, upper, lower):
        """Calculate the broadening gamma of transition between upper ->
        lower transition, return gamma/omega"""
        # TODO
        return 0.1

    def calc_FoM(self, upper, lower):
        """Calculate Figure of Merit.
        This function must be called after solving for wavefunctions

        Parameters
        ----------
        upper : int
        lower : int
            define the transition from upper to lower

        Reterun
        -------
        float: Figure of Merit

        Yield
        ------
        tauLower : float
            the lower state lifetime from LO scattering
        tauUpper : float
            the upper state lifetime from LO scattering
        tauUpperLower : float
            the transition rate from upper to lower due to LO scattering
        FoM : float
            the Figure of Merit in unit angstrom^2 ps
        """
        self.tauLower = self.loLifeTime(lower)
        self.tauUpper = self.loLifeTime(upper)
        self.tauUpperLower = 1/self.loTransition(upper, lower)
        self.FoM = self.z**2 * self.tauUpper * (
            1 - self.tauLower / self.tauUpperLower)
        return self.FoM

    def effective_ridx(self, wl):
        """Return the effective refractive index for TM mode"""
        self.mtrlRIdx = [(m.moleFrac * rIdx[m.A.name](wl) +
                          (1 - m.moleFrac) * rIdx[m.B.name](wl))
                         for m in self.mtrlAlloys]
        self.layerRIdx = np.array([self.mtrlRIdx[n] for n in self.layerMtrls])
        return np.average(1/self.layerRIdx**2,
                          weights=self.layerWidths)**(-1/2)

    def coupleBroadening(self, upper, lower):
        """Boradening of energy difference between upper state and lower state
        due to coupling to other states"""
        # TODO
        return 0

    def ifrBroadening(self, upper, lower):
        """Interface roughness induced broadening"""
        # TODO
        return 0

    def gainCoefficient(self, upper, lower):
        """Calculate the gain coefficient from upper -> lower transition,
        in unit cm/kA (gain/current density)

        Parameters
        ----------
        upper : int
        lower : int
            define the transition from upper to lower

        Return
        -------
        float : gain coefficient

        Yield
        ------
        wl : float
            the wavelength corresponds to upper -> lower transition, unit um
        neff : float
            the effective refractive index for wl
        gaincoef : float
            the gain coefficient on wl from upper -> lower transition,
            unit kA/cm
        """
        self.wl = c0 * h / (e0 *
                            np.abs(self.eigenEs[upper] - self.eigenEs[lower]))
        self.wl *= 1E6  # m->um
        self.neff = self.effective_ridx(self.wl)
        delta = self.dephasing(upper, lower)
        Lp = np.sum(self.layerWidths) * 1E-10  # m
        # 1E-32 angstrom^2 ps -> m^2 s, 1E5 m/A -> cm/kA
        self.gaincoef = e0 * self.FoM * 1E-27 / (
            delta * hbar * c0 * eps0 * self.effective_ridx(self.wl) * Lp)
        return self.gaincoef

    # Optimization
    def optimizeLayer(self, n, upper, lower):
        """Optmize FoM*lorenzian for n-th layer thickness"""
        # TODO: this part need to be improved!!!
        wl = self.wl
        resonancew = self.deltaE/hbar
        FoMnow = self.calc_FoM(upper, lower)/(self.tauUpperLower**2 +
                                              (resonancew-wl)**2)
        width = self.layerWidths[n]
        print(("Start Optimizing Layer NO %d " % n) +
              ("for FoM between state %d and %d.\n" % (upper, lower)) +
              ("\tStart at width=%.1f, FoM=%.1f" % (width, FoMnow)))

        for i in range(20):
            self.layerWidths[n] = width - self.xres
            self.populate_x()
            self.solve_whole()
            self.dipole(upper, lower)
            resonancew = self.deltaE/hbar
            FoMback = self.calc_FoM(upper, lower)/(self.tauUpperLower**2 +
                                                   (resonancew-wl)**2)

            self.layerWidths[n] = width + self.xres
            self.populate_x()
            self.solve_whole()
            self.dipole(upper, lower)
            resonancew = self.deltaE/hbar
            FoMforward = self.calc_FoM(upper, lower)/(self.tauUpperLower**2 +
                                                      (resonancew-wl)**2)

            FoMpp = (FoMforward + FoMback - 2*FoMnow)
            FoMp = (FoMforward - FoMback)/2
            diff = -FoMp/FoMpp
            print("%d-th interation, FoMs=[%f, %f, %f], diff=%.2f" % (
                i+1, FoMback, FoMnow, FoMforward, diff))
            if abs(diff) < 0.5:
                print("Converged at width=%.2f!" % width)
                self.layerWidths[n] = width
                self.populate_x()
                self.solve_whole()
                return
            elif FoMpp > 0:
                print("FoM'' > 0, Newton cannot go maximum. set 5*xres")
                diff = 5 if FoMp > 0 else -5
            elif abs(diff) > width/self.xres * 0.2 and abs(diff) > 5:
                print("Difference too big.. set it to 5*xres")
                diff = 5 if diff > 0 else -5
            else:
                diff = round(diff)

            width += diff * self.xres
            self.layerWidths[n] = width
            self.populate_x()
            self.solve_whole()
            self.dipole(upper, lower)
            FoMnow = self.calc_FoM(upper, lower)
            wl = h * c0 / (e0 * self.deltaE) * 1e6
            print("\twidth=%.2f, FoM=%.2f, lambda=%.1f" % (width, FoMnow, wl))
            FoMnow = FoMnow/(self.tauUpperLower**2 + (resonancew-wl)**2)

        print("Maximum iteration reached! width=%.2f" % width)
        self.layerWidths[n] = width
        self.populate_x()
        self.solve_whole()

# vim: ts=4 sw=4 sts=4 expandtab
