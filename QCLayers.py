"""
This file defines the QCLayer class for simulating QC structure
"""
import numpy as np
from numpy import sqrt
from scipy.constants import (e as e0, epsilon_0 as eps0, h as h,
                             hbar as hbar, electron_mass as m0, c as c0)
import scipy.linalg as slg
import scipy.sparse as sparse
import scipy.sparse.linalg as splg
import OneDQuantum as onedq
import Material
from OptStrata import rIdx
import copy
from typing import List

EUnit = 1e-5    # E field unit from kV/cm to V/Angstrom
# upper and lower bond extension for periodic solver
PeriodU = 1
PeriodL = 0.2

qcMaterial = {
    "InP":  ["InGaAs", "AlInAs"],
    "GaAs": ["AlGaAs"],
    "GaSb": ["InAsSb", "AlGaSb"]
}


class SchrodingerLayer(object):
    """Class for layer structure for Schrodinger solver using different
    Eigen solver.

    This is used as the base class of :class:`.QCLayers` for separation of
    material property and the solver.

    Parameters
    ----------
    xres : float
        Position resolution, in Armstrong
    Eres : float
        Energy resolution, in meV.
        This number being too large may results in miss of some states while
        this being too small will make a long computation time.
        The parameter does not mean the accuracy of the eigen energy. It's
        required for algorithm reasons because of lack of a unversal global
        root finding

    layerWidths : list of float, len = No. of layers
        Width of each layer, in unit angstrom
    layerVc : list of float, len = No. of layers
        The conduction band offset of each layer, in unit eV
    layerMc : list of float, len = No. of layers
        The conduction band effecttive mass of each layer, in unit m0
        the free space electron mass
    layerARs : list of bool, len = No. of layers
        Binaries indicating if the layer is active(True) or not(False),
        only affects basis solver

    EField : float
        External (static) electrical field, in kV/cm = 1e5 V/m
    repeats : int
        Number of repeat times for the given structure

    crystalType : str
        Default being "simple", meaning a simple parabolic effective mass is
        used for the calculation. For setting other than "simple",
        `populate_material` should be implemented.

    basisARonly : bool
        For basis solver if only the Active Region (AR) should be solved.
    basisARInjector : bool
        For basis solver if there should be seperator between AR->Injector
    basisInjectorAR : bool
        For basis solver if there should be seperator between Injector->AR

    solver : str
        The solver used for the eigen problem: 'ODE' or 'matrix'.
        By default 'ODE'. 'matrix' solver is exmperimental.
    matrixEigenCount : int
        The number of eigen pairs to calculate in the 'matrix' solver.
        It would be very expensive to calculate all of them.
    """
    xres: float
    Eres: float
    layerWidths: List[float]
    _layerVc: List[float]
    _layerMc: List[float]
    layerARs: List[bool]
    EField: float
    repeats: int
    crystalType: str
    basisARInjector: bool
    basisInjectorAR: bool
    basisARonly: bool
    matrixEigenCount: int

    def __init__(self, xres: float = 0.5, Eres: float = 0.5,
                 layerWidths: List[float] = [0.0],
                 layerARs: List[bool] = None,
                 layerVc: List[float] = None,
                 layerMc: List[float] = None,
                 EField: float = 0.0, repeats: int = 1):
        self.xres = xres
        self.Eres = Eres
        assert(isinstance(layerWidths, list))
        self.layerWidths = layerWidths
        N = len(layerWidths)
        self._layerVc = layerVc if layerVc is not None else [0.0]*N
        self._layerMc = layerMc if layerMc is not None else [1.0]*N
        self.layerARs = layerARs if layerARs is not None else [True]*N
        self.EField = EField
        self.repeats = repeats

        self.crystalType = 'simple'
        self.solver = 'ODE'
        self.matrixEigenCount = 50

        self.basisARonly = False
        self.basisARInjector = True
        self.basisInjectorAR = True

        # This is an experimental feature for periodic solver
        self.periodic = False
        # self.periodic = True

    def layerVc(self, n: int) -> float:
        """The conduction band offset at n-th layer in eV"""
        return self._layerVc[n]

    def layerMc(self, n: int) -> float:
        """The conduction band effective mass at n-th layer, in m0"""
        return self._layerMc[n]

    def populate_x(self):
        """Calculate the properties in terms of position

        Yield
        -----
        xPoints : np.ndarray of float
            position grid
        xLayerNums : np.ndarray of int
            at xPoints[i] it's xLayerNums[i]-th layer
        xVc : np.ndarray of float
            The band offset energy at each point
        xMc : np.ndarray of float
            The effective mass at different points
        """
        layerCumSum = [0] + np.cumsum(self.layerWidths).tolist()
        periodL = layerCumSum[-1]
        self.xPoints = np.arange(0, periodL * self.repeats, self.xres)
        N = self.xPoints.size
        self.xLayerNums = np.empty(N, dtype=int)
        self.xVc = np.empty(N)
        self.xMc = np.empty(N)
        self.xARs = np.empty(N, dtype=bool)
        # self.xRepeats = np.empty(N, dtype=int)

        for n in range(len(self.layerWidths)):
            indices = np.logical_or.reduce([
                (self.xPoints >= layerCumSum[n] + k * periodL)
                & (self.xPoints < layerCumSum[n+1] + k * periodL)
                for k in range(self.repeats)])
            self.xLayerNums[indices] = n
            self.xVc[indices] = self.layerVc(n)
            self.xMc[indices] = self.layerMc(n)
            self.xARs[indices] = self.layerARs[n]
        # for k in range(self.repeats):
        #     self.xRepeats[(self.xPoints < (k+1)*periodL)
        #                   & (self.xPoints >= k*periodL)] = k
        if self.crystalType != 'simple':
            self.populate_material()

        ExtField = self.xPoints * self.EField * EUnit
        self.xVc -= ExtField

    def xLayerMask(self, n: int) -> np.ndarray:
        """Return the mask for the given layer number `n`.
        A left and right extra point is included for plotting purposes."""
        xSlt = (self.xLayerNums == n)
        xSlt = (xSlt | np.roll(xSlt, 1) | np.roll(xSlt, -1))
        return ~xSlt

    def shiftPeriod(self, ns: List[int], psis0: np.ndarray,
                    eigenEs0: np.ndarray, xPoints: np.ndarray = None):
        """Shift all wave functions in `psis0` for `n` (in ns) period(s) and
        return correlated wave functions and EigenEs. """
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

    def populate_material(self):
        """This should be overrided to yeild bandParams"""
        raise NotImplementedError('Material property is not implemented')

    def solve_whole(self) -> np.ndarray:
        if self.solver == 'ODE':
            self.solve_whole_ode()
        elif self.solver == 'matrix':
            self.solve_whole_matrix()
        else:
            raise NotImplementedError(
                'The {} solver is not implemented'.format(self.solver))
        self.loMatrix = [[None]*len(self.eigenEs)
                         for _ in range(len(self.eigenEs))]
        return self.eigenEs

    def solve_whole_ode(self) -> np.ndarray:
        """
        solve eigen modes for the whole structure

        Yield
        -----
        eigenEs : np.array of float
            the eigenenergy of the layer structure
        psis : np.array of float
            the wave function
        """
        Es = np.arange(np.min(self.xVc), np.max(self.xVc), self.Eres/1E3)
        if self.crystalType == 'simple':
            self.eigenEs = onedq.cSimpleSolve1D(self.xres, Es,
                                                self.xVc, self.xMc)
            self.psis = onedq.cSimpleFillPsi(self.xres, self.eigenEs,
                                             self.xVc, self.xMc)
            return self.eigenEs
        xBand = onedq.Band(self.crystalType, *self.bandParams)
        if not self.periodic:
            self.eigenEs = onedq.cBandSolve1D(
                self.xres, Es, self.xVc, xBand)
            self.psis = onedq.cBandFillPsi(self.xres, self.eigenEs,
                                           self.xVc, xBand)
        else:
            bandOffsets = self.xVc + self.xPoints * self.EField * EUnit
            offset = max(bandOffsets) - min(bandOffsets)
            # mass = self.xMc[np.argmin(self.xVc)]
            # # ground state for triangular well
            # Emin = 2.33810741 * (hbar**2*(self.EField*EUnit)**2/(
            #     2*m0*mass*e0**2))**(1/3)
            # Es = np.linspace(np.min(self.xVc)+Emin, np.max(self.xVc), 1024)
            self.Emin = (min(bandOffsets) - PeriodL*offset)
            self.Emax = (max(bandOffsets) + PeriodU*offset)
            Eshift = self.EField*EUnit*sum(self.layerWidths)
            Es = np.arange(self.Emin - Eshift, self.Emin, self.Eres/1E3)
            eigenEs = onedq.cBandSolve1DBonded(
                self.xres, Es, self.Emin, self.Emax,
                self.EField, self.xVc, xBand)
            psis = onedq.cBandFillPsi(
                self.xres, eigenEs, self.xVc, xBand,
                Elower=self.Emin, Eupper=self.Emax, field=self.EField)
            self.psis, self.eigenEs = self.shiftPeriod(
                                          (-1, 0, 1, 2), psis, eigenEs)
        return self.eigenEs

    def solve_whole_matrix(self) -> np.ndarray:
        """
        Solve eigen modes for the whole structure with matrix eigen-solver.
        This is an experimental feature.

        Yield
        -----
        eigenEs : np.array of float
            the eigenenergy of the layer structure
        psis : np.array of float
            the wave function

        """
        # unit eV/step^2
        unit = hbar**2/(2*e0*m0*(1E-10*self.xres)**2)
        if self.crystalType == 'simple':
            # populate mass half grid self.xMc[i] is m at i+-0.5
            # TODO: reconsider how half step mass should be used
            layerCumSum = [0] + np.cumsum(self.layerWidths).tolist()
            periodL = layerCumSum[-1]
            self.xMcplus = np.empty(self.xPoints.size)
            self.xMcminus = np.empty(self.xPoints.size)
            for n in range(len(self.layerWidths)):
                Indices = np.logical_or.reduce([
                    (self.xPoints >= layerCumSum[n]+k*periodL-self.xres/2)
                    & (self.xPoints < layerCumSum[n+1]+k*periodL-self.xres/2)
                    for k in range(self.repeats)])
                self.xMcplus[Indices] = self.layerMc(n)
            self.xMcplus[-1] = self.xMcplus[-2]
            self.xMcminus[1:] = self.xMcplus[:-1]
            self.xMcminus[0] = self.xMcminus[1]

            # diagonal and sub-diagonal of Hamiltonian
            self.Hdiag = unit*(1/self.xMcplus + 1/self.xMcminus) + self.xVc
            self.Hsubd = -unit / self.xMcplus[:-1]
            self.eigenEs, self.psis = slg.eigh_tridiagonal(
                self.Hdiag, self.Hsubd, select='v',
                select_range=(np.min(self.xVc), np.max(self.xVc)))
            self.psis /= sqrt(self.xres)
            self.psis = self.psis.T
            return self.eigenEs
        if self.crystalType == 'ZincBlende':
            # The 3 band model
            print('matrix solver for ZincBlende')
            N = len(self.xPoints)
            xEg, xF, xEp, xESO = self.bandParams
            xFhalf = np.empty(N)
            xFhalf[1:] = (xF[1:] + xF[:-1])/2
            xFhalf[0] = xF[0]
            xVcHalf = np.empty(N)
            xVcHalf[1:] = (self.xVc[1:] + self.xVc[:-1])/2
            xVcHalf[0] = self.xVc[0]
            self.HBanded = np.zeros((4, 3*N))
            self.Hdiag = self.HBanded[3]

            self.HBanded[3, ::3] = 2 * (1 + 2*xFhalf) * unit + xVcHalf
            self.HBanded[3, 1::3] = self.xVc - xEg
            self.HBanded[3, 2::3] = self.xVc - xEg - xESO

            P = sqrt(xEp * unit)
            self.HBanded[2, 1::3] = sqrt(2/3)*P
            self.HBanded[2, 3::3] = sqrt(1/3)*P[:-1]

            self.HBanded[1, 2::3] = -sqrt(1/3)*P
            self.HBanded[1, 3::3] = -sqrt(2/3)*P[:-1]

            self.HBanded[0, 3::3] = -(1 + 2*xF[:-1]) * unit

            self.Hsparse = sparse.diags(
                [self.HBanded[0, 3:], self.HBanded[1, 2:], self.HBanded[2, 1:],
                 self.HBanded[3, :],
                 self.HBanded[2, 1:], self.HBanded[1, 2:], self.HBanded[0, 3:]
                 ], [-3, -2, -1, 0, 1, 2, 3], shape=(3*N, 3*N)
            )
            Es_low = np.min(self.xVc)
            Es_hi = np.max(self.xVc)
            self.eigen_all, self.psi_all = splg.eigsh(
                self.Hsparse, self.matrixEigenCount, sigma=(Es_low + Es_hi)/2)
            # self.eigen_all, self.psi_all = slg.eig_banded(
            #     self.HBanded, select='v', select_range=(Es_low, Es_hi))
            # normalization should be sum(self.psi_all**2)*self.xres = 1
            self.psi_all /= sqrt(self.xres)
            self.psis = np.zeros((self.eigen_all.shape[0], N))
            psis = self.psi_all[::3, :].T
            self.psis[:, :-1] = (psis[:, 1:] + psis[:, :-1])/2
            for psi in self.psis:
                if psi[np.argmax(np.abs(psi) > 1E-3)] < 0:
                    psi[:] = -psi
            self.eigenEs = self.eigen_all
            return self.eigenEs
        raise NotImplementedError('Matrix solver is not implemented '
                                  'for {}'.format(self.crystalType))

    def solve_basis(self) -> np.ndarray:
        """
        solve basis for the QC device, with each basis being the eigen mode of
        a separate part of the layer structure

        Yield
        -----
        eigenEs : np.array of float
            the eigenenergy of the layer structure
        psis : np.array of float
            the wave function
        """
        StartInd = []
        EndInd = []
        # Get the region of interest
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
            dCL.reset_for_basis(StartInd[n], EndInd[n])
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
        return self.eigenEs

    def reset_for_basis(self, start: int, end: int):
        """Reset the parameters for only solving the layers of [start:end].
        This is a helper method for solve_basis"""
        self.repeats = 1
        self.layerWidths = self.layerWidths[start:end]
        self.layerARs = self.layerARs[start:end]
        self._layerVc = self._layerVc[start:end]
        self._layerMc = self._layerMc[start:end]
        self.periodic = False

    def dipole(self, upper: int, lower: int) -> float:
        """Return Electrical dipole between upper and lower states,
        update self.dipole. Should be called for any other related physics
        quantities."""
        if upper < lower:
            upper, lower = lower, upper
        psi_i = self.psis[upper, :]
        psi_j = self.psis[lower, :]
        Ei = self.eigenEs[upper]
        Ej = self.eigenEs[lower]
        # interpolate half grid wave function and effective mass
        # avgpsi_i = (psi_i[:-1] + psi_i[1:])/2
        # avgpsi_j = (psi_j[:-1] + psi_j[1:])/2
        # avgxMc = (self.xMc[:-1] + self.xMc[1:])/2
        # d_z 1/m + 1/m d_z = [d_z 1/m d_z, z] ~ [P^2, z] ~ [H, z],
        # where d_z means spatial derivative d/d z, P is momentum
        # self.z = np.sum(avgpsi_i * np.diff(psi_j/self.xMc)
        #                 # + 1/avgxMc * (avgpsi_i * np.diff(psi_j)))
        # Eq.(8) in PhysRevB.50.8663
        xInvMc_i = self._xBandMassInv(Ei)
        xInvMc_j = self._xBandMassInv(Ej)
        self.z = np.trapz(psi_i * xInvMc_j * np.gradient(psi_j) -
                          np.gradient(psi_i) * xInvMc_i * psi_j)
        self.z *= hbar**2 / (2 * (Ei - Ej) * e0 * m0) / (1E-10)**2  # Angstrom
        return np.abs(self.z)

    def _xBandMassInv(self, energy):
        """
        Return the energy dependent effective mass in the form of m0/m
        for the given eigen energy.

        The xPoints should have been populated.
        This method is not used in state solver. To completely re-write this
        method it is required to also do it in the C library.
        """
        if self.crystalType == 'simple':
            return self.xMc
        if self.crystalType == 'ZincBlende':
            xEg, xF, xEp, xESO = self.bandParams
            E = energy - self.xVc
            return 1 + 2*xF + 1/3 * xEp/(E+xEg+xESO) + 2/3 * xEp/(E+xEg)
        else:
            raise NotImplementedError(
                 'Material property for {} crystal is not implemented'.format(
                     self.crystalType
                 ))


class QCLayers(SchrodingerLayer):
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
    Position resolution, in Armstrong
Eres : float
    Energy resolution, in meV.
    This number being too large may results in miss of some states while
    this being too small will make a long computation time.
    The parameter does not mean the accuracy of the eigen energy. It's
    required for algorithm reasons because of lack of a unversal global
    root finding.
wl : float
    The wavelength for the design, in unit um, for book keeping and
    optimization, but doesn't go into quantum solver

layerWidths : list of float, len = No. of layers
    Width of each layer, in angstrom
layerMtrls : list of int, len = No. of layers
    Label of materials, depending on substrate
layerDopings : list of float, len = No. of layers
    Doping per volume in unit 1e17 cm-3
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
    """
    def __init__(self, substrate="InP", materials=["InGaAs", "AlInAs"],
                 moleFracs=[0.53, 0.52], xres=0.5, Eres=0.5,
                 layerWidths=[0.0], layerMtrls=None, layerDopings=None,
                 layerARs=None, EField=0, repeats=3, T=300.0, solver="ODE",
                 description="", wl=3.0):
        assert(isinstance(layerWidths, list))
        N = len(layerWidths)
        super().__init__(xres, Eres, layerWidths, layerARs,
                         # layerVc and layerMc is not used for this sub class
                         None, None,
                         EField, repeats)
        self.substrate = substrate
        self.crystalType = Material.MParm[substrate]["Crystal"]
        self.materials = materials
        self.moleFracs = moleFracs
        self.layerMtrls = [0]*N if layerMtrls is None else layerMtrls
        self.layerDopings = [0.0]*N if layerDopings is None else layerDopings
        self.Temperature = T
        self.solver = solver
        self.description = description
        self.wl = 3.0

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
        """Add a new material possibility"""
        self.materials.append(mtrl if mtrl else
                              qcMaterial[self.substrate][0])
        self.moleFracs.append(moleFrac if moleFrac else 0.0)
        self.mtrlAlloys.append(Material.Alloy(
            self.materials[-1], self.moleFracs[-1], self.Temperature))
        self.mtrlAlloys[-1].set_strain(self.a_parallel)

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
        """Return average LO phonon energy in unit eV"""
        if sum(self.layerWidths) <= 1e-5:
            return -1
        sumhwlo = sum(self.mtrlAlloys[self.layerMtrls[n]].parm['hwLO']
                      * self.layerWidths[n]
                      for n in range(len(self.layerWidths)))
        return sumhwlo / sum(self.layerWidths)

    def layerVc(self, n: int):
        return self.mtrlAlloys[self.layerMtrls[n]].parm['EcG']

    def layerMc(self, n: int):
        return self.mtrlAlloys[self.layerMtrls[n]].parm['me0']

    def populate_material(self):
        """
        Update following band structure parameters (with type *np.array
        of float*): **xVc, xVX, xVL, xVLH, xVSO,
        bandParams = (xEg, xF, xEp, xESO)**
        """
        if self.crystalType == 'ZincBlende':
            N = self.xPoints.size
            self.xDopings = np.empty(N)

            self.xVX = np.empty(N)   # X point in the band
            self.xVL = np.empty(N)   # L point in the band
            self.xVLH = np.empty(N)  # The light hole valence band
            self.xVSO = np.empty(N)  # The heavy hole valence band

            # band parameters
            xF = np.empty(N)
            xEg = np.empty(N)
            xESO = np.empty(N)
            xEp = np.empty(N)

            for n in range(len(self.layerWidths)):
                indices = (self.xLayerNums == n)
                self.xDopings[indices] = self.layerDopings[n]
                for (p, key) in ((self.xVLH, 'EvLH'), (self.xVSO, 'EvSO'),
                                 (self.xVX, 'EcX'), (self.xVL, 'EcL'),
                                 (xEg, 'EgLH'), (xESO, 'ESO'),
                                 (xEp, 'Ep'), (xF, 'F')):
                    p[indices] = self.mtrlAlloys[self.layerMtrls[n]].parm[key]
            # xF = -np.ones(N)/2
            self.bandParams = (xEg, xF, xEp, xESO)

            ExtField = self.xPoints * self.EField * EUnit
            for p in (self.xVX, self.xVL, self.xVLH, self.xVSO):
                p -= ExtField
        else:
            raise NotImplementedError(
                 'Material property for {} crystal is not implemented'.format(
                     self.crystalType
                 ))

    def reset_for_basis(self, start, end):
        super().reset_for_basis(start, end)
        self.layerMtrls = self.layerMtrls[start:end]
        self.layerDopings = self.layerDopings[start:end]

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

    def loTransition(self, upper, lower):
        INV_INF = 1e-20  # for infinite small decay rate (ns-1)
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
                # LO phonon scattering doesn't happen
                return INV_INF

            # TODO: improve expression for eff mass
            mass = m0 * sqrt(np.sum(self.xMc * psi_i**2 * self.xres)
                             * np.sum(self.xMc * psi_j**2 * self.xres))
            kl = sqrt(2 * mass / hbar**2 * (Ei-Ej-hwLO) * e0)
            # to improve this by adding the knowledge of zero's of psi
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
        # return 1/sum(self.loTransition(state, q) for q in range(state))
        Ei = self.eigenEs[state]
        psi_i = self.psis[state]
        hwLO = self.avghwLO()
        # return 1/sum(self.loTransition(state, q) for q in range(stat
        #                if self.eigenEs[q] <= Ei - hwLO)
        idxs = self.eigenEs <= Ei - hwLO
        psi_js = self.psis[idxs]
        Ejs = self.eigenEs[idxs]
        masses = m0 * sqrt(np.sum(self.xMc*psi_i**2*self.xres) *
                           np.sum(self.xMc*psi_js**2*self.xres, axis=1))
        kls = sqrt(2 * masses / hbar**2 * (Ei - Ejs - hwLO) * e0)
        epsInf = np.array([a.parm["epsInf"] for a in self.mtrlAlloys])
        epss = np.array([a.parm["epss"] for a in self.mtrlAlloys])
        epsrho = 1 / (1/epsInf - 1/epss)
        epsrho = (np.sum(epsrho[self.layerMtrls] * self.layerWidths)
                  / sum(self.layerWidths))
        fjs = (masses * e0**2 * hwLO * e0 / hbar
               / (4 * hbar**2 * epsrho * eps0 * kls))
        Iijtotal = onedq.OneDSchrodinger.cLOtotal(
            self.xres, kls, psi_i, psi_js, fjs)
        return 1e12 / Iijtotal if Iijtotal > 0 else 1E20

    def dephasing(self, upper, lower):
        """Calculate the broadening gamma of transition between upper ->
        lower transition, return gamma/omega"""
        # TODO
        return 0.1

    def calc_FoM(self, upper, lower):
        """Calculate Figure of Merit.
        This function must be called after solving for wave functions

        Parameters
        ----------
        upper : int
        lower : int
            define the transition from upper to lower

        Return
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
        if sum(self.layerWidths) == 0:
            return 1.0
        self.mtrlRIdx = [(m.moleFrac * rIdx[m.A.name](wl) +
                          (1 - m.moleFrac) * rIdx[m.B.name](wl))
                         for m in self.mtrlAlloys]
        self.layerRIdx = np.array([self.mtrlRIdx[n] for n in self.layerMtrls])
        return np.average(1/self.layerRIdx**2,
                          weights=self.layerWidths)**(-1/2)

    def coupleBroadening(self, upper, lower):
        """Broadening of energy difference between upper state and lower state
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
        """Optimize FoM*Lorenzian for n-th layer thickness"""
        # TODO: this part need to be improved!!!
        wl = self.wl
        Ei = self.eigenEs[upper]
        Ej = self.eigenEs[lower]
        self.deltaE = Ei - Ej
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
            print("%d-th iteration, FoMs=[%f, %f, %f], diff=%.2f" % (
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
