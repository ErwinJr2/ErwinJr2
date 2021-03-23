"""
This file defines the QCLayer class for simulating QC structure
"""
from re import L
import numpy as np
from numpy import sqrt, pi, exp
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

    ifrDelta : List of float, len = No. of layers
        The standard deviation of the interface roughness for the interface at
        layer n and layer n+1, in unit angstrom.
        Default zero
    ifrLambda : List of float, len = No. of layers
        The correlation length of the interface roughness for the interface at
        layer n and layer n+1, in unit angstrom.
        Default zero
    avghwLO : float
        The average LO phonon energy in unit eV.
        This value is used for LO phonon scattering calculation.
    epsrho : float
        The effective relative permitivity for LO phonon.
        1/epsrho = 1/(epsilon for high frequency) - 1/(epsilon for static)
        The value is used for LO phonon scattering calculation.

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
    includeIFR : bool
        Weather to include IFR scattering for performance estimation.
    statePerRepeat : int
        Number of states per repeat, used for calculating matrixEigenCount
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
    ifrDelta: List[float]
    ifrLambda: List[float]
    avghwLO: float
    epsrho: float
    EField: float
    repeats: int
    crystalType: str
    basisARInjector: bool
    basisInjectorAR: bool
    basisARonly: bool
    matrixEigenCount: int

    solver: str
    includeIFR: bool

    def __init__(self, xres: float = 0.5, Eres: float = 0.5,
                 statePerRepeat: int = 20,
                 layerWidths: List[float] = [0.0],
                 layerARs: List[bool] = None,
                 layerVc: List[float] = None,
                 layerMc: List[float] = None,
                 ifrDelta: List[float] = None,
                 ifrLambda: List[float] = None,
                 avghwLO: float = 35E-3,
                 epsrho: float = 69.0,
                 EField: float = 0.0, repeats: int = 1):
        self.xres = xres
        self.Eres = Eres
        self.statePerRepeat = statePerRepeat
        assert(isinstance(layerWidths, list))
        self.layerWidths = layerWidths
        N = len(layerWidths)
        self._layerVc = layerVc if layerVc is not None else [0.0]*N
        self._layerMc = layerMc if layerMc is not None else [1.0]*N
        self.layerARs = layerARs if layerARs is not None else [True]*N
        self.ifrDelta = ifrDelta if ifrDelta is not None else [0.0]*N
        self.ifrLambda = ifrLambda if ifrLambda is not None else [0.0]*N
        self.EField = EField
        self.repeats = repeats
        self.avghwLO = avghwLO
        self.epsrho = epsrho

        self.crystalType = 'simple'
        self.solver = 'ODE'
        self.includeIFR = False
        self.matrixEigenCount = statePerRepeat * repeats

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
        self.xPoints = np.arange(0, periodL*self.repeats, self.xres)
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

        self.offset = max(self.xVc) - min(self.xVc)
        ExtField = self.xPoints * self.EField * EUnit
        self.xVc -= ExtField

        self.Es = np.arange(np.min(self.xVc), np.max(self.xVc), self.Eres/1E3)
        self.matrixEigenCount = self.repeats * self.statePerRepeat
        # solving for eigen states energy close to sigma
        self.matrixSigma = (np.min(self.xVc) + np.max(self.xVc))/2

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
        self.ifrMatrix = [[None]*len(self.eigenEs)
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
        if self.crystalType == 'simple':
            self.eigenEs = onedq.cSimpleSolve1D(self.xres, self.Es,
                                                self.xVc, self.xMc)
            self.psis = onedq.cSimpleFillPsi(self.xres, self.eigenEs,
                                             self.xVc, self.xMc)
            return self.eigenEs
        xBand = onedq.Band(self.crystalType, *self.bandParams)
        if not self.periodic:
            self.eigenEs = onedq.cBandSolve1D(
                self.xres, self.Es, self.xVc, xBand)
            self.psis = onedq.cBandFillPsi(
                self.xres, self.eigenEs, self.xVc, xBand)
        else:
            # This is experimental and may removed in the future if turn out to
            # be useless. So is the corresponding C functions.
            netVc = self.xVc + self.xPoints * self.EField * EUnit
            # mass = self.xMc[np.argmin(self.xVc)]
            # # ground state for triangular well
            # Emin = 2.33810741 * (hbar**2*(self.EField*EUnit)**2/(
            #     2*m0*mass*e0**2))**(1/3)
            # Es = np.linspace(np.min(self.xVc)+Emin, np.max(self.xVc), 1024)
            self.Emin = (min(netVc) - 0.2*self.offset)
            self.Emax = (max(netVc) + 1.1*self.offset)
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

        if self.crystalType == 'ZincBlende':
            # To restore lh and so band to keep consistent with matrix solver
            xEg, xF, xEp, xESO = self.bandParams
            kunit = hbar**2/(2*e0*m0*(1E-10*self.xres)**2)
            xP = sqrt(xEp * kunit)
            dphic = np.zeros(self.psis.shape)
            dphic[:, 1:-1] = (self.psis[:, 2:] - self.psis[:, :-2])/2
            self.philh = np.zeros(self.psis.shape)
            xE = np.broadcast_to(self.eigenEs.reshape(-1, 1), self.psis.shape)
            xE = xE - self.xVc
            self.philh = -sqrt(2/3)*xP/(xE+xEg) * dphic
            self.phiso = np.zeros(self.psis.shape)
            self.phiso = sqrt(1/3)*xP/(xE+xEg+xESO) * dphic
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
        kunit = hbar**2/(2*e0*m0*(1E-10*self.xres)**2)
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
            self.Hdiag = kunit*(1/self.xMcplus + 1/self.xMcminus) + self.xVc
            self.Hsubd = -kunit / self.xMcplus[:-1]
            # self.eigenEs, self.psis = slg.eigh_tridiagonal(
            #     self.Hdiag, self.Hsubd, select='v',
            #     select_range=(np.min(self.xVc), np.max(self.xVc)))
            self.Hsubd = -kunit / self.xMcplus
            N = len(self.xPoints)
            self.Hsparse = sparse.diags([self.Hsubd, self.Hdiag, self.Hsubd],
                                        [-1, 0, 1], shape=(N, N))
            self.eigenEs, self.psis = splg.eigsh(
                self.Hsparse, self.matrixEigenCount, sigma=self.matrixSigma,
                tol=1E-7)
            self.psis /= sqrt(self.xres)
            self.psis = self.psis.T
            return self.eigenEs
        if self.crystalType == 'ZincBlende':
            self.populate_Kane_matrix()
            # self.eigen_all, self.psi_all = slg.eig_banded(
            #     self.HBanded, select='v', select_range=(Es_low, Es_hi))
            self.eigen_all, self.psi_all = splg.eigsh(
                self.Hsparse, self.matrixEigenCount, sigma=self.matrixSigma,
                tol=1E-7)
            # normalization should be sum(self.psi_all**2)*self.xres = 1
            self.psi_all /= sqrt(self.xres)
            for n in range(len(self.eigen_all)):
                phic = self.psi_all[::3, n]
                # for consistency of the phase definition with the ODE solver
                if phic[np.argmax(np.abs(phic) > 1E-3)] < 0:
                    self.psi_all[:, n] = -self.psi_all[:, n]
            self.psis = np.zeros((self.eigen_all.shape[0], self.xPoints.size))
            psis = self.psi_all[::3, :].T
            self.psis[:, :-1] = (psis[:, 1:] + psis[:, :-1])/2
            self.philh = self.psi_all[1::3, :].T
            self.phiso = self.psi_all[2::3, :].T
            self.eigenEs = self.eigen_all
            return self.eigenEs
        raise NotImplementedError('Matrix solver is not implemented '
                                  'for {}'.format(self.crystalType))

    def populate_Kane_matrix(self) -> sparse.spmatrix:
        """
        Populate the finite difference Hamiltonian operation for the Kane
        3 band model.

        Return
        ------
        Hsparse : scipy.sparse.spmatrix
            A sparse matrix for the Hamiltonian

        Yield
        -----
        Hsparse : scipy.sparse.spmatrix
            Same above
        Hbanded : np.ndarray
            Upper banded form of the matrix for the banded solver, as
            used in scipy.linalg.eig_banded
        """
        assert(self.crystalType == 'ZincBlende')
        kunit = hbar**2/(2*e0*m0*(1E-10*self.xres)**2)
        N = len(self.xPoints)
        xEg, xF, xEp, xESO = self.bandParams
        xFhalf = np.empty(N)
        xFhalf[1:] = (xF[1:] + xF[:-1])/2
        xFhalf[0] = xF[0]
        xVcHalf = np.empty(N)
        xVcHalf[1:] = (self.xVc[1:] + self.xVc[:-1])/2
        xVcHalf[0] = self.xVc[0]
        self.HBanded = np.zeros((4, 3*N))
        self.HBanded[3, ::3] = 2*(1 + 2*xFhalf)*kunit + xVcHalf
        self.HBanded[3, 1::3] = self.xVc - xEg  # lh band
        self.HBanded[3, 2::3] = self.xVc - xEg - xESO  # so band
        P = sqrt(xEp * kunit)
        self.HBanded[2, 1::3] = sqrt(2/3)*P
        self.HBanded[2, 3::3] = sqrt(1/3)*P[:-1]
        self.HBanded[1, 2::3] = -sqrt(1/3)*P
        self.HBanded[1, 3::3] = -sqrt(2/3)*P[:-1]
        self.HBanded[0, 3::3] = -(1 + 2*xF[:-1]) * kunit
        if hasattr(self, 'luttinger'):
            gamma1, gamma2, _ = self.luttinger
            tlh = gamma1 + 2*gamma2 - 2*xEp/xEg/3
            tso = gamma1 - xEp/xEg/3
            tlhhalf = (tlh[1:] + tlh[:-1])/2
            tsohalf = (tso[1:] + tso[:-1])/2
            tlh[1:-1] = (tlhhalf[1:] + tlhhalf[:-1])/2
            tso[1:-1] = (tsohalf[1:] + tsohalf[:-1])/2
            self.HBanded[3, 1::3] -= 2*tlh*kunit
            self.HBanded[3, 2::3] -= 2*tso*kunit
            self.HBanded[0, 4::3] = tlhhalf*kunit
            self.HBanded[0, 5::3] = tsohalf*kunit
        self.Hsparse = sparse.diags(
            [self.HBanded[0, 3:], self.HBanded[1, 2:], self.HBanded[2, 1:],
             self.HBanded[3, :],
             self.HBanded[2, 1:], self.HBanded[1, 2:], self.HBanded[0, 3:]
             ], [-3, -2, -1, 0, 1, 2, 3], shape=(3*N, 3*N))
        return self.Hsparse

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
        """Return Electrical dipole between upper and lower states
        in unit angstrom, update self.dipole.
        Should be called for any other related physics quantities."""
        if self.solver == 'ODE':
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
            self.z *= hbar**2 / (2*(Ej-Ei)*e0*m0) / (1E-10)**2  # Angstrom
        elif self.solver == 'matrix' and self.crystalType == 'ZincBlende':
            self.z = self.xres * sum(
                np.trapz(self.psi_all[t::3, upper] * self.xPoints
                         * self.psi_all[t::3, lower]) for t in range(3))
        else:
            raise NotImplementedError(
                '{} not implemented for diple'.format(self.solver))
        return self.z

    def loTransition(self, upper, lower):
        """The LO phonon transition lifetime from upper to lower,
        at zero temperature"""
        # TODO: finite temperature version
        INV_INF = 1e-20  # for infinite small decay rate (ns-1)
        if self.loMatrix[upper][lower] is None:

            psi_i = self.psis[upper, :]
            psi_j = self.psis[lower, :]
            Ei = self.eigenEs[upper]
            Ej = self.eigenEs[lower]

            if Ei - Ej - self.avghwLO < 0:
                # energy difference is smaller than a LO phonon
                # LO phonon scattering doesn't happen
                return INV_INF

            # mi = m0 * np.trapz(self.xMc * psi_i**2) * self.xres
            mj = m0 * np.trapz(self.xMc * psi_j**2) * self.xres
            kl = sqrt(2 * mj / hbar**2 * (Ei - Ej - self.avghwLO) * e0)
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
            self.loMatrix[upper][lower] = (
                mj * e0**2 * self.avghwLO * e0 / hbar * Iij
                / (4 * hbar**2 * self.epsrho * eps0 * kl))
        return self.loMatrix[upper][lower] / 1e12  # unit ps^-1

    def loLifeTime(self, state):
        """ Return the life time due to LO phonon scattering of the
        given state(label)"""
        # return 1/sum(self.loTransition(state, q) for q in range(state))
        Ei = self.eigenEs[state]
        psi_i = self.psis[state]
        # return 1/sum(self.loTransition(state, q) for q in range(stat
        #                if self.eigenEs[q] <= Ei - self.avghwLO)
        idxs = self.eigenEs <= Ei - self.avghwLO
        psi_js = self.psis[idxs]
        Ejs = self.eigenEs[idxs]
        # mi = m0 * np.trapz(self.xMc * psi_i**2) * self.xres
        mjs = m0 * np.trapz(self.xMc * psi_js**2, axis=1) * self.xres
        kls = sqrt(2 * mjs / hbar**2 * (Ei - Ejs - self.avghwLO) * e0)
        fjs = (mjs * e0**2 * self.avghwLO * e0 / hbar
               / (4 * hbar**2 * self.epsrho * eps0 * kls))
        Iijtotal = onedq.OneDSchrodinger.cLOtotal(
            self.xres, kls, psi_i, psi_js, fjs)
        return 1e12 / Iijtotal if Iijtotal > 0 else 1E20

    def ifrTransition(self, upper: int, lower: int) -> float:
        r"""Calculate the interface roughness (IFR) transition rate from
        upper to lower state at zero temperature, in unit ps^-1.

        .. math::
            \frac{1}{\tau_{ij}^\text{IFR}} =
            \frac{\pi m^*_j}{\hbar^3} \sum_n
            \Delta_n^2\Lambda_n^2\delta U_n^2
            \left|\psi_i(z_n)\psi_j^*(z_n)\right|^2
            \mathrm e^{- \Lambda^2 m_j^* (E_i - E_j))/2\hbar^2}
        """
        # TODO: finite temperature
        if upper < lower:
            return 1e-20
        if not self.ifrMatrix[upper][lower] is None:
            return self.ifrMatrix[upper][lower]
        psi_i = self.psis[upper, :]
        psi_j = self.psis[lower, :]
        Ei = self.eigenEs[upper]
        Ej = self.eigenEs[lower]
        mj = m0 * np.trapz(self.xMc * psi_j**2) * self.xres
        kl = sqrt(2 * mj / hbar**2 * (Ei-Ej) * e0)
        tauInv = 0
        zn = 0
        layerN = len(self.layerWidths)
        for _ in range(self.repeats):
            for n in range(layerN):
                lamb = self.ifrLambda[n] * 1E-10  # to m
                delt = self.ifrDelta[n]
                dU = (self.layerVc((n+1) % layerN) - self.layerVc(n))
                dU *= e0  # to J
                # find interface
                zn += self.layerWidths[n]
                zIdx = np.argmax(self.xPoints >= zn)
                if zIdx == 0 or zIdx == len(self.xPoints)-1:
                    continue
                z1 = self.xPoints[zIdx-1]
                z2 = self.xPoints[zIdx]
                psi_iz = (psi_i[zIdx-1]*(zn-z2) - psi_i[zIdx]*(z1-zn))/(z2-z1)
                psi_jz = (psi_j[zIdx-1]*(zn-z2) - psi_j[zIdx]*(z1-zn))/(z2-z1)
                tauInv += (pi * mj / hbar**3 * delt**2 * lamb**2 * dU**2
                           * (psi_iz*psi_jz)**2 * exp(-lamb**2 * kl**2 / 4))
        self.ifrMatrix[upper][lower] = tauInv / 1e12  # unit ps^-1
        return self.ifrMatrix[upper][lower]

    def ifrLifeTime(self, state):
        """ Return to total life time due to IFR scattering"""
        return 1/sum(self.ifrTransition(state, q) for q in range(state))

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

materials : list of str, len >= 2
    Name of alloys
moleFrac : list of float, len = Mp. of materials
    mole fraction for each possible layer material,
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

customIFR : bool
    Wether to use a customized IFR parameter rather than a material decidede
    parameter.
mtrlIFRLambda : list of float, len = No. of materials
    The interface roughness lambda after materials[n]
mtrlIFRDelta : list of float, len = No. of materials
    The interface roughness delta after materials[n]

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
                 customIFR=False, mtrlIFRLambda=None, mtrlIFRDelta=None,
                 ifrDelta=None, ifrLambda=None,
                 layerARs=None, EField=0, repeats=3, T=300.0, solver="ODE",
                 description="", wl=3.0):
        assert(isinstance(layerWidths, list))
        assert(isinstance(materials, list))
        assert(isinstance(moleFracs, list))
        N = len(layerWidths)
        M = len(materials)
        assert(M >= 1)
        assert(len(moleFracs) == M)
        self.substrate = substrate
        self.materials = materials
        self.moleFracs = moleFracs
        self.layerMtrls = [0]*N if layerMtrls is None else layerMtrls
        self.layerDopings = [0.0]*N if layerDopings is None else layerDopings
        self.Temperature = T
        self.customIFR = customIFR
        if not customIFR:
            if isinstance(mtrlIFRDelta, list):
                assert(len(mtrlIFRDelta) == M)
                assert(isinstance(mtrlIFRLambda, list))
                assert(len(mtrlIFRLambda) == M)
                self.mtrlIFRDelta = mtrlIFRDelta
                self.mtrlIFRLambda = mtrlIFRLambda
            else:
                self.mtrlIFRDelta = [mtrlIFRDelta or 0.0] * M
                self.mtrlIFRLambda = [mtrlIFRLambda or 0.0] * M
            ifrDelta, ifrLambda = self._get_IFRList()
        self.description = description
        super().__init__(xres=xres, Eres=Eres,
                         layerWidths=layerWidths, layerARs=layerARs,
                         ifrDelta=ifrDelta, ifrLambda=ifrLambda,
                         EField=EField, repeats=repeats)
        self.crystalType = Material.MParm[substrate]["Crystal"]
        self.subM = Material.Material(self.substrate, self.Temperature)
        self.wl = 3.0
        self.solver = solver
        # self.solver = 'matrix'
        self.update_strain()

    def _get_IFRList(self):
        """Get IFR parameters for SchrodingerLayer. Should be called
        every time the material list changes."""
        assert(not self.customIFR)
        if self.mtrlIFRDelta is not None:
            ifrDelta = [self.mtrlIFRDelta[m] for m in self.layerMtrls]
        else:
            ifrDelta = [0] * len(self.layerMtrls)
        if self.mtrlIFRLambda is not None:
            ifrLambda = [self.mtrlIFRLambda[m] for m in self.layerMtrls]
        else:
            ifrLambda = [0] * len(self.layerMtrls)
        return ifrDelta, ifrLambda

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

    def mtrlOffset(self):
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
        # LO phonon
        if sum(self.layerWidths) <= 1e-5:
            self.avghwLO = -1
            self.epsrho = 69.0
        else:
            sumhwlo = sum(self.mtrlAlloys[self.layerMtrls[n]].parm['hwLO']
                          * self.layerWidths[n]
                          for n in range(len(self.layerWidths)))
            self.avghwLO = sumhwlo / sum(self.layerWidths)
            epsInf = np.array([a.parm["epsInf"] for a in self.mtrlAlloys])
            epss = np.array([a.parm["epss"] for a in self.mtrlAlloys])
            epsrho = 1 / (1/epsInf - 1/epss)
            self.epsrho = (np.sum(epsrho[self.layerMtrls] * self.layerWidths)
                           / sum(self.layerWidths))

        # IFR
        if not self.customIFR:
            self.ifrDelta, self.ifrLambda = self._get_IFRList()

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
        tauLO_l : float
            the lower state lifetime from LO scattering
        tauLO_u : float
            the upper state lifetime from LO scattering
        tauLO_ul : float
            the transition rate from upper to lower due to LO scattering
        tauIFR_l : float
            the lower state lifetime from IFR scattering
        tauIFR_u : float
            the upper state lifetime from IFR scattering
        tauIFR_ul : float
            the transition rate from upper to lower due to IFR scattering
        tau_u : float
            1/(1/tauLO_u + 1/tauIFR_u)
        tau_l : float
            1/(1/tauLO_l + 1/tauIFR_l)
        FoM : float
            the Figure of Merit in unit angstrom^2 ps
        """
        self.tauLO_ul = 1/self.loTransition(upper, lower)
        self.tauLO_l = self.loLifeTime(lower)
        self.tauLO_u = self.loLifeTime(upper)
        if self.includeIFR:
            self.tauIFR_ul = 1/self.ifrTransition(upper, lower)
            self.tauIFR_u = self.ifrLifeTime(upper)
            self.tauIFR_l = self.ifrLifeTime(lower)
            self.tau_u = 1/(1/self.tauLO_u + 1/self.tauIFR_u)
            self.tau_l = 1/(1/self.tauLO_l + 1/self.tauIFR_l)
            self.tau_ul = 1/(1/self.tauLO_ul + 1/self.tauIFR_ul)
        else:
            self.tau_u = self.tauLO_u
            self.tau_l = self.tauLO_l
            self.tau_ul = self.tauLO_ul
        self.FoM = self.z**2 * self.tau_u * (
            1 - self.tau_l / self.tau_ul)
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
        FoMnow = self.calc_FoM(upper, lower)/(self.tauLO_ul**2 +
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
            FoMback = self.calc_FoM(upper, lower)/(self.tauLO_ul**2 +
                                                   (resonancew-wl)**2)

            self.layerWidths[n] = width + self.xres
            self.populate_x()
            self.solve_whole()
            self.dipole(upper, lower)
            resonancew = self.deltaE/hbar
            FoMforward = self.calc_FoM(upper, lower)/(self.tauLO_ul**2 +
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
            FoMnow = FoMnow/(self.tauLO_ul**2 + (resonancew-wl)**2)

        print("Maximum iteration reached! width=%.2f" % width)
        self.layerWidths[n] = width
        self.populate_x()
        self.solve_whole()

# vim: ts=4 sw=4 sts=4 expandtab
