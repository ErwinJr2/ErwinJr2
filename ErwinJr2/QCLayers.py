import numpy as np
from numpy import sqrt, pi, exp
from numpy import fft
from scipy.constants import (e as e0, epsilon_0 as eps0, h as h,
                             hbar as hbar, electron_mass as m0, c as c0)
from scipy.linalg import null_space
import scipy.sparse as sparse
import scipy.sparse.linalg as splg
try:
    from . import OneDQuantum as onedq
except OSError:
    onedq = None
    print('C library is not compiled. Features are limited.')
from . import Material
import copy
from typing import List, Tuple, Union

EUNIT = 1e-5    # E field unit from kV/cm to V/Angstrom
BASISPAD = 100  #: padding barrier for basis solver, unit Angstrom
INV_INF = 1e-20  #: for infinite small decay rate (ns-1)
# used for typing as either an array or a float number
ScalerOrArray = Union[float, np.ndarray]

QCMaterial = {
    "InP":  ["InGaAs", "AlInAs"],
    "GaAs": ["AlGaAs"],
    "GaSb": ["InAsSb", "AlGaSb"]
}  #: supported substrate/material set


class StateRecognizeError(Exception):
    """Raised when QClayers cannot recognize a period of states.
    This typically means the number of periods is too small, or for single
    period structure the barrier at the end of the structure is not long
    enough to bound a state.

    Attributes:
        expression -- input expression in which the error occurred."""
    def __init__(self, expression=''):
        self.expression = expression


class SchrodingerLayer(object):
    """Class for layer structure for Schrodinger solver.

    This is used as the base class of :class:`.QCLayers` for separation of
    material property and the solver.

Parameters
----------
xres :
    Position resolution, in Armstrong
Eres :
    Energy resolution, in meV.
    This number being too large may results in miss of some states while
    this being too small will make a long computation time.
    The parameter does not mean the accuracy of the eigen energy. It's
    required for algorithm reasons because of lack of a universal global
    root finding.
statePerRepeat :
    Number of states per repeat, used for calculating matrixEigenCount

layerWidths :
    Width of each layer, in angstrom. len = No. of layers
layerVc :
    The conduction band offset of each layer, in unit eV, len = No. of layers
layerMc :
    The conduction band effective mass of each layer, in unit m0,
    the free space electron mass, len = No. of layers
layerARs :
    Binaries indicating if the layer is active(True) or not(False),
    only affects basis solver, len = No. of layers

ifrDelta :
    The standard deviation of the interface roughness for the interface at
    layer n and layer n+1, in unit angstrom, len = No. of layers.
    Default zero
ifrLambda :
    The correlation length of the interface roughness for the interface at
    layer n and layer n+1, in unit angstrom, len = No. of layers.
    Default zero
avghwLO :
    The average LO phonon energy in unit eV.
    This value is used for LO phonon scattering calculation.
epsrho :
    The effective relative permittivity for LO phonon.
    1/epsrho = 1/(epsilon for high frequency) - 1/(epsilon for static)
    The value is used for LO phonon scattering calculation.

EField :
    External (static) electrical field, in kV/cm = 1e5 V/m
repeats :
    Number of repeat times for the given structure

crystalType :
    Default being "simple", meaning a simple parabolic effective mass is
    used for the calculation. For setting other than "simple",
    `populate_material` should be implemented.

basisAROnly :
    For basis solver if only the Active Region (AR) should be solved.
basisARInjector :
    For basis solver if there should be separator between AR->Injector
basisInjectorAR :
    For basis solver if there should be separator between Injector->AR

solver :
    The solver used for the eigen problem: 'ODE' or 'matrix'.
    By default 'ODE' if C library exists, 'matrix' is a full back.
includeIFR :
    Weather to include IFR scattering for performance estimation.
matrixEigenCount :
    The number of eigen pairs to calculate in the 'matrix' solver.
    It would be very expensive to calculate all of them.

status :
    - 'unsolved' meaning the structure is not solved yet.
    - 'basis' meaning the eigen problem is solved for basis
    - 'solved' meaning the eigen problem is solved.
    - 'solved-full' meaning the population distribution is known.
    """
    xres: float
    Eres: float
    statePerRepeat: int
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
    basisAROnly: bool
    matrixEigenCount: int

    solver: str
    includeIFR: bool

    status: str

    def __init__(self, xres: float = 0.5, Eres: float = 0.5,
                 statePerRepeat: int = 20,
                 layerWidths: List[float] = [10.0],
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

        self.basisAROnly = False
        self.basisARInjector = True
        self.basisInjectorAR = True

        self.status = 'unsolved'

    def layerVc(self, n: int) -> float:
        """The conduction band offset at n-th layer in eV"""
        return self._layerVc[n]

    def layerMc(self, n: int) -> float:
        """The conduction band effective mass at n-th layer, in m0"""
        return self._layerMc[n]

    def rotate_layer(self):
        for layerList in (self.layerWidths, self.layerARs,
                          self._layerVc, self._layerMc,
                          self.ifrDelta, self.ifrLambda):
            layerList.insert(0, layerList.pop())
        self.status = 'unsolved'

    def del_layer(self, n: int):
        for layerList in (self.layerWidths, self.layerARs,
                          self._layerVc, self._layerMc,
                          self.ifrDelta, self.ifrLambda):
            layerList.pop(n)
        self.status = 'unsolved'

    def add_layer(self, n: int, width: float,
                  Vc: float = 0.0, mc: float = 1.0, AR: bool = True,
                  ifrDelta: float = 0.0, ifrLambda: float = 0.0):
        self.layerWidths.insert(n, width)
        self.layerARs.insert(n, AR)
        self._layerVc.insert(n, Vc)
        self._layerMc.insert(n, mc)
        self.ifrDelta.insert(n, ifrDelta)
        self.ifrLambda.insert(n, ifrLambda)
        self.status = 'unsolved'

    def invert_layer(self):
        self.layerWidths = self.layerWidths[::-1]
        self._layerVc = self._layerVc[::-1]
        self._layerMc = self._layerMc[::-1]
        self.layerARs = self.layerARs[::-1]
        self.ifrDelta = self.ifrDelta[:-1:-1] + [self.ifrDelta[-1]]
        self.ifrLambda = self.ifrLambda[:-1:-1] + [self.ifrLambda[-1]]
        self.status = 'unsolved'

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
        self.periodL = layerCumSum[-1]
        self.xPoints = np.arange(0, self.periodL*self.repeats, self.xres)
        N = self.xPoints.size
        if N == 0:
            self.xPoints = np.array([0])
            N = 1
        self.xLayerNums = np.empty(N, dtype=int)
        self.xVc = np.empty(N)
        self.xMc = np.empty(N)

        for n in range(len(self.layerWidths)):
            indices = np.logical_or.reduce([
                (self.xPoints >= layerCumSum[n] + k * self.periodL)
                & (self.xPoints < layerCumSum[n+1] + k * self.periodL)
                for k in range(self.repeats)])
            self.xLayerNums[indices] = n
            self.xVc[indices] = self.layerVc(n)
            self.xMc[indices] = self.layerMc(n)
        if self.crystalType != 'simple':
            self.populate_material()

        self.offset = np.max(self.xVc) - np.min(self.xVc)
        self.xVField = self.xPoints * self.EField * EUNIT
        self.xVc -= self.xVField

        self.Es = np.arange(np.min(self.xVc), np.max(self.xVc), self.Eres/1E3)
        self.matrixEigenCount = self.repeats * self.statePerRepeat
        # solving for eigen states energy close to sigma
        self.matrixSigma = (np.min(self.xVc) + np.max(self.xVc))/2
        self.Eshift = self.periodL * self.EField * EUNIT
        self.status = 'unsolved'

    def xLayerMask(self, n: int) -> np.ndarray:
        """Return the mask for the given layer number `n`.
        A left and right extra point is included for plotting purposes."""
        xSlt = (self.xLayerNums == n)
        xSlt = (xSlt | np.roll(xSlt, 1) | np.roll(xSlt, -1))
        return ~xSlt

    def _shift_psi(self, psi, n) -> np.ndarray:
        if n == 0:
            return psi
        else:
            return np.interp(self.xPoints-self.periodL*n, self.xPoints, psi,
                             left=0, right=0)

    def _basis_shifter(self, ns: List[int], psis0: np.ndarray,
                       eigenEs0: np.ndarray, xPoints: np.ndarray = None
                       ) -> Tuple[np.ndarray, np.ndarray]:
        if xPoints is None:
            xPoints = self.xPoints
        period = self.periodL
        psis = np.empty((0, self.xPoints.size))
        eigenEs = np.empty(0)
        emax = np.max(self.xVc)
        for n in sorted(ns, reverse=True):
            psisn = np.array([np.interp(self.xPoints-period*n, xPoints, psi,
                                        left=0, right=0) for psi in psis0])
            # filter almost zero solutions)
            if len(psisn) > 0:
                es = eigenEs0 - self.Eshift * n
                idx = (np.sum(psisn**2, axis=1)*self.xres > 0.1) & (es < emax)
                psis = np.concatenate((psis, psisn[idx]))
                eigenEs = np.concatenate((eigenEs, es[idx]))
        return psis, eigenEs

    def populate_material(self):
        """This should be overridden to yield bandParams"""
        raise NotImplementedError('Material property is not implemented')

    def _reset_cache(self):
        # TODO: compatibility with shift
        self._cacheDipole = [[None]*len(self.eigenEs)
                             for _ in range(len(self.eigenEs))]
        self._cacheLO = [[None]*len(self.eigenEs)
                         for _ in range(len(self.eigenEs))]
        self.reset_IFR_cache()

    def reset_IFR_cache(self):
        """Reset cache for IFR scattering. This is exposed s.t. eigen solver
        results can be maintained after changing IFR setting."""
        self._cacheIFR = [[None]*len(self.eigenEs)
                          for _ in range(len(self.eigenEs))]
        self._cacheGamma = [[None]*len(self.eigenEs)
                            for _ in range(len(self.eigenEs))]
        if self.status.startswith('solved'):
            self.status = 'solved'

    def solve_whole(self) -> np.ndarray:
        """Solve the whole structure. Will choose from _solve_whole_ode or
        _solve_whole_matrix according to self.solver.

        Yield
        -----
        eigenEs : np.array of float
            the eigenenergy of the layer structure
        psis : np.array of float
            the wave function
        """
        if self.solver == 'ODE':
            self._solve_whole_ode()
        elif self.solver == 'matrix':
            self._solve_whole_matrix()
        else:
            raise NotImplementedError(
                'The {} solver is not implemented'.format(self.solver))
        self._reset_cache()
        self.status = 'solved'
        return self.eigenEs

    def _solve_whole_ode(self) -> np.ndarray:
        if self.crystalType == 'simple':
            self.eigenEs = onedq.cSimpleSolve1D(self.xres, self.Es,
                                                self.xVc, self.xMc)
            self.psis = onedq.cSimpleFillPsi(self.xres, self.eigenEs,
                                             self.xVc, self.xMc)
            return self.eigenEs
        xBand = onedq.Band(self.crystalType, *self.bandParams)
        self.eigenEs = onedq.cBandSolve1D(self.xres, self.Es, self.xVc, xBand)
        self.psis = onedq.cBandFillPsi(
            self.xres, self.eigenEs, self.xVc, xBand)

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
            self.phiso = sqrt(1/3)*xP/(xE+xEg+xESO) * dphic
        return self.eigenEs

    def _solve_whole_matrix(self) -> np.ndarray:
        # unit eV/step^2
        kunit = hbar**2/(2*e0*m0*(1E-10*self.xres)**2)
        if self.crystalType == 'simple':
            # populate mass half grid self.xMc[i] is m at i+-0.5
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
                tol=1E-8)
            self.psis /= sqrt(self.xres)
            self.psis = self.psis.T
            return self.eigenEs
        if self.crystalType == 'ZincBlende':
            self.populate_Kane_matrix()
            # self.eigen_all, self.psi_all = slg.eig_banded(
            #     self.HBanded, select='v', select_range=(Es_low, Es_hi))
            self.eigen_all, self.psi_all = splg.eigsh(
                self.Hsparse, self.matrixEigenCount, sigma=self.matrixSigma,
                tol=1E-8)
            # filter artifacts
            idx = np.abs(self.psi_all[0, :]) < 1E-1
            self.eigen_all = self.eigen_all[idx]
            self.psi_all = self.psi_all[:, idx]
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

    def psi_overlap(self, upper: int, lower: int, shift=0) -> np.ndarray:
        """Return psi[upper] * psi[lower] with psi[lower] shifted by shift
        number of periods."""
        if self.crystalType == 'ZincBlende':
            return sum(
                phi[upper] * self._shift_psi(phi[lower], shift)
                for phi in (self.psis, self.philh, self.phiso))
        # default fallback and crystalType == 'simple'
        return self.psis[upper] * self._shift_psi(self.psis[lower], shift)

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

    def _isBasisBreak(self, n: int) -> bool:
        if self.basisInjectorAR:
            if not self.layerARs[n-1] and self.layerARs[n]:
                return True
        if self.basisARInjector:
            if self.layerARs[n-1] and not self.layerARs[n]:
                return True
        return False

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
        startIdx = []
        endIdx = []
        # Get the region of interest
        if self.basisAROnly:
            if self.layerARs[0]:
                startIdx.append(0)
            for n in range(1, len(self.layerARs)):
                if not self.layerARs[n-1] and self.layerARs[n]:
                    startIdx.append(n)
                if self.layerARs[n-1] and not self.layerARs[n]:
                    endIdx.append(n)
            if self.layerARs[-1]:
                endIdx.append(len(self.layerARs))
        else:
            startIdx.append(-1)
            for n in range(1, len(self.layerARs)):
                if self._isBasisBreak(n):
                    barrier = n if self.layerVc(n) > self.layerVc(n-1) else n-1
                    if barrier != startIdx[-1]:
                        startIdx.append(barrier)
                        endIdx.append(barrier + 1)
            if self._isBasisBreak(0):
                barrier = 0 if self.layerVc(0) > self.layerVc(-1) else -1
                barrier += len(self.layerWidths)
                startIdx.append(barrier)
                endIdx.append(barrier + 1)
            if len(endIdx) == 0:
                startIdx = [0]
                endIdx = [len(self.layerWidths)]
            else:
                startIdx = startIdx[1:]
                endIdx = endIdx[1:] + [endIdx[0] + len(self.layerWidths)]

        self.eigenEs = np.empty((0))
        self.psis = np.empty((0, self.xPoints.size))
        for n in range(0, len(startIdx)):
            dCL = copy.copy(self)
            dCL._reset_for_basis(startIdx[n], endIdx[n])
            dCL.populate_x()
            dCL.solve_whole()

            # map dCL result back to self
            shift = sum(self.layerWidths[:startIdx[n]]) - BASISPAD
            psis, eigenEs = self._basis_shifter(
                range(-1, self.repeats), dCL.psis,
                dCL.eigenEs - shift*self.EField*EUNIT,
                dCL.xPoints + shift)
            self.eigenEs = np.concatenate((self.eigenEs, eigenEs))
            self.psis = np.concatenate((self.psis, psis))
        self._reset_cache()
        self.status = 'basis'
        return self.eigenEs

    def _reset_for_basis(self, start: int, end: int):
        """Reset the parameters for only solving the layers of [start:end].
        This is a helper method for solve_basis"""
        self.repeats = 1
        self.statePerRepeat = (
            self.statePerRepeat * (end - start) // len(self.layerWidths))
        self.layerWidths = (self.layerWidths*2)[start:end]
        self.layerWidths[0] += BASISPAD
        self.layerWidths[-1] += BASISPAD
        self.layerARs = (self.layerARs*2)[start:end]
        self._layerVc = (self._layerVc*2)[start:end]
        self._layerMc = (self._layerMc*2)[start:end]

    def _dipole(self, upper: int, lower: int, shift: int = 0) -> float:
        """Return Electrical dipole between upper and lower states
        in unit angstrom, update self.dipole.
        shift means lower state is shifted by this number of period.
        Should be called for any other related physics quantities."""
        # TODO: clean up self cache for upper -> lower states
        if self.solver == 'ODE':
            psi_u = self.psis[upper, :]
            psi_l = self._shift_psi(self.psis[lower, :], shift)
            Eu = self.eigenEs[upper]
            El = self.eigenEs[lower] - shift * self.Eshift
            xInvMc_u = self._xBandMassInv(Eu)
            xInvMc_l = self._xBandMassInv(El)
            # Eq.(8) in PhysRevB.50.8663, with more precise mass
            # m_j comes from the density of states of transition final state
            z = np.trapz(psi_u * xInvMc_l * np.gradient(psi_l) -
                         np.gradient(psi_u) * xInvMc_u * psi_l)
            z *= hbar**2 / (2*(El-Eu)*e0*m0) / (1E-10)**2  # Angstrom
        elif self.solver == 'matrix' and self.crystalType == 'ZincBlende':
            # In principle this and above is equivalent, and numerically
            # approximately the same as tested in unit tests.
            # The above is remained for legacy reasons.
            z = self.xres * np.trapz(
                self.psi_overlap(upper, lower, shift) * self.xPoints)
        else:
            raise NotImplementedError(
                '{} not implemented for dipole'.format(self.solver))
        return z

    def dipole(self, upper, lower) -> float:
        if self._cacheDipole[upper][lower] is None:
            self._cacheDipole[upper][lower] = self._dipole(upper, lower)
        return self._cacheDipole[upper][lower]

    def _lo_transition(self, upper: int, lower: int, shift: int = 0) -> float:
        # TODO: finite temperature version
        Eu = self.eigenEs[upper]
        El = self.eigenEs[lower] - shift * self.Eshift
        if Eu - El - self.avghwLO < 0:
            # energy difference is smaller than a LO phonon
            # LO phonon scattering doesn't happen
            return INV_INF
        psi_l_sq = self.psi_overlap(lower, lower)
        ml = m0 * np.trapz(self.xMc * psi_l_sq) * self.xres
        kl = sqrt(2 * ml / hbar**2 * (Eu - El - self.avghwLO) * e0)
        N = self.xPoints.size
        if onedq is None:
            convpsi = fft.irfft(np.abs(fft.rfft(
                self.psi_overlap(upper, lower, shift), 2*N))**2)[:N]
            Iij = 2*self.xres**2*np.trapz(
                exp(-kl*self.xPoints*1E-10)*convpsi)
        # C implementation
        else:
            Iij = onedq.OneDSchrodinger.cLOphononScatter(
                self.xres, kl, self.psi_overlap(upper, lower, shift))
        return (ml * e0**2 * self.avghwLO * e0 / hbar * Iij
                / (4 * hbar**2 * self.epsrho * eps0 * kl)) / 1e12  # unit ps^-1

    def lo_transition(self, upper: int, lower: int) -> float:
        """The LO phonon transition lifetime from upper to lower,
        at zero temperature.
        This is using cached results.
        if status is 'solved-full' and the state is a recognized state of a
        period, it's calculated via translation of the wavefunction, otherwise
        it's calculated based on row self.psis
        """
        if self.status == 'solved-full':
            try:
                pu, ushift = self.periodMap[upper]
                pl, lshift = self.periodMap[lower]
                if lshift - ushift not in (0, 1, -1):
                    return INV_INF
                return self._pLO[lshift-ushift][pl][pu]
            except (TypeError, AttributeError):
                # TypeError is when self.periodMap return None
                # AttributeError is when self.periodMap does not exist
                pass
        if self._cacheLO[upper][lower] is None:
            self._cacheLO[upper][lower] = self._lo_transition(
                upper, lower)
        return self._cacheLO[upper][lower]

    def lo_lifetime(self, state: int) -> float:
        """ Return the life time due to LO phonon scattering of the
        given state(label)
        This is using cached results.
        if status is 'solved-full' and the state is a recognized state of a
        period, it's calculated via translation of the wavefunction, otherwise
        it's calculated based on row self.psis
        """
        if self.status == 'solved-full':
            try:
                return 1/sum(np.sum(self._pLO[n][:, self.periodMap[state][0]])
                             for n in range(3))
            except (TypeError, AttributeError):
                # TypeError is when self.periodMap return None
                # AttributeError is when self.periodMap does not exist
                pass
        Ei = self.eigenEs[state]
        if onedq is None:
            return 1/sum(self.lo_transition(state, q) for q in range(state)
                         if self.eigenEs[q] <= Ei - self.avghwLO)
        idxs = self.eigenEs <= Ei - self.avghwLO
        Ejs = self.eigenEs[idxs]
        idxs, = idxs.nonzero()
        psi_js_sq = np.array([self.psi_overlap(idx, idx) for idx in idxs])
        mjs = m0 * np.trapz(self.xMc * psi_js_sq, axis=1) * self.xres
        kls = sqrt(2 * mjs / hbar**2 * (Ei - Ejs - self.avghwLO) * e0)
        fjs = (mjs * e0**2 * self.avghwLO * e0 / hbar
               / (4 * hbar**2 * self.epsrho * eps0 * kls))
        psi_ijs = np.array([self.psi_overlap(state, idx) for idx in idxs])
        Iijtotal = onedq.OneDSchrodinger.cLOtotal(self.xres, kls, psi_ijs, fjs)
        return 1e12 / Iijtotal if Iijtotal > 0 else 1E20

    def _ifr_transition(self, upper: int, lower: int, shift: int = 0
                        ) -> Tuple[float, float]:
        # TODO: finite temperature
        psi_usq = self.psi_overlap(upper, upper)
        psi_lsq = self.psi_overlap(lower, lower)
        psi_ul = self.psi_overlap(upper, lower, shift)
        Eu = self.eigenEs[upper]
        El = self.eigenEs[lower] - shift * self.Eshift
        if Eu < El:
            return INV_INF, -1
        mu = m0 * np.trapz(self.xMc * psi_usq) * self.xres
        ml = m0 * np.trapz(self.xMc * psi_lsq) * self.xres
        kl = sqrt(2 * ml / hbar**2 * (Eu-El) * e0)
        tauInv = 0
        gamma = 0
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

                def interpZ(psi):
                    return (psi[zIdx-1]*(zn-z2) - psi[zIdx]*(zn-z1))/(z2-z1)
                psi_usqz = interpZ(psi_usq)
                psi_lsqz = interpZ(psi_lsq)
                psi_ulz = interpZ(psi_ul)
                scale = pi / hbar**3 * delt**2 * lamb**2 * dU**2
                tauInv += scale * ml * psi_ulz**2 * exp(-lamb**2 * kl**2 / 4)
                gamma += (scale * (psi_usqz - psi_lsqz)
                          * (psi_usqz * mu - psi_lsqz * ml))/2
        return tauInv/1E12, gamma/1E12  # unit ps^-1

    def ifr_transition(self, upper: int, lower: int) -> float:
        r"""Calculate the interface roughness (IFR) transition rate from
        upper to lower state at zero temperature, in unit ps^-1.

        .. math::
            \frac{1}{\tau_{ij}^\text{IFR}} =
            \frac{\pi m^*_j}{\hbar^3} \sum_n
            \Delta_n^2\Lambda_n^2\delta U_n^2
            \left|\psi_i(z_n)\psi_j^*(z_n)\right|^2
            \mathrm e^{- \Lambda^2 m_j^* (E_i - E_j))/2\hbar^2}

        This is using cached results.
        if status is 'solved-full' and the state is a recognized state of a
        period, it's calculated via translation of the wavefunction, otherwise
        it's calculated based on row self.psis
        """
        assert(self.includeIFR)
        if self.status == 'solved-full':
            try:
                pu, ushift = self.periodMap[upper]
                pl, lshift = self.periodMap[lower]
                if lshift - ushift not in (0, 1, -1):
                    return INV_INF
                return self._pIFR[lshift-ushift][pl][pu]
            except (TypeError, AttributeError):
                # TypeError is when self.periodMap return None
                # AttributeError is when self.periodMap does not exist
                pass
        if self._cacheIFR[upper][lower] is not None:
            return self._cacheIFR[upper][lower]
        tauInv, gamma = self._ifr_transition(upper, lower)
        self._cacheIFR[upper][lower] = tauInv
        if gamma > 0:
            self._cacheGamma[upper][lower] = gamma
            self._cacheGamma[lower][upper] = gamma
        return self._cacheIFR[upper][lower]

    def ifr_lifetime(self, state: int) -> float:
        """Return to total life time due to IFR scattering.
        This is using cached results.
        if status is 'solved-full' and the state is a recognized state of a
        period, it's calculated via translation of the wavefunction, otherwise
        it's calculated based on row self.psis
        """
        assert(self.includeIFR)
        if self.status == 'solved-full':
            try:
                return 1/sum(np.sum(self._pIFR[n][:, self.periodMap[state][0]])
                             for n in range(3))
            except (TypeError, AttributeError):
                # TypeError is when self.periodMap return None
                # AttributeError is when self.periodMap does not exist
                pass
        return 1/sum(self.ifr_transition(state, q) for q in range(state))

    def lifetime(self, state: int) -> float:
        """A convenience wrap of return total lifetime of LO and IFR scattering
        or only LO scattering depending on self.includeIFR."""
        if self.status == 'solved-full':
            try:
                return 1/self.decayRates[self.periodMap[state][0]]
            except (TypeError, AttributeError):
                # TypeError is when self.periodMap return None
                # AttributeError is when self.periodMap does not exist
                pass
        if self.includeIFR:
            return 1/(1/self.ifr_lifetime(state) + 1/self.lo_lifetime(state))
        else:
            return self.lo_lifetime(state)

    def ifr_broadening(self, upper: int, lower: int) -> float:
        """Interface roughness induced broadening"""
        # TODO: shift life time
        if self.status == 'solved-full':
            try:
                pu, ushift = self.periodMap[upper]
                pl, lshift = self.periodMap[lower]
                if lshift - ushift not in (0, 1, -1):
                    return 1/INV_INF
                return self._pGamma[lshift-ushift][pl][pu]
            except (TypeError, AttributeError):
                # TypeError is when self.periodMap return None
                # AttributeError is when self.periodMap does not exist
                pass
        if self._cacheGamma[upper][lower] is None:
            self.ifr_transition(upper, lower)
            self.ifr_transition(lower, upper)
        return self._cacheGamma[upper][lower]

    def dephasing(self, upper: int, lower: int) -> float:
        r"""Calculate the broadening gamma of transition between upper ->
        lower transition, return gamma in unit eV as in Lorentzian:

        .. math::
            \mathcal L(\omega) =
            \frac{1}{\pi} \frac{\gamma}{\gamma^2 + (\omega - \omega_0)^2}

        If IFR scattering is included the broadening is calculated dominantly
        from IFR broadening and finite lifetime of upper and lower states.
        Otherwise 0.1 is returned.
        """
        if not self.includeIFR:
            Eu = self.eigenEs[upper]
            El = self.eigenEs[lower]
            de = np.abs(Eu - El)
            return 0.05 * de
        gamma_u = 1/self.ifr_lifetime(upper) + 1/self.lo_lifetime(upper)
        gamma_l = 1/self.ifr_lifetime(lower) + 1/self.lo_lifetime(lower)
        gamma_parallel = self.ifr_broadening(upper, lower)
        # 1E12: ps^-1 -> Hz
        return (gamma_parallel + (gamma_u + gamma_l)/2) * 1E12 * hbar / e0

    def period_recognize(self, tol: float = 5E-5) -> np.ndarray:
        """Pick a set of eigen states as states in a period.

        Return
        ------
        singlePeriodIdx : np.array of int
            These are indices for the recognized states of a single period.

        Yield
        -----
        unBound : set of int
            includes index of states that are not well bounded.
        """
        # TODO: try look from the low energy side
        periodIdx = self.periodL / self.xres
        psisq = np.abs(self.psis)**2
        self.starts = np.argmax(psisq > tol, axis=1)
        self.ends = np.argmax(psisq[:, ::-1] > tol, axis=1)
        self.periodIdx = np.arange(len(self.eigenEs))[
            (self.starts > periodIdx/3) & (self.starts < 4*periodIdx/3)]
        # check if all states are far away enough from boundary
        self.unBound = set()
        barrierBound = np.max(self.xVc + self.xVField) - self.xVField
        for n in self.periodIdx:
            if self.ends[n] < periodIdx / 3:
                self.unBound.add(n)
                if self.eigenEs[n] < barrierBound[-self.ends[n]]:
                    # TODO: use warning package
                    print('State No.{} is close to the right boundary. '
                          'More repeats may be needed'.format(n))
        # print(self.singlePeriodIdx, self.looselyBounded)
        if len(self.periodIdx) - len(self.unBound) == 0:
            raise StateRecognizeError('No well bounded state found. '
                                      'Try increase repeats.')
        return self.periodIdx

    def period_map_build(self, tol: float = 3E-2, etol: float = 1E-3
                         ) -> List[Tuple[int, int]]:
        """Map states to self.singlePeriodIdx, self.periodMap[n] is a tuple of
        (state index in self.singlePeriodIdx, shift of period(s) or
        None meaning it's not mapped."""
        self.periodMap = [None] * len(self.eigenEs)
        singleE = self.eigenEs[self.periodIdx]
        for n, state in enumerate(self.periodIdx):
            self.periodMap[state] = (n, 0)
        for state, en in enumerate(self.eigenEs):
            if self.periodMap[state] is not None:
                continue
            psi = self.psis[state]
            for shift in range(1, self.repeats):
                en_shifted = en + shift * self.Eshift
                psi_shifted = self._shift_psi(psi, -shift)
                for n, sE in enumerate(singleE):
                    if sE < en_shifted - etol:
                        continue
                    if sE > en_shifted + etol:
                        break
                    sState = self.periodIdx[n]
                    wfDiff = np.trapz((psi_shifted - self.psis[sState])**2)
                    wfDiff *= self.xres
                    if wfDiff < tol:
                        # effective an L2 norm here, tested better than
                        # L1 or L-max norm
                        self.periodMap[state] = (n, shift)
                        break
                if self.periodMap[n] is not None:
                    break
        return self.periodMap

    def full_population(self) -> np.ndarray:
        """Calculate the electron full population on states, assuming the
        result of solve_whole and periodRecognize is valid and no state has
        coupling with states two or more periods away.

        Return
        ------
        population : np.array of float, dim = len(periodIdx)
            The population of electrons in the first recognized period, state
            label self.PeriodIdx[n]

        Yield
        -----
        transitions : np.array of float, dim = len(periodIdx)*len(periodIdx)
            The transition rate from self.PeriodIdx[i] to self.PeriodIdx[j]
        decayRates : np.array of float, dim = len(periodIdx)
            Inverse of lifetimes
        flow : float
            The flow of carrier, in unit ps^-1 (carrier density normalize to 1)
        """
        assert(self.status.startswith('solved'))
        idxPeriod = len(self.periodIdx)
        # p for cache for the periodic version
        self._pLO = [np.zeros((idxPeriod, idxPeriod)) for _ in range(3)]
        self._pIFR = [np.zeros((idxPeriod, idxPeriod)) for _ in range(3)]
        self._pGamma = [np.zeros((idxPeriod, idxPeriod)) for _ in range(3)]
        for i, lower in enumerate(self.periodIdx):
            for j, upper in enumerate(self.periodIdx):
                # keep this consistent with loMatrix etc
                for s in ((1, -1) if lower == upper else (1, 0, -1)):
                    self._pLO[s][j][i] = self._lo_transition(lower, upper, s)
                    self._pIFR[s][j][i], gamma = (
                        self._ifr_transition(lower, upper, s))
                    if gamma > 0:
                        self._pGamma[s][j][i] = self._pGamma[s][i][j] = gamma
        forward = self._pLO[1] + self._pIFR[1]
        backward = self._pLO[-1] + self._pIFR[-1]
        internal = self._pLO[0] + self._pIFR[0]
        self.transitions = internal + forward + backward
        self.decayRates = np.sum(self.transitions, axis=0)
        self.population = null_space(
            np.diag(-self.decayRates) + self.transitions)
        if self.population.shape[1] != 1:
            raise ValueError('More than one steady state found. ',
                             self.population.shape[1])
        self.population = self.population[:, 0].T
        self.population /= np.sum(self.population)
        if self.carrierLeak > 5E-2:
            # TODO: use warning package
            print("The structure seems highly leak or more period needed.")
            print("    Most leaking (%.2f%%) unbounded state: %d" %
                  max((self.population[n], state)
                      for n, state in enumerate(self.periodIdx)
                      if state in self.unBound))
        self.flow = np.sum((forward - backward) @ self.population)
        # print("transition: ", self.transitions)
        # print("lifetimes: ", -1/transfer.diagonal())
        # print("population: ", self.population)
        self.status = 'solved-full'
        return self.population

    @property
    def carrierLeak(self) -> float:
        return sum(
            self.population[n] for n, state in enumerate(self.periodIdx)
            if state in self.unBound)

    def state_population(self, state: int) -> float:
        """This method is only valid after fullPopulation has been called"""
        if self.periodMap[state] is None:
            return None
        return self.population[self.periodMap[state][0]]

    def _xBandMassInv(self, energy: float) -> np.ndarray:
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

materials :
    Name of alloys for the heterostructure materials, len >= 2
moleFracs :
    mole fraction for each possible layer material, len = Mp. of materials
xres :
    Position resolution, in Armstrong
Eres :
    Energy resolution, in meV.
    This number being too large may results in miss of some states while
    this being too small will make a long computation time.
    The parameter does not mean the accuracy of the eigen energy. It's
    required for algorithm reasons because of lack of a universal global
    root finding.
statePerRepeat :
    Number of states per repeat, used for calculating matrixEigenCount
wl :
    The wavelength for the design, in unit um, gain spectrum and optimization,
    but doesn't go into quantum solver

layerWidths :
    Width of each layer, in angstrom. len = No. of layers
layerMtrls :
    Label of materials, depending on substrate. len = No. of layers
layerDopings :
    Doping per volume in unit 1e17 cm-3. len = No. of layers
layerARs :
    Binaries indicating if the layer is active(True) or not(False),
    only affects basis solver. len = No. of layers

customIFR :
    Wether to use a customized IFR parameter rather than a material determined
    parameter.
mtrlIFRLambda :
    The interface roughness lambda after materials[n], len = No. of materials
mtrlIFRDelta :
    The interface roughness delta after materials[n], len = No. of materials

EField :
    External (static) electrical field, in kV/cm = 1e5 V/m
repeats :
    Number of repeat times for the given structure
T :
    Temperature of the device, affecting material property

basisAROnly :
    For basis solver if only the Active Region (AR) should be solved.
basisARInjector :
    For basis solver if there should be separator between AR->Injector
basisInjectorAR :
    For basis solver if there should be separator between Injector->AR

solver :
    The solver used for the eigen problem: 'ODE' or 'matrix'.
    By default 'ODE' if C library exists, 'matrix' is a full back.
includeIFR :
    Weather to include IFR scattering for performance estimation.
matrixEigenCount :
    The number of eigen pairs to calculate in the 'matrix' solver.
    It would be very expensive to calculate all of them.

status :
    - 'unsolved' meaning the structure is not solved yet.
    - 'basis' meaning the eigen problem is solved for basis
    - 'solved' meaning the eigen problem is solved.
    - 'solved-full' meaning the population distribution is known.

description :
    Description of the data. For book-keeping purposes.
    """
    materials: List[str]
    moleFracs: List[float]
    wl: float
    layerMtrls: List[int]
    layerDopings: List[float]
    customIFR: bool
    mtrlIFRLambda: List[float]
    mtrlIFRDelta: List[float]
    T: float
    description: str

    def __init__(self, substrate="InP", materials=["InGaAs", "AlInAs"],
                 moleFracs=[0.53, 0.52], xres=0.5, Eres=0.5, statePerRepeat=20,
                 layerWidths=[10.0], layerMtrls=None, layerDopings=None,
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
        self.temperature = T
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
        super().__init__(xres=xres, Eres=Eres, statePerRepeat=statePerRepeat,
                         layerWidths=layerWidths, layerARs=layerARs,
                         ifrDelta=ifrDelta, ifrLambda=ifrLambda,
                         EField=EField, repeats=repeats)
        self.crystalType = Material.MParam[substrate]["Crystal"]
        self.subM = Material.Material(self.substrate, self.temperature)
        self.wl = wl
        self.solver = solver
        if onedq is None:
            self.solver = 'matrix'
        self.update_mtrls()

    def __copy__(self):
        return QCLayers(
            substrate=self.substrate, materials=copy.copy(self.materials),
            moleFracs=copy.copy(self.moleFracs), xres=self.xres,
            Eres=self.Eres, statePerRepeat=self.statePerRepeat,
            layerWidths=copy.copy(self.layerWidths),
            layerMtrls=copy.copy(self.layerMtrls),
            layerDopings=copy.copy(self.layerDopings),
            customIFR=self.customIFR,
            mtrlIFRLambda=copy.copy(self.mtrlIFRLambda),
            mtrlIFRDelta=copy.copy(self.mtrlIFRDelta),
            ifrDelta=copy.copy(self.ifrDelta),
            ifrLambda=copy.copy(self.ifrLambda),
            layerARs=copy.copy(self.layerARs),
            EField=self.EField, repeats=self.repeats, T=self.temperature,
            solver=self.solver, description=self.description, wl=self.wl
            )

    def _get_IFRList(self) -> Tuple[List[float], List[float]]:
        """Get IFR parameters for SchrodingerLayer. Should be called
        every time the material list changes."""
        assert(not self.customIFR)
        if self.mtrlIFRDelta is not None:
            ifrDelta = [self.mtrlIFRDelta[m] for m in self.layerMtrls]
        else:
            ifrDelta = [0.0] * len(self.layerMtrls)
        if self.mtrlIFRLambda is not None:
            ifrLambda = [self.mtrlIFRLambda[m] for m in self.layerMtrls]
        else:
            ifrLambda = [0.0] * len(self.layerMtrls)
        return ifrDelta, ifrLambda

    def update_mtrls(self):
        """Update properties for the materials.
        This should be called every time the material parameters are changed
        directly via `materials`, `moleFracs` or `temperature`.

        Yield
        -----
        a_parallel : float
            The parallel crystal constant of the structure, determined by the
            substrate material.

        mtrlAlloys : List[Material.Alloy]
            list of the Alloy objects for processing material properties.
        """
        self.a_parallel = self.subM.param['alc']
        self.mtrlAlloys = [Material.Alloy(self.materials[idx],
                                          self.moleFracs[idx],
                                          self.temperature)
                           for idx in range(len(self.materials))]
        for al in self.mtrlAlloys:
            al.set_strain(self.a_parallel)

    def set_mtrl(self, n: int, mtrl: str = None, moleFrac: float = None):
        """Set material[n] to new material (mtrl) and/or moleFrac"""
        if mtrl is None and moleFrac is None:
            raise Exception("Nothing changed")
        if mtrl is None:
            mtrl = self.materials[n]
        if moleFrac is None:
            moleFrac = self.moleFracs[n]
        self.moleFracs[n] = moleFrac
        self.materials[n] = mtrl
        self.mtrlAlloys[n] = Material.Alloy(mtrl, moleFrac, self.temperature)
        self.mtrlAlloys[n].set_strain(self.a_parallel)

    def add_mtrl(self, mtrl: str = None, moleFrac: float = None,
                 IFRLambda: float = None, IFRDelta: float = None):
        """Add a new material possibility"""
        self.materials.append(mtrl if mtrl else
                              QCMaterial[self.substrate][0])
        self.moleFracs.append(moleFrac if moleFrac else 0.0)
        self.mtrlIFRLambda.append(IFRLambda if IFRLambda else
                                  self.mtrlIFRLambda[-1])
        self.mtrlIFRDelta.append(IFRDelta if IFRDelta else
                                 self.mtrlIFRDelta[-1])
        self.mtrlAlloys.append(Material.Alloy(
            self.materials[-1], self.moleFracs[-1], self.temperature))
        self.mtrlAlloys[-1].set_strain(self.a_parallel)

    def del_mtrl(self, n: int):
        """Delete materials labeled n.
        All layers of this material will become previous
        `materials[n-1 if n >0 else 1]`.  There should be at least two
        materials otherwise there will be error."""
        if len(self.materials) <= 2:
            raise ValueError("There should be at least 2 materials")
        self.materials.pop(n)
        self.moleFracs.pop(n)
        for i in range(len(self.layerMtrls)):
            if self.layerMtrls[i] >= n:
                self.layerMtrls[i] = (self.layerMtrls[i] - 1
                                      if self.layerMtrls[i] > 0
                                      else 0)

    def add_layer(self, n: int, width: int, mtrlIdx: int,
                  AR: bool, doping: float):
        self.layerMtrls.insert(n, mtrlIdx)
        self.layerDopings.insert(n, doping)
        super().add_layer(n, width, AR=AR)
        if not self.customIFR:
            self.ifrDelta, self.ifrLambda = self._get_IFRList()

    def rotate_layer(self):
        super().rotate_layer()
        for layerList in (self.layerMtrls, self.layerDopings):
            layerList.insert(0, layerList.pop())

    def del_layer(self, n: int):
        super().del_layer(n)
        for layerList in (self.layerMtrls, self.layerDopings):
            layerList.pop(n)

    def invert_layer(self):
        super().invert_layer()
        self.layerMtrls = self.layerMtrls[::-1]
        self.layerDopings = self.layerDopings[::-1]
        if not self.customIFR:
            self.ifrDelta, self.ifrLambda = self._get_IFRList()

    def set_substrate(self, subs: str):
        if subs in QCMaterial:
            self.substrate = subs
            self.crystalType = Material.MParam[subs]["Crystal"]
            matlN = len(self.materials)
            self.materials = (QCMaterial[subs]*matlN)[0:matlN]
            self.update_mtrls()
        else:
            raise TypeError("Substrate %s not supported" % subs)

    def set_temperature(self, T: float):
        self.temperature = T
        self.subM.set_temperature(T)
        self.update_mtrls()

    @property
    def mtrlOffset(self) -> float:
        """Return the conduction band offset (difference between highest
        conduction band and lowest conduction band energy) of materials,
        in unit eV"""
        ecgs = [alloy.param['EcG'] for alloy in self.mtrlAlloys]
        return max(ecgs) - min(ecgs)

    @property
    def netStrain(self) -> float:
        """Return average strain perpendicular to the layer plane, in
        percentage."""
        if sum(self.layerWidths) <= 1e-5:
            return -1
        totalStrain = sum(self.mtrlAlloys[self.layerMtrls[n]].eps_perp
                          * self.layerWidths[n]
                          for n in range(len(self.layerWidths)))
        return 100 * totalStrain / sum(self.layerWidths)

    def layerVc(self, n: int) -> float:
        return self.mtrlAlloys[self.layerMtrls[n]].param['EcG']

    def layerMc(self, n: int) -> float:
        return self.mtrlAlloys[self.layerMtrls[n]].param['me0']

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
                    p[indices] = self.mtrlAlloys[self.layerMtrls[n]].param[key]
            self.bandParams = (xEg, xF, xEp, xESO)

            ExtField = self.xPoints * self.EField * EUNIT
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
            sumhwlo = sum(self.mtrlAlloys[self.layerMtrls[n]].param['hwLO']
                          * self.layerWidths[n]
                          for n in range(len(self.layerWidths)))
            self.avghwLO = sumhwlo / sum(self.layerWidths)
            epsInf = np.array([a.param["epsInf"] for a in self.mtrlAlloys])
            epss = np.array([a.param["epss"] for a in self.mtrlAlloys])
            epsrho = 1 / (1/epsInf - 1/epss)
            self.epsrho = (np.sum(epsrho[self.layerMtrls] * self.layerWidths)
                           / sum(self.layerWidths))
        # IFR
        if not self.customIFR:
            self.ifrDelta, self.ifrLambda = self._get_IFRList()

    def _reset_for_basis(self, start: int, end: int):
        super()._reset_for_basis(start, end)
        self.layerMtrls = (self.layerMtrls*2)[start:end]
        self.layerDopings = (self.layerDopings*2)[start:end]

    def figure_of_merit(self, upper: int, lower: int) -> float:
        """Calculate Figure of Merit.
        This function must be called after solving for wave functions

        Parameters
        ----------
        upper, lower :
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
        self.tauLO_ul = 1/self.lo_transition(upper, lower)
        self.tauLO_l = self.lo_lifetime(lower)
        self.tauLO_u = self.lo_lifetime(upper)
        if self.includeIFR:
            self.tauIFR_ul = 1/self.ifr_transition(upper, lower)
            self.tauIFR_u = self.ifr_lifetime(upper)
            self.tauIFR_l = self.ifr_lifetime(lower)
            self.tau_u = 1/(1/self.tauLO_u + 1/self.tauIFR_u)
            self.tau_l = 1/(1/self.tauLO_l + 1/self.tauIFR_l)
            self.tau_ul = 1/(1/self.tauLO_ul + 1/self.tauIFR_ul)
        else:
            self.tau_u = self.tauLO_u
            self.tau_l = self.tauLO_l
            self.tau_ul = self.tauLO_ul
        return self.dipole(upper, lower)**2 * self.tau_u * (
            1 - self.tau_l / self.tau_ul)

    def effective_ridx(self, wl: ScalerOrArray) -> ScalerOrArray:
        """Return the effective refractive index for TM mode"""
        if sum(self.layerWidths) == 0:
            return 1.0
        self.mtrlRIdx = [(m.moleFrac * Material.rIdx[m.A.name](wl) +
                          (1 - m.moleFrac) * Material.rIdx[m.B.name](wl))
                         for m in self.mtrlAlloys]
        self.layerRIdx = np.array([self.mtrlRIdx[n] for n in self.layerMtrls])
        neff = np.average(1/self.layerRIdx**2, axis=0,
                          weights=self.layerWidths)**(-1/2)
        return neff

    @property
    def sheet_density(self) -> float:
        """Return the sheet density of doping per period, in unit cm^-2"""
        # 1E9 -> 1E17 cm^-3 * Angstrom -> cm^-2
        return sum(self.layerDopings[n] * self.layerWidths[n]
                   for n in range(0, len(self.layerWidths))) * 1E9

    def full_population(self) -> np.array:
        """Apart from SchrodingerLayer.full_population, this also yield
        current density, with knowledge of doping"""
        res = super().full_population()
        # ps^-1 -> kA/cm^2, where 1E9 = 1E12 * 1E-3, 1E12 is (ps^-1 -> s^-1)
        self.current = self.flow * self.sheet_density * e0 * 1E9
        return res

    def full_gain_spectrum(self, wl: ScalerOrArray = None) -> ScalerOrArray:
        """Perform fully automatic calculation for the gain on wavelength(s).
        """
        if wl is None:
            wl = self.wl
        neff = self.effective_ridx(wl)
        de0 = h * c0 / (wl * 1E-6) / e0
        gain = 0
        # self.gainul = {}
        for i in range(len(self.periodIdx)):
            upper = self.periodIdx[i]
            for j in range(i+1, len(self.periodIdx)):
                lower = self.periodIdx[j]
                for shift in (-1, 0, 1):
                    dipole = self._dipole(upper, lower, shift) * 1E-8  # to cm
                    Eu = self.eigenEs[upper]
                    El = self.eigenEs[lower] - shift * self.Eshift
                    de = Eu - El
                    dpop = self.population[i] - self.population[j]
                    if de < 0:
                        de, dpop = -de, -dpop
                    if self.includeIFR:
                        gamma_para = self._pGamma[shift][j][i]
                    else:
                        gamma_para = self.dephasing(upper, lower)
                        gamma_para /= (1E12 * hbar/e0)  # unit to Hz
                    gamma = (gamma_para +
                             (self.decayRates[i] + self.decayRates[j])/2
                             ) * 1E12 * hbar/e0
                    gainul = dpop * dipole**2 * gamma/(gamma**2 + (de-de0)**2)
                    # if np.max(np.abs(gainul)) > 4E-13:
                    #     self.gainul[(upper, lower)] = gainul
                    #     print(dpop * dipole**2, h*c0/(de*e0)*1E6, gamma/de,
                    #           upper, lower, shift)
                    gain = gain + gainul
        # e0^2 / (hbar * c0 * eps0) is dimension 1.. 1E-8 makes unit cm^-1
        gain *= e0**2*de0*self.sheet_density / (
                hbar*neff*c0*eps0*self.periodL*1E-8)
        return gain


def optimize_layer(qcl: QCLayers, n: int, upper: int, lower: int,
                   iter: int = 50):
    """Optimize FoM*Lorentzian for n-th layer thickness, assuming the state
    index does not change. optimization is performed by searching on the
    position resolution steps.

    Warning: This cannot specify a correct state if there are states
    index crossing.

    TODO: Use period recognizer to improve the algorithm
    """
    Eu = qcl.eigenEs[upper]
    El = qcl.eigenEs[lower]
    if Eu < El:
        upper, lower = lower, upper
        Eu, El = El, Eu
    # 1E-6 um -> m... in unit eV
    w0 = h * c0 / (qcl.wl * 1E-6) / e0

    def reduceFoM():
        wul = Eu - El
        gamma = qcl.dephasing(upper, lower)
        return qcl.figure_of_merit(upper, lower) * gamma * w0/(
            gamma**2 + (wul - w0)**2)
    width = round(qcl.layerWidths[n] / qcl.xres) * qcl.xres
    FoMnow = reduceFoM()
    print(("Start Optimizing Layer NO %d " % n) +
          ("for FoM between state %d and %d.\n" % (upper, lower)) +
          ("\tStart at width=%.1f, FoM=%.5g" % (width, FoMnow)))

    def newFoM(newWidth):
        qcl.layerWidths[n] = newWidth
        qcl.populate_x()
        qcl.solve_whole()
        qcl.dipole(upper, lower)
        return reduceFoM()
    FoMminus = newFoM(width - qcl.xres)
    FoMplus = newFoM(width + qcl.xres)
    for _ in range(iter):
        if FoMnow < FoMplus:
            FoMminus = FoMnow
            FoMnow = FoMplus
            width += qcl.xres
            FoMplus = newFoM(width + qcl.xres)
        elif FoMnow < FoMminus:
            FoMplus = FoMnow
            FoMnow = FoMminus
            width -= qcl.xres
            FoMminus = newFoM(width - qcl.xres)
        else:
            print("Maximum iteration reached.")
            break
        print("\twidth=%.1f, FoM=%.5g" % (width, FoMnow))
    qcl.layerWidths[n] = width
    print("finished, width=%.1f, FoM=%.5g" % (width, FoMnow))
    return FoMnow


def auto_gain(qcl: QCLayers, wls: ScalerOrArray = None):
    """Perform automatic gain calculation from newly loaded a qcl object.

    This is equivalent to

    .. code-block:: python

        qcl.populate_x()
        qcl.solve_whole()
        qcl.period_recognize()
        qcl.full_population()
        result = qcl.full_gain_spectrum(wls)
    """
    qcl.populate_x()
    qcl.solve_whole()
    qcl.period_recognize()
    qcl.full_population()
    return qcl.full_gain_spectrum(wls)


def optimize_global(qcl: QCLayers, iter: int = 50):
    """A global optimization using gradient descent."""
    layerNum = len(qcl.layerWidths)
    now = auto_gain(qcl)
    for n in range(iter):
        changed = False
        print("iteration {}:".format(n))
        print("    initial gain: {}".format(now))
        for n in range(layerNum):
            w = qcl.layerWidths[n]
            maxDiff = min(1, round(0.1*w/qcl.xres))
            qcl.layerWidths[n] = w - qcl.xres
            newMinus = auto_gain(qcl)
            qcl.layerWidths[n] = w + qcl.xres
            newPlus = auto_gain(qcl)
            print(("   layer No.{}.. Width w={}, "
                   "gain (0:{}, +:{}, -:{})").format(
                   n, w, now, newPlus, newMinus))
            if not (now > newMinus and now > newPlus):
                changed = True
                if now < newPlus and now < newMinus:
                    dw = 1 if newPlus > newMinus else -1
                else:
                    dif = (newPlus - newMinus)/2
                    ddif = newPlus + newMinus - 2*now
                    dw = -dif/ddif
                    if ddif > 0 or abs(dw) > maxDiff:
                        dw = maxDiff if dif > 0 else -maxDiff
                    elif round(dw) == 0:
                        dw = 1 if dw > 0 else -1
                    else:
                        dw = round(dw)
                qcl.layerWidths[n] = w + dw*qcl.xres
                print("    new width:", qcl.layerWidths[n])
                if dw == 1:
                    now = newPlus
                elif dw == -1:
                    now = newMinus
                else:
                    now = qcl._auto_gain()
        # TODO: optimize field
        if not changed:
            break
    print('Finished')
