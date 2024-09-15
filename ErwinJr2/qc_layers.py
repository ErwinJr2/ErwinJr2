"""Defines the classes describing layerd structure for quantum cascade lasers."""

import copy
import logging
from typing import List, Tuple, Union

import numpy as np
import scipy.sparse as sparse
import scipy.sparse.linalg as splg
from numpy import exp, fft, pi, sqrt
from scipy.constants import (
    c as c0,
    e as e0,
    electron_mass as m0,
    epsilon_0 as eps0,
    h,
    hbar,
)
from scipy.linalg import null_space

from ErwinJr2 import material
from ErwinJr2.oned_quantum import c_schrodinger as onedq

EUNIT = 1e-5  # E field unit from kV/cm to V/Angstrom
BASISPAD = 100  #: padding barrier for basis solver, unit Angstrom
INV_INF = 1e-20  #: for infinite small decay rate (ns-1)
# used for typing as either an array or a float number
ScalerOrArray = Union[float, np.ndarray]

QC_MATERIAL = {
    "InP": ["InGaAs", "AlInAs"],
    "GaAs": ["AlGaAs"],
    "GaSb": ["InAsSb", "AlGaSb"],
}  #: supported substrate/material set

_logger = logging.getLogger(__name__)


class StateRecognizeError(Exception):
    """Raised when QClayers cannot recognize a period of states.
    This typically means the number of periods is too small, or for single
    period structure the barrier at the end of the structure is not long
    enough to bound a state.

    Attributes:
        expression -- input expression in which the error occurred."""

    def __init__(self, expression=""):
        self.expression = expression


class SchrodingerLayer(object):
    """Class for layer structure for Schrodinger solver.

        This is used as the base class of :class:`.QCLayers` for separation of
        material property and the solver.

    Parameters
    ----------
    x_res :
        Position resolution, in Armstrong
    e_res :
        Energy resolution, in meV.
        This number being too large may results in miss of some states while
        this being too small will make a long computation time.
        The parameter does not mean the accuracy of the eigen energy. It's
        required for algorithm reasons because of lack of a universal global
        root finding.
    state_per_repeat :
        Number of states per repeat, used for calculating matrix_eigen_count

    layer_widths :
        Width of each layer, in angstrom. len = No. of layers
    layerVc :
        The conduction band offset of each layer, in unit eV, len = No. of layers
    layerMc :
        The conduction band effective mass of each layer, in unit m0,
        the free space electron mass, len = No. of layers
    layer_ar :
        Binaries indicating if the layer is active(True) or not(False),
        only affects basis solver, len = No. of layers

    ifr_deltas :
        The standard deviation of the interface roughness for the interface at
        layer n and layer n+1, in unit angstrom, len = No. of layers.
        Default zero
    ifr_lambdas :
        The correlation length of the interface roughness for the interface at
        layer n and layer n+1, in unit angstrom, len = No. of layers.
        Default zero
    avg_hw_lo :
        The average LO phonon energy in unit eV.
        This value is used for LO phonon scattering calculation.
    eps_rho :
        The effective relative permittivity for LO phonon.
        1/eps_rho = 1/(epsilon for high frequency) - 1/(epsilon for static)
        The value is used for LO phonon scattering calculation.

    e_field :
        External (static) electrical field, in kV/cm = 1e5 V/m
    repeats :
        Number of repeat times for the given structure

    crystalType :
        Default being "simple", meaning a simple parabolic effective mass is
        used for the calculation. For setting other than "simple",
        `populate_material` should be implemented.

    basis_ar_only :
        For basis solver if only the Active Region (AR) should be solved.
    basis_ar_injector :
        For basis solver if there should be separator between AR->Injector
    basis_injector_ar :
        For basis solver if there should be separator between Injector->AR

    solver :
        The solver used for the eigen problem: 'ODE' or 'matrix'.
        By default 'ODE' if C library exists, 'matrix' is a full back.
    include_ifr :
        Weather to include IFR scattering for performance estimation.
    matrix_eigen_count :
        The number of eigen pairs to calculate in the 'matrix' solver.
        It would be very expensive to calculate all of them.

    status :
        - 'unsolved' meaning the structure is not solved yet.
        - 'basis' meaning the eigen problem is solved for basis
        - 'solved' meaning the eigen problem is solved.
        - 'solved-full' meaning the population distribution is known.
    """

    x_step: float
    e_step: float
    state_per_repeat: int
    layer_widths: List[float]
    _layer_vc: List[float]
    _layer_mc: List[float]
    layer_ar: List[bool]
    ifr_deltas: List[float]
    ifr_lambdas: List[float]
    avg_hw_lo: float
    eps_rho: float
    e_field: float
    repeats: int
    crystal_type: str
    basis_ar_injector: bool
    basis_injector_ar: bool
    basis_ar_only: bool
    matrix_eigen_count: int

    solver: str
    include_ifr: bool

    status: str

    def __init__(
        self,
        x_res: float = 0.5,
        e_res: float = 0.5,
        state_per_repeat: int = 20,
        layer_widths: List[float] = None,
        layer_ar: List[bool] = None,
        layer_vc: List[float] = None,
        layer_mc: List[float] = None,
        ifr_deltas: List[float] = None,
        ifr_lambdas: List[float] = None,
        avg_hw_lo: float = 35e-3,
        eps_rho: float = 69.0,
        e_field: float = 0.0,
        repeats: int = 1,
    ):
        if layer_widths is None:
            layer_widths = [10.0]
        self.x_step = x_res
        self.e_step = e_res
        self.state_per_repeat = state_per_repeat
        assert isinstance(layer_widths, list)
        self.layer_widths = layer_widths
        length = len(layer_widths)
        self._layer_vc = layer_vc if layer_vc is not None else [0.0] * length
        self._layer_mc = layer_mc if layer_mc is not None else [1.0] * length
        self.layer_ar = layer_ar if layer_ar is not None else [True] * length
        self.ifr_deltas = ifr_deltas if ifr_deltas is not None else [0.0] * length
        self.ifr_lambdas = ifr_lambdas if ifr_lambdas is not None else [0.0] * length
        self.e_field = e_field
        self.repeats = repeats
        self.avg_hw_lo = avg_hw_lo
        self.eps_rho = eps_rho

        self.crystal_type = "simple"
        self.solver = "ODE"
        self.include_ifr = False
        self.matrix_eigen_count = state_per_repeat * repeats

        self.basis_ar_only = False
        self.basis_ar_injector = True
        self.basis_injector_ar = True

        self.status = "unsolved"

    def layer_vc(self, n: int) -> float:
        """The conduction band offset at n-th layer in eV"""
        return self._layer_vc[n]

    def layer_mc(self, n: int) -> float:
        """The conduction band effective mass at n-th layer, in m0"""
        return self._layer_mc[n]

    def rotate_layer(self):
        for layer_list in (
            self.layer_widths,
            self.layer_ar,
            self._layer_vc,
            self._layer_mc,
            self.ifr_deltas,
            self.ifr_lambdas,
        ):
            layer_list.insert(0, layer_list.pop())
        self.status = "unsolved"

    def del_layer(self, n: int):
        for layer_list in (
            self.layer_widths,
            self.layer_ar,
            self._layer_vc,
            self._layer_mc,
            self.ifr_deltas,
            self.ifr_lambdas,
        ):
            layer_list.pop(n)
        self.status = "unsolved"

    def add_layer(
        self,
        n: int,
        width: float,
        vc: float = 0.0,
        mc: float = 1.0,
        is_ar: bool = True,
        ifr_delta: float = 0.0,
        ifr_lambda: float = 0.0,
    ):
        self.layer_widths.insert(n, width)
        self.layer_ar.insert(n, is_ar)
        self._layer_vc.insert(n, vc)
        self._layer_mc.insert(n, mc)
        self.ifr_deltas.insert(n, ifr_delta)
        self.ifr_lambdas.insert(n, ifr_lambda)
        self.status = "unsolved"

    def invert_layer(self):
        self.layer_widths = self.layer_widths[::-1]
        self._layer_vc = self._layer_vc[::-1]
        self._layer_mc = self._layer_mc[::-1]
        self.layer_ar = self.layer_ar[::-1]
        self.ifr_deltas = self.ifr_deltas[:-1:-1] + [self.ifr_deltas[-1]]
        self.ifr_lambdas = self.ifr_lambdas[:-1:-1] + [self.ifr_lambdas[-1]]
        self.status = "unsolved"

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
        layer_cumsum = [0] + np.cumsum(self.layer_widths).tolist()
        self.period_l = layer_cumsum[-1]
        self.x_points = np.arange(0, self.period_l * self.repeats, self.x_step)
        length = self.x_points.size
        if length == 0:
            self.x_points = np.array([0])
            length = 1
        self.x_layer_nums = np.empty(length, dtype=int)
        self.x_vc = np.empty(length)
        self.x_mc = np.empty(length)

        for n in range(len(self.layer_widths)):
            indices = np.logical_or.reduce(
                [
                    (self.x_points >= layer_cumsum[n] + k * self.period_l)
                    & (self.x_points < layer_cumsum[n + 1] + k * self.period_l)
                    for k in range(self.repeats)
                ]
            )
            self.x_layer_nums[indices] = n
            self.x_vc[indices] = self.layer_vc(n)
            self.x_mc[indices] = self.layer_mc(n)
        if self.crystal_type != "simple":
            self.populate_material()

        self.offset = np.max(self.x_vc) - np.min(self.x_vc)
        self.x_v_field = self.x_points * self.e_field * EUNIT
        self.x_vc -= self.x_v_field

        self.es = np.arange(np.min(self.x_vc), np.max(self.x_vc), self.e_step / 1e3)
        self.matrix_eigen_count = self.repeats * self.state_per_repeat
        # solving for eigen states energy close to sigma
        self.matrix_sigma = (np.min(self.x_vc) + np.max(self.x_vc)) / 2
        self.e_shift = self.period_l * self.e_field * EUNIT
        self.status = "unsolved"

    def x_layer_mask(self, n: int) -> np.ndarray:
        """Return the mask for the given layer number `n`.
        A left and right extra point is included for plotting purposes."""
        x_slt = self.x_layer_nums == n
        x_slt = x_slt | np.roll(x_slt, 1) | np.roll(x_slt, -1)
        return ~x_slt

    def _shift_psi(self, psi, n) -> np.ndarray:
        if n == 0:
            return psi
        else:
            return np.interp(
                self.x_points - self.period_l * n, self.x_points, psi, left=0, right=0
            )

    def _basis_shifter(
        self,
        ns: List[int],
        psis0: np.ndarray,
        eigen_es0: np.ndarray,
        x_points: np.ndarray = None,
    ) -> Tuple[np.ndarray, np.ndarray]:
        if x_points is None:
            x_points = self.x_points
        period = self.period_l
        psis = np.empty((0, self.x_points.size))
        eigen_es = np.empty(0)
        emax = np.max(self.x_vc)
        for n in sorted(ns, reverse=True):
            psisn = np.array(
                [
                    np.interp(
                        self.x_points - period * n, x_points, psi, left=0, right=0
                    )
                    for psi in psis0
                ]
            )
            # filter almost zero solutions)
            if len(psisn) > 0:
                es = eigen_es0 - self.e_shift * n
                idx = (np.sum(psisn**2, axis=1) * self.x_step > 0.1) & (es < emax)
                psis = np.concatenate((psis, psisn[idx]))
                eigen_es = np.concatenate((eigen_es, es[idx]))
        return psis, eigen_es

    def populate_material(self):
        """This should be overridden to yield band_params"""
        raise NotImplementedError("Material property is not implemented")

    def _reset_cache(self):
        # TODO: compatibility with shift
        self._cache_dipole = [
            [None] * len(self.eigen_es) for _ in range(len(self.eigen_es))
        ]
        self._cache_lo = [
            [None] * len(self.eigen_es) for _ in range(len(self.eigen_es))
        ]
        self.reset_ifr_cache()

    def reset_ifr_cache(self):
        """Reset cache for IFR scattering. This is exposed s.t. eigen solver
        results can be maintained after changing IFR setting."""
        self._cache_ifr = [
            [None] * len(self.eigen_es) for _ in range(len(self.eigen_es))
        ]
        self._cache_gamma = [
            [None] * len(self.eigen_es) for _ in range(len(self.eigen_es))
        ]
        if self.status.startswith("solved"):
            self.status = "solved"

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
        if self.solver == "ODE":
            self._solve_whole_ode()
        elif self.solver == "matrix":
            self._solve_whole_matrix()
        else:
            raise NotImplementedError(f"The {self.solver} solver is not implemented")
        self._reset_cache()
        self.status = "solved"
        return self.eigen_es

    def _solve_whole_ode(self) -> np.ndarray:
        if self.crystal_type == "simple":
            self.eigen_es = onedq.cSimpleSolve1D(
                self.x_step, self.es, self.x_vc, self.x_mc
            )
            self.psis = onedq.cSimpleFillPsi(
                self.x_step, self.eigen_es, self.x_vc, self.x_mc
            )
            return self.eigen_es
        x_band = onedq.PyBand(self.crystal_type, *self.band_params)
        self.eigen_es = onedq.cBandSolve1D(self.x_step, self.es, self.x_vc, x_band)
        self.psis = onedq.cBandFillPsi(self.x_step, self.eigen_es, self.x_vc, x_band)

        if self.crystal_type == "ZincBlende":
            # To restore lh and so band to keep consistent with matrix solver
            x_eg, _, x_ep, x_eso = self.band_params
            kunit = hbar**2 / (2 * e0 * m0 * (1e-10 * self.x_step) ** 2)
            x_p = sqrt(x_ep * kunit)
            dphic = np.zeros(self.psis.shape)
            dphic[:, 1:-1] = (self.psis[:, 2:] - self.psis[:, :-2]) / 2
            self.philh = np.zeros(self.psis.shape)
            x_e = np.broadcast_to(self.eigen_es.reshape(-1, 1), self.psis.shape)
            x_e = x_e - self.x_vc
            self.philh = -sqrt(2 / 3) * x_p / (x_e + x_eg) * dphic
            self.phiso = sqrt(1 / 3) * x_p / (x_e + x_eg + x_eso) * dphic
        return self.eigen_es

    def _solve_whole_matrix(self) -> np.ndarray:
        # unit eV/step^2
        kunit = hbar**2 / (2 * e0 * m0 * (1e-10 * self.x_step) ** 2)
        if self.crystal_type == "simple":
            # populate mass half grid self.xMc[i] is m at i+-0.5
            layer_cumsum = [0] + np.cumsum(self.layer_widths).tolist()
            period_l = layer_cumsum[-1]
            self.x_mc_plus = np.empty(self.x_points.size)
            self.x_mc_minus = np.empty(self.x_points.size)
            for n in range(len(self.layer_widths)):
                indices = np.logical_or.reduce(
                    [
                        (
                            self.x_points
                            >= layer_cumsum[n] + k * period_l - self.x_step / 2
                        )
                        & (
                            self.x_points
                            < layer_cumsum[n + 1] + k * period_l - self.x_step / 2
                        )
                        for k in range(self.repeats)
                    ]
                )
                self.x_mc_plus[indices] = self.layer_mc(n)
            self.x_mc_plus[-1] = self.x_mc_plus[-2]
            self.x_mc_minus[1:] = self.x_mc_plus[:-1]
            self.x_mc_minus[0] = self.x_mc_minus[1]

            # diagonal and sub-diagonal of Hamiltonian
            self.h_diag = kunit * (1 / self.x_mc_plus + 1 / self.x_mc_minus) + self.x_vc
            self.h_subd = -kunit / self.x_mc_plus[:-1]
            # self.eigenEs, self.psis = slg.eigh_tridiagonal(
            #     self.Hdiag, self.Hsubd, select='v',
            #     select_range=(np.min(self.xVc), np.max(self.xVc)))
            self.h_subd = -kunit / self.x_mc_plus
            length = len(self.x_points)
            self.h_sparse = sparse.diags(
                [self.h_subd, self.h_diag, self.h_subd],
                [-1, 0, 1],
                shape=(length, length),
            )
            self.eigen_es, self.psis = splg.eigsh(
                self.h_sparse,
                self.matrix_eigen_count,
                sigma=self.matrix_sigma,
                tol=1e-8,
            )
            self.psis /= sqrt(self.x_step)
            self.psis = self.psis.T
            return self.eigen_es
        if self.crystal_type == "ZincBlende":
            self.populate_kane_matrix()
            # self.eigen_all, self.psi_all = slg.eig_banded(
            #     self.HBanded, select='v', select_range=(Es_low, Es_hi))
            self.eigen_all, self.psi_all = splg.eigsh(
                self.h_sparse,
                self.matrix_eigen_count,
                sigma=self.matrix_sigma,
                tol=1e-8,
            )
            # filter artifacts
            idx = np.abs(self.psi_all[0, :]) < 1e-1
            self.eigen_all = self.eigen_all[idx]
            self.psi_all = self.psi_all[:, idx]
            # normalization should be sum(self.psi_all**2)*self.x_res = 1
            self.psi_all /= sqrt(self.x_step)
            for n in range(len(self.eigen_all)):
                phic = self.psi_all[::3, n]
                # for consistency of the phase definition with the ODE solver
                if phic[np.argmax(np.abs(phic) > 1e-3)] < 0:
                    self.psi_all[:, n] = -self.psi_all[:, n]
            self.psis = np.zeros((self.eigen_all.shape[0], self.x_points.size))
            psis = self.psi_all[::3, :].T
            self.psis[:, :-1] = (psis[:, 1:] + psis[:, :-1]) / 2
            self.philh = self.psi_all[1::3, :].T
            self.phiso = self.psi_all[2::3, :].T
            self.eigen_es = self.eigen_all
            return self.eigen_es
        raise NotImplementedError(
            f"Matrix solver is not implemented for {self.crystal_type}"
        )

    def psi_overlap(self, upper: int, lower: int, shift=0) -> np.ndarray:
        """Return psi[upper] * psi[lower] with psi[lower] shifted by shift
        number of periods."""
        if self.crystal_type == "ZincBlende":
            return sum(
                phi[upper] * self._shift_psi(phi[lower], shift)
                for phi in (self.psis, self.philh, self.phiso)
            )
        # default fallback and crystalType == 'simple'
        return self.psis[upper] * self._shift_psi(self.psis[lower], shift)

    def populate_kane_matrix(self) -> sparse.spmatrix:
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
        assert self.crystal_type == "ZincBlende"
        kunit = hbar**2 / (2 * e0 * m0 * (1e-10 * self.x_step) ** 2)
        length = len(self.x_points)
        x_eg, x_f, x_ep, x_eso = self.band_params
        x_f_halft = np.empty(length)
        x_f_halft[1:] = (x_f[1:] + x_f[:-1]) / 2
        x_f_halft[0] = x_f[0]
        x_vc_half = np.empty(length)
        x_vc_half[1:] = (self.x_vc[1:] + self.x_vc[:-1]) / 2
        x_vc_half[0] = self.x_vc[0]
        self.h_banded = np.zeros((4, 3 * length))
        self.h_banded[3, ::3] = 2 * (1 + 2 * x_f_halft) * kunit + x_vc_half
        self.h_banded[3, 1::3] = self.x_vc - x_eg  # lh band
        self.h_banded[3, 2::3] = self.x_vc - x_eg - x_eso  # so band
        p = sqrt(x_ep * kunit)
        self.h_banded[2, 1::3] = sqrt(2 / 3) * p
        self.h_banded[2, 3::3] = sqrt(1 / 3) * p[:-1]
        self.h_banded[1, 2::3] = -sqrt(1 / 3) * p
        self.h_banded[1, 3::3] = -sqrt(2 / 3) * p[:-1]
        self.h_banded[0, 3::3] = -(1 + 2 * x_f[:-1]) * kunit
        if hasattr(self, "luttinger"):
            gamma1, gamma2, _ = self.luttinger
            tlh = gamma1 + 2 * gamma2 - 2 * x_ep / x_eg / 3
            tso = gamma1 - x_ep / x_eg / 3
            tlhhalf = (tlh[1:] + tlh[:-1]) / 2
            tsohalf = (tso[1:] + tso[:-1]) / 2
            tlh[1:-1] = (tlhhalf[1:] + tlhhalf[:-1]) / 2
            tso[1:-1] = (tsohalf[1:] + tsohalf[:-1]) / 2
            self.h_banded[3, 1::3] -= 2 * tlh * kunit
            self.h_banded[3, 2::3] -= 2 * tso * kunit
            self.h_banded[0, 4::3] = tlhhalf * kunit
            self.h_banded[0, 5::3] = tsohalf * kunit
        self.h_sparse = sparse.diags(
            [
                self.h_banded[0, 3:],
                self.h_banded[1, 2:],
                self.h_banded[2, 1:],
                self.h_banded[3, :],
                self.h_banded[2, 1:],
                self.h_banded[1, 2:],
                self.h_banded[0, 3:],
            ],
            [-3, -2, -1, 0, 1, 2, 3],
            shape=(3 * length, 3 * length),
        )
        return self.h_sparse

    def _is_basis_break(self, n: int) -> bool:
        if self.basis_injector_ar:
            if not self.layer_ar[n - 1] and self.layer_ar[n]:
                return True
        if self.basis_ar_injector:
            if self.layer_ar[n - 1] and not self.layer_ar[n]:
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
        start_idx = []
        end_idx = []
        # Get the region of interest
        if self.basis_ar_only:
            if self.layer_ar[0]:
                start_idx.append(0)
            for n in range(1, len(self.layer_ar)):
                if not self.layer_ar[n - 1] and self.layer_ar[n]:
                    start_idx.append(n)
                if self.layer_ar[n - 1] and not self.layer_ar[n]:
                    end_idx.append(n)
            if self.layer_ar[-1]:
                end_idx.append(len(self.layer_ar))
        else:
            start_idx.append(-1)
            for n in range(1, len(self.layer_ar)):
                if self._is_basis_break(n):
                    barrier = n if self.layer_vc(n) > self.layer_vc(n - 1) else n - 1
                    if barrier != start_idx[-1]:
                        start_idx.append(barrier)
                        end_idx.append(barrier + 1)
            if self._is_basis_break(0):
                barrier = 0 if self.layer_vc(0) > self.layer_vc(-1) else -1
                barrier += len(self.layer_widths)
                start_idx.append(barrier)
                end_idx.append(barrier + 1)
            if len(end_idx) == 0:
                start_idx = [0]
                end_idx = [len(self.layer_widths)]
            else:
                start_idx = start_idx[1:]
                end_idx = end_idx[1:] + [end_idx[0] + len(self.layer_widths)]

        self.eigen_es = np.empty((0))
        self.psis = np.empty((0, self.x_points.size))
        for n in range(0, len(start_idx)):
            d_cl = copy.copy(self)
            d_cl.reset_for_basis(start_idx[n], end_idx[n])
            d_cl.populate_x()
            d_cl.solve_whole()

            # map dCL result back to self
            shift = sum(self.layer_widths[: start_idx[n]]) - BASISPAD
            psis, eigen_es = self._basis_shifter(
                range(-1, self.repeats),
                d_cl.psis,
                d_cl.eigen_es - shift * self.e_field * EUNIT,
                d_cl.x_points + shift,
            )
            self.eigen_es = np.concatenate((self.eigen_es, eigen_es))
            self.psis = np.concatenate((self.psis, psis))
        self._reset_cache()
        self.status = "basis"
        return self.eigen_es

    def reset_for_basis(self, start: int, end: int):
        """Reset the parameters for only solving the layers of [start:end].
        This is a helper method for solve_basis"""
        self.repeats = 1
        self.state_per_repeat = (
            self.state_per_repeat * (end - start) // len(self.layer_widths)
        )
        self.layer_widths = (self.layer_widths * 2)[start:end]
        self.layer_widths[0] += BASISPAD
        self.layer_widths[-1] += BASISPAD
        self.layer_ar = (self.layer_ar * 2)[start:end]
        self._layer_vc = (self._layer_vc * 2)[start:end]
        self._layer_mc = (self._layer_mc * 2)[start:end]

    def _dipole(self, upper: int, lower: int, shift: int = 0) -> float:
        """Return Electrical dipole between upper and lower states
        in unit angstrom, update self.dipole.
        shift means lower state is shifted by this number of period.
        Should be called for any other related physics quantities."""
        # TODO: clean up self cache for upper -> lower states
        if self.solver == "ODE":
            psi_u = self.psis[upper, :]
            psi_l = self._shift_psi(self.psis[lower, :], shift)
            e_u = self.eigen_es[upper]
            e_l = self.eigen_es[lower] - shift * self.e_shift
            x_inv_mc_u = self._x_band_mass_inv(e_u)
            x_inv_mc_l = self._x_band_mass_inv(e_l)
            # Eq.(8) in PhysRevB.50.8663, with more precise mass
            # m_j comes from the density of states of transition final state
            z = np.trapezoid(
                psi_u * x_inv_mc_l * np.gradient(psi_l)
                - np.gradient(psi_u) * x_inv_mc_u * psi_l
            )
            z *= hbar**2 / (2 * (e_l - e_u) * e0 * m0) / (1e-10) ** 2  # Angstrom
        elif self.solver == "matrix" and self.crystal_type == "ZincBlende":
            # In principle this and above is equivalent, and numerically
            # approximately the same as tested in unit tests.
            # The above is remained for legacy reasons.
            z = self.x_step * np.trapezoid(
                self.psi_overlap(upper, lower, shift) * self.x_points
            )
        else:
            raise NotImplementedError("{self.solver} not implemented for dipole")
        return z

    def dipole(self, upper, lower) -> float:
        if self._cache_dipole[upper][lower] is None:
            self._cache_dipole[upper][lower] = self._dipole(upper, lower)
        return self._cache_dipole[upper][lower]

    def _lo_transition(self, upper: int, lower: int, shift: int = 0) -> float:
        # TODO: finite temperature version
        e_u = self.eigen_es[upper]
        e_l = self.eigen_es[lower] - shift * self.e_shift
        if e_u - e_l - self.avg_hw_lo < 0:
            # energy difference is smaller than a LO phonon
            # LO phonon scattering doesn't happen
            return INV_INF
        psi_l_sq = self.psi_overlap(lower, lower)
        ml = m0 * np.trapezoid(self.x_mc * psi_l_sq) * self.x_step
        kl = sqrt(2 * ml / hbar**2 * (e_u - e_l - self.avg_hw_lo) * e0)
        length = self.x_points.size
        if onedq is None:
            convpsi = fft.irfft(
                np.abs(fft.rfft(self.psi_overlap(upper, lower, shift), 2 * length)) ** 2
            )[:length]
            i_ij = (
                2
                * self.x_step**2
                * np.trapezoid(exp(-kl * self.x_points * 1e-10) * convpsi)
            )
        # C implementation
        else:
            i_ij = onedq.cLOphononScatter(
                self.x_step, kl, self.psi_overlap(upper, lower, shift)
            )
        return (
            ml
            * e0**2
            * self.avg_hw_lo
            * e0
            / hbar
            * i_ij
            / (4 * hbar**2 * self.eps_rho * eps0 * kl)
        ) / 1e12  # unit ps^-1

    def lo_transition(self, upper: int, lower: int) -> float:
        """The LO phonon transition lifetime from upper to lower,
        at zero temperature.
        This is using cached results.
        if status is 'solved-full' and the state is a recognized state of a
        period, it's calculated via translation of the wavefunction, otherwise
        it's calculated based on row self.psis
        """
        if self.status == "solved-full":
            try:
                pu, ushift = self.period_map[upper]
                pl, lshift = self.period_map[lower]
                if lshift - ushift not in (0, 1, -1):
                    return INV_INF
                return self._p_lo[lshift - ushift][pl][pu]
            except (TypeError, AttributeError):
                # TypeError is when self.periodMap return None
                # AttributeError is when self.periodMap does not exist
                pass
        if self._cache_lo[upper][lower] is None:
            self._cache_lo[upper][lower] = self._lo_transition(upper, lower)
        return self._cache_lo[upper][lower]

    def lo_lifetime(self, state: int) -> float:
        """Return the life time due to LO phonon scattering of the
        given state(label)
        This is using cached results.
        if status is 'solved-full' and the state is a recognized state of a
        period, it's calculated via translation of the wavefunction, otherwise
        it's calculated based on row self.psis
        """
        if self.status == "solved-full":
            try:
                return 1 / sum(
                    np.sum(self._p_lo[n][:, self.period_map[state][0]])
                    for n in range(3)
                )
            except (TypeError, AttributeError):
                # TypeError is when self.periodMap return None
                # AttributeError is when self.periodMap does not exist
                pass
        e_i = self.eigen_es[state]
        if onedq is None:
            return 1 / sum(
                self.lo_transition(state, q)
                for q in range(state)
                if self.eigen_es[q] <= e_i - self.avg_hw_lo
            )
        idxs = self.eigen_es <= e_i - self.avg_hw_lo
        e_js = self.eigen_es[idxs]
        (idxs,) = idxs.nonzero()
        psi_js_sq = np.array([self.psi_overlap(idx, idx) for idx in idxs])
        mjs = m0 * np.trapezoid(self.x_mc * psi_js_sq, axis=1) * self.x_step
        kls = sqrt(2 * mjs / hbar**2 * (e_i - e_js - self.avg_hw_lo) * e0)
        fjs = (
            mjs
            * e0**2
            * self.avg_hw_lo
            * e0
            / hbar
            / (4 * hbar**2 * self.eps_rho * eps0 * kls)
        )
        psi_ijs = np.array([self.psi_overlap(state, idx) for idx in idxs])
        i_ij_total = onedq.cLOtotal(self.x_step, kls, psi_ijs, fjs)
        return 1e12 / i_ij_total if i_ij_total > 0 else 1e20

    def _ifr_transition(
        self, upper: int, lower: int, shift: int = 0
    ) -> Tuple[float, float]:
        # TODO: finite temperature
        psi_usq = self.psi_overlap(upper, upper)
        psi_lsq = self.psi_overlap(lower, lower)
        psi_ul = self.psi_overlap(upper, lower, shift)
        e_u = self.eigen_es[upper]
        e_l = self.eigen_es[lower] - shift * self.e_shift
        if e_u < e_l:
            return INV_INF, -1
        mu = m0 * np.trapezoid(self.x_mc * psi_usq) * self.x_step
        ml = m0 * np.trapezoid(self.x_mc * psi_lsq) * self.x_step
        kl = sqrt(2 * ml / hbar**2 * (e_u - e_l) * e0)
        tau_inv = 0
        gamma = 0
        zn = 0
        layer_n = len(self.layer_widths)
        for _ in range(self.repeats):
            for n in range(layer_n):
                lamb = self.ifr_lambdas[n] * 1e-10  # to m
                delt = self.ifr_deltas[n]
                d_u = self.layer_vc((n + 1) % layer_n) - self.layer_vc(n)
                d_u *= e0  # to J
                # find interface
                zn += self.layer_widths[n]
                z_idx = np.argmax(self.x_points >= zn)
                if z_idx == 0 or z_idx == len(self.x_points) - 1:
                    continue
                z1 = self.x_points[z_idx - 1]
                z2 = self.x_points[z_idx]
                interp_a = (z2 - zn) / (z2 - z1)
                interp_b = 1 - interp_a
                psi_usqz = psi_usq[z_idx - 1] * interp_a + psi_usq[z_idx] * interp_b
                psi_lsqz = psi_lsq[z_idx - 1] * interp_a + psi_lsq[z_idx] * interp_b
                psi_ulz = psi_ul[z_idx - 1] * interp_a + psi_ul[z_idx] * interp_b
                scale = pi / hbar**3 * delt**2 * lamb**2 * d_u**2
                tau_inv += scale * ml * psi_ulz**2 * exp(-(lamb**2) * kl**2 / 4)
                gamma += (
                    scale * (psi_usqz - psi_lsqz) * (psi_usqz * mu - psi_lsqz * ml)
                ) / 2
        return tau_inv / 1e12, gamma / 1e12  # unit ps^-1

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
        assert self.include_ifr
        if self.status == "solved-full":
            try:
                pu, ushift = self.period_map[upper]
                pl, lshift = self.period_map[lower]
                if lshift - ushift not in (0, 1, -1):
                    return INV_INF
                return self._p_ifr[lshift - ushift][pl][pu]
            except (TypeError, AttributeError):
                # TypeError is when self.periodMap return None
                # AttributeError is when self.periodMap does not exist
                pass
        if self._cache_ifr[upper][lower] is not None:
            return self._cache_ifr[upper][lower]
        tau_inv, gamma = self._ifr_transition(upper, lower)
        self._cache_ifr[upper][lower] = tau_inv
        if gamma > 0:
            self._cache_gamma[upper][lower] = gamma
            self._cache_gamma[lower][upper] = gamma
        return self._cache_ifr[upper][lower]

    def ifr_lifetime(self, state: int) -> float:
        """Return to total life time due to IFR scattering.
        This is using cached results.
        if status is 'solved-full' and the state is a recognized state of a
        period, it's calculated via translation of the wavefunction, otherwise
        it's calculated based on row self.psis
        """
        assert self.include_ifr
        if self.status == "solved-full":
            try:
                return 1 / sum(
                    np.sum(self._p_ifr[n][:, self.period_map[state][0]])
                    for n in range(3)
                )
            except (TypeError, AttributeError):
                # TypeError is when self.periodMap return None
                # AttributeError is when self.periodMap does not exist
                pass
        return 1 / sum(self.ifr_transition(state, q) for q in range(state))

    def lifetime(self, state: int) -> float:
        """A convenience wrap of return total lifetime of LO and IFR scattering
        or only LO scattering depending on self.include_ifr."""
        if self.status == "solved-full":
            try:
                return 1 / self.decay_rates[self.period_map[state][0]]
            except (TypeError, AttributeError):
                # TypeError is when self.periodMap return None
                # AttributeError is when self.periodMap does not exist
                pass
        if self.include_ifr:
            return 1 / (1 / self.ifr_lifetime(state) + 1 / self.lo_lifetime(state))
        else:
            return self.lo_lifetime(state)

    def ifr_broadening(self, upper: int, lower: int) -> float:
        """Interface roughness induced broadening"""
        # TODO: shift life time
        if self.status == "solved-full":
            try:
                pu, ushift = self.period_map[upper]
                pl, lshift = self.period_map[lower]
                if lshift - ushift not in (0, 1, -1):
                    return 1 / INV_INF
                return self._p_gamma[lshift - ushift][pl][pu]
            except (TypeError, AttributeError):
                # TypeError is when self.periodMap return None
                # AttributeError is when self.periodMap does not exist
                pass
        if self._cache_gamma[upper][lower] is None:
            self.ifr_transition(upper, lower)
            self.ifr_transition(lower, upper)  # pylint: disable=arguments-out-of-order
        return self._cache_gamma[upper][lower]

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
        if not self.include_ifr:
            e_u = self.eigen_es[upper]
            e_l = self.eigen_es[lower]
            de = np.abs(e_u - e_l)
            return 0.05 * de
        gamma_u = 1 / self.ifr_lifetime(upper) + 1 / self.lo_lifetime(upper)
        gamma_l = 1 / self.ifr_lifetime(lower) + 1 / self.lo_lifetime(lower)
        gamma_parallel = self.ifr_broadening(upper, lower)
        # 1E12: ps^-1 -> Hz
        return (gamma_parallel + (gamma_u + gamma_l) / 2) * 1e12 * hbar / e0

    def period_recognize(self, tol: float = 5e-5) -> np.ndarray:
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
        period_idx = self.period_l / self.x_step
        psisq = np.abs(self.psis) ** 2
        self.starts = np.argmax(psisq > tol, axis=1)
        self.ends = np.argmax(psisq[:, ::-1] > tol, axis=1)
        self.period_idx = np.arange(len(self.eigen_es))[
            (self.starts > period_idx / 3) & (self.starts < 4 * period_idx / 3)
        ]
        # check if all states are far away enough from boundary
        self.un_bound = set()
        barrier_bound = np.max(self.x_vc + self.x_v_field) - self.x_v_field
        for n in self.period_idx:
            if self.ends[n] < period_idx / 3:
                self.un_bound.add(n)
                if self.eigen_es[n] < barrier_bound[-self.ends[n]]:
                    _logger.warning(
                        "State No.%d is close to the right boundary. "
                        "More repeats may be needed",
                        n,
                    )
        if len(self.period_idx) - len(self.un_bound) == 0:
            raise StateRecognizeError(
                "No well bounded state found. " "Try increase repeats."
            )
        return self.period_idx

    def period_map_build(
        self, tol: float = 3e-2, etol: float = 1e-3
    ) -> List[Tuple[int, int]]:
        """Map states to self.singlePeriodIdx, self.periodMap[n] is a tuple of
        (state index in self.singlePeriodIdx, shift of period(s) or
        None meaning it's not mapped."""
        self.period_map = [None] * len(self.eigen_es)
        single_e = self.eigen_es[self.period_idx]
        for n, state in enumerate(self.period_idx):
            self.period_map[state] = (n, 0)
        for state, en in enumerate(self.eigen_es):
            if self.period_map[state] is not None:
                continue
            psi = self.psis[state]
            for shift in range(1, self.repeats):
                en_shifted = en + shift * self.e_shift
                psi_shifted = self._shift_psi(psi, -shift)
                for n, s_e in enumerate(single_e):
                    if s_e < en_shifted - etol:
                        continue
                    if s_e > en_shifted + etol:
                        break
                    e_state = self.period_idx[n]
                    wf_diff = np.trapezoid((psi_shifted - self.psis[e_state]) ** 2)
                    wf_diff *= self.x_step
                    if wf_diff < tol:
                        # effective an L2 norm here, tested better than
                        # L1 or L-max norm
                        self.period_map[state] = (n, shift)
                        break
                if self.period_map[n] is not None:
                    break
        return self.period_map

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
        assert self.status.startswith("solved")
        idx_period = len(self.period_idx)
        # p for cache for the periodic version
        self._p_lo = [np.zeros((idx_period, idx_period)) for _ in range(3)]
        self._p_ifr = [np.zeros((idx_period, idx_period)) for _ in range(3)]
        self._p_gamma = [np.zeros((idx_period, idx_period)) for _ in range(3)]
        for i, lower in enumerate(self.period_idx):
            for j, upper in enumerate(self.period_idx):
                # keep this consistent with loMatrix etc
                for s in (1, -1) if lower == upper else (1, 0, -1):
                    self._p_lo[s][j][i] = self._lo_transition(lower, upper, s)
                    self._p_ifr[s][j][i], gamma = self._ifr_transition(lower, upper, s)
                    if gamma > 0:
                        self._p_gamma[s][j][i] = self._p_gamma[s][i][j] = gamma
        forward = self._p_lo[1] + self._p_ifr[1]
        backward = self._p_lo[-1] + self._p_ifr[-1]
        internal = self._p_lo[0] + self._p_ifr[0]
        self.transitions = internal + forward + backward
        self.decay_rates = np.sum(self.transitions, axis=0)
        self.population = null_space(np.diag(-self.decay_rates) + self.transitions)
        if self.population.shape[1] != 1:
            raise ValueError(
                "More than one steady state found. ", self.population.shape[1]
            )
        self.population = self.population[:, 0].T
        self.population /= np.sum(self.population)
        if self.carrier_leak > 5e-2:
            most_leak = max(
                (self.population[n], state)
                for n, state in enumerate(self.period_idx)
                if state in self.un_bound
            )
            _logger.warning(
                "The structure seems highly leak or more period needed. "
                "Most leaking (%.2f%%) unbounded state: %d",
                *most_leak,
            )
        self.flow = np.sum((forward - backward) @ self.population)
        self.status = "solved-full"
        return self.population

    @property
    def carrier_leak(self) -> float:
        return sum(
            self.population[n]
            for n, state in enumerate(self.period_idx)
            if state in self.un_bound
        )

    def state_population(self, state: int) -> float:
        """This method is only valid after fullPopulation has been called"""
        if self.period_map[state] is None:
            return None
        return self.population[self.period_map[state][0]]

    def _x_band_mass_inv(self, energy: float) -> np.ndarray:
        if self.crystal_type == "simple":
            return self.x_mc
        if self.crystal_type == "ZincBlende":
            x_eg, x_f, x_ep, x_eso = self.band_params
            e_mass = energy - self.x_vc
            return (
                1
                + 2 * x_f
                + 1 / 3 * x_ep / (e_mass + x_eg + x_eso)
                + 2 / 3 * x_ep / (e_mass + x_eg)
            )
        else:
            raise NotImplementedError(
                f"Material property for {self.crystal_type} crystal is not implemented"
            )


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
    x_res :
        Position resolution, in Armstrong
    e_res :
        Energy resolution, in meV.
        This number being too large may results in miss of some states while
        this being too small will make a long computation time.
        The parameter does not mean the accuracy of the eigen energy. It's
        required for algorithm reasons because of lack of a universal global
        root finding.
    state_per_repeat :
        Number of states per repeat, used for calculating matrix_eigen_count
    wl :
        The wavelength for the design, in unit um, gain spectrum and optimization,
        but doesn't go into quantum solver

    layer_widths :
        Width of each layer, in angstrom. len = No. of layers
    layer_mtrls :
        Label of materials, depending on substrate. len = No. of layers
    layer_dopings :
        Doping per volume in unit 1e17 cm-3. len = No. of layers
    layer_ar :
        Binaries indicating if the layer is active(True) or not(False),
        only affects basis solver. len = No. of layers

    customIFR :
        Wether to use a customized IFR parameter rather than a material determined
        parameter.
    mtrlIFRLambda :
        The interface roughness lambda after materials[n], len = No. of materials
    mtrlIFRDelta :
        The interface roughness delta after materials[n], len = No. of materials

    e_field :
        External (static) electrical field, in kV/cm = 1e5 V/m
    repeats :
        Number of repeat times for the given structure
    temp :
        Temperature of the device, affecting material property

    basis_ar_only :
        For basis solver if only the Active Region (AR) should be solved.
    basis_ar_injector :
        For basis solver if there should be separator between AR->Injector
    basis_injector_ar :
        For basis solver if there should be separator between Injector->AR

    solver :
        The solver used for the eigen problem: 'ODE' or 'matrix'.
        By default 'ODE' if C library exists, 'matrix' is a full back.
    include_ifr :
        Weather to include IFR scattering for performance estimation.
    matrix_eigen_count :
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
    mole_fracs: List[float]
    wl: float
    layer_mtrls: List[int]
    layer_dopings: List[float]
    custom_ifr: bool
    mtrl_ifr_lambda: List[float]
    mtrl_ifr_delta: List[float]
    temp: float
    description: str

    def __init__(
        self,
        substrate="InP",
        materials=None,
        mole_fracs=None,
        x_res=0.5,
        e_res=0.5,
        state_per_repeat=20,
        layer_widths=None,
        layer_matrls=None,
        layer_dopings=None,
        custom_ifr=False,
        mtrl_ifr_lambda=None,
        mtrl_ifr_delta=None,
        ifr_delta=None,
        ifr_lambda=None,
        layer_ar=None,
        e_field=0,
        repeats=3,
        temp=300.0,
        solver="ODE",
        description="",
        wl=3.0,
    ):
        if materials is None:
            materials = ["InGaAs", "AlInAs"]
        if mole_fracs is None:
            mole_fracs = [0.53, 0.52]
        if layer_widths is None:
            layer_widths = [10.0]
        assert isinstance(layer_widths, list)
        assert isinstance(materials, list)
        assert isinstance(mole_fracs, list)
        layer_len = len(layer_widths)
        mtrl_len = len(materials)
        assert mtrl_len >= 1
        assert len(mole_fracs) == mtrl_len
        self.substrate = substrate
        self.materials = materials
        self.mole_fracs = mole_fracs
        self.layer_mtrls = [0] * layer_len if layer_matrls is None else layer_matrls
        self.layer_dopings = (
            [0.0] * layer_len if layer_dopings is None else layer_dopings
        )
        self.temperature = temp
        self.custom_ifr = custom_ifr
        if not custom_ifr:
            if isinstance(mtrl_ifr_delta, list):
                assert len(mtrl_ifr_delta) == mtrl_len
                assert isinstance(mtrl_ifr_lambda, list)
                assert len(mtrl_ifr_lambda) == mtrl_len
                self.mtrl_ifr_delta = mtrl_ifr_delta
                self.mtrl_ifr_lambda = mtrl_ifr_lambda
            else:
                self.mtrl_ifr_delta = [mtrl_ifr_delta or 0.0] * mtrl_len
                self.mtrl_ifr_lambda = [mtrl_ifr_lambda or 0.0] * mtrl_len
            ifr_delta, ifr_lambda = self._get_ifr_list()
        self.description = description
        super().__init__(
            x_res=x_res,
            e_res=e_res,
            state_per_repeat=state_per_repeat,
            layer_widths=layer_widths,
            layer_ar=layer_ar,
            ifr_deltas=ifr_delta,
            ifr_lambdas=ifr_lambda,
            e_field=e_field,
            repeats=repeats,
        )
        self.crystal_type = material.MTRL_PARAM[substrate]["Crystal"]
        self.sub_mtrl = material.Material(self.substrate, self.temperature)
        self.wl = wl
        self.solver = solver
        self.update_mtrls()

    def __copy__(self):
        return QCLayers(
            substrate=self.substrate,
            materials=copy.copy(self.materials),
            mole_fracs=copy.copy(self.mole_fracs),
            x_res=self.x_step,
            e_res=self.e_step,
            state_per_repeat=self.state_per_repeat,
            layer_widths=copy.copy(self.layer_widths),
            layer_matrls=copy.copy(self.layer_mtrls),
            layer_dopings=copy.copy(self.layer_dopings),
            custom_ifr=self.custom_ifr,
            mtrl_ifr_lambda=copy.copy(self.mtrl_ifr_lambda),
            mtrl_ifr_delta=copy.copy(self.mtrl_ifr_delta),
            ifr_delta=copy.copy(self.ifr_deltas),
            ifr_lambda=copy.copy(self.ifr_lambdas),
            layer_ar=copy.copy(self.layer_ar),
            e_field=self.e_field,
            repeats=self.repeats,
            temp=self.temperature,
            solver=self.solver,
            description=self.description,
            wl=self.wl,
        )

    def _get_ifr_list(self) -> Tuple[List[float], List[float]]:
        """Get IFR parameters for SchrodingerLayer. Should be called
        every time the material list changes."""
        assert not self.custom_ifr
        if self.mtrl_ifr_delta is not None:
            ifr_delta = [self.mtrl_ifr_delta[m] for m in self.layer_mtrls]
        else:
            ifr_delta = [0.0] * len(self.layer_mtrls)
        if self.mtrl_ifr_lambda is not None:
            ifr_lambda = [self.mtrl_ifr_lambda[m] for m in self.layer_mtrls]
        else:
            ifr_lambda = [0.0] * len(self.layer_mtrls)
        return ifr_delta, ifr_lambda

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
        self.a_parallel = self.sub_mtrl.param["alc"]
        self.mtrl_alloys = [
            material.Alloy(self.materials[idx], self.mole_fracs[idx], self.temperature)
            for idx in range(len(self.materials))
        ]
        for al in self.mtrl_alloys:
            al.set_strain(self.a_parallel)

    def set_mtrl(self, n: int, mtrl: str = None, mole_frac: float = None):
        """Set material[n] to new material (mtrl) and/or moleFrac"""
        if mtrl is None and mole_frac is None:
            raise ValueError("Nothing changed")
        if mtrl is None:
            mtrl = self.materials[n]
        if mole_frac is None:
            mole_frac = self.mole_fracs[n]
        self.mole_fracs[n] = mole_frac
        self.materials[n] = mtrl
        self.mtrl_alloys[n] = material.Alloy(mtrl, mole_frac, self.temperature)
        self.mtrl_alloys[n].set_strain(self.a_parallel)

    def add_mtrl(
        self,
        mtrl: str = None,
        mole_frac: float = None,
        ifr_lambda: float = None,
        ifr_delta: float = None,
    ):
        """Add a new material possibility"""
        self.materials.append(mtrl if mtrl else QC_MATERIAL[self.substrate][0])
        self.mole_fracs.append(mole_frac if mole_frac else 0.0)
        self.mtrl_ifr_lambda.append(
            ifr_lambda if ifr_lambda else self.mtrl_ifr_lambda[-1]
        )
        self.mtrl_ifr_delta.append(ifr_delta if ifr_delta else self.mtrl_ifr_delta[-1])
        self.mtrl_alloys.append(
            material.Alloy(self.materials[-1], self.mole_fracs[-1], self.temperature)
        )
        self.mtrl_alloys[-1].set_strain(self.a_parallel)

    def del_mtrl(self, n: int):
        """Delete materials labeled n.
        All layers of this material will become previous
        `materials[n-1 if n >0 else 1]`.  There should be at least two
        materials otherwise there will be error."""
        if len(self.materials) <= 2:
            raise ValueError("There should be at least 2 materials")
        self.materials.pop(n)
        self.mole_fracs.pop(n)
        for i in range(len(self.layer_mtrls)):
            if self.layer_mtrls[i] >= n:
                self.layer_mtrls[i] = (
                    self.layer_mtrls[i] - 1 if self.layer_mtrls[i] > 0 else 0
                )

    def add_layer(self, n: int, width: int, mtrlIdx: int, AR: bool, doping: float):
        self.layer_mtrls.insert(n, mtrlIdx)
        self.layer_dopings.insert(n, doping)
        super().add_layer(n, width, is_ar=AR)
        if not self.custom_ifr:
            self.ifr_deltas, self.ifr_lambdas = self._get_ifr_list()

    def rotate_layer(self):
        super().rotate_layer()
        for layer_list in (self.layer_mtrls, self.layer_dopings):
            layer_list.insert(0, layer_list.pop())

    def del_layer(self, n: int):
        super().del_layer(n)
        for ayer_list in (self.layer_mtrls, self.layer_dopings):
            ayer_list.pop(n)

    def invert_layer(self):
        super().invert_layer()
        self.layer_mtrls = self.layer_mtrls[::-1]
        self.layer_dopings = self.layer_dopings[::-1]
        if not self.custom_ifr:
            self.ifr_deltas, self.ifr_lambdas = self._get_ifr_list()

    def set_substrate(self, subs: str):
        if subs in QC_MATERIAL:
            self.substrate = subs
            self.crystal_type = material.MTRL_PARAM[subs]["Crystal"]
            mtrl_n = len(self.materials)
            self.materials = (QC_MATERIAL[subs] * mtrl_n)[0:mtrl_n]
            self.update_mtrls()
        else:
            raise ValueError(f"Substrate {subs} not supported")

    def set_temperature(self, temp: float):
        self.temperature = temp
        self.sub_mtrl.set_temperature(temp)
        self.update_mtrls()

    @property
    def mtrl_offset(self) -> float:
        """Return the conduction band offset (difference between highest
        conduction band and lowest conduction band energy) of materials,
        in unit eV"""
        ecgs = [alloy.param["EcG"] for alloy in self.mtrl_alloys]
        return max(ecgs) - min(ecgs)

    @property
    def net_strain(self) -> float:
        """Return average strain perpendicular to the layer plane, in
        percentage."""
        if sum(self.layer_widths) <= 1e-5:
            return -1
        total_strain = sum(
            self.mtrl_alloys[self.layer_mtrls[n]].eps_perp * self.layer_widths[n]
            for n in range(len(self.layer_widths))
        )
        return 100 * total_strain / sum(self.layer_widths)

    def layer_vc(self, n: int) -> float:
        return self.mtrl_alloys[self.layer_mtrls[n]].param["EcG"]

    def layer_mc(self, n: int) -> float:
        return self.mtrl_alloys[self.layer_mtrls[n]].param["me0"]

    def populate_material(self):
        """
        Update following band structure parameters (with type *np.array
        of float*): **xVc, xVX, xVL, xVLH, xVSO,
        band_params = (xEg, xF, xEp, xESO)**
        """
        if self.crystal_type == "ZincBlende":
            length = self.x_points.size
            self.x_dopings = np.empty(length)

            self.x_vx = np.empty(length)  # X point in the band
            self.x_vl = np.empty(length)  # L point in the band
            self.x_vlh = np.empty(length)  # The light hole valence band
            self.x_vso = np.empty(length)  # The heavy hole valence band

            # band parameters
            x_f = np.empty(length)
            x_eg = np.empty(length)
            x_eso = np.empty(length)
            x_ep = np.empty(length)

            for n in range(len(self.layer_widths)):
                indices = self.x_layer_nums == n
                self.x_dopings[indices] = self.layer_dopings[n]
                for p, key in (
                    (self.x_vlh, "EvLH"),
                    (self.x_vso, "EvSO"),
                    (self.x_vx, "EcX"),
                    (self.x_vl, "EcL"),
                    (x_eg, "EgLH"),
                    (x_eso, "ESO"),
                    (x_ep, "Ep"),
                    (x_f, "F"),
                ):
                    p[indices] = self.mtrl_alloys[self.layer_mtrls[n]].param[key]
            self.band_params = (x_eg, x_f, x_ep, x_eso)

            ext_field = self.x_points * self.e_field * EUNIT
            for p in (self.x_vx, self.x_vl, self.x_vlh, self.x_vso):
                p -= ext_field
        else:
            raise NotImplementedError(
                f"Material property for {self.crystal_type} crystal is not implemented"
            )
        # LO phonon
        if sum(self.layer_widths) <= 1e-5:
            self.avg_hw_lo = -1
            self.eps_rho = 69.0
        else:
            sumhwlo = sum(
                self.mtrl_alloys[self.layer_mtrls[n]].param["hwLO"]
                * self.layer_widths[n]
                for n in range(len(self.layer_widths))
            )
            self.avg_hw_lo = sumhwlo / sum(self.layer_widths)
            eps_inf = np.array([a.param["epsInf"] for a in self.mtrl_alloys])
            epss = np.array([a.param["epss"] for a in self.mtrl_alloys])
            eps_rho = 1 / (1 / eps_inf - 1 / epss)
            self.eps_rho = np.sum(eps_rho[self.layer_mtrls] * self.layer_widths) / sum(
                self.layer_widths
            )
        # IFR
        if not self.custom_ifr:
            self.ifr_deltas, self.ifr_lambdas = self._get_ifr_list()

    def reset_for_basis(self, start: int, end: int):
        super().reset_for_basis(start, end)
        self.layer_mtrls = (self.layer_mtrls * 2)[start:end]
        self.layer_dopings = (self.layer_dopings * 2)[start:end]

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
        self.tau_lo_ul = 1 / self.lo_transition(upper, lower)
        self.tau_lo_l = self.lo_lifetime(lower)
        self.tau_lo_u = self.lo_lifetime(upper)
        if self.include_ifr:
            self.tau_ifr_ul = 1 / self.ifr_transition(upper, lower)
            self.tau_ifr_u = self.ifr_lifetime(upper)
            self.tau_ifr_l = self.ifr_lifetime(lower)
            self.tau_u = 1 / (1 / self.tau_lo_u + 1 / self.tau_ifr_u)
            self.tau_l = 1 / (1 / self.tau_lo_l + 1 / self.tau_ifr_l)
            self.tau_ul = 1 / (1 / self.tau_lo_ul + 1 / self.tau_ifr_ul)
        else:
            self.tau_u = self.tau_lo_u
            self.tau_l = self.tau_lo_l
            self.tau_ul = self.tau_lo_ul
        return (
            self.dipole(upper, lower) ** 2 * self.tau_u * (1 - self.tau_l / self.tau_ul)
        )

    def effective_ridx(self, wl: ScalerOrArray) -> ScalerOrArray:
        """Return the effective refractive index for TM mode"""
        if sum(self.layer_widths) == 0:
            return 1.0
        self.mtrl_r_idx = [
            (
                m.mole_frac * material.rIdx[m.elm_a.name](wl)
                + (1 - m.mole_frac) * material.rIdx[m.elm_b.name](wl)
            )
            for m in self.mtrl_alloys
        ]
        self.layer_r_idx = np.array([self.mtrl_r_idx[n] for n in self.layer_mtrls])
        neff = np.average(
            1 / self.layer_r_idx**2, axis=0, weights=self.layer_widths
        ) ** (-1 / 2)
        return neff

    @property
    def sheet_density(self) -> float:
        """Return the sheet density of doping per period, in unit cm^-2"""
        # 1E9 -> 1E17 cm^-3 * Angstrom -> cm^-2
        return (
            sum(
                self.layer_dopings[n] * self.layer_widths[n]
                for n in range(0, len(self.layer_widths))
            )
            * 1e9
        )

    def full_population(self) -> np.array:
        """Apart from SchrodingerLayer.full_population, this also yield
        current density, with knowledge of doping"""
        res = super().full_population()
        # ps^-1 -> kA/cm^2, where 1E9 = 1E12 * 1E-3, 1E12 is (ps^-1 -> s^-1)
        self.current = self.flow * self.sheet_density * e0 * 1e9
        return res

    def full_gain_spectrum(self, wl: ScalerOrArray = None) -> ScalerOrArray:
        """Perform fully automatic calculation for the gain on wavelength(s)."""
        if wl is None:
            wl = self.wl
        neff = self.effective_ridx(wl)
        de0 = h * c0 / (wl * 1e-6) / e0
        gain = 0
        # self.gainul = {}
        for i in range(len(self.period_idx)):
            upper = self.period_idx[i]
            for j in range(i + 1, len(self.period_idx)):
                lower = self.period_idx[j]
                for shift in (-1, 0, 1):
                    dipole = self._dipole(upper, lower, shift) * 1e-8  # to cm
                    e_u = self.eigen_es[upper]
                    e_l = self.eigen_es[lower] - shift * self.e_shift
                    de = e_u - e_l
                    dpop = self.population[i] - self.population[j]
                    if de < 0:
                        de, dpop = -de, -dpop
                    if self.include_ifr:
                        gamma_para = self._p_gamma[shift][j][i]
                    else:
                        gamma_para = self.dephasing(upper, lower)
                        gamma_para /= 1e12 * hbar / e0  # unit to Hz
                    gamma = (
                        (gamma_para + (self.decay_rates[i] + self.decay_rates[j]) / 2)
                        * 1e12
                        * hbar
                        / e0
                    )
                    gainul = dpop * dipole**2 * gamma / (gamma**2 + (de - de0) ** 2)
                    gain = gain + gainul
        # e0^2 / (hbar * c0 * eps0) is dimension 1.. 1E-8 makes unit cm^-1
        gain *= (
            e0**2
            * de0
            * self.sheet_density
            / (hbar * neff * c0 * eps0 * self.period_l * 1e-8)
        )
        return gain


def optimize_layer(qcl: QCLayers, n: int, upper: int, lower: int, max_iter: int = 50):
    """Optimize FoM*Lorentzian for n-th layer thickness, assuming the state
    index does not change. optimization is performed by searching on the
    position resolution steps.

    Warning: This cannot specify a correct state if there are states
    index crossing.

    TODO: Use period recognizer to improve the algorithm
    """
    e_u = qcl.eigen_es[upper]
    e_l = qcl.eigen_es[lower]
    if e_u < e_l:
        upper, lower = lower, upper
        e_u, e_l = e_l, e_u
    # 1E-6 um -> m... in unit eV
    w0 = h * c0 / (qcl.wl * 1e-6) / e0

    def reduce_fom():
        wul = e_u - e_l
        gamma = qcl.dephasing(upper, lower)
        return (
            qcl.figure_of_merit(upper, lower)
            * gamma
            * w0
            / (gamma**2 + (wul - w0) ** 2)
        )

    width = round(qcl.layer_widths[n] / qcl.x_step) * qcl.x_step
    fom_now = reduce_fom()
    _logger.info(
        "Start Optimizing Layer NO %d for FoM between state %d and %d.\n"
        "\tStart at width=%.1f, FoM=%.5g",
        n,
        upper,
        lower,
        width,
        fom_now,
    )

    def new_fom(new_width):
        qcl.layer_widths[n] = new_width
        qcl.populate_x()
        qcl.solve_whole()
        qcl.dipole(upper, lower)
        return reduce_fom()

    fom_minus = new_fom(width - qcl.x_step)
    fom_plus = new_fom(width + qcl.x_step)
    for _ in range(max_iter):
        if fom_now < fom_plus:
            fom_minus = fom_now
            fom_now = fom_plus
            width += qcl.x_step
            fom_plus = new_fom(width + qcl.x_step)
        elif fom_now < fom_minus:
            fom_plus = fom_now
            fom_now = fom_minus
            width -= qcl.x_step
            fom_minus = new_fom(width - qcl.x_step)
        else:
            _logger.info("Maximum iteration reached.")
            break
        _logger.info("\twidth=%.1f, FoM=%.5g", width, fom_now)
    qcl.layer_widths[n] = width
    _logger.info("finished, width=%.1f, FoM=%.5g", width, fom_now)
    return fom_now


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
