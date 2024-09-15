"""Module for Maxwell solver using transfer matrix method."""

import typing
from collections import defaultdict
from typing import List

import numpy as np
from numpy import cos, exp, pi, sin, sinc, sqrt

from ErwinJr2.material import MTRL_PARAM, REFRATICE_INDICES

# Effective mass in unit of free electron mass.
EFF_MASS_MAP = {
    mtrl: MTRL_PARAM[mtrl]["me0"] for mtrl in ("GaAs", "InAs", "AlAs", "InP")
}

ALLOY_MAP = {
    "AlxGa1-xAs": ("AlAs", "GaAs"),
    "InxGa1-xAs": ("InAs", "GaAs"),
    "Al1-xInxAs": ("InAs", "AlAs"),
}
ALLOY_NICKNAMES = {
    "AlGaAs": "AlxGa1-xAs",
    "InGaAs": "InxGa1-xAs",
    "AlInAs": "Al1-xInxAs",
}
DOPABLE_MATERIALS = set(["GaAs", "InAs", "AlAs", "InP"] + list(ALLOY_MAP.keys()))


class MaxwellLayer(object):
    """Class for layer structure for Maxwell solver using transfer matrix
    method.

    This is used as the base class of :class:`.OptStrata` for separation of
    material property and the solver.

    Parameters
    ----------
    wl :
        Wavelength in vacuum to guide in the stratum

    indices :
        The Refractive indices of the layers except for the top and the
        substrate. The substrate and the top layer is not part of the list
        because they decides the boundary condition for the solver.

    index_0 :
        The refractive index for the top layer

    index_s :
        The refractive index for the substrate layer

    layer_widths :
        Thickness of the stratum, same unit as wl. The first and last elements
        are for top and substrate and is not used for calculation
    """

    wl: float
    index_0: complex
    index_s: complex
    indices: np.ndarray
    layer_widths: List[complex]

    def __init__(self, wl: float, layer_widths: List = None, indices: List = None):
        if layer_widths is None:
            layer_widths = [1.0, 3.0]
        if indices is None:
            indices = [1.0, 1.0]
        self.wl = wl
        self.index_0 = indices[0]
        self.index_s = indices[-1]
        self.indices = np.array(indices[1:-1], dtype=np.complex128)
        self.layer_widths = (
            np.array(layer_widths)
            if len(layer_widths) == len(indices)
            else np.array([1.0] + list(layer_widths) + [2.0])
        )

    def transfer_tm(self, beta: complex) -> np.ndarray:
        """The transfer matrix for a TM wave

        Calculate the transfer matrix for TM wave with frequency
        :math:`\\omega = c/\\text{wl}` on the stratum structure described
        with the thickness and index list.
        wl and Ls should be same unit

        Parameters
        ----------
        beta: complex
            The normalized wavenumber along z.
            :math:`k_z = 2\\pi/\\lambda *\\beta`

        Returns
        -------
        np.ndarray
            A 2*2 matrix, the transfer matrix from substrate to top
        """
        alpha = self._alpha(beta)
        k = 2 * pi / self.wl
        phi = alpha * k * self.layer_widths[1:-1]
        ms = np.moveaxis(
            np.array(
                [
                    [cos(phi), -1j * sin(phi) * alpha / self.indices**2],
                    [
                        -1j
                        * sinc(phi / pi)
                        * self.indices**2
                        * k
                        * self.layer_widths[1:-1],
                        cos(phi),
                    ],
                ]
            ),
            -1,
            0,
        )
        # np.sinc (x) is defined as sin(pi*x)/(pi*x)
        return np.linalg.multi_dot(ms)

    def chi_m_tm(self, beta: complex) -> complex:
        """Modal-dispersion function for TM wave

        Calculate the with modal-dispersion function for TM wave with frequency
        :math:`\\omega = c/\\text{wl}` on the stratum structure described
        with the thickness and index list;
        top/substrate defined by index and index.
        wl and Ls should be same unit

        Parameters
        ----------
        beta : complex
            The normalized wavenumber along z.
            :math:`k_z = 2\\pi/\\lambda*\\beta`

        Returns
        -------
        complex
            The modal-dispersion function
        """
        gamma0 = np.sqrt((1 + 0j) * self.index_0**2 - beta**2) / self.index_0**2
        if gamma0.imag < 0:
            gamma0 = -gamma0
        gammas = np.sqrt((1 + 0j) * self.index_s**2 - beta**2) / self.index_s**2
        if gammas.imag < 0:
            gammas = -gammas
        m = self.transfer_tm(beta)
        return m[0, 0] * gammas + m[0, 1] + (m[1, 0] * gammas + m[1, 1]) * gamma0

    def bound_mode_tm(self, beta: typing.Optional[complex] = None) -> complex:
        """Solve for TM bounded mode near beta

        Solve for TM bounded mode near beta (as first guess in root finding)
        with frequency :math:`\\omega = c/\\text{wl}` on the stratum structure
        described with the thickness and index list;
        top/substrate defined by index0 and indexS.
        wl and Ls should be same unit

        Parameters
        ----------
        beta : complex
            The effective refractive index traveling along z (initial guess)

        Returns
        -------
        complex
            the effective refractive index for bounded mode

        Raises
        ------
        Exception:
            The solver doesn't converge.
            Most likely a bounded mode doesn't exist.

        """
        if beta is None:
            if len(self.indices) == 0:
                beta = 1.5 * self.index_s
            else:
                beta = max(self.indices.real)
                beta = min(
                    beta,
                    2.0
                    * np.average(self.indices.real, weights=self.layer_widths[1:-1]),
                )
        # # Minimization algorithm. for complex function zero searching
        # beta0 = newton(lambda beta:
        #                chiMTM(beta, wl, Ls, indices, index0, indexS).imag,
        #                x0=max(indices))
        # res = minimize(lambda x: abs(chiMTM(
        #                x[0] + 1j*x[1], wl, Ls, indices, index0, indexS))**2,
        #                x0 = [beta0+0.001, -1E-3], tol=1e-12, method="BFGS")
        # beta = res.x[0] + 1j*res.x[1]
        # residual = abs(chiMTM(beta, wl, Ls, indices, index0, indexS))
        # print(res)
        residual = self.chi_m_tm(beta)
        beta_diff = beta
        t = 0
        dbeta = 1e-6
        while abs(beta_diff) > 1e-5 and abs(residual) > 1e-10:
            # print(beta, residual)
            t += 1
            fp = (self.chi_m_tm(beta + dbeta) - self.chi_m_tm(beta - dbeta)) / (
                2 * dbeta
            )
            beta_diff = residual / fp
            beta = beta - beta_diff
            residual = self.chi_m_tm(beta)
            if t > 200:
                raise TimeoutError("Doesn't converge")
        return beta

    def _alpha(self, beta: complex) -> np.ndarray:
        """return alpha = np.sqrt((1+0j)*self.indices**2 - beta**2)
        for isotropic material"""
        return np.sqrt((1 + 0j) * self.indices**2 - beta**2)

    def populate_mode(
        self, beta: complex, xs: np.ndarray
    ) -> typing.Union[np.ndarray, np.ndarray, np.ndarray]:
        """Generate TM modes (field) on position array xs

        Generate TM modes (field) on position array xs, assuming beta is a
        bounded Mode (i.e. return value of :meth:`~OneDMaxwell.boundModeTM`)
        return Ey, Hx, Ez, normalized to max(abs(Ey)) = 1; unit of Hx is
        :math:`([\\text{nit of beta}]\\times \\sqrt{\\mu_0/\\epsilon_0})^-1`

        Parameters
        ----------
        beta : complex
            The effective refractive index traveling along z

        xs : np.ndarray
            The position coordinate to calculate field on

        Returns
        -------
        np.ndarray:
            E field perpendicular to layers and propagation

        np.ndarray:
            H field parallel to layers and perpendicular to propagation

        np.ndarray:
            E field in propagation
        """
        h_x = np.zeros(xs.shape, dtype=np.complex128)
        e_y = np.zeros(xs.shape, dtype=np.complex128)
        e_z = np.zeros(xs.shape, dtype=np.complex128)
        k = 2 * pi / self.wl

        # top outside medium
        alpha0 = np.sqrt((1 + 0j) * self.index_0**2 - beta**2)
        if alpha0.imag < 0:
            alpha0 = -alpha0
        gamma0 = alpha0 / self.index_0**2
        e_z0 = -gamma0
        h_x0 = 1
        h_x[xs < 0] = exp(-1j * alpha0 * k * xs[xs < 0])
        e_z[xs < 0] = -gamma0 * h_x[xs < 0]
        e_y[xs < 0] = -beta * h_x[xs < 0] / self.index_0**2

        # middle stratum
        lsum = np.cumsum(self.layer_widths[:-1]) - self.layer_widths[0]
        alpha = self._alpha(beta)
        # [[cos(phi), 1j*sin(phi)*alpha/indices**2],
        #  [1j*sinc(phi/pi)*indices**2*k * Ls, cos(phi)]]
        for n in range(len(lsum) - 1):
            idx = (xs >= lsum[n]) & (xs < lsum[n + 1])
            phi = alpha[n] * k * (xs[idx] - lsum[n])
            e_z[idx] = cos(phi) * e_z0 + 1j * (
                alpha[n] / self.indices[n] ** 2 * sin(phi) * h_x0
            )
            h_x[idx] = cos(phi) * h_x0 + 1j * k * (xs[idx] - lsum[n]) * (
                sinc(phi / pi) * self.indices[n] ** 2 * e_z0
            )
            e_y[idx] = -beta * h_x[idx] / self.indices[n] ** 2
            phi_l = alpha[n] * k * self.layer_widths[n + 1]
            e_z0, h_x0 = (
                cos(phi_l) * e_z0
                + 1j * alpha[n] / self.indices[n] ** 2 * (sin(phi_l) * h_x0),
                cos(phi_l) * h_x0
                + 1j
                * k
                * self.layer_widths[n + 1]
                * sinc(phi_l / pi)
                * (self.indices[n] ** 2 * e_z0),
            )

        # last substrate
        alphas = np.sqrt((1 + 0j) * self.index_s**2 - beta**2)
        if alphas.imag < 0:
            alphas = -alphas
        gammas = alphas / self.index_s**2
        h_x[xs >= lsum[-1]] = h_x0 * exp(
            1j * alphas * k * (xs[xs >= lsum[-1]] - lsum[-1])
        )
        e_z[xs >= lsum[-1]] = gammas * h_x[xs >= lsum[-1]]
        e_y[xs >= lsum[-1]] = -beta * h_x[xs >= lsum[-1]] / self.index_s**2

        scale = e_y[np.argmax(np.abs(e_y))]
        self.h_x = h_x / scale
        self.e_y = e_y / scale
        self.e_z = e_z / scale
        return self.e_y, self.h_x, self.e_z

    def populate_indices(self, xs: np.ndarray) -> np.ndarray:
        """Generate indices position array xs

        Parameters
        ----------
        xs : np.ndarray
            The position coordinate to calculate field on

        Returns
        -------
        np.ndarray:
            refractive indices of the material at position xs
        """
        lsum = np.cumsum(self.layer_widths[:-1]) - self.layer_widths[0]
        if len(self.indices) > 0:
            n = np.piecewise(
                xs + 0j,
                [(xs >= lsum[i]) & (xs < lsum[i + 1]) for i in range(len(lsum) - 1)],
                self.indices,
            )
        else:
            n = np.empty(xs.shape, dtype=np.complex128)
        n[xs < 0] = self.index_0
        n[xs >= lsum[-1]] = self.index_s
        return n

    def confinementy(
        self,
        beta: complex,
        actives: List[np.ndarray],
        xs: typing.Optional[np.ndarray] = None,
        e_y: typing.Optional[np.ndarray] = None,
    ) -> float:
        """Return the confinement factor corresponds to mode with effective
        refractive index beta. Assuming active only couple to E_y filed.
        The active region is defined in `actives`. If xs and Ey is None, they
        will be generated.

        Parameters
        ----------
        beta : complex
            The refractive index of the mode to calculate

        actives : list(np.ndarray(bool))
            The list of active region for confinement factor

        xs : np.ndarray
            The array for positions: controls the accuracy of the numerical
            integral for confinement factor calculation

        Ey : np.ndarray(complex)
            The field to integral on

        Returns
        -------
        float:
            The confinement factor of the structure
        """
        # TODO: an analytical version of it.
        if xs is None:
            xs = np.linspace(-3, sum(self.layer_widths[1:]), 5000)
        if e_y is None:
            e_y, _, _ = self.populate_mode(beta, xs)
        confinement = 0
        nx = self.populate_indices(xs).real
        for ar in actives:
            confinement += np.trapezoid(nx[ar] * np.abs(self.e_y[ar]) ** 2, xs[ar])
        confinement = (
            beta.real * confinement / np.trapezoid((nx * np.abs(e_y)) ** 2, xs)
        )
        return confinement


class MaxwellLayerAnisotropic(MaxwellLayer):
    """class for anisotropic Maxwell layers, cannot deal with anisotropy for
    top air and bottom substrate."""

    def __init__(self, wl, Ls=None, indexz=None, indexy=None):
        super().__init__(wl, Ls, indexz)
        if indexy is None:
            self.indexy = np.copy(self.indices)
        else:
            self.indexy = np.array(indexy, dtype=np.complex128)

    def _alpha(self, beta):
        return self.indices / self.indexy * np.sqrt((1 + 0j) * self.indexy**2 - beta**2)

    def populate_indices(self, xs):
        lsum = np.cumsum(self.layer_widths[:-1]) - self.layer_widths[0]
        if len(self.indices) > 0:
            n = np.piecewise(
                xs + 0j,
                [(xs >= lsum[i]) & (xs < lsum[i + 1]) for i in range(len(lsum) - 1)],
                self.indices,
            )
            ny = np.piecewise(
                xs + 0j,
                [(xs >= lsum[i]) & (xs < lsum[i + 1]) for i in range(len(lsum) - 1)],
                self.indexy,
            )
        else:
            n = np.empty(xs.shape, dtype=np.complex128)
            ny = np.empty(xs.shape, dtype=np.complex128)
        n[xs < 0] = self.index_0
        n[xs >= lsum[-1]] = self.index_s
        ny[xs < 0] = self.index_0
        ny[xs >= lsum[-1]] = self.index_s
        return n, ny


class OptStrata(MaxwellLayer):
    """Class for groups of stratum, a wrapper of :class:`.MaxwellLayer`
    with material parsing support

    Parameters
    ----------
    wl : float
        Wavelength in vacuum to guide in the stratum

    materials : list(str)
        Name of the material for each strata, with materials[0] being the top
        (usually air) and materials[-1] the substrate.

    mole_fracs : list(float)
        Mole fractions for each material. The number should be between 0 and 1.
        For strata that the parameter is not applicable, the number doesn't
        influence the result.

    dopings : list(float)
        The doping level in unit 1E17 cm^-3 for the material.
        For strata that the parameter is not applicable, the number doesn't
        influence the result.

    layer_widths : list(float)
        Thickness of  stratum, same unit as wl. The first and last elements
        are for top and substrate and is not used for calculation

    mobilities : None or list(float)
        Mobility influence the relaxation rate for plasmon resonance.
        When it's None, it is assumed to have 1E13 s^1 relaxation.

    cstm_idx : dict
        A dictionary of customized material, with key the name and value
        the complex refractive index

    cstm_prd : dict
        A dictionary of customized material, with key the name and value a
        list of two elements, 0-th being the period length (float in unit
        Angstrom) and 1-st (int) being the number of periods.
        If the 0-th value is 0 or not key not exists, it's not a periodic
        structure. The variable is only for book keeping and is
        not used to validate Ls within the class

    cstm_gain : dict
        A dictionary of customized material, with key the name and value
        the gain coefficient. The variable is only for book keeping and is
        not used for any calculation.
    """

    def __init__(
        self,
        wl=3.0,
        materials=None,
        mole_fracs=None,
        dopings=None,
        layer_widths=None,
        mobilities=None,
        cstm_idx=None,
        cstm_prd=None,
        cstm_gain=None,
    ):
        if materials is None:
            materials = ["Air", "InP"]
        if layer_widths is None:
            layer_widths = [1.0, 1.0]
        if cstm_idx is None:
            cstm_idx = dict()
        if cstm_prd is None:
            cstm_prd = defaultdict(float)
        if cstm_gain is None:
            cstm_gain = defaultdict(float)
        super().__init__(wl)
        mtrl_count = len(materials)
        self.materials = list(
            ALLOY_NICKNAMES[m] if m in ALLOY_NICKNAMES else m for m in materials
        )
        self.mole_fracs = list(mole_fracs) if mole_fracs else [0.0] * mtrl_count
        self.dopings = list(dopings) if dopings else [0.0] * mtrl_count
        self.layer_widths = (
            np.array(layer_widths)
            if len(layer_widths) == len(materials)
            else np.array([1.0] + list(layer_widths) + [2.0])
        )
        self.mobilities = [None] * len(materials) if mobilities is None else mobilities
        assert len(self.mole_fracs) == mtrl_count
        assert len(self.dopings) == mtrl_count
        assert len(self.layer_widths) == mtrl_count
        self.cstm_idx = cstm_idx
        self.cstm_prd = cstm_prd
        self.cstm_gain = cstm_gain
        self.update_indices()

    def __str__(self):
        return "\n".join(
            (
                f"wavelength: {self.wl}",
                f"materials: {self.materials}",
                f"moleFracs: {self.mole_fracs}",
                f"dopings: {self.dopings}",
                f"layers: {self.layer_widths}",
                f"mobilities: {self.mobilities}",
                f"custom: {self.cstm_idx}",
            )
        )

    def set_wl(self, wl):
        self.wl = wl

    def insert(
        self, row, material=None, mole_frac=None, doping=None, width=None, mobility=None
    ):
        """Insert a strata indexed with row (top air and bottom substrate
        included) with parameters listed"""
        if row <= 0 or row >= len(self.materials):
            raise IndexError("Cannot insert strata beyond top and substrate.")
        self.materials.insert(row, material if material else self.materials[row - 1])
        self.mole_fracs.insert(
            row, mole_frac if mole_frac is not None else self.mole_fracs[row - 1]
        )
        self.dopings.insert(row, doping if doping is not None else 0)
        self.layer_widths = np.insert(
            self.layer_widths, row, width if width is not None else 1.0
        )
        self.mobilities.insert(row, mobility if mobility else self.mobilities[row - 1])

    def delete(self, row):
        """Delete the strata indexed with row (top air and bottom substrate
        included)"""
        del self.materials[row]
        del self.mole_fracs[row]
        del self.dopings[row]
        del self.mobilities[row]
        self.layer_widths = np.delete(self.layer_widths, row)

    def index_of(self, n):
        """Return the refractive of the strata indexed with n (top air and
        bottom substrate included). The result is a combination of linear
        interpolation of mole fraction and plasmonic loss by Drude model"""
        mtrl = self.materials[n]
        if mtrl in self.cstm_idx:
            return self.cstm_idx[mtrl]
        if mtrl in ALLOY_MAP:
            # linear interpolation
            layer_r_idx = self.mole_fracs[n] * REFRATICE_INDICES[ALLOY_MAP[mtrl][0]](
                self.wl
            ) + (1 - self.mole_fracs[n]) * REFRATICE_INDICES[ALLOY_MAP[mtrl][1]](
                self.wl
            )
        else:
            layer_r_idx = REFRATICE_INDICES[mtrl](self.wl)

        plasmon_unit = 8.9698e-5 * self.wl**2
        # 8.9698e-05 = mu_0*1E23*e**2*(1E-6)**2/(electron_mass*4*pi**2)
        # omega_P = plasmon_unit * doping / me / (1 + 1j gamma/omega)
        # doping in unit 1E17 cm^-3,
        # electron effective mass me in unit free electron mass
        gamma_unit = 5.3088e-3
        # 5.3088E-3 = 1E7/(2*pi*c)
        # gamma_unit * gamma_c (unit 1E13 s-1) * wl (unit um) = gamma_c/omega
        if mtrl in DOPABLE_MATERIALS:
            if mtrl in ALLOY_MAP:
                me = EFF_MASS_MAP[ALLOY_MAP[mtrl][0]] * self.mole_fracs[
                    n
                ] + EFF_MASS_MAP[ALLOY_MAP[mtrl][1]] * (1 - self.mole_fracs[n])
            else:
                me = EFF_MASS_MAP[mtrl]
            nue = (
                1.0
                if self.mobilities[n] is None
                else 1.7588e11 / (me * self.mobilities[n] * 1e-4)
            )
            # 1.7588E11 = electron charge / free electron mass
            layer_r_idx = sqrt(
                layer_r_idx**2
                - (
                    plasmon_unit
                    * self.dopings[n]
                    / me
                    / (1 + 1j * gamma_unit * self.wl * nue)
                )
            )
        return layer_r_idx

    def update_indices(self):
        """Update indices information according to material info

        Yield
        -----
        indices : list(complex)
            complex refractive index of the stratum

        index0 : complex
            Refractive index of the top (before Ls[0]) strata

        indexS : complex
            refractive index of the substrate (after Ls[-1])
        """
        self.indices = np.array([self.index_of(n) for n in range(len(self.materials))])
        self.index_0 = self.indices[0]
        self.index_s = self.indices[-1]
        self.indices = self.indices[1:-1]

    def populate_mtrl(self, xs, mtrl_list=None):
        """Populate a boolean array for index slicing on xs for material
        in the mtrlList. If mtrlList is None, the active region is labelled by
        anything start with "Active".

        Parameters
        ----------
        mtrlList : list(str) or tuple(str)
            name of materials to be labeled

        xs : np.ndarray
            The position coordinate to label material

        Returns
        -------
        list(np.ndarray(bool))
            list of slicing indices to label materials
        """
        res = []
        lsum = np.empty(len(self.layer_widths) + 1)
        lsum[1:] = np.cumsum(self.layer_widths) - self.layer_widths[0]
        lsum[0] = -np.inf
        lsum[-1] = np.inf
        for i in range(len(self.materials)):
            if (
                self.materials[i] in mtrl_list
                if mtrl_list is not None
                else self.materials[i].startswith("Active")
            ):
                res.append((xs >= lsum[i]) & (xs <= lsum[i + 1]))
        return res

    def confinementy(self, beta, xs=None, Ey=None):
        """Return the confinement factor corresponds to mode with effective
        refractive index beta. If xs and Ey is None, they will be generated.
        The active region is labelled by anything start with "Active".

        Parameters
        ----------
        beta : complex
            The refractive index of the mode to calculate

        xs : np.ndarray
            The array for positions: controls the accuracy of the numerical
            integral for confinement factor calculation

        Ey : np.ndarray(complex)
            The field to integral on
        """
        if xs is None:
            xs = np.linspace(-3, sum(self.layer_widths[1:]), 5000)
        if Ey is None:
            Ey, _, _ = self.populate_mode(beta, xs)
        return super().confinementy(beta, self.populate_mtrl(xs), xs, Ey)


def optimize_opt_strata(
    stratum: OptStrata,
    alpha_m,
    to_optimize: List[int],
    max_length: float,
    max_iter: int = 20,
    tol=0.05,
) -> float:
    """Optimize strata with threshold gain as the minimize objective function
    g_th = (alphaM + alphaW)/confinement, where the waveguide loss alphaW
    and the confinement factor is a character of strata.
    The optimization is done on layer indexed by elements in toOptimize,
    and conditioned on total length to be smaller than maxLength.
    Newton's method with an increasing penalty is used.

    Returns
    -------
    float: The optimized threshold gain in cm^-1
    """
    assert to_optimize
    assert all(n > 0 and n < len(stratum.layer_widths) - 1 for n in to_optimize)

    def objective(penalty: float):
        beta = stratum.bound_mode_tm()
        gamma = stratum.confinementy(beta)
        alpha_w = 4 * pi / (stratum.wl / 1e4) * beta.imag  # cm^-1
        length = sum(stratum.layer_widths[1:-1])
        return (alpha_m + alpha_w) / gamma + penalty * max(0, length - max_length) ** 2

    penalty = 0.5
    now = objective(penalty)
    for n in range(max_iter):
        changed = False
        for n in to_optimize:
            w = stratum.layer_widths[n]
            step = 1e-5 * w
            stratum.layer_widths[n] = w - step
            new_minus = objective(penalty)
            stratum.layer_widths[n] = w + step
            new_plus = objective(penalty)
            dif = (new_plus - new_minus) / 2 / step
            ddif = (new_plus + new_minus - 2 * now) / step**2
            dw = -dif / ddif
            if ddif < 0 or abs(dw) > 0.2 * w:
                dw = 0.2 * w if dif < 0 else -0.2 * w
            elif abs(dw) < tol:
                continue
            changed = True
            res = tol / 5
            stratum.layer_widths[n] = res * round((w + dw) / res)
            now = objective(penalty)
        if not changed and sum(stratum.layer_widths[1:-1]) <= max_length + tol:
            break
        penalty *= 2
    print("Finished with penalty: ", penalty)
