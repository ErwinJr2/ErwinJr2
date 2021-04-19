import numpy as np
from numpy import sqrt, exp, pi, sin, cos, sinc
from .Material import rIdx, MParam
from collections import defaultdict
import typing
from typing import List

# Effective mass in unit of free electron mass.
EffectiveMass = {mtrl: MParam[mtrl]['me0']
                 for mtrl in ('GaAs', 'InAs', 'AlAs', 'InP')}

Alloy = {"AlxGa1-xAs": ("AlAs", "GaAs"),
         "InxGa1-xAs": ("InAs", "GaAs"),
         "Al1-xInxAs": ("InAs", "AlAs")}
AlloyNick = {'AlGaAs': 'AlxGa1-xAs',
             'InGaAs': 'InxGa1-xAs',
             'AlInAs': 'Al1-xInxAs'}
Dopable = set(['GaAs', 'InAs', 'AlAs', 'InP'] + list(Alloy.keys()))


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

    index0 :
        The refractive index for the top layer

    indexS :
        The refractive index for the substrate layer

    Ls :
        Thickness of the stratum, same unit as wl. The first and last elements
        are for top and substrate and is not used for calculation
    """
    wl: float
    index0: complex
    indexS: complex
    indices: np.ndarray
    Ls: List[complex]

    def __init__(self, wl: float, Ls: List = [1.0, 3.0],
                 allIndex: List = [1.0, 1.0]):
        self.wl = wl
        self.index0 = allIndex[0]
        self.indexS = allIndex[-1]
        self.indices = np.array(allIndex[1:-1], dtype=np.complex128)
        self.Ls = (np.array(Ls) if len(Ls) == len(allIndex)
                   else np.array([1.0] + list(Ls) + [2.0]))

    def transferTM(self, beta: complex) -> np.ndarray:
        """Transfer matrix for TM wave

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
        phi = alpha * k * self.Ls[1:-1]
        ms = np.moveaxis(np.array([
            [cos(phi), -1j*sin(phi)*alpha/self.indices**2],
            [-1j*sinc(phi/pi)*self.indices**2*k * self.Ls[1:-1], cos(phi)]]),
            -1, 0)
        # np.sinc (x) is defined as sin(pi*x)/(pi*x)
        return np.linalg.multi_dot(ms)

    def chiMTM(self, beta: complex) -> complex:
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
        gamma0 = np.sqrt((1+0j)*self.index0**2 - beta**2)/self.index0**2
        if gamma0.imag < 0:
            gamma0 = -gamma0
        gammas = np.sqrt((1+0j)*self.indexS**2 - beta**2)/self.indexS**2
        if gammas.imag < 0:
            gammas = -gammas
        m = self.transferTM(beta)
        return m[0, 0]*gammas+m[0, 1]+(m[1, 0]*gammas+m[1, 1])*gamma0

    def boundModeTM(self, beta: typing.Optional[complex] = None) -> complex:
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
                beta = 1.5*self.indexS
            else:
                beta = max(self.indices.real)
                beta = min(beta, 2.0*np.average(self.indices.real,
                           weights=self.Ls[1:-1]))
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
        residual = self.chiMTM(beta)
        betaDiff = beta
        t = 0
        dbeta = 1E-6
        while abs(betaDiff) > 1E-5 and abs(residual) > 1E-10:
            # print(beta, residual)
            t += 1
            fp = (self.chiMTM(beta+dbeta) - self.chiMTM(beta-dbeta))/(2*dbeta)
            betaDiff = residual / fp
            beta = beta - betaDiff
            residual = self.chiMTM(beta)
            if t > 200:
                raise TimeoutError("Doesn't converge")
                break
        return beta

    def _alpha(self, beta: complex) -> np.ndarray:
        """return alpha = np.sqrt((1+0j)*self.indices**2 - beta**2)
        for isotropic material"""
        return np.sqrt((1+0j)*self.indices**2 - beta**2)

    def populateMode(self, beta: complex, xs: np.ndarray
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
        Hx = np.zeros(xs.shape, dtype=np.complex128)
        Ey = np.zeros(xs.shape, dtype=np.complex128)
        Ez = np.zeros(xs.shape, dtype=np.complex128)
        k = 2 * pi / self.wl

        # top outside medium
        alpha0 = np.sqrt((1+0j)*self.index0**2 - beta**2)
        if alpha0.imag < 0:
            alpha0 = -alpha0
        gamma0 = alpha0/self.index0**2
        Ez0 = -gamma0
        Hx0 = 1
        Hx[xs < 0] = exp(-1j*alpha0*k*xs[xs < 0])
        Ez[xs < 0] = -gamma0*Hx[xs < 0]
        Ey[xs < 0] = -beta*Hx[xs < 0]/self.index0**2

        # middle stratum
        lsum = np.cumsum(self.Ls[:-1]) - self.Ls[0]
        alpha = self._alpha(beta)
        # [[cos(phi), 1j*sin(phi)*alpha/indices**2],
        #  [1j*sinc(phi/pi)*indices**2*k * Ls, cos(phi)]]
        for n in range(len(lsum)-1):
            idx = (xs >= lsum[n]) & (xs < lsum[n+1])
            phi = alpha[n]*k*(xs[idx] - lsum[n])
            Ez[idx] = cos(phi)*Ez0 + 1j*(
                alpha[n]/self.indices[n]**2*sin(phi)*Hx0)
            Hx[idx] = cos(phi)*Hx0 + 1j*k*(xs[idx] - lsum[n])*(
                                sinc(phi/pi)*self.indices[n]**2*Ez0)
            Ey[idx] = -beta*Hx[idx]/self.indices[n]**2
            phiL = alpha[n]*k*self.Ls[n+1]
            Ez0, Hx0 = (cos(phiL)*Ez0 + 1j*alpha[n]/self.indices[n]**2*(
                                           sin(phiL)*Hx0),
                        cos(phiL)*Hx0 + 1j*k*self.Ls[n+1]*sinc(phiL/pi)*(
                                           self.indices[n]**2*Ez0))

        # last substrate
        alphas = np.sqrt((1+0j)*self.indexS**2 - beta**2)
        if alphas.imag < 0:
            alphas = -alphas
        gammas = alphas/self.indexS**2
        Hx[xs >= lsum[-1]] = Hx0 * exp(1j*alphas * k *
                                       (xs[xs >= lsum[-1]] - lsum[-1]))
        Ez[xs >= lsum[-1]] = gammas*Hx[xs >= lsum[-1]]
        Ey[xs >= lsum[-1]] = -beta*Hx[xs >= lsum[-1]]/self.indexS**2

        scale = Ey[np.argmax(np.abs(Ey))]
        self.Hx = Hx/scale
        self.Ey = Ey/scale
        self.Ez = Ez/scale
        return self.Ey, self.Hx, self.Ez

    def populateIndices(self, xs: np.ndarray) -> np.ndarray:
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
        lsum = np.cumsum(self.Ls[:-1]) - self.Ls[0]
        if len(self.indices) > 0:
            n = np.piecewise(xs+0j, [(xs >= lsum[i]) & (xs < lsum[i+1])
                                     for i in range(len(lsum)-1)],
                             self.indices)
        else:
            n = np.empty(xs.shape, dtype=np.complex128)
        n[xs < 0] = self.index0
        n[xs >= lsum[-1]] = self.indexS
        return n

    def confinementy(self, beta: complex, actives: List[np.ndarray],
                     xs: typing.Optional[np.ndarray] = None,
                     Ey: typing.Optional[np.ndarray] = None) -> float:
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
            xs = np.linspace(-3, sum(self.Ls[1:]), 5000)
        if Ey is None:
            Ey, _, _ = self.populateMode(beta, xs)
        confinement = 0
        nx = self.populateIndices(xs).real
        for ar in actives:
            confinement += np.trapz(
                nx[ar] * np.abs(self.Ey[ar])**2, xs[ar])
        confinement = beta.real * confinement / np.trapz(
            (nx * np.abs(Ey))**2, xs)
        return confinement


class MaxwellLayer_anisotropic(MaxwellLayer):
    """class for anisotropic Maxwell layers, cannot deal with anisotropy for
    top air and bottom substrate."""
    def __init__(self, wl, Ls=[1.0, 1.0], indexz=[1.0, 1.0], indexy=None):
        super(MaxwellLayer_anisotropic, self).__init__(wl, Ls, indexz)
        if indexy is None:
            self.indexy = np.copy(self.indices)
        else:
            self.indexy = np.array(indexy, dtype=np.complex128)

    def _alpha(self, beta):
        return self.indices/self.indexy*np.sqrt(
                    (1+0j)*self.indexy**2 - beta**2)

    def populateIndices(self, xs):
        lsum = np.cumsum(self.Ls[:-1]) - self.Ls[0]
        if len(self.indices) > 0:
            n = np.piecewise(xs+0j, [(xs >= lsum[i]) & (xs < lsum[i+1])
                                     for i in range(len(lsum)-1)],
                             self.indices)
            ny = np.piecewise(xs+0j, [(xs >= lsum[i]) & (xs < lsum[i+1])
                                      for i in range(len(lsum)-1)],
                              self.indexy)
        else:
            n = np.empty(xs.shape, dtype=np.complex128)
            ny = np.empty(xs.shape, dtype=np.complex128)
        n[xs < 0] = self.index0
        n[xs >= lsum[-1]] = self.indexS
        ny[xs < 0] = self.index0
        ny[xs >= lsum[-1]] = self.indexS
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

    moleFracs : list(float)
        Mole fractions for each material. The number should be between 0 and 1.
        For strata that the parameter is not applicable, the number doesn't
        influence the result.

    dopings : list(float)
        The doping level in unit 1E17 cm^-3 for the material.
        For strata that the parameter is not applicable, the number doesn't
        influence the result.

    Ls : list(float)
        Thickness of  stratum, same unit as wl. The first and last elements
        are for top and substrate and is not used for calculation

    mobilities : None or list(float)
        Mobility influence the relaxation rate for plasmon resonance.
        When it's None, it is assumed to have 1E13 s^1 relaxation.

    cstmIdx : dict
        A dictionary of customized material, with key the name and value
        the complex refractive index

    cstmPrd : dict
        A dictionary of customized material, with key the name and value a
        list of two elements, 0-th being the period length (float in unit
        Angstrom) and 1-st (int) being the number of periods.
        If the 0-th value is 0 or not key not exists, it's not a periodic
        structure. The variable is only for book keeping and is
        not used to validate Ls within the class

    cstmGain : dict
        A dictionary of customized material, with key the name and value
        the gain coefficient. The variable is only for book keeping and is
        not used for any calculation.
    """
    def __init__(self, wl=3.0, materials=['Air', 'InP'], moleFracs=None,
                 dopings=None, Ls=[1.0, 1.0], mobilities=None,
                 cstmIndx=dict(), cstmPrd=defaultdict(float),
                 cstmGain=defaultdict(float)):
        super(OptStrata, self).__init__(wl)
        N = len(materials)
        self.materials = list(AlloyNick[m] if m in AlloyNick else m
                              for m in materials)
        self.moleFracs = list(moleFracs) if moleFracs else [0.0] * N
        self.dopings = list(dopings) if dopings else [0.0] * N
        self.Ls = (np.array(Ls) if len(Ls) == len(materials)
                   else np.array([1.0] + list(Ls) + [2.0]))
        self.mobilities = ([None]*len(materials)
                           if mobilities is None else mobilities)
        assert len(self.moleFracs) == N
        assert len(self.dopings) == N
        assert len(self.Ls) == N
        self.cstmIndx = cstmIndx
        self.cstmPrd = cstmPrd
        self.cstmGain = cstmGain
        self.updateIndices()

    def __str__(self):
        return "\n".join(("wavelength: %f" % self.wl,
                          "materials: %s" % str(self.materials),
                          "moleFracs: %s" % str(self.moleFracs),
                          "dopings: %s" % str(self.dopings),
                          "Ls: %s" % str(self.Ls),
                          "mobilities: %s" % str(self.mobilities),
                          "custom: %s" % str(self.cstmIndx)))

    def setWl(self, wl):
        self.wl = wl

    def insert(self, row, material=None, moleFrac=None, doping=None,
               L=None, mobility=None):
        """Insert a strata indexed with row (top air and bottom substrate
        included) with parameters listed"""
        if row <= 0 or row >= len(self.materials):
            raise IndexError("Cannot insert strata beyond top and substrate.")
        self.materials.insert(row, material if material else
                              self.materials[row-1])
        self.moleFracs.insert(row, moleFrac if moleFrac is not None else
                              self.moleFracs[row-1])
        self.dopings.insert(row, doping if doping is not None else 0)
        self.Ls = np.insert(self.Ls, row, 1.0)
        self.mobilities.insert(row, mobility if mobility else
                               self.mobilities[row-1])

    def delete(self, row):
        """Delete the strata indexed with row (top air and bottom substrate
        included)"""
        del self.materials[row]
        del self.moleFracs[row]
        del self.dopings[row]
        del self.mobilities[row]
        self.Ls = np.delete(self.Ls, row)

    def indexOf(self, n):
        """Return the refractive of the strata indexed with n (top air and
        bottom substrate included). The result is a combination of linear
        interpolation of mole fraction and plasmonic loss by Drude model"""
        mtrl = self.materials[n]
        if mtrl in self.cstmIndx:
            return self.cstmIndx[mtrl]
        if mtrl in Alloy:
            # linear interpolation
            layerRIdx = (
                self.moleFracs[n] * rIdx[Alloy[mtrl][0]](self.wl)
                + (1 - self.moleFracs[n]) * rIdx[Alloy[mtrl][1]](self.wl))
        else:
            layerRIdx = rIdx[mtrl](self.wl)

        plasmonUnit = 8.9698E-5 * self.wl**2
        # 8.9698e-05 = mu_0*1E23*e**2*(1E-6)**2/(electron_mass*4*pi**2)
        # omega_P = plasmonUnit * doping / me / (1 + 1j gamma/omega)
        # doping in unit 1E17 cm^-3,
        # electron effective mass me in unit free electron mass
        gammaUnit = 5.3088E-3
        # 5.3088E-3 = 1E7/(2*pi*c)
        # gammaUnit * gamma_c (unit 1E13 s-1) * wl (unit um) = gamma_c/omega
        if mtrl in Dopable:
            if mtrl in Alloy:
                me = (EffectiveMass[Alloy[mtrl][0]] * self.moleFracs[n] +
                      EffectiveMass[Alloy[mtrl][1]] * (1 - self.moleFracs[n]))
            else:
                me = EffectiveMass[mtrl]
            nue = (1.0 if self.mobilities[n] is None else
                   1.7588E11/(me * self.mobilities[n] * 1E-4))
            # 1.7588E11 = electron charge / free electron mass
            layerRIdx = sqrt(layerRIdx**2 - (
                plasmonUnit * self.dopings[n]/me
                / (1 + 1j*gammaUnit * self.wl * nue)))
        return layerRIdx

    def updateIndices(self):
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
        self.indices = np.array(
            [self.indexOf(n) for n in range(len(self.materials))])
        self.index0 = self.indices[0]
        self.indexS = self.indices[-1]
        self.indices = self.indices[1:-1]

    def populateMtrl(self, xs, mtrlList=None):
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
        lsum = np.empty(len(self.Ls)+1)
        lsum[1:] = np.cumsum(self.Ls) - self.Ls[0]
        lsum[0] = -np.inf
        lsum[-1] = np.inf
        for i in range(len(self.materials)):
            if (self.materials[i] in mtrlList if mtrlList is not None else
                    self.materials[i].startswith('Active')):
                res.append((xs >= lsum[i]) & (xs <= lsum[i+1]))
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
            xs = np.linspace(-3, sum(self.Ls[1:]), 5000)
        if Ey is None:
            Ey, _, _ = self.populateMode(beta, xs)
        return super(OptStrata, self).confinementy(
            beta, self.populateMtrl(xs), xs, Ey)


def optimizeOptStrata(stratum: OptStrata, alphaM,
                      toOptimize: List[int], maxLength: float,
                      iter: int = 20, tol=0.05) -> float:
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
    assert(toOptimize)
    assert(all(n > 0 and n < len(stratum.Ls)-1 for n in toOptimize))

    def objective(penalty: float):
        beta = stratum.boundModeTM()
        gamma = stratum.confinementy(beta)
        alphaW = 4*pi/(stratum.wl/1E4) * beta.imag  # cm^-1
        length = sum(stratum.Ls[1:-1])
        return (alphaM+alphaW)/gamma+penalty*max(0, length-maxLength)**2

    penalty = 0.5
    now = objective(penalty)
    for n in range(iter):
        changed = False
        for n in toOptimize:
            w = stratum.Ls[n]
            step = 1E-5 * w
            stratum.Ls[n] = w - step
            newMinus = objective(penalty)
            stratum.Ls[n] = w + step
            newPlus = objective(penalty)
            dif = (newPlus - newMinus)/2/step
            ddif = (newPlus + newMinus - 2*now)/step**2
            dw = -dif/ddif
            if ddif < 0 or abs(dw) > 0.2*w:
                dw = 0.2*w if dif < 0 else -0.2*w
            elif abs(dw) < tol:
                continue
            changed = True
            res = tol/5
            stratum.Ls[n] = res*round((w + dw)/res)
            now = objective(penalty)
        if not changed and sum(stratum.Ls[1:-1]) <= maxLength+tol:
            break
        penalty *= 2
    print('Finished with penalty: ', penalty)
