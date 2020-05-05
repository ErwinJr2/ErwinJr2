"""
This file defined the Stratum class for confined optical mode in 1D using
transfer matrix method
"""

import numpy as np
from numpy import sqrt, exp, pi, sin, cos, sinc
from rFittings import AlGaAsIndex, SiNxIndex, SiO2Index
from collections import defaultdict
# from OneDQuantum.OneDMaxwell import *

# refractive indices
rIdx = {
    # J. Appl. Phys., 94, 6447-6455 (2003)
    # https://refractiveindex.info/?shelf=main&book=GaAs&page=Skauli
    # 0.97um - 17um
    "GaAs": lambda wl: sqrt(5.372514 + 5.466742/(1-(0.4431307/wl)**2)
                            + 0.02429960/(1-(0.8746453/wl)**2)
                            + 1.957522/(1-(36.9166/wl)**2)),
    # J. Appl. Phys. 36, 1841-1844 (1965)
    # https://refractiveindex.info/?shelf=main&book=InAs&page=Lorimor
    # 3.7um - 31.3um
    "InAs": lambda wl: sqrt(11.1 + 0.71/(1-(2.551/wl)**2)
                            + 2.75/(1-(45.66/wl)**2)),
    # Handbook of Optics, 2nd edition, Vol. 2. McGraw-Hill 1994
    # https://refractiveindex.info/?shelf=main&book=InP&page=Pettit
    # 0.95um - 10um
    "InP": lambda wl: sqrt(7.255 + 2.316/(1-(0.6263/wl)**2)
                           + 2.765/(1-(32.935/wl)**2)),
    # J. Appl. Phys. 42, 3499-3500 (1971)
    # https://refractiveindex.info/?shelf=main&book=AlAs&page=Fern
    # 0.56um - 2.2um
    "AlAs": lambda wl: sqrt(2.0792 + 6.0840/(1-(0.2822/wl)**2)
                            + 1.9/(1-(27.62/wl)**2)),
    "Au": lambda wl: (-0.1933-0.382j + (0.3321+6.8522j)*wl
                      + (0.0938-0.1289j)*wl**2),
    # from old ej
    # "Au": lambda wl: (-0.1933-0.382j + (0.3321+6.8522j)*wl
    #                   + (0.0938**2+0.1289**2*1j)*wl),
    "AlxGa1-xAs": lambda wl, x: AlGaAsIndex(wl, x),
    "SiNx": lambda wl: SiNxIndex(wl),
    "SiO2": lambda wl: SiO2Index(wl),
    "Air": lambda wl: 1
}


# Effective mass in unit of free electron mass.
# Not using Material.py data to reduce dependence
EffectiveMass = {'GaAs': 0.067, 'InAs': 0.026, 'AlAs': 0.15, 'InP': 0.0795}

Alloy = {"AlxGa1-xAs": ("AlAs", "GaAs"),
         "InxGa1-xAs": ("InAs", "GaAs"),
         "Al1-xInxAs": ("InAs", "AlAs")}
AlloyNick = {'AlGaAs': 'AlxGa1-xAs',
             'InGaAs': 'InxGa1-xAs',
             'AlInAs': 'Al1-xInxAs'}
Dopable = set(['GaAs', 'InAs', 'AlAs', 'InP'] + list(Alloy.keys()))


class MaxwellLayer(object):
    """Class for layer structure for Maxwell solver using transfer matrix
    method

    Parameters
    ----------
    wl : float
        Wavelength in vacuum to guide in the stratum
    indices : list(complex)
        refractive indicies of the layers
    Ls : list(float)
        Thickness of  stratum, same unit as wl. The first and last elements
        are for top and substrate and is not used for calculation
    """
    def __init__(self, wl, indices=[1.0, 1.0], Ls=[1.0, 1.0]):
        self.wl = wl
        self.index0 = indices[0]
        self.indexs = indices[-1]
        self.indices = np.array(indices[1:-1], dtype=np.complex128)
        self.Ls = (np.array(Ls) if len(Ls) == len(indices)
                   else np.array([1.0] + list(Ls) + [2.0]))

    def transferTM(self, beta):
        """Tranfer matrix for TM wave

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
        alpha = np.sqrt((1+0j)*self.indices**2 - beta**2)
        k = 2 * pi / self.wl
        phi = alpha * k * self.Ls[1:-1]
        ms = np.moveaxis(np.array([
            [cos(phi), -1j*sin(phi)*alpha/self.indices**2],
            [-1j*sinc(phi/pi)*self.indices**2*k * self.Ls[1:-1], cos(phi)]]),
            -1, 0)
        # np.sinc (x) is defined as sin(pi*x)/(pi*x)
        return np.linalg.multi_dot(ms)

    def chiMTM(self, beta):
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
        gammas = np.sqrt((1+0j)*self.indexs**2 - beta**2)/self.indexs**2
        if gammas.imag < 0:
            gammas = -gammas
        m = self.transferTM(beta)
        return m[0, 0]*gammas+m[0, 1]+(m[1, 0]*gammas+m[1, 1])*gamma0

    def boundModeTM(self, beta=None):
        """Solve for TM bounded mode near beta

        Solve for TM bounded mode near beta (as first guess in root finding)
        with frequency :math:`\\omega = c/\\text{wl}` on the stratum structure
        discribed with the thickness and index list;
        top/substrate defined by index0 and indexs.
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
            beta = max(self.indices.real, default=1.5*self.indexs)
        # # Minimization algo. for complex function zero searching
        # beta0 = newton(lambda beta:
        #                chiMTM(beta, wl, Ls, indices, index0, indexs).imag,
        #                x0=max(indices))
        # res = minimize(lambda x: abs(chiMTM(
        #                x[0] + 1j*x[1], wl, Ls, indices, index0, indexs))**2,
        #                x0 = [beta0+0.001, -1E-3], tol=1e-12, method="BFGS")
        # beta = res.x[0] + 1j*res.x[1]
        # residule = abs(chiMTM(beta, wl, Ls, indices, index0, indexs))
        # print(res)
        residule = self.chiMTM(beta)
        betaDiff = beta
        t = 0
        dbeta = 1E-6
        while abs(betaDiff) > 1E-5 and abs(residule) > 1E-10:
            # print(beta, residule)
            t += 1
            fp = (self.chiMTM(beta+dbeta) - self.chiMTM(beta-dbeta))/(2*dbeta)
            betaDiff = residule / fp
            beta = beta - betaDiff
            residule = self.chiMTM(beta)
            if t > 500:
                raise TimeoutError("Doesn't converge")
                break
        return beta

    def populateMode(self, beta, xs):
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
        alpha = np.sqrt((1+0j)*self.indices**2 - beta**2)
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
        alphas = np.sqrt((1+0j)*self.indexs**2 - beta**2)
        if alphas.imag < 0:
            alphas = -alphas
        gammas = alphas/self.indexs**2
        Hx[xs >= lsum[-1]] = Hx0 * exp(1j*alphas * k *
                                       (xs[xs >= lsum[-1]] - lsum[-1]))
        Ez[xs >= lsum[-1]] = gammas*Hx[xs >= lsum[-1]]
        Ey[xs >= lsum[-1]] = -beta*Hx[xs >= lsum[-1]]/self.indexs**2

        scale = Ey[np.argmax(np.abs(Ey))]
        Hx = Hx/scale
        Ey = Ey/scale
        Ez = Ez/scale
        return Ey, Hx, Ez

    def populateIndices(self, xs):
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
        n[xs >= lsum[-1]] = self.indexs
        return n


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
        Angstrom) and 1-st (int) being the nubmer of periods.
        If the 0-th value is 0 or not key not exists, it's not a peirodic
        structure. The variable is only for book keeping and is
        not used to validate Ls within the class
    cstmGain : dict
        A dictionary of customized material, with key the name and value
        the gain coefficient. The variable is only for book keeping and is
        not used for any calculation.
    """
    def __init__(self, wl=3.0, materials=['Air', 'InP'], moleFracs=[0.0, 0.0],
                 dopings=[0.0, 0.0], Ls=[1.0, 1.0], mobilities=None,
                 cstmIndx=dict(), cstmPrd=defaultdict(float),
                 cstmGain=defaultdict(float)):
        super(OptStrata, self).__init__(wl)
        self.materials = list(AlloyNick[m] if m in AlloyNick else m
                              for m in materials)
        self.moleFracs = list(moleFracs)
        self.dopings = list(dopings)
        self.Ls = (np.array(Ls) if len(Ls) == len(materials)
                   else np.array([1.0] + list(Ls) + [2.0]))
        self.mobilities = ([None]*len(materials)
                           if mobilities is None else mobilities)
        self.cstmIndx = cstmIndx
        self.cstmPrd = cstmPrd
        self.cstmGain = cstmGain

        self.plasmonUnit = 8.9698E-5 * self.wl**2
        # 8.9698e-05 = mu_0*1E23*e**2*(1E-6)**2/(electron_mass*4*pi**2)
        # omega_P = plasmonUnit * doping / me / (1 + 1j gamma/omega)
        # doping in unit 1E17 cm^-3,
        # electron effective mass me in unit free electron mass
        self.gammaUnit = 5.3088E-3
        # 5.3088E-3 = 1E7/(2*pi*c)
        # gammaUnit * gamma_c (unit 1E13 s-1) * wl (unit um) = gamma_c/omega

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
        self.plasmonUnit = 8.9698E-5 * self.wl**2

    def insert(self, row, material=None, moleFrac=None, doping=None,
               L=None, mobility=None):
        """Insert a strata indexed with row (top air and bottom substrate
        included) with parameters listed"""
        if row <= 0 or row >= len(self.materials):
            raise IndexError("Cannot inser strata beyond top and substrate.")
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
        """Return the refractive of the strata indexed with row
        (top air and bottom substrate included)"""
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
                self.plasmonUnit * self.dopings[n]/me
                / (1 + 1j*self.gammaUnit * nue)))
        return layerRIdx

    def updateIndices(self):
        """Update indices information according to material info

        Yield
        -----
        indices : list(complex)
            complex refractive index of the stratum
        index0 : complex
            Refractive index of the top (before Ls[0]) strata
        indexs : complex
            refractive index of the substrate (after Ls[-1])
        """
        self.indices = np.array(
            [self.indexOf(n) for n in range(len(self.materials))])
        self.index0 = self.indices[0]
        self.indexs = self.indices[-1]
        self.indices = self.indices[1:-1]

    def populateMtrl(self, xs, mtrlList=None):
        """Populate a boolean array for index slicing on xs for material
        in the mtrlList

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
            if (self.materials[i] in mtrlList if mtrlList else
                    self.materials[i].startswith('Active')):
                res.append((xs > lsum[i]) & (xs < lsum[i+1]))
        return res

# vim: ts=4 sw=4 sts=4 expandtab
