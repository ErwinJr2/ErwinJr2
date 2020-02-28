#!/usr/bin/env python
# -*- coding:utf-8 -*-
import numpy as np
from numpy import pi, cos, sin, sinc, abs, exp
from scipy.optimize import newton, minimize
from ctypes import *
import os
path = os.path.dirname(__file__)

_clib = np.ctypeslib.load_library('1DMaxwell', path)
_doubleArray = np.ctypeslib.ndpointer(
    dtype=np.float64, ndim=1, flags="C_CONTIGUOUS")
__all__ = ['cCoulombField', 'cCoulombField0', 'chiMTM', 'boundModeTM', 'genModeTM']

_clib.CoulombField.argtypes = [c_double, c_int, _doubleArray, _doubleArray, 
                              _doubleArray]
_clib.CoulombField.restype = c_double

def cCoulombField(step, eDensity, eps, xmin=0, xmax=None): 
    """ 
    from e density to coulomb field 
    """
    if not xmax:
        xmax = eDensity.size
    if not isinstance(eps, np.ndarray):
        eps = eps*np.ones(eDensity.size)
    Vc = np.empty(xmax-xmin)
    _clib.CoulombField(c_double(step), xmax-xmin, eDensity[xmin:xmax],
                      eps[xmin:xmax], Vc)
    return Vc

_clib.CoulombField0.argtypes = [c_double, c_int, _doubleArray, _doubleArray, 
                              _doubleArray]
_clib.CoulombField0.restype = c_double

def cCoulombField0(step, eDensity, eps, xmin=0, xmax=None): 
    """
    from e density to Coulomb field
    """
    if not xmax:
        xmax = eDensity.size
    if not isinstance(eps, np.ndarray):
        eps = eps*np.ones(eDensity.size)
    Vc = np.empty(xmax-xmin)
    _clib.CoulombField0(c_double(step), xmax-xmin, eDensity[xmin:xmax],
                      eps[xmin:xmax], Vc)
    return Vc

def transferTM(beta, wl, Ls, indices):
    """
    Calculate the tranfer matrix for the system described with thickness (Ls) and index, 
    wl and Ls should be same unit
    """
    alpha = np.sqrt((1+0j)*indices**2 - beta**2)
    k = 2 * pi / wl
    phi = alpha * k * Ls
    ms = np.moveaxis(np.array([
        [cos(phi), -1j*sin(phi)*alpha/indices**2],
        [-1j*sinc(phi/pi)*indices**2*k * Ls, cos(phi)]]), -1, 0)
    # np.sinc (x) is defined as sin(pi*x)/(pi*x)
    return np.linalg.multi_dot(ms)

def chiMTM(beta, wl, Ls, indices, index0, indexs):
    """
    Calculate model dispersion function for wave omega = c/wl on stratum structure 
    discribed with thickness and index
    """
    gamma0 = np.sqrt((1+0j)*index0**2 - beta**2)/index0**2
    if gamma0.imag < 0:
        gamma0 = -gamma0
    gammas = np.sqrt((1+0j)*indexs**2 - beta**2)/indexs**2
    if gammas.imag < 0:
        gammas = -gammas
    m = transferTM(beta, wl, Ls, indices)
    return m[0,0]*gammas + m[0,1] + m[1,0]*gammas*gamma0 + m[1,1]*gamma0

def boundModeTM(beta, wl, Ls, indices, index0, indexs):
    """
    Solve for bounded mode near beta, with settings descibed by wl, Ls, indices and 
    substrate index0/indexs
    """
    # # Minization algo. for complex function zero searching
    # beta0 = newton(lambda beta: 
    #                chiMTM(beta, wl, Ls, indices, index0, indexs).imag, 
    #                x0=max(indices))
    # res = minimize(lambda x: 
    #                abs(chiMTM(x[0] + 1j*x[1], wl, Ls, indices, index0, indexs))**2, 
    #                x0 = [beta0+0.001, -1E-3], tol=1e-12, method="BFGS")
    # beta = res.x[0] + 1j*res.x[1]
    # residule = abs(chiMTM(beta, wl, Ls, indices, index0, indexs))
    # print(res)
    residule = chiMTM(beta, wl, Ls, indices, index0, indexs)
    t = 0
    dbeta = 1E-6
    while abs(residule) > 1E-10:
        t += 1
        fp = (chiMTM(beta+dbeta, wl, Ls, indices, index0, indexs) - 
            chiMTM(beta-dbeta, wl, Ls, indices, index0, indexs))/(2*dbeta)
        beta = beta - residule / fp
        residule = chiMTM(beta, wl, Ls, indices, index0, indexs)
        if t > 200:
            raise Exception("Doesn't converge")
            break
    return beta

def genModeTM(beta, xs, wl, Ls, indices, index0, indexs):
    """
    Generate TM modes (field) on position array xs, assuming beta is a bounded Mode
    return Ey, Hx, Ez, normalized to max(abs(Ey)) = 1; unit of Hx is 
    ([unit of beta]*sqrt(mu_0/epsilon_0))^-1
    """
    Hx = np.zeros(xs.shape, dtype=np.complex128)
    Ey = np.zeros(xs.shape, dtype=np.complex128)
    Ez = np.zeros(xs.shape, dtype=np.complex128)
    k = 2 * pi / wl

    # top outside medium
    alpha0 = np.sqrt((1+0j)*index0**2 - beta**2)
    if alpha0.imag < 0:
        alpha0 = -alpha0
    gamma0 = alpha0/index0**2
    Ez0 = -gamma0
    Hx0 = 1
    Hx[xs<0] = exp(-1j*alpha0*k*xs[xs<0])
    Ez[xs<0] = -gamma0*Hx[xs<0]
    Ey[xs<0] = -beta*Hx[xs<0]/index0**2

    # middle stratum
    lsum = np.zeros(len(Ls)+1)
    lsum[1:] = np.cumsum(Ls)
    alpha = np.sqrt((1+0j)*indices**2 - beta**2)
    # [[cos(phi), 1j*sin(phi)*alpha/indices**2],
    #  [1j*sinc(phi/pi)*indices**2*k * Ls, cos(phi)]]
    for n in range(len(Ls)):
        msk = (xs>=lsum[n]) & (xs<lsum[n+1])
        phi = alpha[n]*k*(xs[msk] - lsum[n])
        Ez[msk] = cos(phi)*Ez0 + 1j*alpha[n]/indices[n]**2*sin(phi)*Hx0
        Hx[msk] = cos(phi)*Hx0 + 1j*k*(xs[msk] - lsum[n])*sinc(phi/pi)*indices[n]**2*Ez0
        Ey[msk] = -beta*Hx[msk]/indices[n]**2
        phiL = alpha[n]*k*Ls[n]
        Ez0, Hx0 = (cos(phiL)*Ez0 + 1j*alpha[n]/indices[n]**2*sin(phiL)*Hx0, 
            cos(phiL)*Hx0 + 1j*k*Ls[n]*sinc(phiL/pi)*indices[n]**2*Ez0)

    # last substrate
    alphas = np.sqrt((1+0j)*indexs**2 - beta**2)
    if alphas.imag < 0:
        alphas = -alphas
    gammas = alphas/indexs**2
    Hx[xs>=lsum[-1]] = Hx0 * exp(1j*alphas*k*(xs[xs>=lsum[-1]]-lsum[-1]))
    Ez[xs>=lsum[-1]] = gammas*Hx[xs>=lsum[-1]]
    Ey[xs>=lsum[-1]] = -beta*Hx[xs>=lsum[-1]]/indexs**2

    scale = -np.max(np.abs(Ey))
    Hx = Hx/scale
    Ey = Ey/scale
    Ez = Ez/scale
    return Ey, Hx, Ez


if __name__ == "__main__":
    print(_clib.speedOfLight())

# vim: ts=4 sw=4 sts=4 expandtab
