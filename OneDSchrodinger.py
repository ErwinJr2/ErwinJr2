#!/usr/bin/env python
# -*- coding:utf-8 -*-
import numpy as np
from ctypes import *
_clib = np.ctypeslib.load_library('1DSchrodinger', '.')
_doubleArray = np.ctypeslib.ndpointer(
    dtype=np.float64, ndim=1, flags="C_CONTIGUOUS")

_clib.Numerov.argtypes = [c_double, c_int, c_double, c_double, 
                         c_double, _doubleArray, _doubleArray, _doubleArray]
_clib.Numerov.restype = c_double
def cNumerov(step, y0, y1, E, V, m, xmin=0, xmax=None):
    if not xmax:
        xmax = V.size
    if not isinstance(m, np.ndarray):
        m = m*np.ones(V.size)
    y = np.empty(xmax-xmin)
    yend = _clib.Numerov(c_double(step), xmax-xmin, 
                        c_double(y0), c_double(y1), E, 
                        V[xmin:xmax], m[xmin:xmax], y)
    return yend

_clib.SimpleSolve1D.argtypes = [c_double, c_int, _doubleArray, c_int,
                               _doubleArray, _doubleArray, _doubleArray]
_clib.SimpleSolve1D.restype = c_int
def cSimpleSolve1D(step, Es, V, m, xmin=0, xmax=None): 
    if not xmax:
        xmax = V.size
    if not isinstance(m, np.ndarray):
        m = m*np.ones(V.size)
    EigenE = np.empty(Es.size) 
    EigenEN = _clib.SimpleSolve1D(c_double(step), xmax-xmin, Es, Es.size, 
                                 V[xmin:xmax], m[xmin:xmax], EigenE)
    return EigenE[:EigenEN]

_clib.FillPsi.argtypes = [c_double, c_int, _doubleArray, c_int, 
                         _doubleArray, _doubleArray, _doubleArray]
_clib.FillPsi.restype = None
def cFillPsi(step, EigenEs, V, m, xmin=0, xmax=None): 
    if not xmax:
        xmax = V.size
    if not isinstance(m, np.ndarray):
        m = m*np.ones(V.size)
    psis = np.empty(EigenEs.size*(xmax-xmin))
    _clib.FillPsi(c_double(step), xmax-xmin, EigenEs, EigenEs.size, 
                  V[xmin:xmax], m[xmin:xmax], psis)
    return psis.reshape((EigenEs.size, xmax-xmin))

_clib.CoulombField.argtypes = [c_double, c_int, _doubleArray, _doubleArray, 
                              _doubleArray]
_clib.CoulombField.restype = c_double
def cCoulombField(step, eDensity, eps, xmin=0, xmax=None): 
    if not xmax:
        xmax = eDensity.size
    if not isinstance(eps, np.ndarray):
        eps = eps*np.ones(eDensity.size)
    Vc = np.empty(xmax-xmin)
    _clib.CoulombField(c_double(step), xmax-xmin, eDensity[xmin:xmax],
                      eps[xmin:xmax], Vc)
    return Vc

if __name__ == "__main__":
    print(_clib.invAlpha())

# vim: ts=4 sw=4 sts=4 expandtab
