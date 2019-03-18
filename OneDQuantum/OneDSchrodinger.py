#!/usr/bin/env python
# -*- coding:utf-8 -*-
import numpy as np
from ctypes import *
from . import band as _bd
import os
path = os.path.dirname(__file__)
__all__ = ['cSimpleSolve1D', 'cSimpleFillPsi', 
           'Band', 'cBandFillPsi', 'cBandSolve1D']

_doubleArray = np.ctypeslib.ndpointer(
    dtype=np.float64, ndim=1, flags="C_CONTIGUOUS")

def bindOpenMP(on=True):
    """
    set OpenMP. Default = True
    """
    global _clib
    if(on):
        _clib = np.ctypeslib.load_library('1DSchrodinger_MP', path)
    else:
        _clib = np.ctypeslib.load_library('1DSchrodinger', path)

    _bd.init(_clib)
    global Band
    from .band import cBand, Band
    #  _clib.Numerov.argtypes = [c_double, c_int, c_double, c_double, 
    #                           c_double, _doubleArray, _doubleArray, _doubleArray]
    #  _clib.Numerov.restype = c_double

    _clib.SimpleSolve1D.argtypes = [c_double, c_int, _doubleArray, c_int,
                                   _doubleArray, _doubleArray, _doubleArray]
    _clib.SimpleSolve1D.restype = c_int

    _clib.FillPsi.argtypes = [c_double, c_int, _doubleArray, c_int, 
                             _doubleArray, _doubleArray, _doubleArray, 
                              POINTER(cBand)]
    _clib.FillPsi.restype = None

    _clib.BandSolve1D.argtypes = [c_double, c_int, _doubleArray, c_int,
                                   _doubleArray, POINTER(cBand), _doubleArray]
    _clib.BandSolve1D.restype = c_int

    _clib.LOphononScatter.argtypes = [c_double, c_int, c_double, 
                                      _doubleArray, _doubleArray]
    _clib.LOphononScatter.restype = c_double

bindOpenMP(False)

#  def cNumerov(step, y0, y1, E, V, m, xmin=0, xmax=None):
#      if not xmax:
#          xmax = V.size
#      if not isinstance(m, np.ndarray):
#          m = m*np.ones(V.size)
#      y = np.empty(xmax-xmin)
#      yend = _clib.Numerov(c_double(step), xmax-xmin, 
#                          c_double(y0), c_double(y1), E, 
#                          V[xmin:xmax], m[xmin:xmax], y)
#      return yend

def cSimpleSolve1D(step, Es, V, m, xmin=0, xmax=None): 
    """
    Find eigen energies. Assume mass as given.
    """
    if not xmax:
        xmax = V.size
    if not isinstance(m, np.ndarray):
        m = m*np.ones(V.size)
    EigenE = np.empty(Es.size) 
    EigenEN = _clib.SimpleSolve1D(c_double(step), xmax-xmin, Es, Es.size, 
                                 V[xmin:xmax], m[xmin:xmax], EigenE)
    return EigenE[:EigenEN]

def cSimpleFillPsi(step, EigenEs, V, m, xmin=0, xmax=None): 
    """
    Find wave functions. Assume mass as given.
    """
    if not xmax:
        xmax = V.size
    if not isinstance(m, np.ndarray):
        m = m*np.ones(V.size)
    psis = np.empty(EigenEs.size*(xmax-xmin))
    _clib.FillPsi(c_double(step), xmax-xmin, EigenEs, EigenEs.size, 
                  V[xmin:xmax], m[xmin:xmax], psis, None)
    return psis.reshape((EigenEs.size, xmax-xmin))

def cBandSolve1D(step, Es, V, band, xmin=0, xmax=None): 
    """
    Find eigen energies using band mass.
    """
    if not xmax:
        xmax = V.size
    EigenE = np.empty(Es.size) 
    EigenEN = _clib.BandSolve1D(c_double(step), xmax-xmin, Es, Es.size, 
                                 V[xmin:xmax], band.c, EigenE)
    return EigenE[:EigenEN]

def cBandFillPsi(step, EigenEs, V, band, xmin=0, xmax=None): 
    """
    Find wave functions using band mass.
    """
    if not xmax:
        xmax = V.size
    psis = np.empty(EigenEs.size*(xmax-xmin))
    _clib.FillPsi(c_double(step), xmax-xmin, EigenEs, EigenEs.size, 
                  V[xmin:xmax], None, psis, band.c)
    return psis.reshape((EigenEs.size, xmax-xmin))

def cLOphononScatter(step, kl, psi_i, psi_j, xmin=0, xmax=None):
    if not xmax:
        xmax = psi_i.size
    return  _clib.LOphononScatter(c_double(step), xmax-xmin, kl, 
                                  psi_i, psi_j)

if __name__ == "__main__":
    print(_clib.invAlpha())

# vim: ts=4 sw=4 sts=4 expandtab
