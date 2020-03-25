#!/usr/bin/env python3
# -*- coding:utf-8 -*-
import numpy as np
from ctypes import c_int, c_double, POINTER
from . import band as _bd
import os
path = os.path.dirname(__file__)
__all__ = ['cSimpleSolve1D', 'cSimpleFillPsi',
           'Band', 'cBandFillPsi', 'cBandSolve1D',
           'cLOphononScatter', 'cLOtotal',
           'cBandSolve1DBonded']

_doubleArray = np.ctypeslib.ndpointer(
    dtype=np.float64, ndim=1, flags="C_CONTIGUOUS")
_intArray = np.ctypeslib.ndpointer(
    dtype=np.int64, ndim=1, flags="C_CONTIGUOUS")


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
    #                            c_double, _doubleArray, _doubleArray,
    #                            _doubleArray]
    #  _clib.Numerov.restype = c_double

    _clib.Solve1D.argtypes = [c_double, c_int, _doubleArray, c_int,
                              _doubleArray, _doubleArray,
                              POINTER(cBand), _doubleArray]
    _clib.Solve1D.restype = c_int

    _clib.Solve1DBonded.argtypes = [c_double, c_int, c_double, c_double,
                                    c_double, _doubleArray, c_int,
                                    _doubleArray, _doubleArray,
                                    POINTER(cBand), _doubleArray]
    _clib.Solve1DBonded.restype = c_int

    _clib.FillPsi.argtypes = [c_double, c_int, _doubleArray, c_int,
                              _doubleArray, _doubleArray, _doubleArray,
                              _intArray, _intArray, POINTER(cBand)]
    _clib.FillPsi.restype = None

    _clib.LOphononScatter.argtypes = [c_double, c_int, c_double,
                                      _doubleArray, _doubleArray]
    _clib.LOphononScatter.restype = c_double

    _clib.LOtotal.argtypes = [c_double, c_int, _doubleArray, _doubleArray,
                              _doubleArray, _doubleArray, c_int]
    _clib.LOtotal.restype = c_double


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
    EigenEN = _clib.Solve1D(c_double(step), xmax-xmin, Es, Es.size,
                            V[xmin:xmax], m[xmin:xmax], None, EigenE)
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
                  V[xmin:xmax], m[xmin:xmax], psis,
                  np.zeros(xmax-xmin, dtype=np.int64),
                  (xmax-xmin)*np.ones(xmax-xmin, dtype=np.int64),
                  None)
    return psis.reshape((EigenEs.size, xmax-xmin))


def cBandSolve1D(step, Es, V, band, xmin=0, xmax=None):
    """
    Find eigen energies using band mass.
    """
    if not xmax:
        xmax = V.size
    EigenE = np.empty(Es.size)
    EigenEN = _clib.Solve1D(c_double(step), xmax-xmin, Es, Es.size, 
                            V[xmin:xmax], np.empty(0), band.c, EigenE)
    return EigenE[:EigenEN]


def cBandFillPsi(step, EigenEs, V, band, xmin=0, xmax=None,
                 Elower=None, Eupper=None, field=None):
    """
    Find wave functions using band mass.
    """
    if not xmax:
        xmax = V.size
    psis = np.empty(EigenEs.size*(xmax-xmin))
    if field is not None:
        starts = np.floor((Elower - EigenEs)/(field*step*1E-5)
                          ).astype(np.int64)
        starts[starts < xmin] = xmin
        ends = np.ceil((Eupper - EigenEs)/(field*step*1E-5)
                       ).astype(np.int64)
        ends[ends > xmax] = xmax
    else:
        starts = np.zeros(xmax-xmin, dtype=np.int64)
        ends = (xmax-xmin)*np.ones(xmax-xmin, dtype=np.int64)
    _clib.FillPsi(c_double(step), xmax-xmin, EigenEs, EigenEs.size,
                  V[xmin:xmax], np.empty(0), psis,
                  starts, ends, band.c)
    return psis.reshape((EigenEs.size, xmax-xmin))


def cBandSolve1DBonded(step, Es, Elower, Eupper, field,
                       V, band, xmin=0, xmax=None):
    """
    Find eigen energies using band mass.
    bonded by [Elower-field*x, Eupper-field*x]
    """
    if not xmax:
        xmax = V.size
    EigenE = np.empty(Es.size)
    EigenEN = _clib.Solve1DBonded(c_double(step), xmax-xmin, Elower, Eupper,
                                  field, Es, Es.size, V[xmin:xmax],
                                  np.empty(0), band.c, EigenE)
    return EigenE[:EigenEN]


def cLOphononScatter(step, kl, psi_i, psi_j, xmin=0, xmax=None):
    if not xmax:
        xmax = psi_i.size
    return _clib.LOphononScatter(c_double(step), xmax-xmin, kl,
                                 psi_i, psi_j)


def cLOtotal(step, kls, psi_i, psi_js, fjs):
    return _clib.LOtotal(c_double(step), len(psi_i), kls, psi_i,
                         psi_js.flatten(), fjs, len(psi_js))


if __name__ == "__main__":
    print(_clib.invAlpha())

# vim: ts=4 sw=4 sts=4 expandtab
