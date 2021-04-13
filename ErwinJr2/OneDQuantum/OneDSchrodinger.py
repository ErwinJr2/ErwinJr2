#!/usr/bin/env python3
# -*- coding:utf-8 -*-
import numpy as np
from ctypes import c_int, c_double, POINTER, CDLL
import typing
from . import band as _bd
from .typeDefs import doubleArray, intArray, floatOrArray
import os
path = os.path.dirname(__file__)
__all__ = ['cSimpleSolve1D', 'cSimpleFillPsi',
           'Band', 'cBandFillPsi', 'cBandSolve1D',
           'cLOphononScatter', 'cLOtotal']


def bindOpenMP(on: bool = True) -> typing.Tuple[CDLL, typing.Type[object]]:
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

    _clib.Solve1D.argtypes = [c_double, c_int, doubleArray, c_int,
                              doubleArray, doubleArray,
                              POINTER(cBand), doubleArray]
    _clib.Solve1D.restype = c_int

    _clib.FillPsi.argtypes = [c_double, c_int, doubleArray, c_int,
                              doubleArray, doubleArray, doubleArray,
                              intArray, intArray, POINTER(cBand)]
    _clib.FillPsi.restype = None

    _clib.LOphononScatter.argtypes = [c_double, c_int, c_double, doubleArray]
    _clib.LOphononScatter.restype = c_double

    _clib.LOtotal.argtypes = [c_double, c_int, doubleArray, doubleArray,
                              doubleArray, c_int]
    _clib.LOtotal.restype = c_double
    return _clib, Band


_clib, Band = bindOpenMP(False)

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


def cSimpleSolve1D(step: float, Es: np.ndarray, V: np.ndarray, m: floatOrArray,
                   xmin: int = 0, xmax: typing.Optional[int] = None
                   ) -> np.ndarray:
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


def cSimpleFillPsi(step: float, EigenEs: np.ndarray, V: np.ndarray,
                   m: floatOrArray, xmin: int = 0,
                   xmax: typing.Optional[int] = None) -> np.ndarray:
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
                  np.zeros(xmax-xmin, dtype=np.int32),
                  (xmax-xmin)*np.ones(xmax-xmin, dtype=np.int32),
                  None)
    return psis.reshape((EigenEs.size, xmax-xmin))


def cBandSolve1D(step: float, Es: np.ndarray, V: np.ndarray, band: Band,
                 xmin: int = 0, xmax: typing.Optional[int] = None
                 ) -> np.ndarray:
    """
    Find eigen energies using band mass.
    """
    if not xmax:
        xmax = V.size
    EigenE = np.empty(Es.size)
    EigenEN = _clib.Solve1D(c_double(step), xmax-xmin, Es, Es.size,
                            V[xmin:xmax], np.empty(0), band.c, EigenE)
    return EigenE[:EigenEN]


def cBandFillPsi(step: float, EigenEs: np.ndarray, V: np.ndarray, band: Band,
                 xmin: int = 0, xmax: typing.Optional[int] = None,
                 Elower: typing.Optional[float] = None,
                 Eupper: typing.Optional[float] = None,
                 field: typing.Optional[float] = None) -> np.ndarray:
    """
    Find wave functions using band mass. `field`, `Elower` and `Eupper` is used
    only for bound the energy range of the wave functions: outside the bound
    the wavefunction is promised to be zero.
    """
    if not xmax:
        xmax = V.size
    psis = np.empty(EigenEs.size*(xmax-xmin))
    if field is not None:
        starts = np.floor((Elower - EigenEs)/(field*step*1E-5)
                          ).astype(np.int32)
        starts[starts < xmin] = xmin
        ends = np.ceil((Eupper - EigenEs)/(field*step*1E-5)
                       ).astype(np.int32)
        ends[ends > xmax] = xmax
    else:
        starts = np.zeros(xmax-xmin, dtype=np.int32)
        ends = (xmax-xmin)*np.ones(xmax-xmin, dtype=np.int32)
    _clib.FillPsi(c_double(step), xmax-xmin, EigenEs, EigenEs.size,
                  V[xmin:xmax], np.empty(0), psis,
                  starts, ends, band.c)
    return psis.reshape((EigenEs.size, xmax-xmin))


def cLOphononScatter(step: float, kl: float,
                     psi_ij: np.ndarray, xmin: int = 0,
                     xmax: typing.Optional[int] = None) -> float:
    if not xmax:
        xmax = psi_ij.size
    return _clib.LOphononScatter(c_double(step), xmax-xmin,
                                 c_double(kl), psi_ij)


def cLOtotal(step: float, kls: np.ndarray,
             psi_ijs: np.ndarray, fjs: np.ndarray):
    return _clib.LOtotal(c_double(step), psi_ijs.shape[1], kls,
                         psi_ijs.flatten(), fjs, psi_ijs.shape[0])


def isMP() -> bool:
    return _clib.isMP() == 1


if __name__ == "__main__":
    print(_clib.invAlpha())
    if isMP():
        print("Woo! OpenMP is loaded!")
