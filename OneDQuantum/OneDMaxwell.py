#!/usr/bin/env python3
# -*- coding:utf-8 -*-
import numpy as np
# from scipy.optimize import newton, minimize
from ctypes import c_double, c_int
import os
import typing
from .typeDefs import doubleArray
path = os.path.dirname(__file__)

_clib = np.ctypeslib.load_library('1DMaxwell', path)
__all__ = ['cCoulombField', 'cCoulombField0']

_clib.CoulombField.argtypes = [c_double, c_int, doubleArray, doubleArray,
                               doubleArray]
_clib.CoulombField.restype = c_double


def cCoulombField(step: float, eDensity: np.ndarray, eps: np.ndarray,
                  xmin: int = 0, xmax: typing.Optional[int] = None
                  ) -> np.ndarray:
    """from e density to coulomb field, assuming 0 field on the left"""
    if not xmax:
        xmax = eDensity.size
    if not isinstance(eps, np.ndarray):
        eps = eps*np.ones(eDensity.size)
    Vc = np.empty(xmax-xmin)
    _clib.CoulombField(c_double(step), xmax-xmin, eDensity[xmin:xmax],
                       eps[xmin:xmax], Vc)
    return Vc


_clib.CoulombField0.argtypes = [c_double, c_int, doubleArray, doubleArray,
                                doubleArray]
_clib.CoulombField0.restype = c_double


def cCoulombField0(step: float, eDensity: np.ndarray, eps: np.ndarray,
                   xmin: int = 0, xmax: typing.Optional[int] = None
                   ) -> np.ndarray:
    """from e density to Coulomb field, assuming no external field"""
    if not xmax:
        xmax = eDensity.size
    if not isinstance(eps, np.ndarray):
        eps = eps*np.ones(eDensity.size)
    Vc = np.empty(xmax-xmin)
    _clib.CoulombField0(c_double(step), xmax-xmin, eDensity[xmin:xmax],
                        eps[xmin:xmax], Vc)
    return Vc


if __name__ == "__main__":
    print(_clib.speedOfLight())

# vim: ts=4 sw=4 sts=4 expandtab
