#!/usr/bin/env python
# -*- coding:utf-8 -*-
import numpy as np
from ctypes import *
import os
path = os.path.dirname(__file__)

_clib = np.ctypeslib.load_library('1DMaxwell', path)
_doubleArray = np.ctypeslib.ndpointer(
    dtype=np.float64, ndim=1, flags="C_CONTIGUOUS")
__all__ = ['cCoulombField', 'cCoulombField0']

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

if __name__ == "__main__":
    print(_clib.speedOfLight())

# vim: ts=4 sw=4 sts=4 expandtab
