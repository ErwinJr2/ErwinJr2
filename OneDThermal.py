#!/usr/bin/env python
# -*- coding:utf-8 -*-
import numpy as np
from ctypes import *
_clib = np.ctypeslib.load_library('1DThermal', '.')
_doubleArray = np.ctypeslib.ndpointer(
    dtype=np.float64, ndim=1, flags="C_CONTIGUOUS")
_doubleMatrix = np.ctypeslib.ndpointer(dtype=np.float64, 
                                       ndim=2, flags="C_CONTIGUOUS"), 

_clib.FermiDirac0.argtypes = [c_double, _doubleArray, c_int,
                              _doubleArray, _doubleArray, _doubleMatrix,
                              c_int, c_double, _doubleArray]
_clib.FermiDirac0.restype = c_double
def cFermiDirac0(EF, EigenEs, m, psis, step):
    if not isinstance(m, np.ndarray):
        m = m*np.ones(psis.shape[1])
    eDensity = np.empty(xmax-xmin)
    _clib.FermiDirac0(c_double(EF), EigenEs, EigenEs.size, 
                      m, psis, psis.shape[1], c_double(step), eDensity)
    return eDensity

_clib.Boltzmann.argtypes = [c_double, c_double, _doubleArray, c_int,
                              _doubleArray, _doubleArray, _doubleMatrix,
                              c_int, c_double, _doubleArray]
_clib.Boltzmann.restype = c_double
def cBoltzmann(T, EF, EigenEs, m, psis, step):
    if not isinstance(m, np.ndarray):
        m = m*np.ones(psis.shape[1])
    eDensity = np.empty(xmax-xmin)
    _clib.cBoltzmann(c_double(T), c_double(EF), EigenEs, EigenEs.size, 
                      m, psis, psis.shape[1], c_double(step), eDensity)
    return eDensity

if __name__ == "__main__":
    print(_clib.answer())
# vim: ts=4 sw=4 sts=4 expandtab
