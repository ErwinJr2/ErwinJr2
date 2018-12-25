#!/usr/bin/env python
# -*- coding:utf-8 -*-
import numpy as np
from ctypes import *
import os
path = os.path.dirname(__file__)
_clib = np.ctypeslib.load_library('1DThermal', path)
_doubleArray = np.ctypeslib.ndpointer(
    dtype=np.float64, ndim=1, flags="C_CONTIGUOUS")
_doubleMatrix = np.ctypeslib.ndpointer(dtype=np.float64, 
                                       ndim=2, flags="C_CONTIGUOUS")
__all__ = ['cFermiDirac0', 'cFermiDirac0N', 'cFermiDirac', 'cBoltzmann', 'cBoltzmannN']

_clib.FermiDirac0.argtypes = [c_double, _doubleArray, c_int,
                              _doubleArray, _doubleMatrix, c_int, c_double,
                              _doubleArray]
_clib.FermiDirac0.restype = c_double
def cFermiDirac0(EF, EigenEs, m, psis, step):
    if not isinstance(m, np.ndarray):
        m = m*np.ones(psis.shape[1])
    eDensity = np.empty(psis.shape[1])
    _clib.FermiDirac0(c_double(EF), EigenEs, EigenEs.size, 
                      m, psis, psis.shape[1], c_double(step), eDensity)
    return eDensity

_clib.FermiDirac0N.argtypes = [c_double, _doubleArray, c_int,
                               _doubleArray, _doubleMatrix, c_int, c_double,
                               _doubleArray]
_clib.FermiDirac0N.restype = c_double
def cFermiDirac0N(sheet, EigenEs, m, psis, step):
    if not isinstance(m, np.ndarray):
        m = m*np.ones(psis.shape[1])
    eDensity = np.empty(psis.shape[1])
    EF = _clib.FermiDirac0N(c_double(sheet), EigenEs, EigenEs.size, 
                           m, psis, psis.shape[1], c_double(step), eDensity)
    return eDensity, EF

_clib.FermiDirac.argtypes = [c_double, c_double, _doubleArray, c_int,
                              _doubleArray, _doubleMatrix, c_int, c_double,
                              _doubleArray]
_clib.FermiDirac.restype = c_double
def cFermiDirac(T, EF, EigenEs, m, psis, step):
    if not isinstance(m, np.ndarray):
        m = m*np.ones(psis.shape[1])
    eDensity = np.empty(psis.shape[1])
    _clib.FermiDirac(c_double(T), c_double(EF), EigenEs, EigenEs.size, 
                      m, psis, psis.shape[1], c_double(step), eDensity)
    return eDensity

#_clib.FermiDiracN.argtypes = [c_double, c_double, _doubleArray, c_int, 
#                              _doubleArray, _doubleMatrix, c_int, c_double, 
#                              _doubleArray]
#_clib.FermiDiracN.restype = c_double
#def cFermiDiracN(T, sheet, EigenEs, m, psis, step):
#    if not isinstance(m, np.ndarray):
#        m = m*np.ones(psis.shape[1])
#    eDensity = np.empty(psis.shape[1])
#    EF = _clib.FermiDiracN(c_double(T), c_double(sheet), EigenEs, EigenEs.size, 
#                           m, psis, psis.shape[1], c_double(step), eDensity)
#    return eDensity, EF

_clib.Boltzmann.argtypes = [c_double, c_double, _doubleArray, c_int, 
                            _doubleArray, _doubleMatrix, 
                            c_int, c_double, _doubleArray]
_clib.Boltzmann.restype = c_double
def cBoltzmann(T, EF, EigenEs, m, psis, step):
    if not isinstance(m, np.ndarray):
        m = m*np.ones(psis.shape[1])
    eDensity = np.empty(psis.shape[1])
    _clib.Boltzmann(c_double(T), c_double(EF), EigenEs, 
                    EigenEs.size, m, psis, psis.shape[1], 
                    c_double(step), eDensity)
    return eDensity


def cFermiDirac0(EF, EigenEs, m, psis, step):
    if not isinstance(m, np.ndarray):
        m = m*np.ones(psis.shape[1])
    eDensity = np.empty(psis.shape[1])
    _clib.FermiDirac0(c_double(EF), EigenEs, EigenEs.size,
                      m, psis, psis.shape[1], c_double(step), eDensity)
    return eDensity




_clib.BoltzmannN.argtypes = [c_double, c_double, _doubleArray, c_int,
                              _doubleArray, _doubleMatrix, c_int,
                             c_double, _doubleArray]
_clib.BoltzmannN.restype = c_double

def cBoltzmannN(T, sheet, EigenEs, m, psis, step):
    if not isinstance(m, np.ndarray):
        m = m*np.ones(psis.shape[1])
    eDensity = np.empty(psis.shape[1])
    EF = _clib.BoltzmannN(c_double(T), c_double(sheet), EigenEs,
                           EigenEs.size, m, psis, psis.shape[1], 
                           c_double(step), eDensity)
    return eDensity, EF

if __name__ == "__main__":
    print(_clib.answer())
# vim: ts=4 sw=4 sts=4 expandtab
