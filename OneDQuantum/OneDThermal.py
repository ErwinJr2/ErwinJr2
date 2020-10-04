#!/usr/bin/env python3
# -*- coding:utf-8 -*-
import numpy as np
from ctypes import c_int, c_double
import os
import typing
from .typeDefs import floatOrArray, doubleArray, doubleMatrix
path = os.path.dirname(__file__)
_clib = np.ctypeslib.load_library('1DThermal', path)
__all__ = ['cFermiDirac0', 'cFermiDirac0N',
           'cFermiDirac', 'cFermiDiracN',
           'cBoltzmann', 'cBoltzmannN']

_clib.FermiDirac0.argtypes = [c_double, doubleArray, c_int,
                              doubleArray, doubleMatrix, c_int, c_double,
                              doubleArray]
_clib.FermiDirac0.restype = c_double
_clib.FermiDirac0N.argtypes = [c_double, doubleArray, c_int,
                               doubleArray, doubleMatrix, c_int, c_double,
                               doubleArray]
_clib.FermiDirac0N.restype = c_double
_clib.FermiDirac.argtypes = [c_double, c_double, doubleArray, c_int,
                             doubleArray, doubleMatrix, c_int, c_double,
                             doubleArray]
_clib.FermiDirac.restype = c_double
_clib.FermiDiracN.argtypes = [c_double, c_double, doubleArray, c_int,
                              doubleArray, doubleMatrix, c_int, c_double,
                              doubleArray]
_clib.FermiDiracN.restype = c_double
_clib.Boltzmann.argtypes = [c_double, c_double, doubleArray, c_int,
                            doubleArray, doubleMatrix,
                            c_int, c_double, doubleArray]
_clib.Boltzmann.restype = c_double
_clib.BoltzmannN.argtypes = [c_double, c_double, doubleArray, c_int,
                             doubleArray, doubleMatrix, c_int,
                             c_double, doubleArray]
_clib.BoltzmannN.restype = c_double


def cFermiDirac0(EF: float, EigenEs: np.ndarray, m: floatOrArray,
                 psis: np.ndarray, step: float) -> np.ndarray:
    """0T Fermi-Dirac: Fermi energy to sheet density"""
    if not isinstance(m, np.ndarray):
        m = m*np.ones(psis.shape[1])
    eDensity = np.empty(psis.shape[1])
    _clib.FermiDirac0(c_double(EF), EigenEs, EigenEs.size,
                      m, psis, psis.shape[1], c_double(step), eDensity)
    return eDensity


def cFermiDirac0N(sheet: float, EigenEs: np.ndarray, m: floatOrArray,
                  psis: np.ndarray, step: float
                  ) -> typing.Tuple[np.ndarray, float]:
    """0T Fermi-Dirac: sheet density to eDensity, EF"""
    if not isinstance(m, np.ndarray):
        m = m*np.ones(psis.shape[1])
    eDensity = np.empty(psis.shape[1])
    EF = _clib.FermiDirac0N(c_double(sheet), EigenEs, EigenEs.size,
                            m, psis, psis.shape[1], c_double(step), eDensity)
    return eDensity, EF


def cFermiDirac(T: float, EF: float, EigenEs: np.ndarray, m: floatOrArray,
                psis: np.ndarray, step: float) -> np.ndarray:
    """Finite temperature Fermi-Dirac: Fermi energy to electron density"""
    if not isinstance(m, np.ndarray):
        m = m*np.ones(psis.shape[1])
    eDensity = np.empty(psis.shape[1])
    _clib.FermiDirac(c_double(T), c_double(EF), EigenEs, EigenEs.size,
                     m, psis, psis.shape[1], c_double(step), eDensity)
    return eDensity


def cFermiDiracN(T: float, sheet: float, EigenEs: np.ndarray, m: floatOrArray,
                 psis: np.ndarray, step: float
                 ) -> typing.Union[np.ndarray, float]:
    """Finite-temperature Fermi-Dirac: sheet density to electron density, EF"""
    if not isinstance(m, np.ndarray):
        m = m*np.ones(psis.shape[1])
    eDensity = np.empty(psis.shape[1])
    EF = _clib.FermiDiracN(c_double(T), c_double(sheet), EigenEs, EigenEs.size,
                           m, psis, psis.shape[1], c_double(step), eDensity)
    return eDensity, EF


def cBoltzmann(T: float, EF: float, EigenEs: np.ndarray, m: floatOrArray,
               psis: np.ndarray, step: float) -> np.ndarray:
    """Maxwell-Boltzmann: Fermi energy to electron density"""
    if not isinstance(m, np.ndarray):
        m = m*np.ones(psis.shape[1])
    eDensity = np.empty(psis.shape[1])
    _clib.Boltzmann(c_double(T), c_double(EF), EigenEs,
                    EigenEs.size, m, psis, psis.shape[1],
                    c_double(step), eDensity)
    return eDensity


def cBoltzmannN(T: float, sheet: float, EigenEs: np.ndarray, m: floatOrArray,
                psis: np.ndarray, step: float
                ) -> typing.Union[np.ndarray, float]:
    """Maxwell-Boltzmann: sheet density to Fermi energy"""
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
