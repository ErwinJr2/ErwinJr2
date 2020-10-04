#!/usr/bin/env python3
# -*- coding:utf-8 -*-
import numpy as np
from ctypes import c_void_p, c_int, c_double, POINTER, Structure
__all__ = ['Band', 'cUpdateBand']
from .typeDefs import doubleArray


class cBand(Structure):
    _fields_ = [("update", c_void_p),
                ("N", c_int),
                ("Eg", POINTER(c_double))]


def init(clib):
    global cUpdateBand, cZBband_new, cZBband_free
    clib.UpdateBand.argtypes = [POINTER(cBand), c_double,
                                doubleArray, doubleArray]
    clib.UpdateBand.restype = c_int

    def cUpdateBand(band: Band, E: float, V: np.ndarray, m: np.ndarray) -> int:
        return clib.UpdateBand(band.c, c_double(E), V, m)

    clib.ZBband_new.argtypes = [c_int, doubleArray, doubleArray,
                                doubleArray, doubleArray]
    clib.ZBband_new.restype = POINTER(cBand)

    def cZBband_new(xEg: np.ndarray, xF: np.ndarray, xEp: np.ndarray,
                    xESO: np.ndarray) -> POINTER(cBand):
        return clib.ZBband_new(xEg.size, xEg, xF, xEp, xESO)

    clib.ZBband_free.argtypes = [POINTER(cBand)]
    clib.ZBband_free.restype = None

    def cZBband_free(zbband: POINTER(cBand)) -> None:
        return clib.ZBband_free(zbband)

    return cUpdateBand, cZBband_new, cZBband_free


class Band(object):
    """Python interface for a Band"""
    bandtype: str
    c: POINTER(cBand)

    def __init__(self, bandtype: str, *args, **kwargs):
        super(Band, self).__init__()
        self.bandtype = bandtype
        # Make refs of parameters of band so it's not garbage collected
        self.args = args
        self.kwargs = kwargs
        if bandtype == "ZincBlende":
            self.c = cZBband_new(*args, **kwargs)
        else:
            raise ValueError(
                "crystal structure \"%s\" not implemented" % bandtype)

    def __del__(self):
        if self.bandtype == "ZincBlende":
            cZBband_free(self.c)
        else:
            raise ValueError(
                "crystal structure \"%s\" not implemented" % self.bandtype)


def cUpdateBand(band: Band, E: float, V: np.ndarray, m: np.ndarray) -> int:
    raise NotImplementedError("Not linked to a proper clib.")


def cZBband_new(xEg: np.ndarray, xF: np.ndarray, xEp: np.ndarray,
                xESO: np.ndarray) -> POINTER(cBand):
    raise NotImplementedError("Not linked to a proper clib.")


def cZBband_free(zbband: POINTER(cBand)) -> None:
    raise NotImplementedError("Not linked to a proper clib.")

# vim: ts=4 sw=4 sts=4 expandtab
