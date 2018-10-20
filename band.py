#!/usr/bin/env python
# -*- coding:utf-8 -*-
import numpy as np
from ctypes import *
_doubleArray = np.ctypeslib.ndpointer(
    dtype=np.float64, ndim=1, flags="C_CONTIGUOUS")

class Band(Structure):
    _fields_ = [("update", c_void_p),
                ("N", c_int),
                ("V", POINTER(c_double)),
                ("m", POINTER(c_double)),
                ("Eg", POINTER(c_double))]

def init(clib):
    clib.UpdateBand.argtypes = [POINTER(Band), c_double]
    clib.UpdateBand.restype = c_int
    def cUpdateBand(band, E):
        return clib.UpdateBand(band, c_double(E))
    clib.ZBband_new.argtypes = [c_int, _doubleArray, _doubleArray, 
                                _doubleArray, _doubleArray, _doubleArray]
    clib.ZBband_new.restype = POINTER(Band)
    def cZBband_new(xEg, xVc, xF, xEp, xESO):
        return clib.ZBband_new(xEg.size, xEg, xVc, xF, xEp, xESO)
    clib.ZBband_free.argtypes = [POINTER(Band)]
    clib.ZBband_free.restype = None
    def cZBband_free(zbband):
        return clib.ZBband_free(zbband)
    return cUpdateBand, cZBband_new, cZBband_free

# vim: ts=4 sw=4 sts=4 expandtab
