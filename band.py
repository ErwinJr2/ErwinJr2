#!/usr/bin/env python
# -*- coding:utf-8 -*-
import numpy as np
from ctypes import *
_doubleArray = np.ctypeslib.ndpointer(
    dtype=np.float64, ndim=1, flags="C_CONTIGUOUS")

class Band(Structure):
    _fields_ = [("update", c_void_p),
                ("N", c_int),
                ("V", _doubleArray),
                ("m", _doubleArray),
                ("Eg", _doubleArray)]

class ZBband(Structure):
    _fields_ = [("update", c_void_p),
                ("N", c_int),
                ("V", _doubleArray),
                ("m", _doubleArray),
                ("Eg", _doubleArray),
                ('xF', _doubleArray), 
                ('xEp', _doubleArray), 
                ('xESO', _doubleArray)]
def getZBband(clib):
    clib.ZBband_new.argtypes = [c_int, _doubleArray, _doubleArray, 
                                _doubleArray, _doubleArray, _doubleArray]
    clib.ZBband_new.restype = POINTER(ZBband)
    def cZBband_new(xEg, xVc, xF, xEp, xESO):
        return clib.ZBband_new(xEg.size, xEg, xVc, xF, xEp, xESO)
    clib.ZBband_free.argtypes = [POINTER(ZBband)]
    clib.ZBband_free.restype = None
    def cZBband_free(zbband):
        return clib.ZBband_free(zbband)
    return cZBband_new, cZBband_free

# vim: ts=4 sw=4 sts=4 expandtab
