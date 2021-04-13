import numpy as np
from ctypes import c_void_p, c_int, c_double, POINTER, Structure, CDLL
from typing import Callable, Union
__all__ = ['Band', 'cBandUpdateM', 'cBandNormalize']
from .typeDefs import doubleArray


class cBand(Structure):
    _fields_ = [("updateM", c_void_p),
                ("normalize", c_void_p),
                ("N", c_int),
                ("Eg", POINTER(c_double))]


CBAND_P = POINTER(cBand)


class Band(object):
    """Python interface for a Band"""
    bandtype: str
    c: CBAND_P

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


def cBandUpdateM(band: Band, E: float, V: np.ndarray, m: np.ndarray) -> int:
    raise NotImplementedError("Not linked to a proper clib.")


def cBandNormalize(band: Band, E: float, V: np.ndarray,
                   psi: np.ndarray, xres: float) -> float:
    raise NotImplementedError("Not linked to a proper clib.")


def cZBband_new(xEg: np.ndarray, xF: np.ndarray, xEp: np.ndarray,
                xESO: np.ndarray) -> CBAND_P:
    raise NotImplementedError("Not linked to a proper clib.")


def cZBband_free(zbband: CBAND_P) -> None:
    raise NotImplementedError("Not linked to a proper clib.")


def init(clib: CDLL) -> Union[
        Callable[[Band, float, np.ndarray, np.ndarray], int],
        Callable[[np.ndarray, np.ndarray, np.ndarray, np.ndarray], CBAND_P],
        Callable[[CBAND_P], None]]:
    global cBandUpdateM, cBandNormalize, cZBband_new, cZBband_free
    clib.BandUpdateM.argtypes = [POINTER(cBand), c_double,
                                 doubleArray, doubleArray]
    clib.BandUpdateM.restype = c_int
    clib.BandNormalize.argtypes = [POINTER(cBand), c_double, doubleArray,
                                   doubleArray, c_double]

    def cBandUpdateM(band: Band, E: float, V: np.ndarray, m: np.ndarray
                     ) -> int:
        return clib.BandUpdateM(band.c, c_double(E), V, m)

    def cBandNormalize(band: Band, E: float, V: np.ndarray,
                       psi: np.ndarray, xres: float) -> float:
        return clib.BandNormalize(band.c, c_double(E), V, psi, xres)

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

    return cBandUpdateM, cZBband_new, cZBband_free
