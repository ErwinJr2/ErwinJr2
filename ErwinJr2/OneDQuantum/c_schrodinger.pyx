#cython: language_level=3
cimport numpy as np

import typing

import numpy as np

from libc.stdint cimport int32_t

floatOrArray = typing.Union[float, np.ndarray]

ctypedef int32_t numpyint
cdef extern from "band.h":
    ctypedef struct Band:
        pass

    Band *ZBband_new(numpyint N, const double *xEg, const double *xF,
                     const double *xEp, const double *xESO)
    void ZBband_free(Band * band)
    Band *WZband_new(numpyint N, const double *xEg,
                     const double *xEp, const double *xESO)
    void WZband_free(Band * band)

cdef double* get_c_ptr_f(double[:] array):
    return &array[0]

cdef numpyint* get_c_ptr_i(numpyint[:] array):
    return &array[0]

cdef class PyBand:
    """Python interface for a Band"""
    cdef str _bandtype
    cdef numpyint _N
    cdef const double* _xEg
    cdef const double* _xEp
    cdef const double* _xF
    cdef const double* _xESO
    cdef Band* _c

    def __cinit__(self, bandtype: str, *args):
        self._bandtype = bandtype
        # Make refs of parameters of band so it's not garbage collected
        if bandtype == "ZincBlende":
            arg_xEg, arg_xF, arg_xEp, arg_xESO = args
            self._N = arg_xEg.size
            self._xEg = get_c_ptr_f(arg_xEg)
            self._xEp = get_c_ptr_f(arg_xEp)
            self._xF = get_c_ptr_f(arg_xF)
            self._xESO = get_c_ptr_f(arg_xESO)
            self._c = ZBband_new(self._N, self._xEg, self._xF, self._xEp, self._xESO)
        elif bandtype == "Wurtzite":
            arg_xEg, arg_xEp, arg_xESO = args
            self._N = arg_xEg.size
            self._xEg = get_c_ptr_f(arg_xEg)
            self._xEp = get_c_ptr_f(arg_xEp)
            self._xESO = get_c_ptr_f(arg_xESO)
            self._c = WZband_new(self._N, self._xEg, self._xEp, self._xESO)
        else:
            raise ValueError(
                "crystal structure \"%s\" not implemented" % bandtype)

    @property
    def bandtype(self) -> str:
        return self._bandtype

    def __del__(self):
        if self._bandtype == "ZincBlende":
            ZBband_free(self._c)
        elif self.bandtype == "Wurtzite":
            WZband_free(self._c)
        else:
            raise ValueError(
                "crystal structure \"%s\" not implemented" % self._bandtype)

cdef extern from "1DSchrodinger.h":
    numpyint Solve1D(double step, numpyint N,
                     const double *Es, numpyint EN,
                     const double *V, double *m, const Band *mat,
                     double *EigenE)

    void FillPsi(double step, numpyint N, const double *EigenEs, numpyint EN,
                 const double *V, double *m, double *psis,
                 numpyint *starts, numpyint *ends, const Band *mat)

    double LOphononScatter(double step, numpyint N, double kl,
                           const double *psi_ij)

    double LOtotal(double step, numpyint N, const double *kls,
                   const double *psi_ijs, const double *fjs, numpyint Nj)

    numpyint invAlpha()
    int isMP()

def cSimpleSolve1D(step: float, Es: np.ndarray, V: np.ndarray, m: floatOrArray,
                   xmin: int = 0, xmax: typing.Optional[int] = None
                   ) -> np.ndarray:
    if not xmax:
        xmax = V.size
    if not isinstance(m, np.ndarray):
        m = m*np.ones(V.size)
    EigenE = np.empty(Es.size)
    EigenEN = Solve1D(step, xmax-xmin, get_c_ptr_f(Es), Es.size,
                      get_c_ptr_f(V[xmin:xmax]), get_c_ptr_f(m[xmin:xmax]),
                      NULL, get_c_ptr_f(EigenE))
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
    starts = np.zeros(xmax-xmin, dtype=np.int32)
    ends = (xmax-xmin)*np.ones(xmax-xmin, dtype=np.int32)
    FillPsi(step, xmax-xmin, get_c_ptr_f(EigenEs), EigenEs.size,
            get_c_ptr_f(V[xmin:xmax]), get_c_ptr_f(m[xmin:xmax]), get_c_ptr_f(psis),
            get_c_ptr_i(starts), get_c_ptr_i(ends), NULL)
    return psis.reshape((EigenEs.size, xmax-xmin))


def cBandSolve1D(step: float, Es: np.ndarray, V: np.ndarray, band: PyBand,
                 xmin: int = 0, xmax: typing.Optional[int] = None
                 ) -> np.ndarray:
    """
    Find eigen energies using band mass.
    """
    if not xmax:
        xmax = V.size
    EigenE = np.empty(Es.size)
    EigenEN = Solve1D(step, xmax-xmin, get_c_ptr_f(Es), Es.size,
                      get_c_ptr_f(V[xmin:xmax]), NULL, band._c,
                      get_c_ptr_f(EigenE))
    return EigenE[:EigenEN]


def cBandFillPsi(step: float, EigenEs: np.ndarray, V: np.ndarray, band: PyBand,
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
    FillPsi(step, xmax-xmin, get_c_ptr_f(EigenEs), EigenEs.size,
            get_c_ptr_f(V[xmin:xmax]), NULL, get_c_ptr_f(psis),
            get_c_ptr_i(starts), get_c_ptr_i(ends), band._c)
    return psis.reshape((EigenEs.size, xmax-xmin))


def cLOphononScatter(step: float, kl: float,
                     psi_ij: np.ndarray, xmin: int = 0,
                     xmax: typing.Optional[int] = None) -> float:
    if not xmax:
        xmax = psi_ij.size
    return LOphononScatter(step, xmax-xmin, kl, get_c_ptr_f(psi_ij))


def cLOtotal(step: float, kls: np.ndarray,
             psi_ijs: np.ndarray, fjs: np.ndarray):
    psi_ijs_flatten = np.ascontiguousarray(psi_ijs).flatten()
    return LOtotal(step, psi_ijs.shape[1], get_c_ptr_f(kls),
                   get_c_ptr_f(psi_ijs_flatten), get_c_ptr_f(fjs),
                   psi_ijs.shape[0])


def cisMP() -> bool:
    return isMP() == 1


if __name__ == "__main__":
    print(invAlpha())
    if isMP():
        print("Woo! OpenMP is loaded!")
