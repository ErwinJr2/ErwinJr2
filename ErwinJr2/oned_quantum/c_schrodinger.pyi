import typing

import numpy as np

floatOrArray = typing.Union[float, np.ndarray]


class PyBand:
    bandtype: str
    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def __del__(self, *args, **kwargs) -> None: ...
    def __reduce__(self): ...


def cSimpleSolve1D(step: float, Es: np.ndarray, V: np.ndarray, m: floatOrArray,
                   xmin: int = 0, xmax: typing.Optional[int] = None
                   ) -> np.ndarray: ...


def cSimpleFillPsi(step: float, EigenEs: np.ndarray, V: np.ndarray,
                   m: floatOrArray, xmin: int = 0,
                   xmax: typing.Optional[int] = None) -> np.ndarray:
    """
    Find wave functions. Assume mass as given.
    """


def cBandSolve1D(step: float, Es: np.ndarray, V: np.ndarray, band: PyBand,
                 xmin: int = 0, xmax: typing.Optional[int] = None
                 ) -> np.ndarray:
    """
    Find eigen energies using band mass.
    """


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


def cLOphononScatter(step: float, kl: float,
                     psi_ij: np.ndarray, xmin: int = 0,
                     xmax: typing.Optional[int] = None) -> float: ...


def cLOtotal(step: float, kls: np.ndarray,
             psi_ijs: np.ndarray, fjs: np.ndarray): ...


def cisMP() -> bool: ...
