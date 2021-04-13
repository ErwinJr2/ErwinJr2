import typing
import numpy as np

floatOrArray = typing.Union[float, np.ndarray]
doubleArray = np.ctypeslib.ndpointer(
    dtype=np.float64, ndim=1, flags="C_CONTIGUOUS")
doubleMatrix = np.ctypeslib.ndpointer(dtype=np.float64,
                                      ndim=2, flags="C_CONTIGUOUS")
intArray = np.ctypeslib.ndpointer(
    dtype=np.int32, ndim=1, flags="C_CONTIGUOUS")
