#!/usr/bin/env python3
# -*- coding:utf-8 -*-

import os
import sys

import numpy
from Cython.Build import cythonize
from setuptools import Extension, setup

_CLIB_PREFIX = "ErwinJr2/OneDQuantum"
_CLIB_FILES = [
    "c_schrodinger.pyx",
    "1DSchrodinger.c",
    "band.c",
    "fftautocorr/fftautocorr.c",
]


if os.environ.get("ENABLE_OMP", "false").lower() not in ("false", "0", "f"):
    enable_omp = False
elif sys.platform in ("linux", "win32"):
    enable_omp = True
else:
    enable_omp = False

setup(
    ext_modules=cythonize([
        Extension(
            name="ErwinJr2.OneDQuantum.c_schrodinger",
            sources=[os.path.join(_CLIB_PREFIX, f) for f in _CLIB_FILES],
            include_dirs=[numpy.get_include()],
            define_macros=[
                ("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")
                ] + [("__MP", None)] if enable_omp else [],
            extra_compile_args=[
                "-Ofast",
                # "-Werror",
            ] + ["-fopenmp"] if enable_omp else [],
            extra_link_args=[
                "-lm",
            ] + ["-lgomp"] if enable_omp else [],
        ),
    ], show_all_warnings=True)
)
