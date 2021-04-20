Introduction
===============

ErwinJr2 is a cross-platform software with a combination of:

1. a C library for solving 1D quantum problem and the related thermal and electrical problems;
2. a Python interface for the C library;
3. a set of Python modules for loading, saving, organizing, solving quantum eigen-states in
   semiconductor quantum wells or super lattices, and calculating relevant physics parameters
   for these states, especially for quantum cascade laser (QCL) design purposes;
4. a GUI front-end for the above Python modules.

The C library is based on ANSI-C and is tested using GCC (under Linux and MacOS) and Visual
Studio (on Windows). It also has optional OpenMP support for parallel computing if the
environment supports.

The C library Python interface is based on :py:mod:`ctypes` in standard Python library and
:py:mod:`numpy`.
The Python module for simulation adds :py:mod:`scipy.constants` requirement for scientific constants,
:py:mod:`scipy.sparse` and :py:mod:`scipy.linalg` for matrix solvers,
and uses :py:mod:`json` for saving and loading super-lattices (or quantum wells) information.

The GUI interface is based on :py:mod:`PyQt5` and :py:mod:`matplotlib`.


About QCL
----------
Quantum cascade lasers (QCLs) are semiconductor lasers that emit light
through inter-subband transitions.
These lasers consist of periodic series of thin
layers of various semiconductor materials which creates a one-dimensional
multiple-quantum-well confinement.
Compare to conventional semiconductor lasers which use single material,
QCLs have the advantage of both a higher output efficiency
due to possible quantum cascades across different quantum wells,
and an improved flexibility in tuning the frequencies.


Models and Formulas
--------------------
The physics model and formulas used in the software are discussed in
:doc:`physics`


.. todo::
   List references and validations
