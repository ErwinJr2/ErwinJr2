Introduction
===============

ErwinJr2 is cross-platform software with a combination of:

1. a C library for solving 1D quantum problems and the related thermal and electrical problems;
2. a Python interface for the C library;
3. a set of Python modules for loading, saving, organizing, solving quantum eigenstates in
   semiconductor quantum wells or superlattices, and calculating relevant physics parameters
   for these states, especially for quantum cascade laser (QCL) design purposes;
4. 1D waveguide simulation optimized for QCL design purposes;
5. a GUI front-end for the above Python modules.

The C library is based on ANSI-C and is tested using GCC (on Linux and macOS),
clang (on macOS), and Visual Studio (on Windows).
It also has optional OpenMP support for parallel computing if the environment
supports it.


About QCLs
----------
Quantum cascade lasers (QCLs) are semiconductor lasers that emit light
through inter-subband transitions.
These lasers consist of periodic series of thin
layers of various semiconductor materials which create one-dimensional
multiple-quantum-well confinement.
Compared to conventional semiconductor lasers which use a single material,
QCLs have the advantage of both a higher output efficiency
due to possible quantum cascades across different quantum wells
and improved flexibility in tuning the frequencies.


Package Dependence
-------------------

The C library-Python interface is based on :py:mod:`ctypes` in the standard Python
library and :py:mod:`numpy`.
The Python module for simulation adds :py:mod:`scipy.constants` requirement for
scientific constants,
:py:mod:`scipy.sparse` and :py:mod:`scipy.linalg` for matrix solvers,
and uses :py:mod:`json` for saving and loading super-lattices (or quantum wells) information.

The GUI interface is based on :py:mod:`PyQt5` and :py:mod:`matplotlib`.


Models and Formulas
--------------------
The physics model and formulas used in the software are discussed in
:doc:`physics`


History
--------
The current version is written by Ming from Princeton University, Gmachl group,
as part of his PhD thesis work.

The early versions of the software started as a quantum well eigen-solver using
the shooting algorithm, developed with QBASIC at Bell Labs in the early 1990s.
Ref.~\cite{PhysRevB.50.8663} is one of the early works with the solver.
It was then translated into C with some modern modifications
at Princeton in the 2000s by Dr. Daniel Wasserman,
when the name `ErwinJr` was introduced.
Later it was improved by Dr. Kale Franz
who implemented an early version of the GUI software with Python2 and
`Qwt` as the plotting library, as well as a MATLAB-based version.
It was then collaboratively improved by Dr. Yu Song.
