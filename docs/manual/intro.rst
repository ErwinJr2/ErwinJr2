Introduction
===============

OneDQ is a corss-platform software with a combination of: 

1. A C library for solving 1D quantum problem and correlated thermal and electrical problem
2. A Python interface for the C library
3. A set of Python module for loading, saving, organizing, solving quantum states in 
   semiconductor quantum wells or super lattices, and calculating relevant physics parameters 
   for these states, especially for quantum cascade laser (QCL) design purposes. 
4. A GUI front-end interface for the above Python module

The C library is based on standard C and is tested using gcc(under Linux and MacOS) or 
visual studio(on windows). It also has optional openmp support for parallel computing if the 
environment has openmp support. 

The C library Python interface is based on :py:mod:`ctypes` in standard Python library and 
:py:mod:`numpy`. 
The Python module for simulation added :py:mod:`scipy.constants` requirement for scientific constants, 
and uses :py:mod:`json` for saving and loading super-lattices (or quantum wells) information. 

The GUI interface is based on :py:mod:`PyQt5` and :py:mod:`matplotlib`. 


About QCL
----------
Quantum cascade lasers (QCLs) are semiconductor lasers that emit light
through inter-subband transitions.
These lasers consist of periodic series of thin
layers of various semiconductor materials which creates a one-dimensional
multiple quantum well confinement.
Compare to conventional semiconductor lasers which use single material,
QCLs have the advantage of both a higher output efficiency
due to possible quantum cascades across different quantum wells,
and an improved flexibility in tuning the frequencies.
