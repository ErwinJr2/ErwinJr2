1D Quantum solver for Quantum Cascade Laser simulation
================

master:
[![master Build Status](https://travis-ci.org/PrincetonUniversity/OneDQ.svg?branch=master)](https://travis-ci.org/PrincetonUniversity/OneDQ)
dev:
[![dev Build Status](https://travis-ci.org/PrincetonUniversity/OneDQ.svg?branch=dev)](https://travis-ci.org/PrincetonUniversity/OneDQ)

![Main Window Screenshot](./docs/figures/mainwindow.png)

This is a extendable Python/C program for 1D quantum problem and Quantum Cascade Laser simulation. 

`OneDQuantum` is A C lib for 1D quantum problem, with python interface. 

## TODO list
- [X] OpenMP support
- [X] Material parameters
- [X] Add Wurtzite crystal structure 
- [ ] ?Add a linear algebra solver using BLAS and/or Lapack
- [ ] ?NEGF solver
- [ ] Add scattering: acoustic phonon, interface roughness, ionized impurity
- [X] Finite temperature Fermi-Dirac distribution
- [X] QCLayer class
- [X] Save/Load using json
- [ ] setup.py update
- [ ] Test case improve
- [X] Documents
- [X] Profile
- [X] ?Travis CI automatic testing
- [ ] ?coveralls.io
- [ ] ?codacy
- [ ] ?CFFI or SWIG
