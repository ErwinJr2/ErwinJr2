1D Quantum solver for Quantum Cascade Laser simulation
================

master:
[![master Build Status](https://travis-ci.com/CareF/ErwinJr2.svg?branch=master)](https://travis-ci.org/CareF/ErwinJr2)
dev:
[![dev Build Status](https://travis-ci.com/CareF/ErwinJr2.svg?branch=dev)](https://travis-ci.org/CareF/ErwinJr2)

![Main Window Screenshot](./docs/figures/mainwindow.png)

This is a extendable Python/C program for 1D quantum problem and Quantum Cascade Laser simulation. 

`OneDQuantum` is a C lib for 1D quantum problem, with python interface. 

In the following a simple installation guide is included. A more comprehensive 
documents can be found [here](https://erwinjr2.readthedocs.io/)


Installation
---------------
The software is based on Python3, and dependent packages are listed in 
`requirements.txt`: 
`doxygen` and `OpenMP` are optional. For Windows users Visual Studio as compiler
is also supported. 

The software requires some compilation to be installed. 
To compile, run the Python script `install.py`. 

### Windows ###
Anaconda with `numpy`, `scipy` and `matplotlib` is recommended. 
To get best performance, using Visual Studio as compiler is recommended: 
go to the directory of ErwinJr2, run the following command in Anaconda Prompt

```
pip install -r requirements.txt
python install.py --msbuild=[PATH to MSBuild.exe]
```
where `[PATH to MSBuild.exe]` will be Visual Studio version dependent, 
for example `'C:\Program Files (x86)\Microsoft Visual Studio\2017\Community\MSBuild\15.0\Bin\MSBuild.exe'`

### Linux and MacOS ###
```
pip3 install -r requirements.txt
python install.py
```
For MacOS specifically, to enable multi-processing with `OpenMP`, `gcc` is 
required. By default `gcc` in MacOS command line is actually `clang`, and to 
call real `gcc` we need to specify the version number, so 
the second command needs to be, 
```
CC=gcc-9 python install.py
```
or to modify the `OneDQuantum/Makefile` to include the `OpenMP` 
library for `clang`. 

### Build local documentation ###
During the installation you will be asked if you want to build local 
documentation. It's optional but if you choose yes, dependencies for building 
`sphinx` documentation is required: run the following command before
`install.py` script.
```
pip install -r docs/requirements.txt
```

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
