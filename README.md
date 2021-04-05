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
Anaconda with `numpy`, `scipy` and `matplotlib` is the recommended Python
distribution. To create shortcut on windows, `winshell` is also necessary.
Compiling with Visual Studio and gcc are both supported. Visual Studio project
file and standard `Makefile` are both provided.
To install, start Anaconda Prompt and
go to the directory of ErwinJr2, run the following command:

```
pip install -r requirements.txt
pip install winshell
python install.py --msbuild=[PATH to MSBuild.exe]
```
where `--msbuild=[PATH to MSBuild.exe]` is when using Visual Studio to compile
the C code, `[PATH to MSBuild.exe]` will be Visual Studio version dependent,
for example `'C:\Program Files (x86)\Microsoft Visual Studio\2017\Community\MSBuild\15.0\Bin\MSBuild.exe'`.
For gcc user it's not necessary.

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
CC=gcc-10 python install.py
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
- [X] Add a linear algebra solver
- [ ] ?NEGF solver
- [X] install.py update
- [X] remove unnecessary C lib
- [X] Add IFR scattering
- [ ] Add impurity scattering (may be important for transport)
- [ ] Add finite temperature  (to improve population distribution)
- [X] Add gain spectrum
- [ ] global optimizer for QCLayers
- [ ] optimizer for optical stratum
- [ ] Test case improve:
    - [ ] LO and IFR scattering results
    - [ ] Consistency with and without C lib
    - [ ] Electron population check
- [X] Documents
- [ ] GUI indication of running computation
- [ ] plot style to global settings
- [X] Profile
- [X] Travis CI automatic testing
- [ ] upload to pip as a library
- [ ] ?coveralls.io
- [ ] ?codacy
- [ ] ?CFFI or SWIG
- [X] Provide binary
