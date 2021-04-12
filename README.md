A software for Quantum Cascade Laser design and simulation
================

master:
[![master Build Status](https://travis-ci.com/CareF/ErwinJr2.svg?branch=master)](https://travis-ci.com/CareF/ErwinJr2)
dev:
[![dev Build Status](https://travis-ci.com/CareF/ErwinJr2.svg?branch=dev)](https://travis-ci.com/CareF/ErwinJr2)

![Main Window Screenshot](./docs/figures/qtab.png)

In the following a simple installation guide is included. A more comprehensive
documents can be found [here](https://erwinjr2.readthedocs.io/)


Installation
---------------
The software is based on Python (>=3.6) and uses `setuptools` for installation.
Make sure you have the latest version of `setuptools` installed
(`python` may be `python3` depending on your platform.):

```bash
python -m pip install setuptools --upgrade
```

If you are not installing the software on a remote server, you may also want to
install `PyQt5` for GUI support:

```bash
python -m pip install pyqt5
```

To install the software, at the code directory, run

```bash
python setup.py install
```

If you don't have full control of your system, add `--user` by the end to
install the software in user directory.

After installation, you can start the software by `ErwinJr` command or

```bash
python -m ErwinJr2
```

To create a shortcut on desktop, run

```
ErwinJr-genshortcut
```

### Windows ###
The software requires some compiled components for best performance. The
default compiler depends on GNU `gcc` and `make`, which may not be supported
on Windows. To use [Visual Studio](https://visualstudio.microsoft.com/),
set the `MSBUILD` environment variable to the corresponding directory for
example:
```
set MSBUILD=C:\Program Files (x86)\Microsoft Visual Studio\2017\Community\MSBuild\15.0\Bin\MSBuild.exe
```
Note that the `C:\Program Files (x86)\...` path depends on where you install
Visual Studio on your computer. With the environment variable `setup.py` should
be able to call Visual Studio for the compilation.

### Linux and MacOS ###
By default most Linux distributions have necessary dev-tools installed but for
MacOS if you haven't installed xcode, run the following command"

```bash
xcode-select —install
```

It is recommended to have `openMP` installed for the best performance.
For MacOS specifically, the default `gcc` is a alias to the native `clang`,
we recommend install via [homebrew](https://brew.sh/) for `gcc` and `openMP`
before install

```bash
brew install gcc
```

and use `gcc` as the compiler by

```
CC=gcc-10 python setup.py install
```

where `gcc-10` is the current latest version by depending on your install
it may be other number.

### Build local documentation ###
The software will look for the [online document](https://erwinjr2.readthedocs.io/)
but if you want to build your local version, you need doxygen and:

```bash
cd docs
python -m install -r requirements.txt
make html
```

## TODO list
- [X] OpenMP support
- [X] Add a linear algebra solver
- [ ] ?NEGF solver
- [ ] upload to pip
- [ ] register to OS
- [X] remove unnecessary C lib
- [X] Add IFR scattering
- [ ] Add impurity scattering (may be important for transport)
- [ ] Add finite temperature (to improve population distribution)
- [X] Add gain spectrum
- [X] global optimizer for QCLayers
- [X] optimizer for optical stratum
- [ ] save to excel (to growth sheet)
- [ ] Test case improve:
    - [ ] LO and IFR scattering results
    - [ ] Consistency with and without C lib
    - [ ] Electron population check
- [X] Documents
- [ ] GUI indication of running computation
- [ ] EJcanvas.config to qt setting
- [X] plot style to global settings
- [X] Profile
- [X] Travis CI automatic testing
- [ ] upload to pip as a library
- [ ] ?coveralls.io
- [ ] ?codacy
- [ ] ?CFFI or SWIG
- [X] Provide binary
