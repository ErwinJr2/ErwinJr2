Installation Guide for OneDQ
=============================

The software is a portable software, but it still requires an 
installation procedure to build C library and documents before 
use. All installation process is included in ``install.py`` file. 

For Windows users without GCC and GNU make support, there is also 
a Visual Studio project file provided as ``OneDQuantum/OneDQuantum.sln``. 
Automatic script in ``install.py`` requires ``msbuild=<path to MSBuild>`` 
arguments to specify |MSBuild|_. 

.. |MSBuild| replace:: ``MSBuild`` 
.. _MSBuild: https://docs.microsoft.com/en-us/visualstudio/msbuild/msbuild


The file ``install.py`` also tries to install OpenMP support library. If the 
environment does not support OpenMP, depending on the platform, there may be
error messages but it can be ignored. The software is fully functional 
without OpenMP. 

Automatic dependence check is not implemented yet, due to limited cross-platform 
support for ``pip``. Make sure you have :py:mod:`numpy`, :py:mod:`scipy` 
installed for CLI and :py:mod:`PyQt5` and :py:mod:`matplotlib` for GUI. 
For Windows users, |anaconda|_ is recommended as the Python distribution.  

.. |anaconda| replace:: ``anaconda``
.. _anaconda: https://www.anaconda.com/

.. todo::
   improve install.py and comments on anaconda
