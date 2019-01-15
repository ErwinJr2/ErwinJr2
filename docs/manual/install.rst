Installation Guide for OneDQ
=============================

The software is a portable software, but it still requires an 
installation procedure to build C library and documents before 
use. All installation process in included in ``install.py`` file. 

For Windows users without gcc and GNU make support, there's also 
a visual studio project file provided as ``OneDQuantum/OneDQuantum.sln``. 
Automatic script in ``install.py`` require ``msbuild=<path to MSBuild>`` 
arguments to specify |MSBuild|_. 

.. |MSBuild| replace:: ``MSBuild`` 
.. _MSBuild: https://docs.microsoft.com/en-us/visualstudio/msbuild/msbuild


``install.py`` also tries to install openMP support library. If the 
environment doesn't support openMP, depending on the platform, there may be
error messages but it can be ignored. The software is fully functional 
without openMP. 

Automatic dependence check is not implemented yet, due to limited cross-platform 
support for ``pip``. Make sure you have :py:mod:`numpy`, :py:mod:`scipy` 
installed for CLI and :py:mod:`PyQt5` and :py:mod:`matplotlib` for GUI. 
For Windows users, |anaconda|_ is recommended as the Python distribution.  

.. |anaconda| replace:: ``anaconda``
.. _anaconda: https://www.anaconda.com/

.. todo::
   improve install.py and comments on anaconda
