Installation Guide for OneDQ
=============================

The software is a portable software, but it still require build 
C library and documents before use. All installation process in 
included in ``install.py`` file. 
For Windows users without gcc and GNU make support, there's also 
a visual studio project file provided as ``OneDQuantum/OneDQuantum.sln``. 
Automatic script in ``install.py`` require ``msbuild=<path to MSBuild>`` 
arguments to specify |MSBuild|_.

.. |MSBuild| replace:: ``MSBuild`` 
.. _MSBuild: https://docs.microsoft.com/en-us/visualstudio/msbuild/msbuild


``install.py`` also try to install openMP support library. 

Automatic dependence check is not implemented yet. Make sure you have 
:py:mod:`numpy`, :py:mod:`scipy` installed 
(:py:mod:`PyQt5` and :py:mod:`matplotlib` for GUI). 

.. todo::
   improve install.py and comments on anaconda
