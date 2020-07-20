Installation Guide for ErwinJr2
===============================

Tge software is available as source code in Github: 
https://github.com/CareF/ErwinJr2 .
You can directly donwload as a zip file, or using `git` to clone the 
repository, which helps to get further updates. 

The software is a portable software, but it still requires an 
installation procedure to build C library and documents before 
use. All installation process is included in ``install.py`` file. 

The script ``install.py`` do the following things: 

1. Compile the C library for :mod:`OneDQuantum` 
2. Tries the OpenMP version. The software is fully functional without OpenMP. 
3. (Optional) Build the documentation locally
4. (Optional) Create a shortcut on desktop

The dependence packages are listed in `requirements.txt`. For windows users 
specifically `winshell` is also required to create a desktop shortcut. 

Install under Windows
-----------------------

To install under Windows, Python is required, as well as a C compiler. 
For Python environment, |anaconda|_ is recommended; for C compiler, 
|vs|_ is supported. Using GNU gcc and make is sometimes tricky on Windows, 
but let me know if you have a neat solution. 

If you don't have a C compiler, you may try the binary release_, and extract
the dll files in the zip file to `OneDQuantum/`

.. _release: https://github.com/CareF/ErwinJr2/releases

For |vs|_ project file is provided as ``OneDQuantum/OneDQuantum.sln``. 
Automatic script in ``install.py`` requires ``msbuild=<path to MSBuild>`` 
arguments to specify |MSBuild|_. 

The precedure to install the software is in Anaconda Prompt, go to the 
directory of ErwinJr2, run the following command

.. code-block:: bash

   pip install -r requirements.txt
   pip install winshell
   python install.py --msbuild=[PATH to MSBuild.exe]

Install under MacOS and Linux
-------------------------------

.. code-block:: bash

   pip3 install -r requirements.txt
   python install.py

For MacOS specifically, to enable multi-processing with ``OpenMP``, ``gcc`` is 
recommended but not provided in MacOS by default. If you have gcc installed via
for example homebrew, than you may replace the second command with the 
following:

.. code-block:: bash

   CC=gcc-9 python install.py

.. |MSBuild| replace:: ``MSBuild`` 
.. _MSBuild: https://docs.microsoft.com/en-us/visualstudio/msbuild/msbuild

.. |anaconda| replace:: ``Anaconda``
.. _anaconda: https://www.anaconda.com/

.. |vs| replace:: ``Visual Studio``
.. _vs: https://visualstudio.microsoft.com/

.. |MinGW| replace:: ``MinGW``
.. _MinGW: https://www.mingw.org/
