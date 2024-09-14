Installation Guide for ErwinJr2
================================

The software is based on Python (>=3.6).
For the Python distributions, it is recommended to use |anaconda|_
for Windows, and |homebrew|_ installed Python3 for macOS.
On this page, the `python` command may be `python3` depending on your platform.

Simple Installation
---------------------

In most cases you can have the software out-of-box from ``pip``, assuming you
already have ``python`` installed.
First to make sure you have the latest ``pip`` installed:

.. code-block:: bash

   python -m pip install pip --upgrade

Then install ``ErwinJr2`` by:

.. code-block:: bash

   python -m pip install ErwinJr2

Now you can start the software via

.. code-block:: bash

   ErwinJr

or if you want to, you can create a shortcut on the desktop via

.. code-block:: bash

   ErwinJr-genshortcut

For reasonably new versions of Windows or macOS with x86_64 CPU, the compiled
library is provided via ``pip`` installation. For Linux or other platforms, it
should be compiled automatically if ``gcc`` and ``make`` are available (which is
the case for almost all Linux distributions).
Otherwise, if the compiled library may not work,
you may see ``C library is not compiled. Features are limited.`` when starting
the software from the command line.

If you are a developer and want to compile the software locally, the following
is a beginner's guide.
The software is available as source code in Github:
https://github.com/ErwinJr2/ErwinJr2 .
You can directly download as a zip file, or using ``git clone`` to clone the
repository, which helps to get further updates.


Run the software
------------------

After installation, you can start the software by ``ErwinJr`` command or

.. code-block:: bash

   python -m ErwinJr2

To create a shortcut on the desktop, run

.. code-block:: bash

   ErwinJr-genshortcut


The following is for advanced users and developers, or for debugging if the
installation fails.


Prepare the environment
------------------------

The software uses ``setuptools`` for installation.
Make sure you have the latest version of `setuptools` installed:

.. code-block:: bash

   python -m pip install setuptools, wheel --upgrade

If you are not installing the software on a remote server, you may also want to
install ``PyQt5`` for GUI support:

.. code-block:: bash

   python -m pip install pyqt5

The software requires some compiled components for best performance.
The default compiler depends on GNU ``gcc`` and ``make``, and optionally ``openMP``.
This compiling environment depends on the Operating System.


Windows
<<<<<<<<

For Windows, the support for GNU compilers may not be easy, but we can use |vs|_
instead. To do so, set the ``MSBUILD`` environment variable to the corresponding
directory for example:

.. code-block::

   set MSBUILD=C:\Program Files (x86)\Microsoft Visual Studio\2017\Community\MSBuild\15.0\Bin\MSBuild.exe

Note that the ``C:\Program Files (x86)\...`` path depends on where you install
Visual Studio on your computer (See more in |MSBuild|_).
With the environment variable, the following command should be able to call
Visual Studio for the compilation.

.. code-block:: bash

   python setup.py install

MacOS
<<<<<<<<

If you haven't installed xcode, run the following command:

.. code-block:: bash

   xcode-select —-install

It is recommended to have ``openMP`` installed for the best performance.
For macOS specifically, the default ``gcc`` is an alias to the native ``clang``,
which does not support ``openMP``.
We recommend install via |homebrew|_ for ``gcc`` before install

.. code-block:: bash

   brew install gcc

and use ``gcc`` as the compiler by

.. code-block:: bash

   CC=gcc-10 python setup.py install

where ``gcc-10`` is the current latest version by depending on your install
it may be another number (like ``gcc-11``).


Linux
<<<<<<<

By default most Linux distributions have necessary dev-tools installed, but
``openMP`` is not necessarily so. If you are using Linux I'm sure you will be
able to install packages via corresponding package management tools :)

To install the software, at the code directory, run

.. code-block:: bash

   python setup.py install

If you don't have full control of your system, add ``--user`` by the end to
install the software in the user directory.



Run the software without installation
--------------------------------------

Sometimes you may want to run the software without installing it to the Python
package directory, especially if you want to change the source code. You can
manually build the C library by

.. code-block:: bash

   cd [PATH_TO_THE_CODE]/ErwinJr2/oned_quantum
   make
   make MP

And run the software via

.. code-block:: bash

   cd [PATH_TO_THE_CODE]
   PYTHONPATH=[PATH_TO_THE_CODE] python ErwinJr2

The installation is basically moving the code set to the Python install path,
so without installation the ``PYTHONPATH`` environment variable is required to
manually add the path to the code so Python can import it.


Build local documentation
--------------------------

The software will look for this online document.
but if you want to build your local version, you need doxygen and:

.. code-block:: bash

   cd docs
   python -m install -r requirements.txt
   make html



.. |MSBuild| replace:: ``MSBuild``
.. _MSBuild: https://docs.microsoft.com/en-us/visualstudio/msbuild/msbuild

.. |homebrew| replace:: ``homebrew``
.. _homebrew: https://brew.sh/

.. |anaconda| replace:: ``Anaconda``
.. _anaconda: https://www.anaconda.com/

.. |vs| replace:: ``Visual Studio``
.. _vs: https://visualstudio.microsoft.com/

.. |MinGW| replace:: ``MinGW``
.. _MinGW: https://www.mingw.org/
