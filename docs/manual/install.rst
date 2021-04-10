Installation Guide for ErwinJr2
================================

Tge software is available as source code in Github:
https://github.com/CareF/ErwinJr2 .
You can directly donwload as a zip file, or using ``git clone`` to clone the
repository, which helps to get further updates.

The software is designed to be portable, but it recommended to install so that
it compiles the C library and register with the Python environment.

Prepare the environment
------------------------

The software is based on Python (>=3.6) and uses ``setuptools`` for installation.
For the Python distributions, it is recommended to use |anaconda|_
for Windows and |homebrew|_ installed Python for MacOS. After installing,
Make sure you have the latest version of `setuptools` installed
(`python` may be `python3` depending on your platform):

.. code-block:: bash

   python -m pip install setuptools, wheel --upgrade

If you are not installing the software on a remote server, you may also want to
install ``PyQt5`` for GUI support:

.. code-block:: bash

   python -m pip install pyqt5

The software requires some compiled components for best performance.
The default compiler depends on GNU ``gcc`` and ``make``, and optionally ``openMP``.
These compiling environment depends on the Operating System.


Windows
<<<<<<<<

For Windows the support for GNU compilers may not be easy, but we can use |vs|_
instead. To do so, set the ``MSBUILD`` environment variable to the corresponding
directory for example:

.. code-block::

   set MSBUILD=C:\Program Files (x86)\Microsoft Visual Studio\2017\Community\MSBuild\15.0\Bin\MSBuild.exe

Note that the ``C:\Program Files (x86)\...`` path depends on where you install
Visual Studio on your computer (See more in |MSBuild|_).
With the environment variable the following command should be able to call
Visual Studio for the compilation.

.. code-block:: bash

   python setup.py install

MacOS
<<<<<<<<

If you haven't installed xcode, run the following command:

.. code-block:: bash

   xcode-select —-install

It is recommended to have ``openMP`` installed for the best performance.
For MacOS specifically, the default ``gcc`` is a alias to the native ``clang``,
we recommend install via |homebrew|_ for ``gcc`` before install

.. code-block:: bash

   brew install gcc

and use ``gcc`` as the compiler by

.. code-block:: bash

   CC=gcc-10 python setup.py install

where ``gcc-10`` is the current latest version by depending on your install
it may be other number.


Linux
<<<<<<<

By default most Linux distributions have necessary dev-tools installed, but
``openMP`` is not necessarily so. If you are using Linux I'm sure you will be
able to install packages via corresponding package management tools :)

To install the software, at the code directory, run

.. code-block:: bash

   python setup.py install

If you don't have full control of your system, add ``--user`` by the end to
install the software in user directory.


Run the software
------------------

After installation, you can start the software by ``ErwinJr`` command or

.. code-block:: bash

   python -m ErwinJr2

To create a shortcut on desktop, run

.. code-block:: bash

   ErwinJr-genshortcut


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


.. todo::
   pip install and developer mode install