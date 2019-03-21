.. OneDQ documentation master file, created by
   sphinx-quickstart on Thu Dec  6 20:04:16 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to ErwinJr2's documentation!
=================================

ErwinJr2 is a cross-platform software for modeling quantum problems in 1D semiconductor
super lattices. It is specifically optimized for Quantum Cascade Laser (QCL) simulation, 
but coded in a way that's easy to use for general problems about static state in finite
and periodic quantum wells. 

See `Github <https://github.com/PrincetonUniversity/ErwinJr2/>`_ for the source code. 

.. figure:: figures/mainwindow.png

   Screenshot of ErwinJr2


.. toctree::
   :caption: User Manual

   manual/intro.rst
   manual/install.rst
   manual/gui.rst
   manual/cli.rst
   ack.rst

.. toctree::
   :maxdepth: 2
   :caption: For Developers

   OneDQuantum_description
   code_struct
   other
   Profiling
   QCLayers
   Material
   SaveLoad
   OneDQuantum
   clib/filelist

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
