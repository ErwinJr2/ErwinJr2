QCLayers module
===============

This module contains `QCLayers` class.

QCLayers.QCLayers class
-----------------------
.. currentmodule:: QCLayers
.. autoclass:: SchrodingerLayer
   :members:
.. autoclass:: QCLayers
   :members:


Example
-------

Here is an example of how to construct a QCLayers class and solve for eigen
states and wave functions. A sample *json* file is can be found
:ref:`here <sample_json_file>`.

.. code-block:: python
   :linenos:

   with open("path/to/file.json") as f:
       qcl = SaveLoad.qclLoad(f)

   qcl.layerSelected = 3
   qcl.NonParabolic = False
   qcl.populate_x()
   qcl.solve_whole()
   qcl.dipole(19, 15)
   qcl.FoM(19, 15)
