qclayers module
================

.. currentmodule:: ErwinJr2.qclayers

This module contains :class:`QCLayers` and :class:`SchrodingerLayer` classes.
The physics model can be seen in :doc:`manual/physics_quantum`

:class:`SchrodingerLayer` is a Schrodinger solver together with transition
properties evaluation,
:class:`QCLayers` is a wrapper of :class:`SchrodingerLayer` with material
information coded.

QCLayers class
---------------
.. autoclass:: QCLayers
   :members:

SchrodingerLayer class
-----------------------
.. autoclass:: SchrodingerLayer
   :members:


Others
---------------
.. automodule:: ErwinJr2.qclayers
   :members:
   :show-inheritance:
   :exclude-members: QCLayers, SchrodingerLayer

Example
-------

Here is an example of how to construct a QCLayers class and solve for eigen
states and wave functions. A sample *json* file is can be found
:ref:`here <sample_json_file>`.

.. code-block:: python
   :linenos:

   with open("path/to/file.json") as f:
       qcl = save_load.qclLoad(f)

   qcl.layerSelected = 3
   qcl.populate_x()
   qcl.solve_whole()
   qcl.dipole(19, 15)
   qcl.FoM(19, 15)
