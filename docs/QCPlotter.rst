QCPlotter module
=================

This module contains functions to plot :class:`ErwinJr2.QCLayers`

.. automodule:: ErwinJr2.QCPlotter
    :members:
    :show-inheritance:



Example
-------

Here is an example of how to construct a QCLayers class and solve for eigen
states and wave functions. A sample *json* file is can be found
:ref:`here <sample_json_file>`.

.. code-block:: python
    :linenos:

    from ErwinJr2 import SaveLoad
    from ErwinJr2.QCPlotter import plotPotential, scaleWF, plotW
    from pylab import *

    with open("path/to/file.json") as f:
        qcl = SaveLoad.qclLoad(f)

    qcl.populate_x()
    qcl.solve_whole()
    axes = gca()
    plotPotential(axes, qcl)
    plotWF(axes, qcl)
