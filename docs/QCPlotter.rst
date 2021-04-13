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

.. plot::
    :include-source:

    from ErwinJr2 import SaveLoad
    from ErwinJr2.QCPlotter import plotPotential, plotWF
    import matplotlib.pyplot as plt

    with open("../ErwinJr2/example/PQLiu.json", 'r') as f:
        qcl = SaveLoad.qclLoad(f)

    qcl.populate_x()
    qcl.solve_whole()
    plotPotential(qcl)
    plotWF(qcl)
    plt.xlabel('Position (Å)')
    plt.ylabel('Energy (eV)')
    plt.show()