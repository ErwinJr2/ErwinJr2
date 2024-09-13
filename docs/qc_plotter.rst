qc_plotter module
=================

This module contains functions to plot :class:`ErwinJr2.qc_layers`

.. automodule:: ErwinJr2.qc_plotter
    :members:
    :show-inheritance:



Example
-------

Here is an example of how to construct a QCLayers class and solve for eigen
states and wave functions. A sample *json* file is can be found
:ref:`here <sample_json_file>`.

.. plot::
    :include-source:

    from ErwinJr2 import save_load
    from ErwinJr2.qc_plotter import plotPotential, plotWF
    import matplotlib.pyplot as plt

    with open("../ErwinJr2/example/PQLiu.json", 'r') as f:
        qcl = save_load.qclLoad(f)

    qcl.populate_x()
    qcl.solve_whole()
    plotPotential(qcl)
    plotWF(qcl)
    plt.xlabel('Position (Å)')
    plt.ylabel('Energy (eV)')
    plt.show()
