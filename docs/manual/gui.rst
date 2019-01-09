Graphical User Interface Guide
===================================
.. currentmodule:: QCLayers

The GUI starting program is defined in ``ErwinJr.pyw``. The software can be start by 
``python ErwinJr.pyw [filename]``. It's mostly a GUI wrapper of :py:class:`QCLayers` 
with a plotting canvas. 

.. figure:: ../figures/mainwindow.png

   Screenshot of OneDQ

The interface includes 4 columns: 

===================== =====================
                 settingBox          
===========================================
 Description           Used as a comment, not for calculation 
--------------------- --------------------- 
 Substrate             Decide substrate, can influence strain and material set
--------------------- ---------------------
 E-field               Global electrical field
--------------------- ---------------------
 Position resolution   Finite-element grid size
--------------------- ---------------------
 Energy resolution     Scan size for eigen solve root finder. 
                       This should be smaller than smallest energy difference.
                       If this is too small it's possible to lose some states
--------------------- ---------------------
 Repeats               Number of the whole structure 
--------------------- ---------------------
 Basis Divisions       Defined for basis solver. See :py:meth:`QCLayers.solve_basis`
--------------------- ---------------------
 Period info           Calculate total length, doping density and sheet density
===================== =====================

===================== =====================
                  layerBox
===========================================
 Layer Buttons         insert above layer and delete selected layer
--------------------- ---------------------
 Layer Table           Show the table that defines the layer structure
===================== =====================

===================== =====================
                  solveBox
===========================================
 Solve basis button    Call :py:meth:`QCLayers.solve_basis` and update plot
--------------------- ---------------------
 Solve whole button    Call :py:meth:`QCLayers.solve_whole` and update plot
--------------------- ---------------------
 Material Table        Define the material used in the structure and 
                       give some information about the material
--------------------- ---------------------
 Plot control          - Layer Select: select layer in Layer Table by clicking
                       - Zoom: in with right button of mouse, out with left
                       - Pan: move with mouse
                       - Reset: make the plot back to default position
                       - Clear: remove all wavefunctions
--------------------- ---------------------
 Calculatin box        - Pair Select: start to pick states by mouse 
                       - FoM: calculating full information of the picked state pair
                       - Filter: remove unbounded state smartly (to improve)
===================== =====================

Export of the figure and data, save and load actions are embaded into `File` menu; 
temperature setting and advanced table settings are in `Edit` menu;
options to choose what is included in the plot is listed in `View` menu. 

.. todo::
   desktop shortcut and register to system

