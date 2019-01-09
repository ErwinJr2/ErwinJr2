SaveLoad module
===============


.. automodule:: SaveLoad
   :members:


.. _sample_json_file:
      
Here is a sample json file

.. code-block:: json
   :linenos:

   {
       "FileType": "ErwinJr2 Data File", 
       "Version": "181107", 
       "Description": "dx.doi.org/10.1038/nphoton.2009.262", 
       "Substrate": "InP", 
       "EField": 102.0, 
       "x resolution": 1.0, 
       "E resolution": 0.5, 
       "Solver": "ODE", 
       "Temperature": 300.0, 
       "Repeats": 3, 
       "Materials": {
           "Compostion": ["InGaAs", "AlInAs"], 
           "Mole Fraction": [0.66, 0.31]
       }, 
       "QC Layers": {
           "Material": [0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1], 
           "Width": [28.0, 26.0, 22.0, 21.0, 18.0, 18.0, 15.0, 13.0, 12.0, 10.0, 42.0, 12.0, 39.0, 14.0, 33.0, 23.0], 
           "Doping": [0.0, 1.5, 1.5, 1.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 
           "Active Region": [true, false, false, false, false, false, true, true, true, true, true, true, true, true, true, false]
       }
   }
