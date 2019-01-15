C code 
==================

General Structure
-----------------

All C library is intended to be compiled into dynamic link library (
``.so`` for Linux, ``.dylib`` for MacOS and ``.dll`` for Windows). 

:doc:`file/1DSchrodinger_8c` implements quantum physics related functions, 
:doc:`file/1DMaxwell_8c`, is for Maxwell solvers (from electron density to 
potentials) with zero in one side and symmetric boundary conditions. 
:doc:`file/1DThermal_8c` solves thermal distributions. 


Band Structure 
-------------------------------

:doc:file/band_8h.rst includes interfaces for band structure and 
:doc:file/band_8c.rst has two implementation of two different bands: 
Zinc-blende structure and Wurtzite structure. 

``struct BAND`` includes a function pointer for updating parameters 
depending on band structure, an int number for the size of the band and 
a pointer to a double array for bandgap (at different positions). 

Both bands should be regarded as `sub-class` of ``struct BAND``, except 
that there's no syntax in C for inheritance. Therefore if new band 
structure is to be added, it has to be promised that the first three components 
are aligned. At the same time add new interface items in :doc:`../bandpy` .


Memory Management 
------------------------

Dynamic memory allocation happened within C functions will be freed within 
the function. Memory space for return values should be provided and managed 
by callers. And in our cases, by Python modules. 

Memory for band parameters is managed by Python, which promises to call 
allocate and free functions in the C library during construction and 
destruction. 

Users are not expected to use the C code directly. 

File List
-----------------------
.. toctree::
   :glob:

   file/*
