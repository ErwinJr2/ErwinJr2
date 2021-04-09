Profiling
=========

> This page is very out-of-date and has some inconsistencies with current versions

The following block of code is profiled by :py:mod:`cProfile`:

.. literalinclude:: ../tool/profQCLayers.py
   :language: python



where *PQLiu.json* is the following file:

.. literalinclude:: ../ErwinJr2/example/PQLiu.json
   :language: json


Current version with population and gain
-----------------------------------------

=========  =======  =======  =======  ======= ===============================================================
   ncalls  tottime  percall  cumtime  percall filename:lineno(function)
=========  =======  =======  =======  ======= ===============================================================
      1    0.001    0.001    5.955    5.955   ./ErwinJr2/QCLayers.py:1493(auto_gain)
      1    0.000    0.000    4.753    4.753   ./ErwinJr2/QCLayers.py:1386(full_population)
      1    0.024    0.024    4.753    4.753   ./ErwinJr2/QCLayers.py:921(full_population)
   3640    1.358    0.000    3.361    0.001   ./ErwinJr2/QCLayers.py:701(_ifr_transition)
  14428    0.014    0.000    1.545    0.000   ./ErwinJr2/QCLayers.py:433(psi_overlap)
  14436    0.316    0.000    1.531    0.000   {built-in method builtins.sum}
   3640    0.038    0.000    1.367    0.000   ./ErwinJr2/QCLayers.py:621(_lo_transition)
  57712    0.292    0.000    1.216    0.000   ./ErwinJr2/QCLayers.py:437(<genexpr>)
  45069    0.090    0.000    1.024    0.000   ./ErwinJr2/QCLayers.py:273(_shift_psi)
  12164    0.012    0.000    0.935    0.000   <__array_function__ internals>:2(interp)
  12164    0.015    0.000    0.914    0.000   [.python-packages]/numpy/lib/function_base.py:1289(interp)
  12164    0.871    0.000    0.871    0.000   {built-in method numpy.core._multiarray_umath.interp}
   1754    0.783    0.000    0.800    0.000   ./ErwinJr2/OneDQuantum/OneDSchrodinger.py:144(cLOphononScatter)
      1    0.017    0.017    0.678    0.678   ./ErwinJr2/QCLayers.py:1394(full_gain_spectrum)
   1785    0.078    0.000    0.656    0.000   ./ErwinJr2/QCLayers.py:587(_dipole)
      1    0.002    0.002    0.513    0.513   ./ErwinJr2/QCLayers.py:322(solve_whole)
      1    0.036    0.036    0.510    0.510   ./ErwinJr2/QCLayers.py:344(_solve_whole_ode)
      1    0.469    0.469    0.469    0.469   ./ErwinJr2/OneDQuantum/OneDSchrodinger.py:101(cBandSolve1D)
 343980    0.372    0.000    0.372    0.000   ./ErwinJr2/QCLayers.py:732(interpZ)
 116482    0.062    0.000    0.355    0.000   <__array_function__ internals>:2(argmax)
   3570    0.281    0.000    0.281    0.000   ./ErwinJr2/QCLayers.py:993(_xBandMassInv)
   7179    0.005    0.000    0.248    0.000   <__array_function__ internals>:2(trapz)
   7179    0.169    0.000    0.237    0.000   [.python-packages]/numpy/lib/function_base.py:4006(trapz)
 116482    0.053    0.000    0.234    0.000   [.python-packages]/numpy/core/fromnumeric.py:1114(argmax)
 116483    0.059    0.000    0.182    0.000   [.python-packages]/numpy/core/fromnumeric.py:52(_wrapfunc)
   3570    0.002    0.000    0.138    0.000   <__array_function__ internals>:2(gradient)
   3570    0.101    0.000    0.129    0.000   [.python-packages]/numpy/lib/function_base.py:802(gradient)
 116482    0.103    0.000    0.103    0.000   {method 'argmax' of 'numpy.ndarray' objects}
   7181    0.005    0.000    0.063    0.000   {method 'sum' of 'numpy.ndarray' objects}
 232976    0.058    0.000    0.058    0.000   ./ErwinJr2/QCLayers.py:1247(layerVc)
   7181    0.003    0.000    0.058    0.000   [.python-packages]/numpy/core/_methods.py:45(_sum)
   7210    0.056    0.000    0.056    0.000   {method 'reduce' of 'numpy.ufunc' objects}
  12164    0.006    0.000    0.021    0.000   <__array_function__ internals>:2(iscomplexobj)
=========  =======  =======  =======  ======= ===============================================================



1.0 and earlier versions
-------------------------

Before optimizing code performamce, there are 814603 function calls taking 10.835
seconds in total. The output of :py:mod:`cProfile` sorted by cumulative time is
listed in the following table.


=========  =======  =======  =======  ======= ===============================================================
   ncalls  tottime  percall  cumtime  percall filename:lineno(function)
=========  =======  =======  =======  ======= ===============================================================
        1    0.000    0.000   10.803   10.803 QCLayers.py:425(calc_FoM)
       35    9.742    0.278   10.803    0.309 QCLayers.py:383(loTransition)
        2    0.000    0.000   10.449    5.224 QCLayers.py:419(loLifeTime)
        2    0.000    0.000   10.449    5.224 QCLayers.py:422(<listcomp>)
   135301    0.132    0.000    0.797    0.000 fromnumeric.py:1821(sum)
   135303    0.155    0.000    0.640    0.000 fromnumeric.py:64(_wrapreduction)
   135331    0.467    0.000    0.467    0.000 {method 'reduce' of 'numpy.ufunc' objects}
   135180    0.262    0.000    0.262    0.000 {built-in method builtins.abs}
        1    0.000    0.000    0.037    0.037 QCLayers.py:260(solve_whole)
        1    0.036    0.036    0.036    0.036 OneDSchrodinger.py:61(cSimpleSolve1D)
   135313    0.025    0.000    0.025    0.000 {built-in method builtins.isinstance}
   135303    0.019    0.000    0.019    0.000 {method 'items' of 'dict' objects}
        1    0.001    0.001    0.002    0.002 QCLayers.py:178(populate_x)
       65    0.000    0.000    0.001    0.000 QCLayers.py:169(avghwLO)
        1    0.001    0.001    0.001    0.001 OneDSchrodinger.py:74(cSimpleFillPsi)
      227    0.000    0.000    0.001    0.000 {built-in method builtins.sum}
       28    0.001    0.000    0.001    0.000 QCLayers.py:227(<listcomp>)
     1885    0.001    0.000    0.001    0.000 QCLayers.py:173(<genexpr>)
=========  =======  =======  =======  ======= ===============================================================

.. currentmodule:: QCLayers



After optimizing code performance by replace most time consuming function call
(:py:meth:`QCLayers.loLifeTime`), by C code and add openMP support,
there are 3613 function calls taking 1.255 seconds in total.
The output of :py:mod:`cProfile` sorted by cumulative time is listed in the
following table.

=========  =======  =======  =======  ======= ===============================================================
   ncalls  tottime  percall  cumtime  percall filename:lineno(function)
=========  =======  =======  =======  ======= ===============================================================
        1    0.000    0.000    1.205    1.205 QCLayers.py:427(calc_FoM)
       35    0.003    0.000    1.205    0.034 QCLayers.py:383(loTransition)
       30    1.199    0.040    1.199    0.040 OneDSchrodinger.py:113(cLOphononScatter)
        2    0.000    0.000    1.166    0.583 QCLayers.py:421(loLifeTime)
        2    0.000    0.000    1.165    0.583 QCLayers.py:424(<listcomp>)
        1    0.000    0.000    0.048    0.048 QCLayers.py:260(solve_whole)
        1    0.046    0.046    0.046    0.046 OneDSchrodinger.py:65(cSimpleSolve1D)
        1    0.001    0.001    0.002    0.002 QCLayers.py:178(populate_x)
       91    0.000    0.000    0.001    0.000 fromnumeric.py:1821(sum)
       65    0.000    0.000    0.001    0.000 QCLayers.py:169(avghwLO)
      227    0.000    0.000    0.001    0.000 {built-in method builtins.sum}
        1    0.001    0.001    0.001    0.001 OneDSchrodinger.py:78(cSimpleFillPsi)
       93    0.000    0.000    0.001    0.000 fromnumeric.py:64(_wrapreduction)
      121    0.001    0.000    0.001    0.000 {method 'reduce' of 'numpy.ufunc' objects}
     1885    0.001    0.000    0.001    0.000 QCLayers.py:173(<genexpr>)
       68    0.000    0.000    0.001    0.000 ctypeslib.py:196(from_param)
       28    0.000    0.000    0.000    0.000 QCLayers.py:227(<listcomp>)
       65    0.000    0.000    0.000    0.000 {built-in method numpy.core.multiarray.array}
=========  =======  =======  =======  ======= ===============================================================
