Profiling
=========

The following block of code is profiled by cProfile.

.. code-block:: python
   :linenos:

   with open("../example/16um.json") as f:
       qcl = SaveLoad.qclLoad(f)

   qcl.layerSelected = 3
   qcl.NonParabolic = False
   qcl.populate_x()
   qcl.solve_whole()
   qcl.dipole(19, 15)
   qcl.FoM(19, 15)




Before optimizing code performamce, there are 814603 function calls in 52.102
seconds. The output of cProfile sorted by cumulative time is listed in the
following table.



=========  =======  =======  =======  ======= ===============================================================
   ncalls  tottime  percall  cumtime  percall filename:lineno(function)
=========  =======  =======  =======  ======= ===============================================================
        1    0.000    0.000   52.037   52.037 QCLayers.py:425(FoM)
       35   50.011    1.429   52.037    1.487 QCLayers.py:383(loTransition)
        2    0.000    0.000   50.303   25.151 QCLayers.py:419(loLifeTime)
        2    0.000    0.000   50.303   25.151 QCLayers.py:422(<listcomp>)
   135301    0.279    0.000    1.608    0.000 fromnumeric.py:1821(sum)
   135303    0.271    0.000    1.276    0.000 fromnumeric.py:64(_wrapreduction)
   135331    0.989    0.000    0.989    0.000 {method 'reduce' of 'numpy.ufunc' objects}
   135180    0.416    0.000    0.416    0.000 {built-in method builtins.abs}
        1    0.000    0.000    0.061    0.061 QCLayers.py:260(solve_whole)
        1    0.058    0.058    0.058    0.058 OneDSchrodinger.py:57(cSimpleSolve1D)
   135313    0.053    0.000    0.053    0.000 {built-in method builtins.isinstance}
   135303    0.017    0.000    0.017    0.000 {method 'items' of 'dict' objects}
        1    0.003    0.003    0.003    0.003 OneDSchrodinger.py:67(cSimpleFillPsi)
        1    0.001    0.001    0.003    0.003 QCLayers.py:178(populate_x)
       65    0.000    0.000    0.001    0.000 QCLayers.py:169(avghwLO)
      227    0.000    0.000    0.001    0.000 {built-in method builtins.sum}
       28    0.001    0.000    0.001    0.000 QCLayers.py:227(<listcomp>)
     1885    0.001    0.000    0.001    0.000 QCLayers.py:173(<genexpr>)
        1    0.000    0.000    0.000    0.000 SaveLoad.py:8(qclLoad)
       65    0.000    0.000    0.000    0.000 {built-in method numpy.core.multiarray.array}
        1    0.000    0.000    0.000    0.000 QCLayers.py:79(__init__)
        1    0.000    0.000    0.000    0.000 QCLayers.py:104(update_strain)
        1    0.000    0.000    0.000    0.000 QCLayers.py:365(dipole)
        1    0.000    0.000    0.000    0.000 QCLayers.py:106(<listcomp>)
        2    0.000    0.000    0.000    0.000 Material.py:171(__init__)
       50    0.000    0.000    0.000    0.000 {built-in method numpy.core.multiarray.empty}
        9    0.000    0.000    0.000    0.000 Material.py:65(set_temperature)
       30    0.000    0.000    0.000    0.000 QCLayers.py:410(<listcomp>)
        5    0.000    0.000    0.000    0.000 Material.py:51(__init__)
        1    0.000    0.000    0.000    0.000 __init__.py:243(load)
        8    0.000    0.000    0.000    0.000 ctypeslib.py:196(from_param)
        1    0.000    0.000    0.000    0.000 function_base.py:25(linspace)
        2    0.000    0.000    0.000    0.000 fromnumeric.py:49(_wrapfunc)
        2    0.000    0.000    0.000    0.000 Material.py:182(set_temperature)
        2    0.000    0.000    0.000    0.000 Material.py:198(set_molefrac)
        2    0.000    0.000    0.000    0.000 function_base.py:1079(diff)
        1    0.000    0.000    0.000    0.000 __init__.py:271(loads)
        1    0.000    0.000    0.000    0.000 fromnumeric.py:2092(cumsum)
        1    0.000    0.000    0.000    0.000 decoder.py:334(decode)
        1    0.000    0.000    0.000    0.000 {built-in method io.open}
        2    0.000    0.000    0.000    0.000 Material.py:95(set_strain)
        1    0.000    0.000    0.000    0.000 fromnumeric.py:36(_wrapit)
        2    0.000    0.000    0.000    0.000 {built-in method numpy.core.multiarray.arange}
        1    0.000    0.000    0.000    0.000 decoder.py:345(raw_decode)
        3    0.000    0.000    0.000    0.000 {built-in method numpy.core.multiarray.concatenate}
       68    0.000    0.000    0.000    0.000 {built-in method builtins.len}
      261    0.000    0.000    0.000    0.000 {method 'endswith' of 'str' objects}
       30    0.000    0.000    0.000    0.000 QCLayers.py:411(<listcomp>)
        3    0.000    0.000    0.000    0.000 numeric.py:156(ones)
        1    0.000    0.000    0.000    0.000 fromnumeric.py:2337(amin)
        1    0.000    0.000    0.000    0.000 numeric.py:433(asarray)
        1    0.000    0.000    0.000    0.000 {method 'read' of '_io.TextIOWrapper' objects}
        1    0.000    0.000    0.000    0.000 fromnumeric.py:1040(argmin)
        4    0.000    0.000    0.000    0.000 numeric.py:504(asanyarray)
        1    0.000    0.000    0.000    0.000 fromnumeric.py:2227(amax)
        8    0.000    0.000    0.000    0.000 _internal.py:247(__init__)
        3    0.000    0.000    0.000    0.000 {built-in method numpy.core.multiarray.copyto}
        1    0.000    0.000    0.000    0.000 {method 'cumsum' of 'numpy.ndarray' objects}
        1    0.000    0.000    0.000    0.000 {method 'argmin' of 'numpy.ndarray' objects}
        5    0.000    0.000    0.000    0.000 {method 'copy' of 'dict' objects}
        1    0.000    0.000    0.000    0.000 _bootlocale.py:23(getpreferredencoding)
        1    0.000    0.000    0.000    0.000 codecs.py:318(decode)
        1    0.000    0.000    0.000    0.000 {method 'reshape' of 'numpy.ndarray' objects}
        1    0.000    0.000    0.000    0.000 codecs.py:308(__init__)
        1    0.000    0.000    0.000    0.000 function_base.py:13(_index_deprecate)
        2    0.000    0.000    0.000    0.000 {method 'match' of '_sre.SRE_Pattern' objects}
        8    0.000    0.000    0.000    0.000 _internal.py:281(get_as_parameter)
        3    0.000    0.000    0.000    0.000 {built-in method builtins.getattr}
        1    0.000    0.000    0.000    0.000 {built-in method numpy.core.multiarray.result_type}
        2    0.000    0.000    0.000    0.000 {built-in method numpy.core.multiarray.normalize_axis_index}
        1    0.000    0.000    0.000    0.000 numeric.py:1927(isscalar)
        1    0.000    0.000    0.000    0.000 {built-in method _locale.nl_langinfo}
        1    0.000    0.000    0.000    0.000 {method 'astype' of 'numpy.ndarray' objects}
        1    0.000    0.000    0.000    0.000 {method 'tolist' of 'numpy.ndarray' objects}
        1    0.000    0.000    0.000    0.000 {built-in method _codecs.utf_8_decode}
        2    0.000    0.000    0.000    0.000 {method 'end' of '_sre.SRE_Match' objects}
        1    0.000    0.000    0.000    0.000 codecs.py:259(__init__)
        5    0.000    0.000    0.000    0.000 {method 'pop' of 'dict' objects}
        1    0.000    0.000    0.000    0.000 {method 'startswith' of 'str' objects}
        1    0.000    0.000    0.000    0.000 {method 'disable' of '_lsprof.Profiler' objects}
        1    0.000    0.000    0.000    0.000 {built-in method _operator.index}
=========  =======  =======  =======  ======= ===============================================================



After optimizing code performance, there are 3613 function calls in 31.973
seconds. The output of cProfile sorted by cumulative time is listed in the
following table.


=========  =======  =======  =======  ======= ===============================================================
   ncalls  tottime  percall  cumtime  percall filename:lineno(function)
=========  =======  =======  =======  ======= ===============================================================
        1    0.000    0.000   31.907   31.907 QCLayers.py:427(FoM)
       35    0.004    0.000   31.907    0.912 QCLayers.py:383(loTransition)
       30   31.899    1.063   31.900    1.063 OneDSchrodinger.py:97(cLOphononScatter)
        2    0.000    0.000   30.845   15.422 QCLayers.py:421(loLifeTime)
        2    0.000    0.000   30.845   15.422 QCLayers.py:424(<listcomp>)
        1    0.000    0.000    0.062    0.062 QCLayers.py:260(solve_whole)
        1    0.059    0.059    0.059    0.059 OneDSchrodinger.py:61(cSimpleSolve1D)
        1    0.003    0.003    0.003    0.003 OneDSchrodinger.py:71(cSimpleFillPsi)
        1    0.001    0.001    0.003    0.003 QCLayers.py:178(populate_x)
       91    0.000    0.000    0.002    0.000 fromnumeric.py:1821(sum)
       93    0.000    0.000    0.001    0.000 fromnumeric.py:64(_wrapreduction)
      121    0.001    0.000    0.001    0.000 {method 'reduce' of 'numpy.ufunc' objects}
       65    0.000    0.000    0.001    0.000 QCLayers.py:169(avghwLO)
      227    0.000    0.000    0.001    0.000 {built-in method builtins.sum}
       68    0.001    0.000    0.001    0.000 ctypeslib.py:196(from_param)
       28    0.001    0.000    0.001    0.000 QCLayers.py:227(<listcomp>)
     1885    0.001    0.000    0.001    0.000 QCLayers.py:173(<genexpr>)
       65    0.000    0.000    0.000    0.000 {built-in method numpy.core.multiarray.array}
        1    0.000    0.000    0.000    0.000 SaveLoad.py:8(qclLoad)
        1    0.000    0.000    0.000    0.000 QCLayers.py:79(__init__)
        1    0.000    0.000    0.000    0.000 QCLayers.py:104(update_strain)
        1    0.000    0.000    0.000    0.000 QCLayers.py:365(dipole)
        1    0.000    0.000    0.000    0.000 QCLayers.py:106(<listcomp>)
        2    0.000    0.000    0.000    0.000 Material.py:171(__init__)
       30    0.000    0.000    0.000    0.000 QCLayers.py:412(<listcomp>)
        9    0.000    0.000    0.000    0.000 Material.py:65(set_temperature)
       68    0.000    0.000    0.000    0.000 _internal.py:247(__init__)
        5    0.000    0.000    0.000    0.000 Material.py:51(__init__)
      163    0.000    0.000    0.000    0.000 {built-in method builtins.isinstance}
        1    0.000    0.000    0.000    0.000 __init__.py:243(load)
       68    0.000    0.000    0.000    0.000 _internal.py:281(get_as_parameter)
       20    0.000    0.000    0.000    0.000 {built-in method numpy.core.multiarray.empty}
        2    0.000    0.000    0.000    0.000 fromnumeric.py:49(_wrapfunc)
        1    0.000    0.000    0.000    0.000 function_base.py:25(linspace)
        1    0.000    0.000    0.000    0.000 __init__.py:271(loads)
        2    0.000    0.000    0.000    0.000 Material.py:198(set_molefrac)
        1    0.000    0.000    0.000    0.000 fromnumeric.py:2092(cumsum)
        2    0.000    0.000    0.000    0.000 Material.py:182(set_temperature)
        2    0.000    0.000    0.000    0.000 function_base.py:1079(diff)
        1    0.000    0.000    0.000    0.000 decoder.py:334(decode)
        1    0.000    0.000    0.000    0.000 {built-in method io.open}
        2    0.000    0.000    0.000    0.000 Material.py:95(set_strain)
        1    0.000    0.000    0.000    0.000 fromnumeric.py:36(_wrapit)
        2    0.000    0.000    0.000    0.000 {built-in method numpy.core.multiarray.arange}
      261    0.000    0.000    0.000    0.000 {method 'endswith' of 'str' objects}
        1    0.000    0.000    0.000    0.000 decoder.py:345(raw_decode)
        3    0.000    0.000    0.000    0.000 {built-in method numpy.core.multiarray.concatenate}
       30    0.000    0.000    0.000    0.000 QCLayers.py:413(<listcomp>)
       93    0.000    0.000    0.000    0.000 {method 'items' of 'dict' objects}
        3    0.000    0.000    0.000    0.000 numeric.py:156(ones)
       68    0.000    0.000    0.000    0.000 {built-in method builtins.len}
        1    0.000    0.000    0.000    0.000 fromnumeric.py:2337(amin)
        1    0.000    0.000    0.000    0.000 {method 'read' of '_io.TextIOWrapper' objects}
        1    0.000    0.000    0.000    0.000 numeric.py:433(asarray)
        1    0.000    0.000    0.000    0.000 fromnumeric.py:1040(argmin)
        1    0.000    0.000    0.000    0.000 {method 'cumsum' of 'numpy.ndarray' objects}
        1    0.000    0.000    0.000    0.000 fromnumeric.py:2227(amax)
        3    0.000    0.000    0.000    0.000 {built-in method numpy.core.multiarray.copyto}
        4    0.000    0.000    0.000    0.000 numeric.py:504(asanyarray)
        1    0.000    0.000    0.000    0.000 {method 'argmin' of 'numpy.ndarray' objects}
        5    0.000    0.000    0.000    0.000 {method 'copy' of 'dict' objects}
        1    0.000    0.000    0.000    0.000 codecs.py:318(decode)
        1    0.000    0.000    0.000    0.000 _bootlocale.py:23(getpreferredencoding)
        1    0.000    0.000    0.000    0.000 codecs.py:308(__init__)
        1    0.000    0.000    0.000    0.000 function_base.py:13(_index_deprecate)
        2    0.000    0.000    0.000    0.000 {method 'match' of '_sre.SRE_Pattern' objects}
        1    0.000    0.000    0.000    0.000 {method 'reshape' of 'numpy.ndarray' objects}
        1    0.000    0.000    0.000    0.000 numeric.py:1927(isscalar)
        1    0.000    0.000    0.000    0.000 {built-in method _locale.nl_langinfo}
        1    0.000    0.000    0.000    0.000 {built-in method numpy.core.multiarray.result_type}
        1    0.000    0.000    0.000    0.000 {method 'tolist' of 'numpy.ndarray' objects}
        1    0.000    0.000    0.000    0.000 {method 'astype' of 'numpy.ndarray' objects}
        1    0.000    0.000    0.000    0.000 {built-in method _codecs.utf_8_decode}
        5    0.000    0.000    0.000    0.000 {method 'pop' of 'dict' objects}
        1    0.000    0.000    0.000    0.000 {method 'startswith' of 'str' objects}
        3    0.000    0.000    0.000    0.000 {built-in method builtins.getattr}
        2    0.000    0.000    0.000    0.000 {built-in method numpy.core.multiarray.normalize_axis_index}
        1    0.000    0.000    0.000    0.000 codecs.py:259(__init__)
        1    0.000    0.000    0.000    0.000 {method 'disable' of '_lsprof.Profiler' objects}
        2    0.000    0.000    0.000    0.000 {method 'end' of '_sre.SRE_Match' objects}
        1    0.000    0.000    0.000    0.000 {built-in method _operator.index}
=========  =======  =======  =======  ======= ===============================================================
