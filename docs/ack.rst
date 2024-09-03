Acknowledgment
===================

The current version is written by Ming Lyu from Princeton University,
`Gmachl Group <http://ee.princeton.edu/people/faculty/claire-gmachl/>`_,
as part of his PhD thesis work.

We would also like to give our appreciation for the contribution from
Kevin Oresick (University of Wisconsin).


Version 1.0 and earlier
-----------------------

The author contributions:

- Ming Lyu (CareF): structure design, C prototype, GUI and general integration
- Yaofeng (Desmond) Zhong: Python module, documentation and profiling
- Xiaowen Chen: C lib testing, physics modeling review/documentation and advanced feature in C code

Sara Kacmoli, Yasin Kaya, Hao Deng also provide ideas about improving the software.



ErwinJr before ErwinJr2
-----------------------

The early versions of the software started as a quantum well eigen-solver using
the shooting algorithm, developed with QBASIC at Bell Labs in the early 1990s.
Ref.~\cite{PhysRevB.50.8663} is one of the early works with the solver.
It was then translated into C with some modern modifications
at Princeton in the 2000s by Dr. Daniel Wasserman,
when the name `ErwinJr` was introduced.
Later it was improved by Dr. Kale Franz
who implemented an early version of the GUI software with Python2 and
`Qwt` as the plotting library, as well as a MATLAB-based version.
It was then collaboratively improved by Dr. Yu Song.
