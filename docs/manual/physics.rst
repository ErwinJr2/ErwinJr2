Physics Model and Formulas
==========================================

The C library solves 1D quantum problem and the related thermal
and electrical problems. A Python interface for the C library is
also included in this package.


Numerically Solving 1D Schrodinger's Equation
---------------------------------------------

The time-independent Schrodinger equation has the form

.. math::
   \left[-\frac{\hbar^2}{2}\frac{\partial}{\partial x} \frac{1}{m} \frac{\partial}{\partial x}
   + V(x)\right]\psi_i(x) = E_i\psi_i(x)

The inputs of the Schrodinger equation solver include: a finite 1D array 
with position :math:`x`, the corresponding potential :math:`V(x)` with the same size, the
effective mass :math:`m`, and an eigenstate range specified by the user,
:math:`\left[E_\text{min}, E_\text{max}\right]`. The outputs are the eigenfunction,
:math:`\psi`, and the eigenvalue, :math:`E`.
The difference from standard form of the mass between spatial derivative is the requirement 
of Hermiticity for spatial dependent mass. 

We solve the 1D Schrodinger's equation numerically. Our
method combines the Newton's method that searches for eigenvalues :math:`E`
and the Numerov's method that solves for the corresponding eigenfunction
:math:`\psi(x)` given any specific :math:`E`. 

Shooting Algorithm for Eigen-problem
------------------------------------

#. Initialize with a range for eigen energies. 
#. For each possible eigen energy, solve for the wavefunction using the
   Numerov's method for second order differential equations (if the mass is constant, 
   otherwise Euler's method), and check whether the solution satisfies boundary condition.
   If so, the energy is an eigen energy.
#. Use secant's method to find eigen energy. Newton's method is not chosen
   because the root finding converges usually ~10 interaction, and scant's
   method (O(n^0.618)) compare to Newton's (O(n^0.5)) doesn't worth the 
   extra wavefunction evaluation for numerical derivative. 
#. One noticeable problem for shooting algorithm is that it can miss state pairs that are 
   almost degenerate. When running the software and seeing missing of state, it is 
   recommended to change the global field slightly and/or rotate the layer design and 
   try again. 

An example of solving simple Schrodinger equation can be found 
:ref:`here<example_schrodinger>`.

Effective Mass in Band 
----------------------

Band theory predicts that movement of a particle in a potential over long
distance can be very different from the movement of the same particle in
vacuum. Usually, the movement is complicated; however, when the electron is
in the highest energies of the valence band or the lowest energies of the
conduction band, it can be shown that electrons behave as free electrons
except with a different mass, the effective mass :math:`m_\text{eff}`.

A particle's effective mass in each band can be approximated by Taylor
expanding the band structure and ignoring higher-than-second-order terms, as
the band structure can be expanded locally as

.. math:: 

   E(k) \approx E_0 + \frac{\hbar^2 k^2}{2 m_\text{eff}}

where :math:`k` is the wave vector, and :math:`E_0` is the edge energy of the band. 

For QCL simulations, because of small layer thickness, constant effective
mass approximation is sometimes not enough. This can be corrected by
including non-parabolic effective mass, or effective mass with energy
dependence.  In this package, we will offer constant effective mass as simple
solver and also non-parabolic effective mass computed using 
:math:`k\cdot p` method. 

Specifically in Zincblende crystal, the effective mass has the form including 
non-parabolic dispersion effect: 

.. math::
    \frac{m_0}{m_\text{eff}} = 1 + 2F 
    + \frac 13 \frac{E_P}{\Delta E_c + E_g + \Delta_{\text{SO}}}
    + \frac 23 \frac{E_P}{\Delta E_c + E_g}

Where :math:`m_0` is electron mass in vacuum, :math:`E_g` is the bandgap 
at :math:`\Gamma` point, :math:`\Delta_{\text{SO}}, E_P, F` are parameters
describing near-:math:`\Gamma` behavior of the conduction band and valence 
band, :math:`\Delta E_c` is the energy of electrons above conduction band, 
or effective kinetic energy. When :math:`\Delta E_c=0` the model reduces to 
standard effective mass model without non-parabolic dispersion. 
See :cite:`vurgaftman2001band`, :cite:`PhysRevB.50.8663`. 

For Wurtzite crystal, :math:`F=0`. 

The computation of effective mass is implemented in
:doc:`../clib/file/band_8c` (also see :doc:`../clib/file/band_8h`).
The code structure is also capable of adding new crystal structures. 
See the material sections for details. 

Self-consistency Solver for Electron Coulomb Potential
------------------------------------------------------

Electron-electron Coulomb interaction can be a determinant part of electron
motion in semiconductors. To first order this interaction is included by
adding a Maxwell-Poisson equation to correct the potential and solve the
equations self-consistently. 

.. math::

   &V = V_0 + V_c\\
   &\nabla^2 V_c = \frac{\rho(x)}{\epsilon} = \sum_i 
   \frac{e n_i}{\epsilon} \mid\psi_i(x)\mid^2

which means that the potential depends on the 
eigenstates as well as the corresponding occupation number :math:`n_i`.

An example comparing the results from solving the simple Schrodinger equation 
and from solving the equation with the electron-electron interaction correction
can be found :ref:`here<example_maxwell>`.

Electron Thermal Distributions
------------------------------

The 1D Schrodinger's equation solver provides the energy bands, which are
useful for calculations of physical properties of the material. Here, we
consider the electron density and the mean energy, predicted by the
Fermi-Dirac statistics, where the occupation frequency for each eigenstate is

.. math:: 
   n_i = \frac{1}{\exp\big[(E_i- \mu)/k_BT\big]+1}.

At zero temperature, Fermi-Dirac statistics becomes

.. math::
   n_i \stackrel{k_BT\to 0}{=} \begin{cases}
   0, & \text{ if } { E_i > \mu, } \\
   1, & \text{ if } { E_i < \mu. }
   \end{cases}

At high temperature, Fermi-Dirac statistics approaches Maxwell-Boltzmann distribution

.. math:: 
    n_i \stackrel{k_BT\gg E-\mu}{=} \exp\left(-\frac{E-\mu}{k_BT}\right).


In this package, we provide the zero- and finite-temperature computation of
the Fermi-Dirac statistics, and the high-temperature approximation with the 
Maxwell-Boltzmann distribution. All distributions will have two methods: given
constant chemical potential :math:`\mu` distribution and return total number of
particles :math:`\sum n_i`, and given total number of particles :math:`\sum n_i` and
return chemical potential :math:`\mu`.

An example of finding the thermal distribution of electrons, 
given eigen energies and wavefunctions,
can be found :ref:`here<example_thermal>`.


.. bibliography:: refs.bib
