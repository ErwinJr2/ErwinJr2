OneDQuantum C Library and Python Interface
==========================================

In the C part of the software we are going to implement numerical computation
of the following problems. 


Numerically solving 1D Schrodinger's Equation
-----------------------------------------------

The time-independent Schrodinger equation has the form

.. math::
   \left[-\frac{\hbar^2}{2m}\frac{\partial}{\partial x^2} 
   + V(x)\right]\psi_i(x) = E_i\psi_i(x)

The inputs of the Schrodinger equation solver include: a finite 1D array 
with position :math:`x`, the corresponding potential :math:`V(x)` with the same size, the
effective mass :math:`m`, and an eigenstate range specified by the user,
:math:`\left[E_\text{min}, E_\text{max}\right]`. The outputs are the eigenfunction,
:math:`\psi`, and the eigenvalue, :math:`E`.

We will implement two methods of solving the 1D Schrodinger's equation. The
first method combines the Newton's method that searches for eigenvalues :math:`E`
and the Numerov's method that solving for the corresponding eigenfunction
:math:`\psi(x)` given any specific :math:`E`. The second method diagonalizes the
Hamiltonian directly using linear algebra.

Effective mass in band 
----------------------
Band theory predicts that movement of particle in a potential over long
distance can be very different from the movement of the same particle in
vacuum. Usually, the movement is complicated; however, when the electron is
in the highest energies of the valence band or the lowest energies of the
conduction band, it can be shown that electrons behave as free electrons
except with a different mass, the effective mass :math:`m_\text{eff}`.

A particle's effective mass in each band can be approximated by Taylor
expanding the band structure and ignoring higher-than-second-order terms, as
the band structure can be expanded locally as

.. math:: 

   E(k) \approx E_0 + \frac{\hbar^2 k^2}{2 m_\text{eff}},
   \frac{1}{m_\text{eff}} \equiv \frac{2}{\hbar^2}\frac{\partial^2 E}{\partial k^2}

where :math:`k` is the wave vector, and :math:`E_0` is the edge energy of the band. 

For QCL simulations, because of small layer thickness, constant effective
mass approximation is sometimes not enough. This can be corrected by
including non-parabolic effective mass, or effective mass with energy
dependence.  In this package, we will offer constant effective mass as simple
solver and also non-parabolic effective mass computed using 
:math:`k\cdot p` method. 

Self-consistency solver
-----------------------

Electron-electron Coulomb interaction can be a determinant part of electron
motion in semiconductors. To the first order this interaction is included by
adding a Maxwell-Poisson equation to correct the potential and solve the
equations self-consistently. 

.. math::

   V = V_0 + V_c,
   \nabla^2 V_c = \frac{\rho(x)}{\epsilon} = \sum_i 
   \frac{e n_i}{\epsilon} |\psi_i(x)|^2

which means that the potential depends on the 
eigenstates as well as the corresponding occupation number :math:`n_i`.

Electron thermal distributions
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
the Fermi-Dirac distribution (Eq.~(\ref{eq:zero_temp}) and
(\ref{eq:finite_temp}) respectively),
and the high-temperature approximation with the Boltzmann distribution
(Eq.~(\ref{eq:boltzmann})). All distributions will have two methods, giving
constant chemical potential :math:`\mu` distribution and return total number of
particles :math:`\sum n_i`, and given total number of particles :math:`\sum n_i` and
return chemical potential :math:`\mu`.