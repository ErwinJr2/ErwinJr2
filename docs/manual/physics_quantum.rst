Quantum Tab
==========================================

Eigenstates in the QC structure
---------------------------------------------

The time-independent Schrodinger equation has the form

.. math::
   \left[T(k) + V(x)\right]\psi_i(x) = E_i\psi_i(x)

where :math:`T(k)` the kinetic energy term is often written as
:math:`-\frac{\hbar^2}{2}\frac{\partial}{\partial x} \frac{1}{m} \frac{\partial}{\partial x}`,
a quadratic kinetic energy. But for thin quantum wells of mid-IR resonance,
it has been shown not accurate enough :cite:`PhysRevB.50.8663`.

The kinetic energy used in ErwinJr2 is the 3 band model, where couplings from
the conduction band, the light hole (LH) band and the split off (SO) band are
included to get a more accurate kinetic energy (dispersion relation).

.. math::
   T(k) =
   \begin{pmatrix}
      E_c + \hbar^2k^2/2m^*_c & \mathrm i\sqrt{\frac 23} kP & -\mathrm i\sqrt{\frac 13} kP\\
      -\mathrm i\sqrt{\frac 23}Pk & E_{\mathrm{LH}} + \hbar^2k^2/2m^*_\mathrm{LH} \\
      \mathrm i\sqrt{\frac 13} Pk & & E_{\mathrm{so}} + \hbar^2k^2/2m^*_\mathrm{SO}
   \end{pmatrix}

Based on the 3 band model :cite:`PhysRevB.50.8663` introduced an energy-dependent
effective mass method to reduce the problem back to a single band problem
(for Zincblende crystals):

.. math::
    \frac{m_0}{m^*} = 1 + 2F
    + \frac 13 \frac{E_P}{\Delta E_c + E_g + \Delta_{\text{SO}}}
    + \frac 23 \frac{E_P}{\Delta E_c + E_g}

where :math:`m_0` is electron mass in vacuum, :math:`E_g` is the bandgap
at :math:`\Gamma` point, :math:`\Delta_{\text{SO}}, E_P, F` are parameters
describing near-:math:`\Gamma` behavior of the conduction band and valence
band, specifically :math:`E_P = \frac{2m_0}{\hbar^2}P^2`,
:math:`\Delta E_c` is the energy of electrons above conduction band,
or effective kinetic energy. When :math:`\Delta E_c=0` the model reduces to
standard effective mass model without non-parabolic dispersion.

This method is used for the ODE-based (see the following shooting algorithm)
solver, while the original 3 band matrix is used for the matrix-based solver.

It is worth pointing out that both method used in ErwinJr2 assumes the interface
condition to be continuous conduction band wavefunction :math:`\phi_c` and
its first derivative times inverse mass :math:`\frac{1}{m} \frac{\partial}{\partial x}\phi_c`.

The band parameters we used in ErwinJr2 mostly comes from :cite:`vurgaftman2001band`.


Shooting Algorithm for the Eigen-problem
-----------------------------------------

The inputs of the Schrodinger equation solver include: a finite 1D array
with position :math:`x`, the corresponding potential :math:`V(x)` with the same size, the
effective mass :math:`m`, and an eigenstate range specified by the user,
:math:`\left[E_\text{min}, E_\text{max}\right]`. The outputs are the eigenfunction,
:math:`\psi`, and the eigenvalue, :math:`E`.
The difference from standard form of the mass between spatial derivative is the requirement
of Hermiticity for spatial dependent mass.

#. Initialize with a range for eigen energies.
#. For each possible eigen energy, solve for the wavefunction using the
   RK4 method for second order differential equations (if the mass is constant,
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

The computation of effective mass is implemented in
:doc:`../clib/file/band_8c` (also see :doc:`../clib/file/band_8h`).
The code structure is also capable of adding new crystal structures.
See the material sections for details.


Matrix Algorithm for the Eigen-problem
--------------------------------------

By properly discretize the :math:`T(k)` operator (:math:`k = \partial/\partial z`),
the Hamiltonian becomes a sparse matrix. With :meth:`scipy.sparse.linalg.eigs`
the eigenstates can be solved efficiently.


Scattering mechanism: LO phonon
--------------------------------

The dominant scattering mechanism for inter-subband transition is Longitudinal
optical phonon transition :cite:`PhysRevB.40.1074`.
The scattering rate between state :math:`\psi_u` and :math:`\psi_l` is:

.. math::
    &\frac{1}{\tau_{ul}} =
    \frac{m_l^* e^2 \omega_{\text{LO}}}{8\pi\hbar^2\epsilon_\rho}
    \int_0^{2\pi} I_{ul}(Q_\theta) \mathrm{d}\theta\\
    &I_{ul}(Q_\theta) = \frac{1}{Q_\theta}\iint \mathrm{d}z\mathrm{d}z' \psi_u(z)\psi_l(z)
    \mathrm{e}^{-Q_\theta\mid z-z'\mid}\psi_u(z')\psi_l(z') \\
    &Q_\theta = \sqrt{k_u^2 + k_l^2 - 2k_u k_l \cos\theta} \\
    &\frac{\hbar^2k_u^2}{2m_u^*} = \frac{\hbar^2k_l^2}{2m_l^*}
    + E_u - E_l - \hbar\omega_{\text{LO}} \\
    &\epsilon_\rho^{-1} = \epsilon_\infty^{-1} - \epsilon_{\text{static}}^{-1}

where :math:`k_u` and :math:`k_l` are upper and lower state electron momentum
in the epitaxy layer plain, and :math:`Q_\theta` is the phonon momentum.
With the assumption that :math:`k_u = 0`, the formula reduces to:

.. math::
    \frac{1}{\tau_{ij}} = \frac{m_l^* e^2 \omega_{\text{LO}}}
    {4\hbar^2 \epsilon_\rho} I_{ij}(k_l)



Scattering mechanism: interface roughness (IFR)
-----------------------------------------------

The IFR is described by the standard deviation :math:`\Delta_n`, the correlation
length :math:`\Lambda_n` and the potential change :math:`\delta U_n` at the
interface :math:`n`, whose scattering rate at zero temperature is (for
:math:`E_i \ge E_j`, otherwise it is 0):

.. math::
    \frac{1}{\tau_{ij}^\mathrm{IFR}} =
    \frac{\pi m^*_j}{\hbar^3} \sum_n \Delta_n^2\Lambda_n^2\delta U_n^2
    \left|\psi_j^*(z_n)\psi_i(z_n)\right|^2
    \mathrm e^{- \Lambda^2 m_j^* (E_i - E_j))/2\hbar^2}


The IFR scattering is critical for electron transport in the injectors, and
currently the only mechanism included in ErwinJr2 for intersubband scattering
between states with energy close to each other. Without IFR scattering in the
model, the electron population calculation and therefore the full gain spectrum
from the software may not be very physical.

The GUI supports constant IFR (:math:`\Delta`, :math:`\Lambda` are independent
of :math:`n`) and material dependent IFR (:math:`\Delta`, :math:`\Lambda` are
determined by the material of the layer BEFORE the interface, meaning for the
interface of layer-n and layer-n+1, the material of layer-n decides the IFR
parameters).

For CLI, users can define IFR parameters for each individual interfaces.
See the :doc:`cli` for details.

.. _quantum_gain:

Optical gain and threshold current
------------------------------------

Using Maxwell-Bloch equation the optical gain from intersubband transition is

.. math::
    &g_{ul} = \frac{N_ee^2|d_{ul}|^2\omega(\rho_{uu} - \rho_{ll})}{\hbar \varepsilon_0 nc}
    \pi\mathcal L(\omega - \omega_{ul})\\
    &\mathcal L(\omega - \omega_{ul}) \equiv \frac 1\pi\frac{\gamma_{ul}}
    {\gamma_{ul}^2 + (\omega - \omega_{ul})^2}

where :math:`N_e` is the electron sheet density of a single period,
:math:`\rho_{uu}` and :math:`\rho_{ll}` is the electron percentage population
density of the upper state and lower state, respectively.

This result can be approximated to two level systems by introducing
:math:`\eta_{\text{inf}}` the injection efficiency:

.. math::
    &g = -2\alpha = \frac{\pi\omega \eta_{\text{inj}} e J}
    {\hbar c\epsilon_0 nL_p}
    \,\text{FoM}\,\mathcal L(\omega) \\
    &\text{FoM}\equiv |d_{ul}|^2\tau_u\left(1-\frac{\tau_l}{\tau_{ul}}\right)\\
    &\mathcal L(\omega) \equiv \frac 1\pi\frac{\gamma_\parallel}
    {\gamma_\parallel^2 + (\omega - \omega_{ul})^2}

:math:`\eta_{\text{inf}}` depends on transitions between all other states but
is assumed to be approximately constant,
The Figure of Merit (FoM) is used to characterize the performance of a
structure. :math:`J` is the current density into the device, and with
information of the loss of the optical cavity we can estimate a threshold
current, assuming an reasonable :math:`\eta` or just put it 1.
This estimation is much underestimated.

To couple the design of quantum wells and waveguide, we define the gain
coefficient as the ratio of gain the current density, and also assume
:math:`\eta_{\text{inj}} = 1`, :math:`\omega = \omega_{ul}`,
:math:`\gamma_\parallel = 0.1\omega`.

.. math::
    g_j = \frac
        {\omega e\, \text{FoM}}
        {\gamma_\parallel\hbar c \varepsilon_0 n L_p}


Without the two level approximation, the optical gain can also be calculated
by summing the two level gain/loss for all state pairs
from the steady state solution of the rate equation, which is what happens when
calculating the electron population (the ``Pupulation`` button) and the gain
spectrum (the ``Gain spec`` button).


References
-----------

.. bibliography:: quantum_refs.bib
