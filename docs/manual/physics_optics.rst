Physics Model and Formulas: Optic tab
==========================================

Transfer Matrix Methods for Optical stratum
--------------------------------------------

The confined optical mode in a 1D waveguide can be analytically solved 
using transfer matrix methods, cand calculating the
root of the modal dispersion function :math:`\chi_M` 

.. math::
   	&\begin{pmatrix}
		E_z(0) \\
		H_x(0)
	\end{pmatrix} = 
	\begin{pmatrix}
		\cos\alpha L & -i\gamma\sin\alpha L \\
		-i\gamma^{-1}\sin\alpha L & \cos\alpha L
	\end{pmatrix}
	\begin{pmatrix}
		E_z (L) \\
		H_x (L)
	\end{pmatrix} \\
	&\alpha \equiv \sqrt{n^2k^2-\beta^2}
	\qquad
	\gamma\equiv \frac{\alpha}{kn^2} \sqrt{\frac{\mu_0}{\epsilon_0}}\\
    &M_n \equiv
	\begin{pmatrix}
		\cos\alpha_n L_n & -i\gamma_n\sin\alpha_n L_n \\
		-i\gamma_n^{-1}\sin\alpha_n L_n & \cos\alpha_n L_n
	\end{pmatrix} \qquad M \equiv M_1M_2\cdots M_N \\
    &\chi_M \equiv M_{11}\gamma_s + M_{12} + M_{21}\gamma_s\gamma_0 + 
    M_{22}\gamma_0 

:cite:`Chilwell84`. 

Although generally it depends on effective refractive index as a complex 
function, it's analytically within the domain of interest, therefore Newton's 
method is applicable. 

The roots of :math:`\chi_M(\beta)=0` is the effective refractive index of 
guided modes. In the software we only capture the the mode with largest 
real part of :math:`\beta`, i.e. the foundmental mode. The imaginary part 
of :math:`\beta` yeilds the waveguide loss 

.. math::
	\alpha_w = \frac{4\pi}{\lambda}\mathrm{imag}\beta

Confinement factor
------------------------------

There are multiple different definitions of confinement factor in different 
literatures, for example :cite:`modeling` and :cite:`yariv2006photonics`. 
Here we use the following (the optcial version of :cite:`PhysRevC.6.114`)

.. math::
    \Gamma = \frac{\beta\int_{\text{AR}}n_z E_z^2\mathrm d y}
    {\int n_z^2E_z^2\mathrm d z}


Mirror loss
-------------

The mirror loss in the software can be chosen from cleaved surface (
refractivity :math:`R = |(\beta - 1)/(\beta + 1)|^2`), perfect high-refraction
coating (perfect HR, :math:`R = 1`), perfect anti-reflection coating 
(perfect AR, :math:`R=10^{-1}`) and customized refractivity. 

The mirror refractivity leads to a effective mirror lose per waveguide length
:math:`\alpha_M = -\ln(R_1 R_2)/2L` where :math:`L` is the waveguide length. 


Threshold current
-----------------

The threshold current is calculated assuming :math:`\eta_{\text{inj}} = 1`. 

.. math::
	J_{th} = \frac{0.1\hbar c\epsilon_0 enL_p\alpha}
	{\eta_{\text{inj}} \,\text{FoM}}

See :ref:`quantum_gain` for detail. 


Effective medium theory for QC layers
--------------------------------------

The QC layers have superlattice with width much less than wavelength, which is
when efficiency medium theory works best. The effective refractivity index of
the active region is not isotropic, for TM mode that we are interested in, 

.. math::
    n = \left\langle\frac{1}{n^2}\right\rangle^{-1/2} 

.. bibliography:: optic_refs.bib
