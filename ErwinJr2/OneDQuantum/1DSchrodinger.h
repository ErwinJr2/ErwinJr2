#include "science.h"
#include "band.h"

#ifdef _WINDLL
__declspec(dllexport)
#endif
/**
 * Fill in wavefunctions in \f$ \psi \f$'s accroding to eigen energy in EigenEs.
 * \f$ \psi + i N \times sizeof(double) \f$ is the wavefunction with Energy EigenEs[i]
 * The result is normalized to 1 (so psi is unit sqrt(Angstrom^-1))
 * if mat == NULL  or the normalization defined in mat.
 *
 * @param[in] step step size
 * @param[in] N number of steps
 * @param[in] EigenEs list of eigen energies
 * @param[in] EN number of eigen energies we consider
 * @param[in] V V[n] is the potential at \f$ x = x_0 + n \times step \f$
 * @param[in] m m[n] is the effective mass at \f$ x = x_0 + n \times step \f$,
 *                in unit \f$ m_0 \f$ (free electron mass), used only when
 *                mat=Null
 * @param[in] starts
 * @param[in] ends wavefuntion limited to psi[starts[i]:ends[i]]
 * @param[in] mat is a pointer to band structure, for updating
 *                effective mass according to energy and perform normalization.
 *                When it's NULL it means using constant mass with parabolic
 *                kinetic energy.
 * @param[out] psis \f$ \psi + i N \times sizeof(double) \f$ is the
 *                  wavefunction with energy EigenEs[i].
 */
void FillPsi(double step, numpyint N, const double *EigenEs,
        numpyint EN, const double *V, double *m, double *psis,
        numpyint *starts, numpyint *ends,
        const Band *mat);

#ifdef _WINDLL
__declspec(dllexport)
#endif
/**
 * Solve 1D Schrodinger's equation with potential \f$ V \f$ and effective
 * mass \f$ m \f$
 * in the region \f$ x_0 \leq x < x_0 + step \times N \f$.
 * Boundary condition: wavefunction is zero at boundaries.
 *
 * Method: First scan in energy Es[0:EN] and look for zeros(EigenE) by either
 * simple linear interpolation if SIMPLE is defined;
 * or calculate zero using secant method if SIMPLE is not defined.
 *
 * Es should be in small to large order.
 *
 * @param[in] step step size
 * @param[in] N number of steps
 * @param[in] Es initial search range of eigen energy
 * @param[in] EN number of eigen energy to find
 * @param[in] V potential
 * @param[in] m effective mass
 * @param[in] mat is a pointer to band structure
 * @param[out] EigenE eigen energy
 *
 * @return total number of eigen states found.
 */
numpyint Solve1D(double step, numpyint N,
        const double *Es, numpyint EN,
        const double *V, double *m, const Band *mat,
        double *EigenE);

#ifdef _WINDLL
__declspec(dllexport)
#endif
/**
 * Calculate the LO phonon scattering rate
 *
 * @param[in] step step size in unit Angstrom
 * @param[in] N number of steps
 * @param[in] kl wavevector of LO phonon in unit m^-1. This is DIFFERENT than
 *            the step unit for convience to use kl.
 * @param[in] psi_ij \f$\psi_i \psi_j\f$ wavefunction overlap
 *
 * @return    \f$I_{ij} = \int\mathrm dx\mathrm dy\, \psi_i(x)\psi_j(x)
 *             \exp\left[-k_l|x-y|\right]\psi_i(y)\psi_j(y) \f$
 */
double LOphononScatter(double step, numpyint N, double kl,
                       const double *psi_ij);

#ifdef _WINDLL
__declspec(dllexport)
#endif
/**
 * Calculate sum LO phonon scattering rate from psi_i to all psi_j's
 *
 * @param[in] step step size in unit Angstrom
 * @param[in] N number of steps
 * @param[in] kls wavevector of LO phonon between psi_i to psi_j's in unit m^-1
 * @param[in] psi_ijs psi_j = psi_js[n*N] \f$\psi_i\psi_j\f$ overlap between
 *            \f$\psi_i\f$ and \f$\psi_j\f$ wavefunction
 * @param[in] fjs the factor \f$f_j\f$ before \f$I_{ij}\f$ before sum
 * @param[in] Nj number of psi_j
 *
 * @return    \f$\sum_j f_j I_{ij} =
 *             \sum_j f_j \int\mathrm dx\mathrm dy\, \psi_i(x)\psi_j(x)
 *             \exp\left[-k_l|x-y|\right]\psi_i(y)\psi_j(y) \f$
 */
double LOtotal(double step, numpyint N, const double *kls,
               const double *psi_ijs, const double *fjs, numpyint Nj);

#ifdef _WINDLL
__declspec(dllexport)
#endif
/**
 * Checkpoint for python-C interface. Output 137.
 */
numpyint invAlpha();

#ifdef _WINDLL
__declspec(dllexport)
#endif
/**
 * Check if openMP is loaded.
 */
int isMP();
