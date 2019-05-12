/**
 * \file
 * 
 * \brief Solve 1D Schrodinger equation.
 *
 *
 */


#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#ifdef __MP /*openmp support*/
#include <omp.h>
#endif 
#include "science.h"
#include "band.h"

#define Y_EPS 0.1 /**< Starting point for ode solver. Default value = 0.1 */

#ifdef __cplusplus
extern "C" {
#endif


/**
 * 
 * An ODE solver for 
 * \f$ -\frac{d}{dx}(\frac{\hbar^2}{2m(x)} \frac{dy(x)}{dx}) + V(x) y(x) = E y(x) \f$
 * with starting \f$ x_0 \f$ and \f$ y_0 \f$, \f$ y_1 \f$, 
 * ends at \f$ x_0 + N \times step \f$.
 * No normalization imposed.
 *
 * \param[in] step \f$ \Delta x \f$, stepsize
 * \param[in] N number of steps
 * \param[in] y0 value of y at \f$ x_0 \f$
 * \param[in] y1 value of y at \f$ x_0 + step \f$
 * \param[in] E energy, unit eV
 * \param[in] *V V[n] is the potential at \f$ x = x_0 + n \times step \f$
 * \param[in] *m m[n] is the effective mass at \f$ x = x_0 + n \times step \f$. 
 *               m is in unit m0 (free electron mass)
 * \param[out] *y (output) value of y at \f$ x = x_0 + n \times step \f$.
 */
#ifndef MACOS
inline 
#endif
double ode(double step, numpyint N, double y0, double y1, 
        double E, const double *V, const double *m, double *y) {
    int n; 
    y[0] = y0;
    y[1] = y1;
    const double unit = 2*m0/sq(hbar)*e0*sq(ANG*step);
    for (n = 1; n < N-1; n++) {
        if(fabs(m[n+1]-m[n])/step < 1E-5*m[n] &&
                fabs(m[n] - m[n-1])/step < 1E-5*m[n-1] ) {
            /* Numerov's method, step error O(step^6) */
            /* is bad for m is in the middle of derivative TODO*/ 
            y[n+1] = (2 * y[n] * (1.0 - 5.0/12 * ( E - V[n]) * unit * m[n]) 
                    - y[n-1] * (1.0 + 1.0/12 * (E - V[n-1]) * unit * m[n-1])) 
                / (1.0 + 1.0/12 * (E - V[n+1]) * unit * m[n+1]);
        }
        else {
            double mmp = (m[n]/m[n+1] - m[n]/m[n-1])/4; // m*(1/m)' to O(^3)
            /* Simple Euler's method, setp error O(step^4), 
             * TODO: try RK4 */
            y[n+1] = (-(E-V[n])*unit*m[n]*y[n] +
                    2*y[n] - (1-mmp)*y[n-1])/(1 + mmp);
        }
    }
    return y[N-1];
}


#ifdef _WINDLL
__declspec(dllexport)
#endif
 
/** 
 * Fill in wavefunctions in \f$ \psi \f$'s accroding to eigen energy in EigenEs. 
 * \f$ \psi + i N \times sizeof(double) \f$ is the wavefunction with Energy EigenEs[i] 
 * The result is normalized to 1 (so psi is unit sqrt(Angstrom^-1)
 * 
 * \param[in] step step size
 * \param[in] N number of steps
 * \param[in] *EigenEs list of eigen energies
 * \param[in] EN number of eigen energies we consider
 * \param[in] *V V[n] is the potential at \f$ x = x_0 + n \times step \f$
 * \param[in] *m m[n] is the effective mass at \f$ x = x_0 + n \times step \f$, 
 *                in unit \f$ m_0 \f$ (free electron mass), used only when 
 *                mat=Null
 * \param[in] *mat is a pointer to band structure, for updating
 *                effective mass according to energy. When it's NULL it means
 *                using constant mass without non-parabolic effective mass. 
 * \param[out] *psis (output) 
 *                   \f$ \psi + i N \times sizeof(double) \f$ is the 
 *                   wavefunction with energy EigenEs[i].
 */
void FillPsi(double step, numpyint N, const double *EigenEs,
        numpyint EN, const double *V, double *m, double* psis, 
        Band * const mat) {
    int i; 
#ifdef _DEBUG
    if(mat != NULL) {
        assert(N == mat->N);
    }
#endif
#ifdef __MP
    #ifdef _DEBUG
    printf("Start a FillPsi with openMP.\n");
    #endif
    #pragma omp parallel
#endif
    {
        double *bandm;
        if (mat != NULL) {
            bandm = (double*)malloc(N*sizeof(double));
        }
        else {
            bandm = m;
        }
#ifdef __MP
    #ifdef _DEBUG
            printf("    From thread %d out of %d\n", 
                    omp_get_thread_num(), omp_get_num_threads());
    #endif
        #pragma omp for
#endif
        for(i=0; i<EN; i++) {
            int j;
            double* psi = psis + i*N;
            double modsq = 0;
            if (mat != NULL) {
                /*MP assume only mass is updated*/
                UpdateBand(mat, EigenEs[i], V, bandm);
            }
            ode(step, N, 0.0, Y_EPS, EigenEs[i], V, bandm, psi);
            /* Normalization */
            for(j=0; j<N; j++) {
                modsq += sq(psi[j]);
            }
            modsq = sqrt(modsq * step);
            for(j=0; j<N; j++) {
                psi[j] /= modsq;
            }
        }
        if (mat != NULL) {
            free(bandm);
            bandm = NULL;
        }
    }
    return; 
}

/** 
 * Find zero between x1 and x2 of function y(x), 
 * s.t. \f$ y(x_n) = 0 \f$
 * using linear interpolation, 
 * i.e assuming y1 and y2 are of opposite signs and
 * returning \f$ x_n \f$.
 */
#define findZero(x1, y1, x2, y2) ((x1*y2 - x2*y1)/(y2-y1))


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
 * \param[in] step step size
 * \param[in] N number of steps
 * \param[in] *Es initial search range of eigen energy
 * \param[in] EN number of eigen energy to find
 * \param[in] *V potential 
 * \param[in] *m effective mass
 * \param[in] *mat is a pointer to band structure
 * \param[out] *EigenE (output) eigen energy
 *
 */
numpyint Solve1D(double step, numpyint N, 
        const double *Es, numpyint EN, 
        const double *V, double *m, Band * const mat,
        double *EigenE) {
    double *yend; 
    int NofZeros=0;
    int i;
#ifdef _DEBUG
    if(mat != NULL) {
        assert(N == mat->N);
    }
#endif
    yend = (double *)malloc(EN * sizeof(double));
#ifdef __MP
    #ifdef _DEBUG
    printf("Start a Solve1D with openMP.\n");
    #endif
    #pragma omp parallel
#endif
    {
        double *y = (double *)malloc(N * sizeof(double));
        double *mband;
        if(mat != NULL) { 
            mband = (double *)malloc(N * sizeof(double));
        }
        else {
            mband = m;
        }
#ifdef __MP
    #ifdef _DEBUG
        printf("    From thread %d out of %d\n", 
                omp_get_thread_num(), omp_get_num_threads());
    #endif
        #pragma omp for
#endif
        for(i=0; i<EN; i++) {
            if(mat != NULL) { 
                UpdateBand(mat, Es[i], V, mband);
            }
            yend[i] = ode(step, N, 0.0, Y_EPS, Es[i], V, mband, y);
        }

#ifdef __MP
        #pragma omp barrier
        #pragma omp for ordered
#endif
        for(i=1; i<EN; i++) {
            double E0, E1, E2;
            double y0, y1, y2;
            if(yend[i] == 0) {
                E0 = Es[i-1];
            }
            else if(yend[i]*yend[i-1] < 0) {
                /* Here secant method is used instead of Newton's 
                 * * because secant method is more stable, and since the 
                 * * extra yend evaluation introduced in discret derivative 
                 * * for Newton's method, secant method is still more efficient
                 * * */
                E0 = findZero(Es[i], yend[i], Es[i-1], yend[i-1]);
#ifndef SIMPLE
                int count=0; 
                if(mat != NULL) { 
                    UpdateBand(mat, E0, V, mband);
                }
                y0 = ode(step, N, 0.0, Y_EPS, E0, V, mband, y);
                y1 = yend[i-1];
                E1 = Es[i-1];
                y2 = yend[i];
                E2 = Es[i];
                /* Filter singular case */
                if(fabs(y0) > fabs(yend[i]) || fabs(y0) > fabs(yend[i-1])){
                    continue;
                }
                while(count < 40 && fabs(y0) > 1e-14 
                        && fmin(fabs(E2-E0), fabs(E1-E0)) > 1e-14){
    #ifdef _DEBUG
                    printf("    Iter No. %d, E0=%.8f, E1=%.8f, E2=%.8f, "
                            "Delta=%g\n", 
                            count, E0, E1, E2, fabs(E2-E1));
    #endif
                    if(y0 * y1 < 0) {
                        y2 = y0;
                        E2 = E0; 
                    }
                    else {
                        y1 = y0;
                        E1 = E0;
                    }
                    E0 = findZero(E1, y1, E2, y2);
                    if(mat != NULL) { 
                        UpdateBand(mat, E0, V, mband);
                    }
                    y0 = ode(step, N, 0.0, Y_EPS, E0, V, mband, y);
                    count++;
                }
    #ifdef _DEBUG
                printf("After %d times secant iteration, E=%.8f Err=%e\n", 
                        count, E0, fabs(y0));
    #endif
#endif
                }
                else {
                    continue;
                }
#ifdef __MP
                #pragma omp ordered
#endif
                EigenE[NofZeros++] = E0;
        }
        free(y);
        y = NULL;
        if (mat != NULL) {
            free(mband);
            mband = NULL;
        }
    }

    free(yend);
    yend = NULL;
    return NofZeros;
}


#ifdef _WINDLL
__declspec(dllexport)
#endif
/**
 * Calculate the LO phonon scattering rate
 */
double LOphononScatter(double step, numpyint N, double kl,
        const double *psi_i, const double *psi_j) {
    double Iij = 0;
    int i;
#ifdef __MP
#pragma omp parallel for reduction(+:Iij)
#endif
    for(i=0; i<N; i++) {
        int j;
        for(j=0; j<N; j++) {
            Iij += psi_i[i] * psi_j[i] * exp(-kl*abs(i-j)*step*ANG) * 
                psi_i[j] * psi_j[j];
        }
    }
    return Iij * sq(step);
}


#ifdef _WINDLL
__declspec(dllexport)
#endif 
/**
 * Checkpoint for python-C interface. Output 137.
 */
numpyint invAlpha()
{return 137;}


#ifdef __cplusplus
}
#endif

