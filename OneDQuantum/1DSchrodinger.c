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
#include "pocketfft/pocketfft.h"

#define Y_EPS 0.1 /**< Starting point for ode solver. Default value = 0.1 */

#ifdef __cplusplus
extern "C" {
#endif


/**
 * 
 * An ODE solver using Numerov's method for 
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
double numerov(double step, numpyint N, double y0, double y1, 
        double E, const double *V, const double *m, double *y) {
    int n; 
    y[0] = y0;
    y[1] = y1;
    const double unit = 2*m0/sq(hbar)*e0*sq(ANG*step);
    for (n = 1; n < N-1; n++) {
        if(fabs(m[n+1]-m[n])/step < 1E-5*m[n] &&
                fabs(m[n] - m[n-1])/step < 1E-5*m[n-1] ) {
        /* if(0){ */
            /* Numerov's method, step error O(step^6) */
            /* is bad for m is in the middle of derivative TODO*/ 
            y[n+1] = (2 * y[n] * (1.0 - 5.0/12 * ( E - V[n]) * unit * m[n]) 
                    - y[n-1] * (1.0 + 1.0/12 * (E - V[n-1]) * unit * m[n-1])) 
                / (1.0 + 1.0/12 * (E - V[n+1]) * unit * m[n+1]);
        }
        else {
            /* m*(1/m)'*step/2 to O(step^3) */
            /* double mmp = (m[n]/m[n+1] - m[n]/m[n-1])/4;  */
            /* double mmp = -(m[n+1]/m[n] - m[n-1]/m[n])/4;  */
            /* Simple Euler's method, setp error O(step^4) */
            /* y[n+1] = (-(E-V[n])*unit*m[n]*y[n] + 2*y[n]  */
            /*           - (1-mmp)*y[n-1])/(1 + mmp); */
            double mplus = (m[n+1] + m[n])/2;
            double mminus = (m[n] + m[n-1])/2;
            y[n+1] = (-(E-V[n])*unit*y[n] + (1/mplus + 1/mminus)*y[n] 
                      - y[n-1]/mminus)*mplus;
        }
    }
    return y[N-1];
}

/**
 * 
 * An ODE solver using RK4 method for 
 * \f$ -\frac{d}{dx}(\frac{\hbar^2}{2m(x)} \frac{dy(x)}{dx}) + V(x) y(x) = E y(x) \f$
 * with starting \f$ x_0 \f$ and \f$ y_0 \f$, \f$ y_1 \f$, 
 * ends at \f$ x_0 + N \times step \f$.
 * No normalization imposed.
 * See numerov() for arguments explaination 
 *
 */
#ifndef MACOS
inline 
#endif
double rk4(double step, numpyint N, double y0, double y1, 
        double E, const double *V, const double *m, double *y) {
    int n; 
    double yp;
    const double unit = 2*m0/sq(hbar)*e0*sq(ANG*step);
    double k11, k12, k21, k22, k31, k32, k41, k42;
    y[0] = y0;
    yp = (y1 - y0)/m[0];
    for (n = 0; n < N-1; n++) {
        /* RK4 ODE solver */
        double mplus = (m[n] + m[n+1])/2;
        double Vplus = (V[n] + V[n+1])/2;

        k11 = m[n] * yp; 
        k12 = (V[n] - E) * y[n] * unit; 

        k21 = mplus * (yp + k12/2);
        k22 = (Vplus - E) * (y[n] + k11/2) * unit;

        k31 = mplus * (yp + k22/2);
        k32 = (Vplus - E) * (y[n] + k21/2) * unit;

        k41 = m[n+1] * (yp + k32);
        k42 = (V[n+1] - E) * (y[n] + k31) * unit;

        y[n+1] = y[n] + (k11 + 2*(k21+k31) + k41)/6;
        yp += (k12 + 2*(k22+k32) + k42)/6;
        // Euler 
        // y[n+1] = y[n] + k11;
        // yp += k12;
    }
    return y[N-1];
}

#define ode rk4 /**< Use RK4 as the default ODE solver */


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
 * \param[in] *starts
 * \param[in] *ends wavefuntion limited to psi[starts[i]:ends[i]]
 * \param[in] *mat is a pointer to band structure, for updating
 *                effective mass according to energy. When it's NULL it means
 *                using constant mass without non-parabolic effective mass. 
 * \param[out] *psis (output) 
 *                   \f$ \psi + i N \times sizeof(double) \f$ is the 
 *                   wavefunction with energy EigenEs[i].
 */
void FillPsi(double step, numpyint N, const double *EigenEs,
        numpyint EN, const double *V, double *m, double *psis, 
        numpyint *starts, numpyint *ends,
        Band *const mat) {
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
            for(j = 0; j<starts[i]; j++) {
                psi[j] = 0;
            }
            for(j = ends[i]; j<N; j++) {
                psi[j] = 0;
            }
            psi += starts[i];
            int len = (N<ends[i]? N : ends[i]) - starts[i];
            double modsq = 0;
            if (mat != NULL) {
                /*MP assume only mass is updated*/
                UpdateBand(mat, EigenEs[i], V, bandm);
            }
            ode(step, len, 0.0, Y_EPS, EigenEs[i], 
                    V+starts[i], bandm+starts[i], psi);
            /* Normalization */
            for(j=0; j<len; j++) {
                modsq += sq(psi[j]);
            }
            modsq = sqrt(modsq * step);
            for(j=0; j<len; j++) {
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
#define findZero(x1, y1, x2, y2) (((x1)*(y2) - (x2)*(y1))/((y2)-(y1)))


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
                E0 = Es[i];
            }
            else if(yend[i]*yend[i-1] < 0) {
                /* Here secant method is used instead of Newton's 
                 * * because secant method is more stable, and since the 
                 * * extra yend evaluation introduced in discret derivative 
                 * * for Newton's method, secant method is still more efficient
                 * * */
                E1 = Es[i-1];
                E2 = Es[i];
                E0 = findZero(E2, yend[i], E1, yend[i-1]);
#ifndef SIMPLE
                int count=0; 
                if(mat != NULL) { 
                    UpdateBand(mat, E0, V, mband);
                }
                y0 = ode(step, N, 0.0, Y_EPS, E0, V, mband, y);
                y1 = yend[i-1];
                y2 = yend[i];
                /* Filter singular case */
                if(fabs(y0) > fabs(yend[i]) || fabs(y0) > fabs(yend[i-1])){
                    continue;
                }
                while(fabs(y0) > 1e-20 && fabs(E2-E1) > 1e-14 && count < 40){
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
 * Bonded version of Solve1D()
 *
 * \param[in] step step size
 * \param[in] N number of steps: as a limit
 * \param[in] Elower: starting point of ode solver, Elower - x*field
 *                    should be smaller than V, and solve start with
 * \param[in] Eupper: ending point of ode solver, Eupper - x*field
 *                    should be larger than V, and solve start with
 * \param[in] field: field in unit kV/cm, for determing energy shift
 * \param[in] *Es initial search range of eigen energy 
 * \param[in] EN number of eigen energy to find
 * \param[in] *V potential 
 * \param[in] *m effective mass
 * \param[in] *mat is a pointer to band structure
 * \param[out] *EigenE (output) eigen energy
 *
 */
numpyint Solve1DBonded(double step, numpyint N, 
        double Elower, double Eupper, double field,
        const double *Es, numpyint EN, 
        const double *V, double *m, Band * const mat, double *EigenE) {
    double *yend; 
    int NofZeros=0;
    int length;
    /* convert kV/cm to eV/pixal */
    field *= ANG*step*1E5;
    length = (int) ceil((Eupper - Elower)/field);
#ifdef _DEBUG
    if(mat != NULL) {
        assert(N == mat->N);
    }
#endif
    yend = (double *)malloc(EN * sizeof(double));
#ifdef __MP
    #ifdef _DEBUG
    printf("Start a Solve1DPeriod with openMP.\n");
    #endif
    #pragma omp parallel
#endif
    {
        int i;
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
        for(i=0; i < EN; i++) {
            double E = Es[i];
            int start = (int) floor((Elower - E)/field);
            if(start < 0)
                start = 0;
            if(mat != NULL) { 
                UpdateBand(mat, E, V, mband);
            }
            yend[i] = ode(step, start+length<=N ? length : N-start, 
                    0.0, Y_EPS, E, V+start, mband+start, y+start);
        }

#ifdef __MP
        #pragma omp barrier
        #pragma omp for ordered
#endif
        for(i=1; i<EN; i++) {
            double E0, E1, E2;
            double y0, y1, y2;
            if(yend[i] == 0) {
                E0 = Es[i];
            }
            else if(yend[i]*yend[i-1] < 0) {
                /* Here secant method is used instead of Newton's 
                 * * because secant method is more stable, and since the 
                 * * extra yend evaluation introduced in discret derivative 
                 * * for Newton's method, secant method is still more efficient
                 * * */
                E1 = Es[i-1];
                E2 = Es[i];
                E0 = findZero(E2, yend[i], E1, yend[i-1]);
#ifndef SIMPLE
                int count=0; 
                int start = (int) floor((Elower - E0)/field);
                if(start < 0)
                    start = 0;
                if(mat != NULL) { 
                    UpdateBand(mat, E0, V, mband);
                }
                y0 = ode(step, start+length<=N ? length : N-start,
                        0.0, Y_EPS, E0, V+start, mband+start, y+start);
                y1 = yend[i-1];
                y2 = yend[i];
                /* Filter singular case */
                if(fabs(y0) > fabs(yend[i]) || fabs(y0) > fabs(yend[i-1])){
                    continue;
                }
                while(fabs(y0) > 1e-20 && fabs(E2-E1) > 1e-14 && count < 40){
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
                    start = (int) floor((Elower - E0)/field);
                    if(start < 0)
                        start = 0;
                    y0 = ode(step, start+length<=N ? length : N-start,
                        0.0, Y_EPS, E0, V+start, mband+start, y+start);
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


# define MINPSI 1E-5 /**< The min cutoff for integral of wavefunctions */
#ifdef _WINDLL
__declspec(dllexport)
#endif
/**
 * Calculate the LO phonon scattering rate
 *
 * \param[in] step step size
 * \param[in] N number of steps
 * \param[in] kl wavevector of LO phonon
 * \param[in] *psi_i \f$\psi_i\f$ wavefunction i
 * \param[in] *psi_j \f$\psi_j\f$ wavefunction j
 * \return    \f$I_{ij} = \int\mathrm dx\mathrm dy\, \psi_i(x)\psi_j(x)
 *             \exp\left[-k_l|x-y|\right]\psi_i(y)\psi_j(y) \f$
 */
double LOphononScatter(double step, numpyint N, double kl,
        const double *psi_i, const double *psi_j) {
    double Iij = 0;
    int i;
    double powerUnit = -kl*step*ANG;
    double *psiij;
    rfft_plan p;
    int start, end;
    for(start = 0; start < N && 
            (fabs(psi_i[start]) < MINPSI || fabs(psi_j[start]) < MINPSI);
            start++);
    for(end = N-1; end >= start  &&
            (fabs(psi_i[end]) < MINPSI || fabs(psi_j[end]) < MINPSI);
            end--);
    if(start == end)
        return 0.0;
    end += 1;
    /* end - start is at least 1 */ 
    psiij = (double *)malloc(sizeof(double)*(end-start)*2);
    p = make_rfft_plan((end-start)*2);
    for(i=start; i<end; i++) {
        psiij[i-start] = psi_i[i] * psi_j[i];
    }
    for(i=0; i<end-start; i++) {
        psiij[end+i-start] = 0;
    }
    end = (end - start);
    if(rfft_forward(p, psiij, 1.0) != 0) {
        #ifdef _DEBUG
        printf("FFT failed!\n");
        #endif
        return -1.0;
    }
    psiij[0] = sq(psiij[0]);
    for(i=1; i<end; i++) {
        psiij[2*i-1] = sq(psiij[2*i]) + sq(psiij[2*i-1]);
        psiij[2*i] = 0;
    }
    psiij[2*end-1] = sq(psiij[2*end-1]);
    rfft_backward(p, psiij, 0.5/end);
    Iij += psiij[0]/2;
    for(i=1; i<end; i++) {
        Iij += psiij[i]*exp(powerUnit*i);
    }
    destroy_rfft_plan(p);
    free(psiij);
    return 2*Iij * sq(step);
}


#ifdef _WINDLL
__declspec(dllexport)
#endif
/**
 * Calculate sum LO phonon scattering rate from psi_i to all psi_j's
 *
 * \param[in] step step size
 * \param[in] N number of steps
 * \param[in] *kls wavevector of LO phonon between psi_i to psi_j's
 * \param[in] *psi_i \f$\psi_i\f$ wavefunction i
 * \param[in] *psi_js psi_j = psi_js[n*N] \f$\psi_j\f$ wavefunction j
 * \param[in] *factor_js the factor \f$f_j\f$ before \f$I_{ij}\f$ before sum
 * \param[in] Nj number of psi_j
 * \return    \f$\sum_j f_j I_{ij} = 
 *             \sum_j f_j \int\mathrm dx\mathrm dy\, \psi_i(x)\psi_j(x)
 *             \exp\left[-k_l|x-y|\right]\psi_i(y)\psi_j(y) \f$
 */
double LOtotal(double step, numpyint N, const double *kls,
        const double *psi_i, const double *psi_js, const double *fjs,
        int Nj) {
    double Iij = 0;
    int j,n;
    int starti, endi;
    double *psiij;
    rfft_plan p;
    for(starti = 0; starti < N && fabs(psi_i[starti]) < MINPSI; starti++);
    for(endi = N-1; endi >= starti && fabs(psi_i[endi]) < MINPSI; endi--);
    if(starti == endi)
        return 0.0;
    endi = endi - starti + 1;
    psi_i += starti;
    p = make_rfft_plan(endi*2);
#ifdef __MP
#pragma omp parallel
#endif
    {
        psiij = (double *)malloc(sizeof(double) * endi * 2);
        #ifdef __MP
        #pragma omp for reduction(+:Iij)
        #endif
        for(n=0; n<Nj; n++) {
            const double *psi_j = psi_js + N*n + starti;
            for(j=0; j<endi; j++) {
                psiij[j] = psi_i[j] * psi_j[j];
            }
            for(j=0; j<endi; j++) {
                psiij[endi+j] = 0;
            }
            if(rfft_forward(p, psiij, 1.0) != 0) {
                #ifdef _DEBUG
                printf("FFT failed!\n");
                #endif
                return -1.0;
            }
            psiij[0] = sq(psiij[0]);
            for(i=1; i<endi; i++) {
                psiij[2*i-1] = sq(psiij[2*i]) + sq(psiij[2*i-1]);
                psiij[2*i] = 0;
            }
            psiij[2*end-1] = sq(psiij[2*endi-1]);
            rfft_backward(p, psiij, 0.5/endi);
            Iij += psiij[0]/2;
            for(i=1; i<end; i++) {
                Iij += psiij[i]*exp(powerUnit*i);
            }
        }
    }
    destroy_rfft_plan(p);
    free(psiij);
    return 2*Iij * sq(step);
}


#ifdef _WINDLL
__declspec(dllexport)
#endif 
/**
 *
 * \param[in] step step size
 * \param[in] N number of steps
 * \param[in] *EigenEs list of eigen energies
 * \param[in] EN number of eigen energies we consider
 * \param[in] *psis psis[n*N:(n+1)*N] is the wavefunction for EigenEs[n]
 * \param[in] *masses masses[n] is the x-y effective mass for psis[n]
 * \param[in] Eshift the energy shift (period*field) between periods
 *                   should >= max(EigenEs) - min(EigenEs)
 * \param[in] xShift position translation in pixal between periods
 * \param[out] *loMatrix gamma[n,m] == loMatrix[EN*n+m] is the LO transition
 *                       rate between psis[n] and psis[m], 
 *                       EigenEs[n] > EigenEs[m] (mod Eshift)
 */
void LOMatrix(double step, numpyint N, const double *EigenEs, numpyint EN, 
        double *psis, const double *masses, double Eshift, numpyint xShift,
        double *loMatrix){
    /*TODO*/
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

