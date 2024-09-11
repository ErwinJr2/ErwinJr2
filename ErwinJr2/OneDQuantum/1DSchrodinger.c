/**
 * @file
 *
 * @brief Solve 1D Schrodinger equation.
 *
 *
 */

#include "1DSchrodinger.h"

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
#include "fftautocorr/fftautocorr.h"

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
 * @param[in] step \f$ \Delta x \f$, stepsize
 * @param[in] N number of steps
 * @param[in] y0 value of y at \f$ x_0 \f$
 * @param[in] y1 value of y at \f$ x_0 + step \f$
 * @param[in] E energy, unit eV
 * @param[in] V V[n] is the potential at \f$ x = x_0 + n \times step \f$
 * @param[in] m m[n] is the effective mass at \f$ x = x_0 + n \times step \f$.
 *               m is in unit m0 (free electron mass)
 * @param[out] y value of y at \f$ x = x_0 + n \times step \f$.
 *
 * @return psiend the last element of y
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
            /* is bad for m is in the middle of derivative */
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

void FillPsi(double step, numpyint N, const double *EigenEs,
             numpyint EN, const double *V, double *m, double *psis,
             numpyint *starts, numpyint *ends,
             const Band *mat)
{
    int i;
#ifdef _DEBUG
    if (mat != NULL)
    {
        assert N == mat->N;
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
        if (mat != NULL)
        {
            bandm = (double *)malloc(N * sizeof(double));
        }
        else
        {
            bandm = m;
        }
#ifdef __MP
#ifdef _DEBUG
        printf("    From thread %d out of %d\n",
               omp_get_thread_num(), omp_get_num_threads());
#endif
#pragma omp for
#endif
        for (i = 0; i < EN; i++)
        {
            int j;
            double *psi = psis + i * N;
            for (j = 0; j < starts[i]; j++)
            {
                psi[j] = 0;
            }
            for (j = ends[i]; j < N; j++)
            {
                psi[j] = 0;
            }
            psi += starts[i];
            int len = (N < ends[i] ? N : ends[i]) - starts[i];
            double modsq = 0;
            if (mat != NULL)
            {
                /*MP assume only mass is updated*/
                BandUpdateM(mat, EigenEs[i], V, bandm);
            }
            ode(step, len, 0.0, Y_EPS, EigenEs[i],
                V + starts[i], bandm + starts[i], psi);
            /* Normalization */
            if (mat != NULL)
            {
                BandNormalize(mat, EigenEs[i], V, psi, step);
            }
            else
            {
                /* default normalization */
                for (j = 0; j < len; j++)
                {
                    modsq += sq(psi[j]);
                }
                modsq = sqrt(modsq * step);
                for (j = 0; j < len; j++)
                {
                    psi[j] /= modsq;
                }
            }
        }
        if (mat != NULL)
        {
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
#define findZero(x1, y1, x2, y2) (((x1) * (y2) - (x2) * (y1)) / ((y2) - (y1)))

    numpyint Solve1D(double step, numpyint N,
                     const double *Es, numpyint EN,
                     const double *V, double *m, const Band *mat,
                     double *EigenE)
    {
        double *yend;
        int NofZeros = 0;
        int i;
#ifdef _DEBUG
    if(mat != NULL) {
        assert N == mat->N;
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
                BandUpdateM(mat, Es[i], V, mband);
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
                    BandUpdateM(mat, E0, V, mband);
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
                        BandUpdateM(mat, E0, V, mband);
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

#define MINPSI 1E-8 /**< The min cutoff for integral of wavefunctions */
double LOphononScatter(double step, numpyint N, double kl,
        const double *psi_ij) {
    double Iij = 0;
    int i;
    double powerUnit = -kl*step*ANG;
    double *psiij;
    int start, end;
    for(start = 0; start < N && fabs(psi_ij[start]) < MINPSI; start++);
    for(end = N-1; end >= start  && fabs(psi_ij[end]) < MINPSI; end--);
    #ifdef _DEBUG
    printf("LO phonon scattering with start at %d, end at %d, length %d.\n",
           start, end, N);
    #endif
    if(start >= end)
        return 0.0;
    end += 1;
    /* end - start is at least 1 */
    psiij = (double *)malloc(sizeof(double)*(end-start));
    for(i=start; i<end; i++) {
        psiij[i-start] = psi_ij[i];
    }
    end = end - start;
    if(autocorr(psiij, end) != 0) {
        #ifdef _DEBUG
        printf("FFT failed!\n");
        #endif
        free(psiij);
        return -1.0;
    }
    Iij += psiij[0] / 2;
    for(i=1; i<end; i++) {
        Iij += psiij[i]*exp(powerUnit*i);
    }
    free(psiij);
    return 2*Iij * sq(step);
}

double LOtotal(double step, numpyint N, const double *kls,
        const double *psi_ijs, const double *fjs, numpyint Nj) {
    double Iij = 0;
    autocorr_plan plan;
    plan = make_autocorr_plan(N);
#ifdef __MP
    int failed = 0;
#pragma omp parallel
#endif
    {
        double *psiij = (double *)malloc(sizeof(double) * mem_len(plan));
        double *mempool = (double *)malloc(sizeof(double) * mem_len(plan));
        int i,n;
        #ifdef __MP
        #pragma omp for reduction(+:Iij)
        #endif
        for(n = 0; n < Nj; n++) {
            #ifdef __MP
            if(failed)
                continue;
            #endif
            const double *psi_ij = psi_ijs + N*n;
            const double powerUnit = -kls[n]*step*ANG;
            for(i = 0; i < N; i++) {
                psiij[i] = psi_ij[i];
            }
            for(; i < mem_len(plan); i++) {
                psiij[i] = 0.0;
            }
            if(autocorr_mem(plan, psiij, mempool) != 0) {
                #ifdef _DEBUG
                printf("FFT failed!\n");
                #endif
                #ifdef __MP
                failed = 1;
                #ifdef _WINDLL
                break;
                #endif
                #else
                free(psiij);
                free(mempool);
                destroy_autocorr_plan(plan);
                return -1.0;
                #endif
            }
            Iij += fjs[n] * psiij[0]/2;
            for(i = 1; i < N; i++) {
                Iij += fjs[n] * psiij[i] * exp(powerUnit*i);
            }
        }
        free(psiij);
        free(mempool);
    }
    destroy_autocorr_plan(plan);
#ifdef __MP
    if (failed)
        return -1;
#endif
    return 2*Iij * sq(step);
}

numpyint invAlpha() {
    return 137;
}

int isMP() {
# ifdef __MP
    return 1;
#else
    return 0;
# endif
}


#ifdef __cplusplus
}
#endif
