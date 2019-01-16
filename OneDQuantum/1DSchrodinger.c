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

/** Newton's method stopping limit. Default value = \f$ 10^{-5} \f$ */
#define NEWTON_S 1E-5 
/** 
 * NSTEP (unit eV) is the step size to numerically calculate 
 * the derivative for Newton's method. Default value = \f$ 10^{-10} \f$.
 *
 * Optimum step size = \f$ (6*\epsilon/M^3)^{1/3} \f$,  
 * where \f$ \epsilon \f$ is numerical error of yend(E) and \f$ M^3 \f$ is maximum yend(E). 
 * Because of the exponantial behavior when \f$ V>E \f$, yend is very sensitive 
 * to \f$ E \f$ near EigenE. \f$ M^3 \f$ is very large. 
 * Still, it's recommanded to make high \f$ V \f$ side the starting point.
 */
#define NSTEP 1E-10
#define Y_EPS 0.1 /**< Starting point for ode solver. Default value = 0.1 */

#ifdef __cplusplus
extern "C" {
#endif


#ifdef _WINDLL
__declspec(dllexport)
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
inline double ode(double step, numpyint N, double y0, double y1, 
		double E, const double *V, const double *m, double *y) {
	int n; 
	y[0] = y0;
	y[1] = y1;
	const double unit = 2*m0/sq(hbar)*e0*sq(ANG*step);
	for (n = 1; n < N-1; n++) {
		/* Numerov's method, step error O(step^6) */
		/* is bad for m is in the middle of derivative TODO*/
		y[n+1] = (2 * y[n] * (1.0 - 5.0/12 * ( E - V[n]) * unit * m[n]) 
				- y[n-1] * (1.0 + 1.0/12 * (E - V[n-1]) * unit * m[n-1])) 
			/ (1.0 + 1.0/12 * (E - V[n+1]) * unit * m[n+1]);
		/* double mmp = (m[n]/m[n+1] - m[n]/m[n-1])/4; [> m*(1/m)' to O(^3)<] */
		/* Simple Euler's method, setp error O(step^4) */
		/* y[n+1] = (-(E-V[n])*unit*m[n]*y[n] +
		         2*y[n] - (1-mmp)*y[n-1])/(1 + mmp); */
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
 *                in unit \f$ m_0 \f$ (free electron mass)
 * \param[out] *psis (output) 
 *                   \f$ \psi + i N \times sizeof(double) \f$ is the 
 *                   wavefunction with energy EigenEs[i].
 */
void SimpleFillPsi(double step, numpyint N, const double *EigenEs,
		numpyint EN, const double *V, const double *m, double* psis) {
	int i; 
#ifdef __MP
	#ifdef _DEBUG
	printf("Start a SimpleFillPsi with openMP.\n");
	#endif
	#pragma omp parallel for
#endif
	for(i=0; i<EN; i++) {
#ifdef __MP
	#ifdef _DEBUG
		printf("    From thread %d out of %d\n", 
				omp_get_thread_num(), omp_get_num_threads());
	#endif
#endif
		int j;
		double* psi = psis + i*N;
		double modsq = 0;
		ode(step, N, 0.0, Y_EPS, EigenEs[i], V, m, psi);
		/* Normalization */
		for(j=0; j<N; j++) {
			modsq += sq(psi[j]);
		}
		modsq = sqrt(modsq * step);
		for(j=0; j<N; j++) {
			psi[j] /= modsq;
		}
	}
	return; 
}

/** 
 * Find zero between x[n-1] and x[n] of function y(x), 
 * s.t. \f$ y(x_n) = 0 \f$
 * using linear interpolation, 
 * i.e assuming y[n] and y[n-1] are of opposite signs and
 * returning \f$ x_n \f$.
 */
#define findZero(x, y, n) ((y[n]*x[n-1] - y[n-1]*x[n])/ (y[n] - y[n-1]))


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
 * or calculate zero using Newton's Method if SIMPLE is not defined.
 *
 * Es should be in small to large order.
 *
 * \param[in] step step size
 * \param[in] N number of steps
 * \param[in] *Es initial search range of eigen energy
 * \param[in] EN number of eigen energy to find
 * \param[in] *V potential 
 * \param[in] *m effective mass
 * \param[out] *EigenE (output) eigen energy
 *
 */
numpyint SimpleSolve1D(double step, numpyint N, 
		const double *Es, numpyint EN, 
		const double *V, const double *m, 
		double *EigenE) {
	double *yend; 
	int NofZeros=0;
	int i;

	yend = (double *)malloc(EN * sizeof(double));
#ifdef __MP
	#ifdef _DEBUG
	printf("Start a simpleSolve1D with openMP.\n");
	#endif
	#pragma omp parallel
#endif
	{
		double *y = (double *)malloc(N * sizeof(double));
#ifdef __MP
	#ifdef _DEBUG
		printf("    From thread %d out of %d\n", 
				omp_get_thread_num(), omp_get_num_threads());
	#endif
		#pragma omp for
#endif
		for(i=0; i<EN; i++) {
			yend[i] = ode(step, N, 0.0, Y_EPS, Es[i], V, m, y);
		}

#ifdef __MP
		#pragma omp barrier
		#pragma omp for ordered
#endif
		for(i=1; i<EN; i++) {
			double E0;
			if(yend[i] == 0) {
				E0 = Es[i-1];
			}
			else if(yend[i]*yend[i-1] < 0) {
				E0 = findZero(Es, yend, i);
#ifndef SIMPLE
			int count=0; 
			double y0;
			y0 = ode(step, N, 0.0, Y_EPS, E0, V, m, y);
			/* Filter singular case */
			if(fabs(y0) > fabs(yend[i]) || fabs(y0) > fabs(yend[i-1])){
				continue;
			}
			while(fabs(y0) > NEWTON_S && count < 20){
				double y1 = ode(step, N, 0.0, Y_EPS, E0+NSTEP, V, m, y);
				double y2 = ode(step, N, 0.0, Y_EPS, E0-NSTEP, V, m, y);
				double dy = (y1 - y2)/(2*NSTEP);
				if(y1*y2 < 0) {
	#ifdef _DEBUG
					printf("  solution error smaller than step at E=%e,"
							" (count=%d).\n", E0, count);
	#endif
					E0 = (E0 + NSTEP)*y0*y2 / ( (y1 - y2)*(y1 - y0) )
						+ E0*y1*y2 / ( (y0 - y2)*(y0 - y1) )
						+ (E0 - NSTEP)*y1*y0 / ( (y2 - y0)*(y2 - y1) );
					break;
				}
				E0 -= y0/dy;
				y0 = ode(step, N, 0.0, Y_EPS, E0, V, m, y);
				count++;
			}
	#ifdef _DEBUG
			printf("After %d times Newton, E=%f Err=%e\n", 
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
	}

	free(yend);
	yend = NULL;
	return NofZeros;
}


#ifdef _WINDLL
__declspec(dllexport)
#endif 
/** 
 * Same as SimpleFillPsi except for using band related mass 
 */
void BandFillPsi(double step, numpyint N, const double *EigenEs,
		numpyint EN, double* psis, const double* V, Band* mat) {
	int i; 
#ifdef _DEBUG
	assert(N == mat->N);
#endif
#ifdef __MP
	#ifdef _DEBUG
	printf("Start a BandFillPsi with openMP.\n");
	#endif
	#pragma omp parallel 
#endif
	{
		double* m = (double*)malloc(N*sizeof(double));
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
			/*MP assume only mass is updated*/
			UpdateBand(mat, EigenEs[i], V, m);
			ode(step, N, 0.0, Y_EPS, EigenEs[i], V, m, psi);
			/* Normalization */
			for(j=0; j<N; j++) {
				modsq += sq(psi[j]);
			}
			modsq = sqrt(modsq * step);
			for(j=0; j<N; j++) {
				psi[j] /= modsq;
			}
		}
		free(m);
		m = NULL;
	}
	return; 
}


#ifdef _WINDLL
__declspec(dllexport)
#endif 
/**
 * Same as SimpleSolve1D except for using band related mass
 */
numpyint BandSolve1D(double step, numpyint N, 
		const double *Es, numpyint EN, const double *V, Band *mat,
		double *EigenE) {
       	double *yend; 
	int NofZeros=0;
	int i;

#ifdef _DEBUG
	assert(N == mat->N);
#endif
	yend = (double *)malloc(EN * sizeof(double));
#ifdef __MP
	#ifdef _DEBUG
	printf("Start a BandSolve1D with openMP.\n");
	#endif
	#pragma omp parallel 
#endif
	{
		double *y = (double *)malloc(N * sizeof(double));
		double *m = (double *)malloc(N * sizeof(double));
#ifdef __MP
	#ifdef _DEBUG
		printf("    From thread %d out of %d\n", 
				omp_get_thread_num(), omp_get_num_threads());
	#endif
		#pragma omp for
#endif
		for(i=0; i<EN; i++) {
			UpdateBand(mat, Es[i], V, m);
			yend[i] = ode(step, N, 0.0, Y_EPS, Es[i], V, m, y);
		}

#ifdef __MP
		#pragma omp barrier
		#pragma omp for ordered
#endif
		for(i=1; i<EN; i++) {
			double E0;
			if(yend[i] == 0) {
				E0 = Es[i-1];
			}
			else if(yend[i]*yend[i-1] < 0) {
				E0 = findZero(Es, yend, i);
#ifndef SIMPLE
				int count=0; 
				double y0;
				UpdateBand(mat, E0, V, m);
				y0 = ode(step, N, 0.0, Y_EPS, E0, V, m, y);
				/* Filter singular case */
				if(fabs(y0) > fabs(yend[i]) || fabs(y0) > fabs(yend[i-1])){
					continue;
				}
				while(fabs(y0) > NEWTON_S && count < 20){
					UpdateBand(mat, E0+NSTEP, V, m);
					double y1 = ode(step, N, 0.0, Y_EPS, E0+NSTEP, 
							V, m, y);
					UpdateBand(mat, E0-NSTEP, V, m);
					double y2 = ode(step, N, 0.0, Y_EPS, E0-NSTEP, 
							V, m, y);
					double dy = (y1 - y2)/(2*NSTEP);
					if(y1*y2 < 0) {
	#ifdef _DEBUG
						printf("  solution error smaller than step at E=%e,"
								" (count=%d).\n", E0, count);
	#endif
						break;
					}
					E0 -= y0/dy;
					UpdateBand(mat, E0, V, m);
					y0 = ode(step, N, 0.0, Y_EPS, E0, V, m, y);
					count++;
				}
	#ifdef _DEBUG
				printf("After %d times Newton, E=%f Err=%e\n", 
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
		free(m);
		m = NULL;
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

