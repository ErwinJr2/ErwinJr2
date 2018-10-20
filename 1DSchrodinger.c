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

#define NEWTON 1E-3 
/* Newton's method stopping limit */
#define NSTEP 1E-10
/* NSTEP (unit eV) is the step size to numerically calculate  
 * the derivative for Newton's method. 
 * Optimum step size = (6*epsilon/M3)^(1/3) 
 * where epsilon is numerical error of yend(E) and M3 is maximum yend'''(E)
 * Because of the exponantial behavior when V>E, yend is very sensitive 
 * to E near EigenE. M3 is very large. 
 * Still, it's recommanded to make high V side the starting point
 */
#define Y_EPS 0.1 /* Default starting point for Numerov */

#ifdef __cplusplus
extern "C" {
#endif


#ifdef _WINDLL
__declspec(dllexport)
#endif 
double Numerov(double step, numpyint N, double y0, double y1, 
		double E, const double *V, const double *m, double *y) {
	/* An ODE solver for -hbar^2/(2*m(x)) * y''(x) + V(x) * y = E * y(x)
	 * with starting x0 and y0, y1, ends at x0 + step*N (x0 label 0)
	 * using Numerov algorithm
	 * E and V are in unit eV, m are in unit m0 (free electron mass)
	 * Don't normalize
	 * V[n] and m[n] means potential and effectice mass at x = x0 + n*step
	 * return y(x+N*step), put y result in *y
	 */
	int n; 
	y[0] = y0;
	y[1] = y1;
	const double unit = 2*m0/sq(hbar)*e0*sq(ANG*step);
	for (n = 1; n < N-1; n++) {
		y[n+1] = (2 * y[n] * (1.0 - 5.0/12 * ( E - V[n]) * unit * m[n]) 
			 - y[n-1] * (1.0 + 1.0/12 * (E - V[n-1]) * unit * m[n])) 
			/ (1.0 + 1.0/12 * (E - V[n+1]) * unit * m[n]);
	}
	return y[N-1];
}


#ifdef _WINDLL
__declspec(dllexport)
#endif 
void SimpleFillPsi(double step, numpyint N, const double *EigenEs,
		numpyint EN, const double *V, const double *m, double* psis) {
	/* Fill in wavefunctions in psis accroding to eigen energy in EigenEs. 
	 * psi + i*N*sizeof(double) is the wavefunction with Energy EigenEs[i] 
	 * The result is normalized to 1
	 */
	int i; 
#ifdef __MP
	#pragma omp parallel for
#endif
	for(i=0; i<EN; i++) {
		int j;
		double* psi = psis + i*N;
		double modsq = 0;
		Numerov(step, N, 0.0, Y_EPS, EigenEs[i], V, m, psi);
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


#define findZero(x, y, n) ((y[n]*x[n-1] - y[n-1]*x[n])/ (y[n] - y[n-1]))
/* Find zero x_n between x[n-1] and x[n] of function y(x), 
 * s.t. y(x_n) = 0
 * using linear interpolation (to be improved?)
 * assuming y[n] and y[n-1] is opposite sign
 * return x_n
 */


#ifdef _WINDLL
__declspec(dllexport)
#endif 
numpyint SimpleSolve1D(double step, numpyint N, 
		const double *Es, numpyint EN, 
		const double *V, const double *m, 
		double *EigenE) {
	/* Solve 1D shrodinger's equation with potential V, effective mass m and 
	 * in the region x0 <= x < x0+step*N with zero boundary. 
	 * First scan in energy Es[0:EN] and look for zeros(EigenE) by either 
	 *    simple linear interpolation if SIMPLE is not defined;
	 * or calculate zero using newton's method if SIMPLE is defined
	 * Es should be in small to large order
	 */
	double *y; 
	double *yend; 
	int NofZeros=0;
	int i;

	yend = (double *)malloc(EN * sizeof(double));
#ifdef __MP
	#pragma omp parallel for private(y)
	for(i=0; i<EN; i++) {
		y = (double *)malloc(N * sizeof(double));
#else 
	y = (double *)malloc(N * sizeof(double));
	for(i=0; i<EN; i++) {
#endif
		yend[i] = Numerov(step, N, 0.0, Y_EPS, Es[i], V, m, y);
#ifdef __MP
		free(y);
		y = NULL;
#endif
	}

#ifdef __MP
	#pragma omp parallel for private(y) ordered
#endif
	for(i=1; i<EN; i++) {
		double E0;
		if(yend[i] == 0) {
			E0 = Es[i-1];
		}
		else if(yend[i] == Es[i-1]) {
			continue;
		}
		else if(yend[i]*yend[i-1] < 0) {
			E0 = findZero(Es, yend, i);
#ifndef SIMPLE
			int count=0; 
			double y0;
	#ifdef __MP
			y = (double *)malloc(N * sizeof(double));
	#endif
			y0 = Numerov(step, N, 0.0, Y_EPS, E0, V, m, y);
			while(fabs(y0) > NEWTON && count < 20){
				double y1 = Numerov(step, N, 0.0, Y_EPS, E0+NSTEP, V, m, y);
				double y2 = Numerov(step, N, 0.0, Y_EPS, E0-NSTEP, V, m, y);
				double dy = (y1 - y2)/(2*NSTEP);
				if(y1*y2 < 0) {
	#ifdef __DEBUG
					printf("  solution error smaller than step at E=%e,"
							" (count=%d).\n", E0, count);
	#endif
					break;
				}
				E0 -= y0/dy;
				y0 = Numerov(step, N, 0.0, Y_EPS, E0, V, m, y);
				count++;
			}
	#ifdef __DEBUG
			printf("After %d times Newton, E=%f Err=%e\n", 
					count, E0, fabs(y0));
	#endif
	#ifdef __MP
			free(y);
			y = NULL;
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

#ifndef __MP
	free(y);
#endif
	free(yend);
	return NofZeros;
}


#ifdef _WINDLL
__declspec(dllexport)
#endif 
void BandFillPsi(double step, numpyint N, const double *EigenEs,
		numpyint EN, double* psis, Band* mat) {
	/* Same as SimpleFillPsi except for using band related mass and V */
	int i; 
#ifdef __DEBUG
	assert(N == mat->N);
#endif
#ifdef __MP
	#pragma omp parallel for
#endif
	for(i=0; i<EN; i++) {
		int j;
		double* psi = psis + i*N;
		double modsq = 0;
		mat->update(mat, EigenEs[i]);
		Numerov(step, N, 0.0, Y_EPS, EigenEs[i], mat->V, mat->m, psi);
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


#ifdef _WINDLL
__declspec(dllexport)
#endif 
numpyint BandSolve1D(double step, numpyint N, 
		const double *Es, numpyint EN, Band* mat,
		double *EigenE) {
	/* Solve 1D shrodinger's equation with potential V, effective mass m and 
	 * in the region x0 <= x < x0+step*N with zero boundary. 
	 * First scan in energy Es[0:EN] and look for zeros(EigenE) by either 
	 *    simple linear interpolation if SIMPLE is not defined;
	 * or calculate zero using newton's method if SIMPLE is defined
	 * Es should be in small to large order
	 */
	double *y; 
	double *yend; 
	int NofZeros=0;
	int i;

#ifdef __DEBUG
	assert(N == mat->N);
#endif
	yend = (double *)malloc(EN * sizeof(double));
#ifdef __MP
	#pragma omp parallel for private(y)
	for(i=0; i<EN; i++) {
		y = (double *)malloc(N * sizeof(double));
#else 
	y = (double *)malloc(N * sizeof(double));
	for(i=0; i<EN; i++) {
#endif
		mat->update(mat, Es[i]);
		yend[i] = Numerov(step, N, 0.0, Y_EPS, Es[i], mat->V, mat->m, y);
#ifdef __MP
		free(y);
		y = NULL;
#endif
	}

#ifdef __MP
	#pragma omp parallel for private(y) ordered
#endif
	for(i=1; i<EN; i++) {
		double E0;
		if(yend[i] == 0) {
			E0 = Es[i-1];
		}
		else if(yend[i] == Es[i-1]) {
			continue;
		}
		else if(yend[i]*yend[i-1] < 0) {
			E0 = findZero(Es, yend, i);
#ifndef SIMPLE
			int count=0; 
			double y0;
	#ifdef __MP
			y = (double *)malloc(N * sizeof(double));
	#endif
			mat->update(mat, E0);
			y0 = Numerov(step, N, 0.0, Y_EPS, E0, mat->V, mat->m, y);
			while(fabs(y0) > NEWTON && count < 20){
				mat->update(mat, E0+NSTEP);
				double y1 = Numerov(step, N, 0.0, Y_EPS, E0+NSTEP, 
						mat->V, mat->m, y);
				mat->update(mat, E0-NSTEP);
				double y2 = Numerov(step, N, 0.0, Y_EPS, E0-NSTEP, 
						mat->V, mat->m, y);
				double dy = (y1 - y2)/(2*NSTEP);
				if(y1*y2 < 0) {
	#ifdef __DEBUG
					printf("  solution error smaller than step at E=%e,"
							" (count=%d).\n", E0, count);
	#endif
					break;
				}
				E0 -= y0/dy;
				mat->update(mat, E0);
				y0 = Numerov(step, N, 0.0, Y_EPS, E0, mat->V, mat->m, y);
				count++;
			}
	#ifdef __DEBUG
			printf("After %d times Newton, E=%f Err=%e\n", 
					count, E0, fabs(y0));
	#endif
	#ifdef __MP
			free(y);
			y = NULL;
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

#ifndef __MP
	free(y);
#endif
	free(yend);
	return NofZeros;
}

#ifdef _WINDLL
__declspec(dllexport)
#endif 
numpyint invAlpha()
{return 137;}


#ifdef __cplusplus
}
#endif

