#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#ifdef __MP /*openmp support*/
#include <omp.h>
#endif 

#define sq(X) ((X)*(X))
#define NEWTON 1E-3 /* Newton's method stopping limit */  
#define NSTEP 1E-8
/* NSTEP (unit eV) is the step size to numerically calculate  
 * the derivative for Newton's method */

const double hbar = 1.0545718e-34; /*J.s*/
const double m0 = 9.10938356e-31;  /*kg*/
const double e0 = 1.60217662e-19;  /*C*/
const double pi = 3.1415926535897932385;
const double ANG = 1E-10;  
/* Angstrom in meter: all unit for length is Angstrom in the program */

#ifdef _WINDLL
typedef int32_t numpyint;
#else
typedef int64_t numpyint;
#endif // _WINDLL

#ifdef __cplusplus
extern "C" {
#endif


#ifdef _WINDLL
__declspec(dllexport)
#endif // _WINDLL
double Numerov(numpyint x0, double step, numpyint N, double y0, double y1, 
		double E, const double *V, const double *m, double *y) {
	/* 
	 * An ODE solver for -hbar^2/(2*m(x)) * y''(x) + V(x) * y = E * y(x)
	 * with starting x0 and y0, y1, ends at x0 + step*N
	 * using Numerov algorithm
	 * E and V are in unit eV, m are in unit m0 (free electron mass)
	 * Don't normalize
	 * V[n] and m[n] means potential and effectice mass at x = x0 + n*step
	 * return y(x+N*step), put y result in *y
	 */
	int n; 
	y[x0] = y0;
	y[x0+1] = y1;
	const double unit = 2*m0/sq(hbar)*e0*sq(ANG*step);
	for (n = x0+1; n < N-1; n++) {
		y[n+1] = (2 * y[n] * (1.0 - 5.0/12 * ( E - V[n]) * unit * m[n]) 
			 - y[n-1] * (1.0 + 1.0/12 * (E - V[n-1]) * unit * m[n])) 
			/ (1.0 + 1.0/12 * (E - V[n+1]) * unit * m[n]);
	}
	return y[N-1];
}

double findZero(const double *x, const double *y, int n) {
	/*
	 * Find zero x_n between x[n-1] and x[n] of function y(x), 
	 * s.t. y(x_n) = 0
	 * using linear interpolation (to be improved?)
	 * assuming y[n] and y[n-1] is opposite sign
	 * return x_n
	 */
	return (y[n]*x[n-1] - y[n-1]*x[n])/ (y[n] - y[n-1]); 
}

#ifdef _WINDLL
__declspec(dllexport)
#endif // _WINDLL
numpyint SimpleSolve1D(numpyint x0, double step, numpyint N, 
		const double *Es, numpyint EN, 
		const double *V, const double *m, 
		double *EigenE) {
	/* Solve 1D shrodinger's equation with potential V, effective mass m and 
	 * in the region x0 <= x < x0+step*N with zero boundary. 
	 * First scan in energy Es[0:EN] and look for zeros(EigenE) by either 
	 *    simple linear interpolation if NEWTON is not defined;
	 * or calculate zero using newton's method if NEWTON is defined
	 * Es should be in small to large order
	 */
	double *y; 
	double *yend; 
	int NofZeros=0;
	int i;

	y = (double *)malloc(N * sizeof(double));
	yend = (double *)malloc(EN * sizeof(double));
	for(i=0; i<EN; i++) {
		yend[i] = Numerov(x0, step, N, 0.0, 1.0, Es[i], V, m, y);
	}

	for(i=1; i<N; i++) {
		if(yend[i] == 0) {
			EigenE[NofZeros++] = Es[i-1]; 
			continue; 
		}
		else if(yend[i] == Es[i-1]) {
			continue;
		}
		else if(yend[i]*yend[i-1] < 0) {
#ifndef NEWTON
			EigenE[NofZeros++] = findZero(Es, yend, i);
#else
			double E0 = findZero(Es, yend, i);
			double E_step = NSTEP;
			double y0 = Numerov(x0, step, N, 0.0, 1.0, E0, V, m, y);
			while(abs(y0) > NEWTON){
				double dy = (Numerov(x0, step, N, 0.0, 1.0, E0+E_step, 
							V, m, y) - y0)/E_step;
				E0 -= y0/dy;
				y0 = Numerov(x0, step, N, 0.0, 1.0, E0, V, m, y);
			}
			EigenE[NofZeros++] = E0;
#endif
		}
	}

	free(y);
	free(yend);
	return NofZeros;
}

#ifdef _WINDLL
__declspec(dllexport)
#endif // _WINDLL
numpyint answer()
{return 42;}

#ifdef __cplusplus
}
#endif
