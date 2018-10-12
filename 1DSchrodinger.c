#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#ifdef __MP /*openmp support*/
#include <omp.h>
#endif 

#define sq(X) ((X)*(X))
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

const double hbar = 1.0545718e-34; /*J.s*/
const double m0 = 9.10938356e-31;  /*kg*/
const double e0 = 1.60217662e-19;  /*C*/
/* const double pi = 3.1415926535897932385; M_PI in math.h*/
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
double Numerov(double step, numpyint N, double y0, double y1, 
		double E, const double *V, const double *m, double *y) {
	/* 
	 * An ODE solver for -hbar^2/(2*m(x)) * y''(x) + V(x) * y = E * y(x)
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
#endif // _WINDLL
void FillPsi(double step, numpyint N, const double *EigenEs, numpyint EN, 
		const double *V, const double *m, double* psis) {
	/* 
	 * Fill in wavefunctions in psis accroding to eigen energy in EigenEs. 
	 * psi + i*N*sizeof(double) is the wavefunction with Energy EigenEs[i] 
	 * The result is normalized to 1
	 */
	int i; 
	double modsq;
	for(i=0; i<EN; i++) {
		int j;
		double* psi = psis + i*N;
		Numerov(step, N, 0.0, Y_EPS, EigenEs[i], V, m, psi);
		modsq = 0; 
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

	y = (double *)malloc(N * sizeof(double));
	yend = (double *)malloc(EN * sizeof(double));
	for(i=0; i<EN; i++) {
		yend[i] = Numerov(step, N, 0.0, Y_EPS, Es[i], V, m, y);
	}

	for(i=1; i<EN; i++) {
		if(yend[i] == 0) {
			EigenE[NofZeros++] = Es[i-1]; 
			continue; 
		}
		else if(yend[i] == Es[i-1]) {
			continue;
		}
		else if(yend[i]*yend[i-1] < 0) {
#ifdef SIMPLE
			EigenE[NofZeros++] = findZero(Es, yend, i);
#else
			int count=0; 
#ifdef __DEBUG
			printf("Looking for zeros near %d\n", i);
#endif
			double E0 = findZero(Es, yend, i);
			double y0 = Numerov(step, N, 0.0, Y_EPS, E0, V, m, y);
			while(fabs(y0) > NEWTON){
				double y1 = Numerov(step, N, 0.0, Y_EPS, E0+NSTEP, V, m, y);
				double y2 = Numerov(step, N, 0.0, Y_EPS, E0-NSTEP, V, m, y);
				double dy = (y1 - y2)/(2*NSTEP);
				if(y1*y2 < 0) {
#ifdef __DEBUG
					printf("solution error smaller than step.\n");
#endif
					break;
				}
				E0 -= y0/dy;
				y0 = Numerov(step, N, 0.0, Y_EPS, E0, V, m, y);
				count++;
#ifdef __DEBUG
				printf("  The %d-th try for newton, E=%.20f Err=%f\n", 
						count, E0, fabs(y0));
#endif
				if(count > 20) {
#ifdef __DEBUG
					printf("Time out for Newton's method. %d\n", count);
#endif
					break;
				}
			}
#ifdef __DEBUG
			printf("After Newton, E=%f Err=%f\n", E0, fabs(y0));
#endif
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
