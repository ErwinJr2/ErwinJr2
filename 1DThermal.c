#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#ifdef __MP /*openmp support*/
#include <omp.h>
#endif 
#include "science.h"

#ifdef __cplusplus
extern "C" {
#endif


#ifdef _WINDLL
__declspec(dllexport)
#endif // _WINDLL
double FermiDirac0(double EF, const double *EigenEs, numpyint EN, 
		const double *m, const double* psis, numpyint N, double step, 
		double* eDensity) {
	/* T=0 Fermi Driac distribution, 
	 * output electron density in eDensity (unit Angstrom^-3)
	 * and return sheet density (unit Angstrom^-2)
	 */
	int i;
	double sheetDensity = 0;
	/* double beta = 1/(kb*T); */
	for(i=0; i<N; i++) {
		eDensity[i] = 0;
	}
	for(i=0; i<EN; i++) {
		int j;
		/* double nbar = 1/(exp(beta * (EigenEs[i]-EF)) + 1); */
		const double* psi = psis + i*N;
		double invM = 0;  
		/* Density of states effective mass should be the Harmonic mean of
		 * effective mass */
		double DoS;
		if(EigenEs[i] > EF)
			break;
		for(j=0; j<N; j++) {
			invM += sq(psi[j])/m[j];
		}
		DoS = 1/(M_PI*sq(hbar)*invM); 
		/* DoS_2D = m/(pi*\hbar^2), spin included. */
		for(j=0; j<N; j++) {
			eDensity[j] += sq(psi[j])*DoS*(EF - EN);
			sheetDensity += eDensity[j]*step;
		}
	}
	return sheetDensity;
}

#ifdef _WINDLL
__declspec(dllexport)
#endif // _WINDLL
double Boltzmann(double T, double EF, const double *EigenEs, numpyint EN, 
		const double *m, const double* psis, numpyint N, double step,
		double* eDensity) {
	/* T=0 Maxwell Boltzmann distribution, 
	 * output electron density in eDensity (unit Angstrom^-3)
	 * and return sheet density (unit Angstrom^-2)
	 */
	const double MAXKT = 50;
	int i;
	double sheetDensity = 0;
	double kt = kb*T;
	/* if(EigenEs[0] < EF + 2*kt) { */
	/*     printf("Warning: " */
	/*             "Ground state lower than Fermi energy + Thermal scale\n"); */
	/*     printf("  E_F=%e; kbT=%e; E_0=%e\n", EF, kt, EigenEs[0]); */
	/* } */
	for(i=0; i<N; i++) {
		eDensity[i] = 0;
	}
	for(i=0; i<EN; i++) {
		int j;
		const double* psi = psis + i*N;
		double invM = 0;  
		/* Density of states effective mass should be the Harmonic mean of
		 * effective mass */
		double DoS;
		double DeltaE = EigenEs[i] - EF;
		if(DeltaE > MAXKT*kt)
			break;
		for(j=0; j<N; j++) {
			invM += sq(psi[j])/m[j];
		}
		DoS = 1/(M_PI*sq(hbar)*invM); 
		/* DoS_2D = m/(pi*\hbar^2), spin included. */
		for(j=0; j<N; j++) {
			eDensity[j] += sq(psi[j])*DoS*kt*exp(-DeltaE/kt);
			sheetDensity += eDensity[j]*step;
		}
	}
	return sheetDensity;
}

#ifdef _WINDLL
__declspec(dllexport)
#endif // _WINDLL
numpyint answer()
{return 42;}

#ifdef __cplusplus
}
#endif

