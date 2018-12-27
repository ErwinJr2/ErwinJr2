/**
 *
 * \file
 *
 * \brief Compute thermodynamic statistics
 *
 *
 */

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
#endif 
 
/**
 * \brief T=0 Fermi-Dirac distribution. Fermi energy -> electron density
 *
 * T=0 Fermi-Dirac distribution, given Fermi energy
 * output electron density in eDensity (unit Angstrom^-3)	 
 * and return sheet density (unit Angstrom^-2)
 */
double FermiDirac0(double EF, const double *EigenEs, numpyint EN, 
		const double *m, const double* psis, numpyint N, double step, 
		double* eDensity) {
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
		double DoS2D;
		if(EigenEs[i] > EF)
			break;
		for(j=0; j<N; j++) {
			invM += sq(psi[j])/m[j];
		}
		invM *= step;
		DoS2D = m0/(M_PI*sq(hbar)*invM)*sq(ANG)*e0; 
		/* DoS2D_2D = m/(pi*\hbar^2), spin included, Unit Angstrom^-2eV^-1 */
		for(j=0; j<N; j++) {
		        eDensity[j] += sq(psi[j])*DoS2D*(EF-EigenEs[i]);
			/* sheetDensity += eDensity[j]*step; --->X */
		}
		sheetDensity += DoS2D*(EF-EigenEs[i]);
		/* psi is normalized.. It's equivlent with X */
	}
	return sheetDensity;
}


#ifdef _WINDLL
__declspec(dllexport)
#endif
 
/**
 * \brief T=0 Fermi-Dirac distribution. Sheet density -> Fermi energy
 *
 * T=0 Fermi-Dirac distribution, given sheet density
 * output electron density in eDensity (unit Angstrom^-3)
 * and return EF
 */
double FermiDirac0N(double sheet, const double *EigenEs, numpyint EN, 
		const double *m, const double* psis, numpyint N, double step, 
		double* eDensity) {
       	int i;
	double sheetDensity = 0;
	double DoS2Dsum = 0;
	for(i=0; i<N; i++) {
		eDensity[i] = 0;
	}
	for(i=0; i<EN; i++) {
		int j;
		const double* psi = psis + i*N;
		double invM = 0;  
		/* Density of states effective mass should be the Harmonic mean of
		 * effective mass */
		double Emax;
		for(j=0; j<N; j++) {
			invM += sq(psi[j])/m[j];
		}
		invM *= step;
		DoS2Dsum += m0/(M_PI*sq(hbar)*invM)*sq(ANG)*e0; 
		/* DoS2D_2D = m/(pi*\hbar^2), spin included, Unit Angstrom^-2eV^-1 */
#ifdef _DEBUG 
		printf("Energy at %e, Dossum=%e, m* = %e\n", 
				EigenEs[i], DoS2Dsum, 1/invM);
#endif
		Emax = (sheet - sheetDensity)/(DoS2Dsum); 
		if(i == EN-1 || EigenEs[i] + Emax <= EigenEs[i+1]) {
			for(j=0; j<N; j++) {
				eDensity[j] += sq(psi[j])*DoS2Dsum*Emax;
			}
			return EigenEs[i] + Emax;
		}
		for(j=0; j<N; j++) {
		  eDensity[j] += sq(psi[j])*DoS2Dsum*(EigenEs[i+1]-EigenEs[i]);
		}
		sheetDensity += DoS2Dsum * (EigenEs[i+1]-EigenEs[i]);
	}
	printf("Error: You shouldn't reach here!\n");
	return NAN;
}

#ifdef _WINDLL
  __declspec(dllexport)
#endif
/**                                                                                     
 * \brief Finite Fermi-Dirac distribution. Fermi energy -> electron density
 *                                                                                      
 * \param EF Fermi energy
 * \param eDensity output electron density(unit Angstrom^-3)
 * \output sheet density (unit Angstrom^-2)
 */
double FermiDirac(double T, double EF, const double *EigenEs, numpyint EN,
		     const double *m, const double* psis, numpyint N, double step,
		     double* eDensity) {
    int i;
    double sheetDensity = 0;
    double beta = 1/(kb*T);
    for(i=0; i<N; i++) {
      eDensity[i] = 0;
    }
    for(i=0; i<EN; i++) {
      int j;
      double nbar = 1/(exp(beta * (EigenEs[i]-EF)) + 1);
      const double* psi = psis + i*N;
      double invM = 0;
      /* Density of states effective mass should be the Harmonic mean of      
       * effective mass */
      double DoS2D;
      for(j=0; j<N; j++) {
	invM += sq(psi[j])/m[j];
      }
      invM *= step;
      DoS2D = m0/(M_PI*sq(hbar)*invM)*sq(ANG)*e0;
      /* DoS2D_2D = m/(pi*\hbar^2), spin included, Unit Angstrom^-2eV^-1 */
      for(j=0; j<N; j++) {
	eDensity[j] += sq(psi[j])*DoS2D*nbar;
	/* sheetDensity += eDensity[j]*step; --->X */
      }
      sheetDensity += DoS2D*nbar;
      /* psi is normalized.. It's equivlent with X */
    }
    return sheetDensity;
}


#ifdef _WINDLL
__declspec(dllexport)
#endif 

/**
 * \brief Maxwell-Boltzmann distribution at T. Fermi level -> electron density
 *
 * Temperature T Maxwell Boltzmann distribution, given Fermi level EF
 * output electron density in eDensity (unit Angstrom^-3)
 * and return sheet density (unit Angstrom^-2)
 */
double Boltzmann(double T, double EF, const double *EigenEs, numpyint EN, 
		const double *m, const double* psis, numpyint N, double step,
		double* eDensity) {
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
		double DoS2D;
		double DeltaE = EigenEs[i] - EF;
		if(DeltaE > MAXKT*kt)
			break;
		for(j=0; j<N; j++) {
			invM += sq(psi[j])/m[j];
		}
		invM *= step;
		DoS2D = m0/(M_PI*sq(hbar)*invM)*sq(ANG)*e0; 
		/* DoS2D_2D = m/(pi*\hbar^2), spin included, Unit Angstrom^-2eV^-1 */
		for(j=0; j<N; j++) {
			eDensity[j] += sq(psi[j])*DoS2D*kt*exp(-DeltaE/kt);
			/* sheetDensity += eDensity[j]*step; --> X */
		}
		sheetDensity += DoS2D*kt*exp(-DeltaE/kt); 
		/* psi is normalized.. It's equivlent with X */
	}
	return sheetDensity;
}


#ifdef _WINDLL
__declspec(dllexport)
#endif 

/**
 * \brief Maxwell-Boltzmann distribution at T. sheet density -> EF
 *
 * Temperature T Maxwell Boltzmann distribution, given sheet density
 * output electron density in eDensity (unit Angstrom^-3)
 * and return EF (unit eV)
 */
double BoltzmannN(double T, double sheet, const double *EigenEs, numpyint EN, 
		const double *m, const double* psis, numpyint N, double step,
		double* eDensity) {
	double sheet0 = Boltzmann(T, 0, EigenEs, EN, m, psis, N, step, eDensity); 
	double EF = kb * T * log(sheet / sheet0);
#ifdef _DEBUG
	printf("Sheet density = %e at EF = 0, so EF = %e\n", sheet0, EF);
#endif
	int i;
	for(i=0; i<N; i++) {
		eDensity[i] *= sheet/sheet0;
	}
	return EF;
}


#ifdef _WINDLL
__declspec(dllexport)
#endif 
numpyint answer()
{return 42;}


#ifdef __cplusplus
}
#endif

