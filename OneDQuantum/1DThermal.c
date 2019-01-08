/**
 *
 * \file
 *
 * Compute thermodynamic statistics given eigen energy and wavefunctions
 * for ??? systems.
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

#define NEWTON 1E-5 /**< Newton's method stopping limit. Default = 1E-5 */

#ifdef _WINDLL
__declspec(dllexport)
#endif 
 
/**
 *
 * Given Fermi energy EF, this function
 * outputs electron density in eDensity (unit Angstrom^-3)
 * and return sheet density (unit Angstrom^-2).
 *
 * Assumes Fermi-Dirac distribution at zero temperature.
 *
 * \param[in] EF Fermi energy
 * \param[in] *EigenEs Eigen energy
 * \param[in] EN number of eigen energies provided
 * \param[in] *m effective mass
 * \param[in] *psis wavefunctions for the given eigen energies
 * \param[in] N number of steps
 * \param[in] step step size
 * \param[in] *eDensity (output) electron density (unit Angstrom^-3)
 */
double FermiDirac0(double EF, const double *EigenEs, numpyint EN, 
		const double *m, const double* psis, numpyint N, double step, 
		double *eDensity) {
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
 *
 * Given sheet density (unit Angstrom^-2), 
 * outputs electron density in eDensity (unit Angstrom^-3)
 * and returns Fermi energy EF.
 *
 * Assumes Fermi-Dirac distribution at zero temperature.
 *
 * \param[in] sheet sheet density (unit Angstrom^-2)
 * \param[in] *EigenEs Eigen energy
 * \param[in] EN number of eigen energies provided
 * \param[in] *m effective mass
 * \param[in] *psis wavefunctions for the given eigen energies
 * \param[in] N number of steps
 * \param[in] step step size
 * \param[in] *eDensity (output) electron density (unit Angstrom^-3)
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
 * Given Fermi energy EF, this function
 * outputs electron density in eDensity (unit Angstrom^-3)
 * and return sheet density (unit Angstrom^-2)
 *
 * Assumes Fermi-Dirac distribution at finite temperature
 *
 * \param[in] T temperature
 *
 * Other parameters are the same as in FermiDirac0.
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
    for(i=0; i<1; i++) {
      int j;
      //double nbar = 1/(exp(beta * (EigenEs[i]-EF)) + 1);
      double nbarIntegral = EF - EigenEs[i] + log(1+exp(beta*(EigenEs[i]-EF)))/beta;
      if (isinf(nbarIntegral)) {
	  nbarIntegral = 0.0;
	} // dirty way to get around log(1+exp)
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
	eDensity[j] += sq(psi[j])*DoS2D*nbarIntegral;
	/* sheetDensity += eDensity[j]*step; --->X */
      }
      sheetDensity += DoS2D*nbarIntegral;
      /* psi is normalized.. It's equivlent with X */
    }
    return sheetDensity;
}

#ifdef _WINDLL
__declspec(dllexport)
#endif
/** 
 * Derivative of sheet density against Fermi energy at finite temperature,
 * using Fermi-Dirac distribution. 
 * This derivative is used to solve Fermi energy given sheet density.
 *
 */
double DFermiDirac(double T, double EF, const double *EigenEs, numpyint EN, 
		   const double *m, const double* psis, numpyint N, double step) {
    int i;
    double dSheetDensity = 0.0;
    double beta = 1.0/(kb*T);
    //for(i=0; i<N; i++) {
    //  eDensity[i] = 0.0;
    //}
    for(i=0; i<1; i++) {
      int j;
      double dnbarIntegral = 1.0 - exp(beta*(EigenEs[i]-EF))/(1.0+exp(beta*(EigenEs[i]-EF)));
      // if isinf(nbarIntegral) nbarIntegral = 0.0; // dirty way to get around log(1+exp) 
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
      dSheetDensity += DoS2D*dnbarIntegral;
      /* psi is normalized.. It's equivlent with X */
    }
    return dSheetDensity;
  }


#ifdef _WINDLL
__declspec(dllexport)
#endif
/** 
 * Given temperature and sheet density (unit Angstrom^-2),
 * outputs electron density in eDensity (unit Angstrom^-3)
 * and returns Fermi energy using Newton's method.
 *
 * Assumes Fermi-Dirac distribution at finite temperature.
 *
 * Other parameters are the same as in FermiDirac0N.
 *
 * \param[in] T temperature.
 */
double FermiDiracN(double T, double sheet, const double *EigenEs, numpyint EN, 
		    const double *m, const double* psis, numpyint N, double step, 
		    double* eDensity) {
    int i;
    /* double sheetDensity; */
    for(i=0; i<N; i++) {
        eDensity[i] = 0;
    }
    double EF_min = EigenEs[0];
    double EF_max = EigenEs[EN-1]+(EigenEs[EN-1]-EigenEs[0])*3;
    double EF0 = (EF_min + EF_max)/2.0;
    FermiDirac(T, EF_max, EigenEs, EN, m, psis, N, step, eDensity);
    FermiDirac(T, EF_min, EigenEs, EN, m, psis, N, step, eDensity);
    double dsheet = FermiDirac(T, EF0, EigenEs, EN, m, 
            psis, N, step, eDensity)-sheet;
    /** Use Newton's method to find EF */
    int count=0; double EF;
    EF = EF0;
    while (count < 20 && fabs(dsheet) > NEWTON) {
        dsheet = FermiDirac(T, EF0, EigenEs, EN, m, psis, 
                N, step, eDensity) - sheet;
        EF = EF0 - dsheet / DFermiDirac(T, EF0, EigenEs, EN, 
                m, psis, N, step);
        EF0 = EF;
        count++;
    }

    FermiDirac(T, EF, EigenEs, EN, m, psis, N, step, eDensity);
    return EF;
}

#ifdef _WINDLL
__declspec(dllexport)
#endif 

/**
 *
 * Given tempereature T and Fermi level EF,
 * output electron density in eDensity (unit Angstrom^-3)
 * and return sheet density (unit Angstrom^-2).
 *
 * Assumes Maxwell-Boltzmann distribution. Good approximation to Fermi-Dirac
 * distribution at high temperature.
 *
 * Parameters are the same as in FermiDirac.
 *
 * Note: assumes zero filling when (eigen energy - EF) > MAXKT * k_B T
 * Default value for MAXKT = 50
 * 
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
	for(i=0; i<1; i++) {
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
 *
 * Given temperature T and sheet density (unit Angstrom^-2),
 * output electron density in eDensity (unit Angstrom^-3)
 * and return EF (unit eV).
 *
 * Assumes Maxwell-Boltzmann distribution. 
 * 
 * Parameters are the same as in FermiDiracN.
 *
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
/**
 * Checkpoint for python-C interface. Output = 42.
 */
numpyint answer()
{return 42;}


#ifdef __cplusplus
}
#endif

