#include <stdlib.h>
#ifndef BAND_H
#define BAND_H
typedef struct BAND Band;

/* an update parameter function pointer that accepts a pointer to 
 * to band struct (with its parameters as members) and according to 
 * the Band parameters and energy, update para (usually effective mass)
 */
typedef int (*UpdateFunc)(Band *, double); 

typedef struct BAND{
	/* Base class for band structure */
	const UpdateFunc update;
	int N;           /* Size of datas */
	double *V;       /* Band bottom in eV */
	double *m;       /* Effective mass in electron mass m0 */
	double *Eg;      /* Band gap in eV */
} Band; 

/* Zincblende structure band, compatiable with structure BAND */
typedef struct ZBBAND ZBBand; 
ZBBand *ZBband_new(int N, const double *xEg, const double *xVc, 
		const double *xF, const double *xEp, const double *xESO, double *x);
/* Update effective mass of a Zincblende band semiconductor */
int ZBupdateM(Band *, double Eq); 

#endif /* ifndef BAND_H */
