#ifndef BAND_H
#define BAND_H
#include <stdlib.h>
#include "science.h"
typedef struct BAND Band;

/* an update parameter function pointer that accepts a pointer to 
 * to band struct (with its parameters as members) and according to 
 * the Band parameters and energy, update para (usually effective mass)
 */
typedef numpyint (*UpdateFunc)(Band *, double); 

typedef struct BAND{
	/* Base class for band structure */
	const UpdateFunc update;
	numpyint N;           /* Size of datas */
	double *V;       /* Band bottom in eV */
	double *m;       /* Effective mass in electron mass m0 */
	double *Eg;      /* Band gap in eV */
} Band; 
numpyint UpdateBand(Band *, double);

Band *ZBband_new(numpyint N, const double *xEg, const double *xVc, 
		const double *xF, const double *xEp, const double *xESO);
void ZBband_free(Band *);

#endif /* ifndef BAND_H */
