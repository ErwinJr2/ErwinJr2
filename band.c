#include "band.h"

struct ZBBAND {
	/* Zincblende structure band */
	UpdateFunc updateM;
	int N;           /* Size of datas */
	const double *xVc;
	double *m;
	const double *xEg;      /* Band gap in eV */
	const double *xF;
	const double *xEp; 
	const double *xESO;
};


int ZBupdateM(Band *mat, double Eq) {
	ZBBand *zbmat = (ZBBand *) mat;
	int q; 
	for(q=0; q<zbmat->N; q++) {
		zbmat->m[q] = 1 / (1 + 2*zbmat->xF[q] + zbmat->xEp[q]/3 * 
				( 2 / (Eq - zbmat->xVc[q] + zbmat->xEg[q]) + 
				  1 / (Eq - zbmat->xVc[q] + zbmat->xEg[q] + zbmat->xESO[q]) 
				));
	}
	return zbmat->N;
}

ZBBand *ZBband_new(int N, const double *xEg, const double *xVc, 
		const double *xF, const double *xEp, const double *xESO, double *x) {
	ZBBand *zbband = (ZBBand *) malloc( sizeof(ZBBand) );
	zbband->updateM = ZBupdateM;
	zbband->N = N; 
	zbband->xEg = xEg; 
	zbband->xVc = xVc; 
	zbband->xF = xF; 
	zbband->xEp = xEp; 
	zbband->xESO = xESO;
	zbband->m = (double *)malloc( N*sizeof(double) );
	return zbband; 
}

