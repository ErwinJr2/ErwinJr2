#include "band.h"
/* Update effective mass of a Zincblende band semiconductor */
numpyint ZBupdateM(Band *mat, double Eq) {
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

ZBBand *ZBband_new(numpyint N, const double *xEg, const double *xVc, 
		const double *xF, const double *xEp, const double *xESO) {
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

void ZBband_free(ZBBand * zbband) {
	free(zbband->m);
	free(zbband);
	return;
}

