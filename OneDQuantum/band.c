/**
 * \file
 *
 * \brief Zincblende and Wurtzite structure band
 *
 * Zincblende and Wurtzite structure band, compatible
 * with structure BAND.
 */


#include "band.h"

#ifdef __cplusplus
extern "C" {
#endif

/** \brief Update effective mass in band */
numpyint UpdateBand(Band *band, double E, const double *xVc, double *m) {
	return band->update(band, E, xVc, m);
}
typedef struct ZBBAND {
        UpdateFunc updateM;    /**< Update effective mass */
        numpyint N;            /**< Number of finite x positions */
        const double *xEg;     /**< Direct energy gap  */
        const double *xF;      /**< Kane parameter  */
	const double *xEp; 
        const double *xESO;    /**< Spin-orbit splitting */
}ZBBand; 

/** \brief  Update effective mass of a Zincblende band semiconductor */
numpyint ZBupdateM(Band *mat, double Eq, const double *xVc, double *m) {
	ZBBand *zbmat = (ZBBand *) mat;
	int q; 
	for(q=0; q<zbmat->N; q++) {
		m[q] = 1 / (1 + 2*zbmat->xF[q] + zbmat->xEp[q]/3 * ( 
					2 / (Eq - xVc[q] + zbmat->xEg[q]) + 
					1 / (Eq - xVc[q] + zbmat->xEg[q] + zbmat->xESO[q])) );
	}
	return zbmat->N;
}

Band *ZBband_new(numpyint N, const double *xEg, const double *xF,
		const double *xEp, const double *xESO) {
	ZBBand *zbband = (ZBBand *) malloc( sizeof(ZBBand) );
	zbband->updateM = ZBupdateM;
	zbband->N = N; 
	zbband->xEg = xEg; 
	zbband->xF = xF; 
	zbband->xEp = xEp; 
	zbband->xESO = xESO;
	return (Band *) zbband; 
}

void ZBband_free(Band *zbband) {
	free( (ZBBand *) zbband );
	return;
}

#ifdef _DEBUG
#include <stdio.h>
void ZBband_check(const Band *band, numpyint N, const double *xEg, 
		const double *xF, const double *xEp, const double *xESO) {
	printf("Checking ZBband\n");
	const ZBBand *zbband = (const ZBBand *) band;
	if(zbband->updateM != ZBupdateM)
		printf("ZBupdateM checkfail\n");
	if(zbband->N != N)
		printf("N checkfail\n");
	if(zbband->xEg != xEg) 
		printf("xEg checkfail\n");
	if(zbband->xF != xF)
		printf("xF checkfail\n");
	if(zbband->xEp != xEp)
		printf("xEp checkfail\n");
	if(zbband->xESO != xESO)
		printf("xESO checkfail\n");
	printf("Checking Finished\n");
	return;
}
#endif

typedef struct WZBAND {
  UpdateFunc updateM;
  numpyint N;
  const double *xEg;
  const double *xEp;
  const double *xESO;
}WZBAND;

/** \brief  Update effective mass of a Wurtzite band semiconductor */
numpyint ZBupdateM(Band *mat, double Eq, const double *xVc, double *m) {
    ZBBand *zbmat = (ZBBand *) mat;
    int q;
    for(q=0; q<zbmat->N; q++) {
      m[q] = 1 / (1 + zbmat->xEp[q]/3 * (
				2 / (Eq - xVc[q] + zbmat->xEg[q]) +
				1 / (Eq - xVc[q] + zbmat->xEg[q] + zbmat->xESO[q])) );
    }
    return zbmat->N;
}


Band *WZband_new(numpyint N, const double *xEg,
		   const double *xEp, const double *xESO) {
    WZBand *wzband = (WZBand *) malloc( sizeof(WZBand) );
    wzband->updateM = WZupdateM;
    wzband->N = N;
    wzband->xEg = xEg;
    wzband->xEp = xEp;
    wzband->xESO = xESO;
    return (Band *) wzband;
}

void WZband_free(Band *wzband) {
    free( (WZBand *) wzband );
    return;
}

#ifdef _DEBUG
#include <stdio.h>
void WZband_check(const Band *band, numpyint N, const double *xEg,
		    const double *xEp, const double *xESO) {
    printf("Checking WZband\n");
    const WZBand *wzband = (const WZBand *) band;
    if(wzband->updateM != WZupdateM)
      printf("WZupdateM checkfail\n");
    if(wzband->N != N)
      printf("N checkfail\n");
    if(wzband->xEg != xEg)
      printf("xEg checkfail\n");
    if(wzband->xEp != xEp)
      printf("xEp checkfail\n");
    if(wzband->xESO != xESO)
      printf("xESO checkfail\n");
    printf("Checking Finished\n");
    return;
  }
#endif


#ifdef __cplusplus
}
#endif
