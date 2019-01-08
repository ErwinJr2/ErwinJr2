/**
 * \file
 *
 * Zincblende and Wurtzite structure band, compatible
 * with structure BAND.
 */


#include "band.h"

#ifdef __cplusplus
extern "C" {
#endif

/** \cond IMPL
 * implemtation of functions in header files should be excluded in doxygen */
numpyint UpdateBand(Band *band, double E, const double *xVc, double *m) {
	return band->update(band, E, xVc, m);
}
/** \endcond */

/**
 * structure for Zinc-blende band
 */
typedef struct ZBBAND {
        UpdateFunc updateM;    /**< Update effective mass */
        numpyint N;            /**< Number of finite x positions */
        const double *xEg;     /**< Direct energy gap  */
        const double *xF;      /**< Kane parameter  */
        const double *xEp;     /**< Matrix element */
        const double *xESO;    /**< Spin-orbit splitting */
}ZBBand; 

/** Update effective mass of a Zincblende band semiconductor */
numpyint ZBupdateM(Band *mat, double Eq, const double *xVc, double *m) {
	ZBBand *zbmat = (ZBBand *) mat;
	int q; 
	for(q=0; q<zbmat->N; q++) {
		double E = Eq - xVc[q];
		if(E < -zbmat->xEg[q])
			E = -zbmat->xEg[q] + 1e-7; /* Avoid singularity */
		m[q] = 1 / (1 + 2*zbmat->xF[q] + zbmat->xEp[q]/3 * ( 
					2 / (E + zbmat->xEg[q]) + 
					1 / (E + zbmat->xEg[q] + zbmat->xESO[q])) );
	}
	return zbmat->N;
}

/** \cond IMPL
 * implemtation of functions in header files should be excluded in doxygen */
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
/** \endcond */

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

/** \brief struct for Wurtzite band */
typedef struct WZBand {
  UpdateFunc updateM;   /**< Update effective mass */
  numpyint N;           /**< Number of finite x positions */
  const double *xEg;    /**< Direct energy gap */
  const double *xEp;    /**< Matrix parameter */
  const double *xESO;   /**< Spin-orbit splitting */
}WZBand;

/** \brief  Update effective mass of a Wurtzite band semiconductor. Output N */
numpyint WZupdateM(Band *mat, double Eq, const double *xVc, double *m) {
    WZBand *wzmat = (WZBand *) mat;
    int q;
    for(q=0; q<wzmat->N; q++) {
		double E = Eq - xVc[q];
		if(E < -wzmat->xEg[q])
			E = -wzmat->xEg[q] + 1e-7; /* Avoid singularity */
      m[q] = 1 / (1 + wzmat->xEp[q]/3 * (
				2 / (E + wzmat->xEg[q]) +
				1 / (E + wzmat->xEg[q] + wzmat->xESO[q])) );
    }
    return wzmat->N;
}

/** \cond IMPL 
 * implemtation of functions in header files should be excluded in doxygen */
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
/** \endcond */

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
