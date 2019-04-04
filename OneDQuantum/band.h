/**
 * \file
 *
 * \brief  header file for band.c
 *
 */


#ifndef BAND_H
#define BAND_H
#include <stdlib.h>
#include "science.h"
typedef struct BAND Band;

/**
 * An update parameter function pointer that accepts a pointer to 
 * to band struct (with its parameters as members) and according to 
 * the Band parameters and energy, update para (usually effective mass). 
 * The latter two double * parameters are potential (including band offset) 
 * and effective mass respectively.
 */
typedef numpyint (*UpdateFunc)(Band *, double, const double *, double *); 

/** \brief Base class for band structure */
typedef struct BAND{
    const UpdateFunc update;  /**< Update band structure */
    numpyint N;               /**< Size of datas */
    double *Eg;               /**< Band gap in eV */
} Band; 

#ifdef _WINDLL
__declspec(dllexport)
#endif 
/** \brief Update effective mass in band */
numpyint UpdateBand(Band *band, double E, const double *xVc, double *m);

#ifdef _WINDLL
__declspec(dllexport)
#endif 
/** \brief Initiate new Zinc-blende semiconductor */
Band *ZBband_new(numpyint N, const double *xEg, const double *xF,
        const double *xEp, const double *xESO);

#ifdef _WINDLL
__declspec(dllexport)
#endif 
/** \brief Free a Zinc-blende semiconductor */
void ZBband_free(Band *);

#ifdef _WINDLL
__declspec(dllexport)
#endif
/** \brief Initiate new Wurtzite semiconductor */
Band *WZband_new(numpyint N, const double *xEg, 
        const double *xEp, const double *xESO);

#ifdef _WINDLL
__declspec(dllexport)
#endif
/** \brief Free a Wurtzite semiconductor */
void WZband_free(Band *);

#endif /* ifndef BAND_H */
