/*
 *  This file is part of fftautocorr
 *
 *  Copyright (C) 2020 CareF
 *  \author CareF 
 *  Licensed under a 3-clause BSD style license - see LICENSE.md
 */

#include <math.h>
#include <string.h>  /* provide memcpy */

#include "fftautocorr.h"
#include "factortable.h"
#define SQ(x) (x)*(x)
#define SWAP(a,b,type) \
  do { type tmp_=(a); (a)=(b); (b)=tmp_; } while(0)

#ifdef __GNUC__
#define NOINLINE __attribute__((noinline))
#define WARN_UNUSED_RESULT __attribute__ ((warn_unused_result))
#else
#define NOINLINE
#define WARN_UNUSED_RESULT
#endif

#ifdef _WINDOWS
/* Make compatible with visual studio's limited support for C99 */
#define restrict __restrict
#define M_PI 3.14159265358979323846
#endif

/* sin/cos related */

/* We avoid using cos here to decrease float truncation error. */
#define cosm1(x) (-2 * SQ(sin((x)/2))) /**< cos(x) - 1 */

NOINLINE static void calc_first_octant(int den, double * restrict res) {
    const int n = (den+4) >> 3, l1 = (int)sqrt(n); 
    int i, start;
    if(n == 0) 
        return;
    res[0] = 1.0;
    res[1] = 0.0;
    if(n == 1)
        return;
    for(i = 1; i < l1; ++i) {
        const double theta = M_PI * 2 * i / den;
        res[2*i] = cosm1(theta);
        res[2*i + 1] = sin(theta);
    }
    for(start = l1; start < n; start += l1) {
        const double theta = M_PI * 2.0 * start / den;
        const double c = cosm1(theta), s = sin(theta);
        int end;
        res[2*start] = c + 1;
        res[2*start + 1] = s;
        end = start + l1 > n ? n - start : l1;
        for(i = 1; i < end; ++i) {
            const double cx = res[2*i], sx = res[2*i+1];
            res[2*(start + i)] = c*cx - s*sx + c + cx + 1;
            res[2*(start + i) + 1] = c*sx + s*cx + s + sx;
        }
    }
    for(i = 1; i < l1; ++i) {
        res[2*i] += 1.0;
    }
}

NOINLINE static void calc_first_quadrant(int n, double * restrict res) {
    double * restrict p = res + n;
    const int ndone = (n + 2) >> 2;
    int i, idx1, idx2;
    calc_first_octant(n<<1, p);
    for(i =  0, idx1 =  0, idx2 =  2 * ndone - 2; i + 1 < ndone; 
        i += 2, idx1 += 2, idx2 -= 2) {
        res[idx1]   = p[2*i];
        res[idx1+1] = p[2*i+1];
        res[idx2]   = p[2*i+3];
        res[idx2+1] = p[2*i+2];
    }
    if(i != ndone) {
        res[idx1]   = p[2*i];
        res[idx1+1] = p[2*i+1];
    }
}

NOINLINE static void calc_first_half(int n, double * restrict res) {
    /* already know n is odd */
    double * p = res + n - 1;
    int i;
    calc_first_octant(n<<2, p);
    for(i = 0; 4*i <= n - 4*i; ++i) {
        /* octant 0 */
        int xm = 8*i;
        res[2*i]   = p[xm];
        res[2*i+1] = p[xm+1];
    }
    for(; 4*i - n <= 0; ++i) {
        /* octant 1 */
        int xm = 2*n - 8*i;
        res[2*i]   = p[xm+1];
        res[2*i+1] = p[xm];
    }
    for(; 4*i <= 3*n - 4*i; ++i) {
        /* octant 2 */
        int xm = 8*i - 2*n;
        res[2*i]   = -p[xm+1];
        res[2*i+1] = p[xm];
    }
    for(; 2 * i < n; ++i) {
        /* octant 3 */
        int xm = 4*n - 8*i;
        res[2*i] = -p[xm];
        res[2*i+1] = p[xm+1];
    }
}

#define hsqrt2 0.707106781186547524400844362104849  /**< sqrt(2)/2 */
NOINLINE static void fill_first_quadrant(int n, double * restrict res) {
    int quart = n >> 2, i, j;
    if(n % 8 == 0)
        res[quart] = res[quart+1] = hsqrt2;
    for(i=2, j=2*quart-2; i<quart; i+=2, j-=2) {
        res[j]   = res[i+1];
        res[j+1] = res[i];
    }
}

NOINLINE static void fill_first_half(int n, double * restrict res) {
    int half = n>>1;
    int i, j;
    if(n % 4 == 0)
        for(i=0; i<half; i+=2) {
            res[i+half]   = -res[i+1];
            res[i+half+1] =  res[i  ];
        }
    else
        for(i=2, j=2*half-2; i<half; i+=2, j-=2) {
        res[j]   = -res[i];
        res[j+1] =  res[i+1];
        }
}

/**
 * @brief calculate sin / cos FFT coefficient
 * put res[2*i] = cos(2 * M_PI * i / n ) and
 * res[2*i+1] = sin(2 * M_PI * i / n ) for i in range(0, nn)
 * when n % 4 == 0, nn = n/2; when 
 * 
 * @param[in] n the nubmer of grid on 2 pi
 * @param[out] res the address to save the result
 */
NOINLINE static void sincos_2pibyn_half(int n, double * restrict res) {
    if(n % 4 == 0) {
        calc_first_octant(n, res);
        fill_first_quadrant(n, res);
        fill_first_half(n, res);
    }
    else if(n % 2 == 0) {
        calc_first_quadrant(n, res);
        fill_first_half(n, res);
    }
    else
        calc_first_half(n, res);
}


/* rfft related functions */
#define NFCT 25
typedef struct fft_fctdata {
    int fct;
    double *tw;
} fft_fctdata;

/**
 * @brief The struct for a rfft plan
 */
typedef struct autocorr_plan_i {
    int memlen;            /**< the length for the FFT algorithm, 
                            *   >= 2*autocorr_plan_i#datalen */
    int nfct;              /**< number of seperation factors */
    size_t datalen;        /**< the logical length, length of the original
                            *   array for autocorrelation */
    double *mem;           /**< the memory pool for fct[i].tw */
    fft_fctdata fct[NFCT]; /**< the factors for each FFT seperation */
} autocorr_plan_i;

#define WA(x,i) wa[(i)+(x)*(ido-1)]
#define PM(a,b,c,d) { a=c+d; b=c-d; }
/* (a+ib) = conj(c+id) * (e+if) */
#define MULPM(a,b,c,d,e,f) { a=c*e+d*f; b=c*f-d*e; }

#define CC(a,b,c) cc[(a)+ido*((b)+l1*(c))]
#define CH(a,b,c) ch[(a)+ido*((b)+cdim*(c))]

NOINLINE static void radf2 (int ido, int l1, const double * restrict cc,
                            double * restrict ch, const double * restrict wa) {
    const int cdim=2;
    int k, i;
    for(k = 0; k < l1; k++)
        PM (CH(0,0,k),CH(ido-1,1,k),CC(0,k,0),CC(0,k,1))
    if(ido % 2 == 0)
        for(k = 0; k < l1; k++) {
            CH(    0,1,k) = -CC(ido-1,k,1);
            CH(ido-1,0,k) =  CC(ido-1,k,0);
        }
    if(ido<=2)
        return;
    for(k = 0; k < l1; k++)
        for(i = 2; i < ido; i += 2) {
            int ic = ido - i;
            double tr2, ti2;
            MULPM(tr2, ti2, WA(0, i-2), WA(0, i-1), CC(i-1, k, 1), CC(i, k, 1))
            PM(CH(i-1, 0, k), CH(ic-1, 1, k), CC(i - 1, k, 0), tr2)
            PM(CH(i  , 0, k), CH(ic  , 1, k), ti2, CC(i  , k, 0))
        }
}

#define hsqrt3 0.866025403784438646763723170752936  /**< sqrt(3)/2 */
NOINLINE static void radf3(int ido, int l1, const double * restrict cc,
                           double * restrict ch, const double * restrict wa) {
    const int cdim = 3;
    int k, i;
    for(k = 0; k < l1; k++) {
        double cr2 = CC(0, k, 1) + CC(0, k, 2);
        CH(0,     0, k) = CC(0, k, 0) + cr2;
        CH(0,     2, k) = hsqrt3 * (CC(0, k, 2) - CC(0, k, 1));
        CH(ido-1, 1, k) = CC(0, k, 0) - 0.5 * cr2;
    }
    if(ido == 1)
        return;
    for(k = 0; k < l1; k++)
        for(i = 2; i < ido; i += 2) {
            int ic = ido - i;
            double di2, di3, dr2, dr3;
            double tr2, tr3, ti2, ti3;
            double cr2, ci2;
            /* d2=conj(WA0)*CC1 */
            MULPM(dr2, di2, WA(0, i-2), WA(0, i-1), CC(i-1, k, 1), CC(i, k, 1));
            /* d3=conj(WA1)*CC2 */
            MULPM(dr3, di3, WA(1, i-2), WA(1, i-1), CC(i-1, k, 2), CC(i, k, 2));
            /* c add */
            cr2 = dr2 + dr3;
            ci2 = di2 + di3;
            /* c add */
            CH(i-1, 0, k) = CC(i-1, k, 0) + cr2;
            CH(i  , 0, k) = CC(i  , k, 0) + ci2;
            /* c add */
            tr2 = CC(i-1, k, 0) - 0.5 * cr2;
            ti2 = CC(i  , k, 0) - 0.5 * ci2;
            /* t3 = hsqrt3*i*(d3-d2)? */
            tr3 = hsqrt3*(di2-di3);
            ti3 = hsqrt3*(dr3-dr2);
            /* PM(i) = t2+t3 */
            PM(CH(i-1,2,k),CH(ic-1,1,k),tr2,tr3);
            /* PM(ic) = conj(t2-t3) */
            PM(CH(i  ,2,k),CH(ic  ,1,k),ti3,ti2);
        }
}

NOINLINE static void radf4(int ido, int l1, const double * restrict cc,
                           double * restrict ch, const double * restrict wa) {
    const int cdim=4;
    int k, i;
    for(k=0; k<l1; k++)
        {
        double tr1,tr2;
        PM (tr1,CH(0,2,k),CC(0,k,3),CC(0,k,1))
        PM (tr2,CH(ido-1,1,k),CC(0,k,0),CC(0,k,2))
        PM (CH(0,0,k),CH(ido-1,3,k),tr2,tr1)
        }
    if(ido % 2 == 0)
        for(k=0; k<l1; k++) {
            double ti1 = -hsqrt2 * (CC(ido-1,k,1)+CC(ido-1,k,3));
            double tr1= hsqrt2*(CC(ido-1,k,1)-CC(ido-1,k,3));
            PM (CH(ido-1,0,k),CH(ido-1,2,k),CC(ido-1,k,0),tr1);
            PM (CH(    0,3,k),CH(    0,1,k),ti1,CC(ido-1,k,2));
        }
    if(ido <= 2)
        return;
    for(k = 0; k < l1; k++)
        for(i = 2; i < ido; i+=2) {
            int ic = ido - i;
            double ci2, ci3, ci4, cr2, cr3, cr4;
            double ti1, ti2, ti3, ti4, tr1, tr2, tr3, tr4;
            MULPM(cr2, ci2, WA(0, i-2), WA(0, i-1), CC(i-1, k, 1), CC(i, k, 1));
            MULPM(cr3, ci3, WA(1, i-2), WA(1, i-1), CC(i-1, k, 2), CC(i, k, 2));
            MULPM(cr4, ci4, WA(2, i-2), WA(2, i-1), CC(i-1, k, 3), CC(i, k, 3));
            PM(tr1, tr4, cr4, cr2);
            PM(ti1, ti4, ci2, ci4);
            PM(tr2, tr3, CC(i-1,k, 0), cr3);
            PM(ti2, ti3, CC(i  ,k, 0), ci3);
            PM(CH(i-1, 0, k), CH(ic-1, 3, k), tr2, tr1);
            PM(CH(i  , 0, k), CH(ic  , 3, k), ti1, ti2);
            PM(CH(i-1, 2, k), CH(ic-1, 1, k), tr3, ti4);
            PM(CH(i  , 2, k), CH(ic  , 1, k), tr4, ti3);
        }
}

#define tr11  0.309016994374947424102293417182819  /**< sin(pi / 10) */
#define ti11  0.951056516295153572116439333379382  /**< sin(2 pi / 5) */
#define tr12 -0.809016994374947424102293417182819  /**< sin(3 pi / 10) */
#define ti12  0.587785252292473129168705954639073  /**< sin(pi / 5) */
NOINLINE static void radf5(int ido, int l1, const double * restrict cc,
                           double * restrict ch, const double * restrict wa) {
    const int cdim=5;
    int k, i;
    for(k = 0; k < l1; k++) {
        double cr2, cr3, ci4, ci5;
        PM(cr2, ci5, CC(0, k, 4), CC(0, k, 1));
        PM(cr3, ci4, CC(0, k, 3), CC(0, k, 2));
        CH(0, 0, k)     = CC(0, k, 0) + cr2      + cr3;
        CH(ido-1, 1, k) = CC(0, k, 0) + tr11*cr2 + tr12*cr3;
        CH(0, 2, k)     = ti11*ci5 + ti12*ci4;
        CH(ido-1, 3, k) = CC(0, k, 0) + tr12*cr2 + tr11*cr3;
        CH(0, 4, k)     = ti12*ci5 - ti11*ci4;
    }
    if(ido == 1) 
        return;
    for(k = 0; k < l1; ++k)
        for(i = 2; i < ido; i += 2) {
            double ci2, ci3, ci4, ci5, di2, di3, di4, di5, ti2, ti3, ti4, ti5;
            double cr2, cr3, cr4, cr5, dr2, dr3, dr4, dr5, tr2, tr3, tr4, tr5;
            int ic = ido - i;
            MULPM(dr2, di2, WA(0, i-2), WA(0, i-1), CC(i-1, k, 1), CC(i, k, 1));
            MULPM(dr3, di3, WA(1, i-2), WA(1, i-1), CC(i-1, k, 2), CC(i, k, 2));
            MULPM(dr4, di4, WA(2, i-2), WA(2, i-1), CC(i-1, k, 3), CC(i, k, 3));
            MULPM(dr5, di5, WA(3, i-2), WA(3, i-1), CC(i-1, k, 4), CC(i, k, 4));
            PM(cr2, ci5, dr5, dr2);
            PM(ci2, cr5, di2, di5);
            PM(cr3, ci4, dr4, dr3);
            PM(ci3, cr4, di3, di4);
            CH(i-1, 0, k) = CC(i-1, k, 0) + cr2 + cr3;
            CH(i  , 0, k) = CC(i  , k, 0) + ci2 + ci3;
            tr2 = CC(i-1, k, 0) + tr11*cr2 + tr12*cr3;
            ti2 = CC(i  , k, 0) + tr11*ci2 + tr12*ci3;
            tr3 = CC(i-1, k, 0) + tr12*cr2 + tr11*cr3;
            ti3 = CC(i  , k, 0) + tr12*ci2 + tr11*ci3;
            MULPM(tr5, tr4, cr5, cr4, ti11, ti12);
            MULPM(ti5, ti4, ci5, ci4, ti11, ti12);
            PM(CH(i-1, 2, k), CH(ic-1, 1, k), tr2, tr5);
            PM(CH(i  , 2, k), CH(ic  , 1, k), ti5, ti2);
            PM(CH(i-1, 4, k), CH(ic-1, 3, k), tr3, tr4);
            PM(CH(i  , 4, k), CH(ic  , 3, k), ti4, ti3);
        }
}

#undef CC
#undef CH
#define CH(a,b,c) ch[(a)+ido*((b)+l1*(c))]
#define CC(a,b,c) cc[(a)+ido*((b)+cdim*(c))]

NOINLINE static void radb2(int ido, int l1, const double * restrict cc,
                           double * restrict ch, const double * restrict wa) {
    const int cdim=2;
    int k, i;
    for(k = 0; k < l1; k++)
        PM(CH(0, k, 0), CH(0, k, 1), CC(0, 0, k), CC(ido-1, 1, k));
    if(ido % 2 == 0)
        for(k = 0; k < l1; k++) {
            CH(ido-1, k, 0) =  2. * CC(ido-1, 0, k);
            CH(ido-1, k, 1) = -2. * CC(0    , 1, k);
        }
    if(ido <= 2)
        return;
    for(k = 0; k < l1; ++k)
        for(i = 2; i < ido; i += 2) {
            int ic = ido - i;
            double ti2, tr2;
            PM(CH(i-1, k, 0), tr2,         CC(i-1, 0, k), CC(ic-1, 1, k));
            PM(ti2,           CH(i, k, 0), CC(i,   0, k), CC(ic,   1, k));
            MULPM(CH(i, k, 1), CH(i-1, k, 1), WA(0, i-2), WA(0, i-1), ti2, tr2);
        }
}

NOINLINE static void radb3(int ido, int l1, const double * restrict cc,
                           double * restrict ch, const double * restrict wa) {
    const int cdim=3;
    int k, i;
    for(k = 0; k < l1; k++) {
        double tr2  = 2. * CC(ido-1, 1, k);
        double cr2  = CC(0, 0, k) - 0.5 * tr2;
        double ci3;
        CH(0, k, 0) = CC(0, 0, k) + tr2;
        ci3 = 2. * hsqrt3 * CC(0, 2, k);
        PM(CH(0, k, 2), CH(0, k, 1), cr2, ci3);
    }
    if(ido == 1)
        return;
    for(k = 0; k < l1; k++)
        for(i = 2; i < ido; i += 2) {
            int ic = ido - i;
            /* t2 = CC(I) + conj(CC(ic)) */
            double tr2 = CC(i-1, 2, k) + CC(ic-1, 1, k);
            double ti2 = CC(i  , 2, k) - CC(ic  , 1, k);
            /* c2 = CC + taur * t2 */
            double cr2 = CC(i-1, 0, k) - 0.5 * tr2;
            double ci2 = CC(i  , 0, k) - 0.5 * ti2;
            double cr3, ci3, di2, di3, dr2, dr3;
            /* CH = CC + t2 */
            CH(i-1, k, 0) = CC(i-1, 0, k) + tr2;
            CH(i  , k, 0) = CC(i  , 0, k) + ti2;
            /* c3 = hsqrt3 * (CC(i) - conj(CC(ic))) */
            cr3 = hsqrt3 * (CC(i-1, 2, k) - CC(ic-1, 1, k));
            ci3 = hsqrt3 * (CC(i  , 2, k) + CC(ic  , 1, k));
            /* d2 = (cr2 - ci3, ci2 + cr3) = c2 + i*c3 */
            PM(dr3, dr2, cr2, ci3);
            /* d3 = (cr2+ci3, ci2-cr3) = c2-i*c3 */
            PM(di2, di3, ci2, cr3);
            /* ch = WA*d2 */
            MULPM(CH(i, k, 1), CH(i-1, k, 1), WA(0, i-2), WA(0, i-1), di2, dr2);
            MULPM(CH(i, k, 2), CH(i-1, k, 2), WA(1, i-2), WA(1, i-1), di3, dr3);
        }
}

#define sqrt2 1.414213562373095048801688724209698
NOINLINE static void radb4(int ido, int l1, const double * restrict cc,
                           double * restrict ch, const double * restrict wa) {
    const int cdim=4;
    int k, i;
    for(k=0; k<l1; k++) {
        double tr1, tr2, tr3, tr4;
        PM(tr2, tr1, CC(0, 0, k), CC(ido-1, 3, k));
        tr3 = 2. * CC(ido-1, 1, k);
        tr4 = 2. * CC(0,     2, k);
        PM (CH(0,k,0),CH(0,k,2),tr2,tr3)
        PM (CH(0,k,3),CH(0,k,1),tr1,tr4)
    }
    if(ido % 2 == 0)
        for(int k=0; k<l1; k++) {
            double tr1, tr2, ti1, ti2;
            PM(ti1, ti2, CC(0    , 3, k), CC(0    , 1, k));
            PM(tr2, tr1, CC(ido-1, 0, k), CC(ido-1, 2, k));
            CH(ido-1, k, 0) = tr2 + tr2;
            CH(ido-1, k, 1) = sqrt2 * (tr1 - ti1);
            CH(ido-1, k, 2) = ti2 + ti2;
            CH(ido-1, k, 3) = -sqrt2 * (tr1 + ti1);
        }
    if(ido <= 2)
        return;
    for(k = 0; k < l1; ++k)
        for(i = 2; i < ido; i += 2) {
            double ci2, ci3, ci4, cr2, cr3, cr4;
            double ti1, ti2, ti3, ti4, tr1, tr2, tr3, tr4;
            int ic = ido - i;
            PM(tr2, tr1, CC(i-1, 0, k), CC(ic-1, 3, k));
            PM(ti1, ti2, CC(i  , 0, k), CC(ic  , 3, k));
            PM(tr4, ti3, CC(i  , 2, k), CC(ic  , 1, k));
            PM(tr3, ti4, CC(i-1, 2, k), CC(ic-1, 1, k));
            PM(CH(i-1, k, 0), cr3, tr2, tr3);
            PM(CH(i  , k, 0), ci3, ti2, ti3);
            PM(cr4, cr2, tr1, tr4);
            PM(ci2, ci4, ti1, ti4);
            MULPM(CH(i, k, 1), CH(i-1, k, 1), WA(0, i-2), WA(0, i-1), ci2, cr2);
            MULPM(CH(i, k, 2), CH(i-1, k, 2), WA(1, i-2), WA(1, i-1), ci3, cr3);
            MULPM(CH(i, k, 3), CH(i-1, k, 3), WA(2, i-2), WA(2, i-1), ci4, cr4);
        }
}

NOINLINE static void radb5(int ido, int l1, const double * restrict cc,
                           double * restrict ch, const double * restrict wa) {
    const int cdim=5;
    int k, i;
    for(k = 0; k < l1; k++) {
        double ci4, ci5, cr2, cr3;
        double ti5 = CC(0, 2, k) + CC(0, 2, k);
        double ti4 = CC(0, 4, k) + CC(0, 4, k);
        double tr2 = CC(ido-1, 1, k) + CC(ido-1, 1, k);
        double tr3 = CC(ido-1, 3, k) + CC(ido-1, 3, k);
        CH(0, k, 0) = CC(0, 0, k) + tr2 + tr3;
        cr2 = CC(0, 0, k) + tr11*tr2 + tr12*tr3;
        cr3 = CC(0, 0, k) + tr12*tr2 + tr11*tr3;
        MULPM(ci5, ci4, ti5, ti4, ti11, ti12);
        PM(CH(0, k, 4), CH(0, k, 1), cr2, ci5);
        PM(CH(0, k, 3), CH(0, k, 2), cr3, ci4);
    }
    if(ido == 1)
        return;
    for(k = 0; k < l1; ++k)
        for(i = 2; i < ido; i += 2) {
            int ic = ido-i;
            double tr2, tr3, tr4, tr5, ti2, ti3, ti4, ti5;
            double cr2, ci2, cr3, ci3;
            double ci4, ci5, cr5, cr4;
            double dr2, dr3, dr4, dr5, di2, di3, di4, di5;
            PM(tr2, tr5, CC(i-1, 2, k), CC(ic-1, 1, k));
            PM(ti5, ti2, CC(i  , 2, k), CC(ic  , 1, k));
            PM(tr3, tr4, CC(i-1, 4, k), CC(ic-1, 3, k));
            PM(ti4, ti3, CC(i  , 4, k), CC(ic  , 3, k));
            CH(i-1,k,0)=CC(i-1,0,k)+tr2+tr3;
            CH(i  ,k,0)=CC(i  ,0,k)+ti2+ti3;
            cr2 = CC(i-1, 0, k) + tr11*tr2 + tr12*tr3;
            ci2 = CC(i  , 0, k) + tr11*ti2 + tr12*ti3;
            cr3 = CC(i-1, 0, k) + tr12*tr2 + tr11*tr3;
            ci3 = CC(i  , 0, k) + tr12*ti2 + tr11*ti3;
            MULPM(cr5,cr4,tr5,tr4,ti11,ti12)
            MULPM(ci5,ci4,ti5,ti4,ti11,ti12)
            PM(dr4, dr3, cr3, ci4);
            PM(di3, di4, ci3, cr4);
            PM(dr5, dr2, cr2, ci5);
            PM(di2, di5, ci2, cr5);
            MULPM(CH(i, k, 1), CH(i-1, k, 1), WA(0, i-2), WA(0, i-1), di2, dr2);
            MULPM(CH(i, k, 2), CH(i-1, k, 2), WA(1, i-2), WA(1, i-1), di3, dr3);
            MULPM(CH(i, k, 3), CH(i-1, k, 3), WA(2, i-2), WA(2, i-1), di4, dr4);
            MULPM(CH(i, k, 4), CH(i-1, k, 4), WA(3, i-2), WA(3, i-1), di5, dr5);
        }
}

#undef CC
#undef CH
#undef PM
#undef MULPM
#undef WA

static void copy_and_norm(double *c, double *p1, int n, double fct) {
    if(p1 != c) {
        if(fct != 1.) {
            for(int i=0; i<n; ++i)
                c[i] = fct*p1[i];
        }
        else {
            memcpy (c,p1,n*sizeof(double));
        }
    }
    else if(fct != 1.) {
        for(int i=0; i<n; ++i)
            c[i] *= fct;
    }
}

/**
 * @brief Calculate forward rFFT
 * 
 * @param plan 
 * @param c 
 * @param fct 
 * @param mem 
 * @return WARN_UNUSED_RESULT 
 */
WARN_UNUSED_RESULT
static int rfftp_forward(autocorr_plan plan, double c[], 
                         double fct, double *mem) {
    if(plan->memlen == 1)
        return 0;
    int n = plan->memlen;
    int l1 = n, nf = plan->nfct;
    double *p1 = c, *p2 = mem;
    int k1;

    for(k1 = 0; k1 < nf; ++k1) {
        int k = nf-k1-1;
        int ip = plan->fct[k].fct;
        int ido = n / l1;
        l1 /= ip;
        switch(ip) {
            case 4:
                radf4(ido, l1, p1, p2, plan->fct[k].tw);
                break;
            case 2:
                radf2(ido, l1, p1, p2, plan->fct[k].tw);
                break;
            case 3:
                radf3(ido, l1, p1, p2, plan->fct[k].tw);
                break;
            case 5:
                radf5(ido, l1, p1, p2, plan->fct[k].tw);
                break;
            default:
                return -1;
        }
        SWAP (p1, p2, double *);
    }
    copy_and_norm(c,p1,n,fct);
    return 0;
}

/**
 * @brief Calculate backward rFFT
 * 
 * @param plan 
 * @param c 
 * @param fct 
 * @param mem 
 * @return WARN_UNUSED_RESULT 
 */
WARN_UNUSED_RESULT
static int rfftp_backward(autocorr_plan plan, double c[], 
                          double fct, double *mem) {
    if(plan->memlen == 1)
        return 0;
    int n = plan->memlen;
    int l1 = 1, nf = plan->nfct;
    double *p1 = c, *p2 = mem;

    for(int k = 0; k < nf; k++) {
        int ip = plan->fct[k].fct, ido = n / (ip*l1);
        switch(ip) {
            case 4:
                radb4(ido, l1, p1, p2, plan->fct[k].tw);
                break;
            case 2:
                radb2(ido, l1, p1, p2, plan->fct[k].tw);
                break;
            case 3:
                radb3(ido, l1, p1, p2, plan->fct[k].tw);
                break;
            case 5:
                radb5(ido, l1, p1, p2, plan->fct[k].tw);
                break;
            default:
                return -1;
        }
        SWAP(p1, p2, double *);
        l1 *= ip;
    }
    copy_and_norm(c,p1,n,fct);
    return 0;
}

#ifdef NO_TABLE

/**
 * @brief Factorize for plan
 * so that prod(plan->fct[:plan->nfct]) = plan->memlen >= 2 * plan->datalen
 * 
 * @param[in] plan 
 * @return int 0 if success, -1 if failed
 */
WARN_UNUSED_RESULT
static int rfftp_factorize (autocorr_plan plan) { 
    int length=1;
    int nfct=0;
    while(length < 2 * plan->datalen) {
        if(nfct >= NFCT) 
            return -1; 
        plan->fct[nfct++].fct=4;
        length <<= 2;
    }
    if(length >= 4 * plan->datalen) {
        plan->fct[0].fct = 2;
        length >>= 1;
    }
    int n3 = nfct-1;
    while(n3 >= 1 && length >= 8.0 * plan->datalen / 3) {
        plan->fct[n3--].fct = 3;
        length = length / 4 * 3;
    }
    plan->memlen = length;
    plan->nfct=nfct;
    return 0;
}

#else

/**
 * @brief Factorize for plan using a built-in table
 * so that prod(plan->fct[:plan->nfct]) = plan->memlen >= 2 * plan->datalen
 * and plan->memlen is the minimum of such that is a composite of 2, 3, 4, 5
 * 
 * @param[in] plan 
 * @return 0 if success, -1 if failed
 */
WARN_UNUSED_RESULT
static int rfftp_factorize (autocorr_plan plan) { 
    int length = find_factor(plan->datalen*2);
    int nfct = 0;
    if(length == -1)
        return -1;
    plan->memlen = length;
    while (length % 4 == 0) {
        plan->fct[nfct++].fct = 4;
        length >>= 2;
    }
    if (length % 2 == 0) {
        plan->fct[nfct++].fct = plan->fct[0].fct;
        plan->fct[0].fct = 2;
        length >>= 1;
    }
    while (length % 3 == 0) {
        plan->fct[nfct++].fct = 3;
        length /= 3;
    }
    while (length % 5 == 0) {
        plan->fct[nfct++].fct = 5;
        length /= 5;
    }
    if (length != 1) 
        return -1;
    plan->nfct=nfct;
    return 0;
}

#endif

/**
 * @brief The total length of the memory needed to store the twiddle factors
 *        for the real FFT calculation.
 * 
 * @param[in] plan The plan to compute the twiddle factors for.
 * @return Number of double float twiddle factors for the plan
 */
static size_t rfftp_twsize(autocorr_plan plan) {
    int twsize = 0, l1 = 1;
    for(int k = 0; k < plan->nfct - 1; ++k) {
        int ip = plan->fct[k].fct, ido = plan->memlen / (l1*ip);
        twsize += (ip - 1) * (ido - 1);
        l1 *= ip;
    }
    return twsize;
}

/**
 * @brief Compute the twiddle factors for a plan
 * 
 * @param[in,out] plan The plan to compute the twiddle factors for
 * @return 0 for success, -1 for fail.
 */
WARN_UNUSED_RESULT NOINLINE static int rfftp_comp_twiddle (autocorr_plan plan) {
    int length=plan->memlen;
    double *twid = (double *)malloc(2*length * sizeof(double)), *ptr;
    int l1=1, k;
    if(twid == NULL) 
        return -1;
    sincos_2pibyn_half(length, twid);
    ptr = plan->mem;
    for(k = 0; k < plan->nfct; ++k) {
        int fct = plan->fct[k].fct;
        int ido = length / (l1 * fct);
        if(k < plan->nfct - 1) {
            /* last factor doesn't need twiddles */
            int i, j;
            plan->fct[k].tw = ptr;
            ptr += (fct - 1) * (ido - 1);
            for(j = 0; j < fct-1; ++j)
            for(i = 0; i <= (ido - 1) / 2 - 1; ++i) {
                int fctidx = j * (ido - 1) + 2*i;
                int twididx = 2 * (j+1) * l1 * (i+1);
                plan->fct[k].tw[fctidx] = twid[twididx];
                plan->fct[k].tw[fctidx + 1] = twid[twididx + 1];
            }
        }
        l1 *= fct;
    }
    free(twid);
    return 0;
}


/* AUTOCORR IMPLEMENTATION */

autocorr_plan make_autocorr_plan(size_t length) {
    autocorr_plan plan;
    if(length==0) 
        return NULL;
    plan = (autocorr_plan)malloc(sizeof(autocorr_plan_i));
    if(!plan) 
        return NULL;
    plan->datalen = length;
    plan->nfct=0;
    plan->mem=NULL;
    for (int i=0; i<NFCT; ++i)
        plan->fct[i]=(fft_fctdata){0,0};
    if(rfftp_factorize(plan)!=0) {
        free(plan); 
        return NULL; 
    }
    size_t tws = rfftp_twsize(plan);
    plan->mem = (double *)malloc(tws * sizeof(double));
    if(!plan->mem){
        free(plan); 
        return NULL;
    }
    if(rfftp_comp_twiddle(plan) != 0){
        free(plan->mem);
        free(plan);
        return NULL;
    }
    return plan;
}

void destroy_autocorr_plan(autocorr_plan plan) {
    free(plan->mem);
    plan->mem = NULL;
    free(plan);
    plan = NULL;
}

size_t mem_len(autocorr_plan plan) {
    return plan->memlen;
}

size_t data_len(autocorr_plan plan) {
    return plan->datalen;
}

int autocorr_mem(autocorr_plan plan, double data[], double *mempool) {
    if(rfftp_forward(plan, data, 1, mempool) != 0)
        return -1;
    data[0] = SQ(data[0]);
    for (int i = 1; 2 * i < mem_len(plan); i++) {
        data[2*i - 1] = SQ(data[2*i - 1]) + SQ(data[2*i]);
        data[2*i] = 0;
    }
    if(mem_len(plan) % 2 == 0)
        data[mem_len(plan) - 1] = SQ(data[mem_len(plan) - 1]);
    if(rfftp_backward(plan, data, 1.0/mem_len(plan), mempool) != 0)
        return -1;
    return 0;
}

int autocorr_p(autocorr_plan plan, double data[]) {
    int result;
    double *mempool = (double *) malloc(mem_len(plan) * sizeof(double));
    if(!mempool) {
        return -1;
    }
    result = autocorr_mem(plan, data, mempool);
    free(mempool);
    return result;
}

int autocorr(double data[], size_t length) {
    autocorr_plan plan = make_autocorr_plan(length);
    double *fftauto = malloc(mem_len(plan) * sizeof(double));
    int i;
    for(i = 0; i < length; i++) {
        fftauto[i] = data[i];
    }
    for(; i < mem_len(plan); i++) {
        fftauto[i] = 0;
    }
    if (autocorr_p(plan, fftauto) != 0) {
        destroy_autocorr_plan(plan);
        return -1;
    }
    for(i = 0; i < length; i++) {
        data[i] = fftauto[i];
    }
    destroy_autocorr_plan(plan);
    return 0;
}