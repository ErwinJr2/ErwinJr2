/*
 * Intigration test for fftautocorr, compared with FFTW and pocketfft
 * 
 *  Copyright (C) 2020 CareF
 *  \author CareF 
 *  Licensed under a 3-clause BSD style license - see LICENSE.md
 */

#include "fftautocorr.h"
#include "pocketfft/pocketfft.h"
#include <stdio.h>
#include <fftw3.h>
#include <time.h>
#include <math.h>
#include <limits.h>

#define ERR 1E-8
#define SQ(x) (x)*(x)

int testAgainst(int N) {
    fftw_plan pr2c, pc2r;
    double *data = malloc(N * sizeof(double));
    double error;
    for(int i=0; i<N; i++) {
        data[i] = (double)rand() / (double)(RAND_MAX);
        /* printf("%f\n", in[i]); */
    }
    clock_t t0, t;
    printf("Testing auto-correlation for length %d:\n", N);

    t0 = clock();
    /* Calculate auto-correlation using FFTW */
    double *infftw = fftw_malloc(2 * N * sizeof(double));
    fftw_complex *freq = fftw_malloc((N+1) * sizeof(fftw_complex));
    pr2c = fftw_plan_dft_r2c_1d(2*N, infftw, freq, FFTW_ESTIMATE);
    pc2r = fftw_plan_dft_c2r_1d(2*N, freq, infftw, FFTW_ESTIMATE);
    t = clock();
    printf("FFTW planning time: %g s\n", (t-t0)/(double)CLOCKS_PER_SEC);
    t0 = clock();
    for(int i=0; i<N; i++) {
        infftw[i] = data[i];
    }
    for(int i=N; i<2*N; i++) {
        infftw[i] = 0;
    }
    fftw_execute(pr2c);
    for(int i=0; i<N+1; i++) {
        freq[i][0] = SQ(freq[i][1]) + SQ(freq[i][0]);
        freq[i][1] = 0;
    }
    fftw_execute(pc2r);
    for(int i=0; i<N; i++) {
        infftw[i] /= (2*N);
    }
    t = clock();
    printf("FFTW execution time: %g s\n", (t-t0)/(double)CLOCKS_PER_SEC);

    t0 = clock();
    /* Calculate auto-correlation using pocketfft */
    rfft_plan pocketp = make_rfft_plan(2 * N);
    double *pocketc = malloc(2 * N * sizeof(double));
    t = clock();
    printf("pocketfft planning time: %g s\n", (t-t0)/(double)CLOCKS_PER_SEC);
    t0 = clock();
    for(int i=0; i<N; i++) {
        pocketc[i] = data[i];
    }
    for(int i=N; i<2 * N; i++) {
        pocketc[i] = 0;
    }
    if(rfft_forward(pocketp, pocketc, 1) != 0) {
        printf("pocketfft executing failed.\n");
        return -1;
    }
    pocketc[0] = SQ(pocketc[0]);
    for(int i=1; i<N; i++) {
        pocketc[2*i - 1] = SQ(pocketc[2*i - 1]) + SQ(pocketc[2*i]);
        pocketc[2*i] = 0;
    }
    pocketc[2*N-1] = SQ(pocketc[2*N-1]);
    if(rfft_backward(pocketp, pocketc, 0.5/N) != 0) {
        printf("pocketfft executing failed.\n");
        return -1;
    }
    t = clock();
    printf("pocketfft execution time: %g s\n", (t-t0)/(double)CLOCKS_PER_SEC);

    t0 = clock();
    /* Calculate auto-correlation using fftautocorr */
    autocorr_plan fftauto_plan = make_autocorr_plan(N);
    double *fftauto = malloc(mem_len(fftauto_plan) * sizeof(double));
    t = clock();
    printf("autocorr planning time: %g s\n", (t-t0)/(double)CLOCKS_PER_SEC);
    t0 = clock();
    for(int i=0; i<N; i++) {
        fftauto[i] = data[i];
    }
    for(int i=N; i<mem_len(fftauto_plan); i++) {
        fftauto[i] = 0;
    }
    if(autocorr_p(fftauto_plan, fftauto) != 0) {
        printf("Autocorr executing failed.\n");
        return -1;
    }
    t = clock();
    printf("autocorr execution time: %g s\n", (t-t0)/(double)CLOCKS_PER_SEC);

    /* Estimate error */
    error = 0;
    for(int i=0; i < N; i++) {
        if(fabs(infftw[i] - fftauto[i]) > ERR) {
            printf("Test failed %d, diff = %g\n", i,
                   infftw[i] - fftauto[i]);
            return -1;
        }
        error += SQ(infftw[i] - fftauto[i]);
    }
    printf("Succeed. total avg err = %g\n", sqrt(error) / N);

    fftw_destroy_plan(pc2r);
    fftw_destroy_plan(pr2c);
    fftw_free(infftw);
    fftw_free(freq);

    destroy_rfft_plan(pocketp);
    free(pocketc);

    destroy_autocorr_plan(fftauto_plan);
    free(fftauto);
    free(data);
    return 0;
}

int main(int argc, char *argv[]) {
    int Ls[] = {
        7,
        1<<22, 
        14348907, /* = 3^15 */
        9765625,  /* = 5^10 */
        5764801,  /* = 7^8 */
        4561727,  /* prime */
        295245/2, /* padded to odd 295245 */
        // INT_MAX/4, 
    };
    size_t n = sizeof(Ls)/sizeof(Ls[0]);
    for(int i = 0; i < n; ++i) {
        printf("Testing for L = %d\n", Ls[i]);
        if(testAgainst(Ls[i]) != 0) {
            return -1;
        }
        putchar('\n');
    }
    return 0;
}
