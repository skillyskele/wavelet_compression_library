#ifndef DWT_H
#define DWT_H

#include "wavedec.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Discrete Wavelet Transform (symmetric extension) stride-based API.
    Matches the signature from dwt.c */
void dwt_sym_stride(double *inp, int N, double *lpd, double *hpd, int lpd_len,
                          double *cA, int len_cA, double *cD, int istride, int ostride);

void dwt(wt_object wt, const double *input);

#ifdef __cplusplus
}
#endif

#endif /* DWT_H */
