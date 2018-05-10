/**
 * @file blas.h
 * @brief      F77 BLAS header file
 * @author     Christian Kirches
 * @date       2009 Jun 1
 */

#pragma once
#ifndef _GLUE_BLAS_H_INCLUDED_
#define _GLUE_BLAS_H_INCLUDED_

#ifdef __cplusplus
extern "C" {
#endif

typedef long blas_int_t;   // 32-bit on x86_32, 64-bit on x86_64
// typedef int blas_int_t;   // always 32-bit even on x86_64

extern blas_int_t idamax_ (
    const blas_int_t   *n,
    const double *x,
    const blas_int_t   *ix
);

extern blas_int_t dcopy_ (
    const blas_int_t   *n,
    const double *x,
    const blas_int_t   *ix,
    double       *y,
    const blas_int_t   *iy
);

extern blas_int_t dswap_ (
    const blas_int_t   *n,
    double       *x,
    const blas_int_t   *ix,
    double       *y,
    const blas_int_t   *iy
);

extern blas_int_t daxpy_ (
    const blas_int_t   *n,
    const double *alpha,
    const double *x,
    const blas_int_t   *ix,
    double       *y,
    const blas_int_t   *iy
);

extern blas_int_t dscal_ (
    const blas_int_t   *n,
    const double *alpha,
    double       *x,
    const blas_int_t   *ix
);

extern double ddot_ (
    const blas_int_t   *n,
    const double *x,
    const blas_int_t   *ix,
    const double *y,
    const blas_int_t   *iy
);

extern double dnrm2_ (
    const blas_int_t   *n,
    const double *x,
    const blas_int_t   *ix
);

extern double dasum_ (
    const blas_int_t   *n,
    const double *x,
    const blas_int_t   *ix
);

extern void drotg_ (
    double *a,
    double *b,
    double *c,
    double *s
);

extern double drot_ (
    const blas_int_t   *n,
    double       *x,
    const blas_int_t   *ix,
    double       *y,
    const blas_int_t   *iy,
    const double *c,
    const double *s
);


/* Level 2 BLAS */

extern blas_int_t dgemv_ (
    const char   *trans,
    const blas_int_t   *m,
    const blas_int_t   *n,
    const double *alpha,
    const double *a,
    const blas_int_t   *lda,
    const double *x,
    const blas_int_t   *ix,
    const double *beta,
    double       *y,
    const blas_int_t   *iy,
    const blas_int_t    trans_len
);

extern blas_int_t dtrmv_ (
    const char   *uplo,
    const char   *trans,
    const char   *diag,
    const blas_int_t   *n,
    const double *a,
    const blas_int_t   *lda,
    const double *x,
    const blas_int_t   *incx,
    const blas_int_t    uplo_len,
    const blas_int_t    trans_len,
    const blas_int_t    diag_len
);

extern blas_int_t dsymv_ (
    const char   *uplo,
    const blas_int_t   *n,
    const double *alpha,
    const double *a,
    const blas_int_t   *lda,
    const double *x,
    const blas_int_t   *ix,
    const double *beta,
    double       *y,
    const blas_int_t   *iy,
    const blas_int_t    uplo_len
);

extern blas_int_t dtrsv_ (
    const char   *uplo,
    const char   *transa,
    const char   *diag,
    const blas_int_t   *n,
    const double *a,
    const blas_int_t   *lda,
    double       *x,
    const blas_int_t   *incx,
    const int    uplo_len,
    const int    trans_len,
    const int    diag_len
);

extern blas_int_t dger_ (
    const blas_int_t   *m,
    const blas_int_t   *n,
    const double *alpha,
    const double *x,
    const blas_int_t   *ix,
    const double *y,
    const blas_int_t   *iy,
    double       *a,
    const blas_int_t   *lda
);

extern blas_int_t dsyr_ (
    const char   *uplo,
    const blas_int_t   *n,
    const double *alpha,
    const double *x,
    const blas_int_t   *ix,
    double       *a,
    const blas_int_t   *lda,
    const blas_int_t    uplo_len
);


/* Level 3 BLAS */

extern blas_int_t dgemm_ (
    const char   *transa,
    const char   *transb,
    const blas_int_t   *m,
    const blas_int_t   *n,
    const blas_int_t   *k,
    const double *alpha,
    const double *a,
    const blas_int_t   *lda,
    const double *b,
    const blas_int_t   *ldb,
    const double *beta,
    const double *c,
    const blas_int_t   *ldc,
    const blas_int_t    transa_len,
    const blas_int_t    transb_len
);

extern blas_int_t dsymm_ (
    const char   *side,
    const char   *uplo,
    const blas_int_t   *m,
    const blas_int_t   *n,
    const double *alpha,
    const double *a,
    const blas_int_t   *lda,
    const double *b,
    const blas_int_t   *ldb,
    const double *beta,
    const double *c,
    const blas_int_t   *ldc,
    const blas_int_t    side_len,
    const blas_int_t    uplo_len
);

extern blas_int_t dsyrk_ (
    const char   *uplo,
    const char   *trans,
    const blas_int_t   *n,
    const blas_int_t   *k,
    const double *alpha,
    const double *a,
    const blas_int_t   *lda,
    const double *beta,
    const double *c,
    const blas_int_t   *ldc,
    const blas_int_t    uplo_len,
    const blas_int_t    trans_len
);

extern blas_int_t dtrsm_ (
    const char   *side,
    const char   *uplo,
    const char   *transa,
    const char   *diag,
    const blas_int_t   *m,
    const blas_int_t   *n,
    const double *alpha,
    const double *a,
    const blas_int_t   *lda,
    const double *b,
    const blas_int_t   *ldb,
    const blas_int_t    side_len,
    const blas_int_t    uplo_len,
    const blas_int_t    trans_len,
    const blas_int_t    diag_len
);

#ifdef __cplusplus
}
#endif

#endif // _GLUE_BLAS_H_INCLUDED_
