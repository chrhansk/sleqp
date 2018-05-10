/** 
 * @file cblas.h
 * @brief      C BLAS header file
 * @author     ATLAS developers, based on cblas.h from http://sourceforge.net/projects/math-atlas/
 * @date       2009 Jun 1
 */

#pragma once
#ifndef _GLUE_CBLAS_H_INCLUDED_
#define _GLUE_CBLAS_H_INCLUDED_

#ifdef __cplusplus
extern "C" {
#endif

typedef long cblas_int_t;   // 32-bit on x86_32, 64-bit on x86_64
// typedef int cblas_int_t;   // always 32-bit even on x86_64

enum CBLAS_ORDER {
    CblasRowMajor = 101, 
    CblasColMajor = 102
};

enum CBLAS_TRANSPOSE {
    CblasNoTrans   = 111, 
    CblasTrans     = 112, 
    CblasConjTrans = 113,
    AtlasConj      = 114
};

enum CBLAS_UPLO {
    CblasUpper = 121, 
    CblasLower = 122
};

enum CBLAS_DIAG {
    CblasNonUnit = 131, 
    CblasUnit    = 132
};

enum CBLAS_SIDE {
    CblasLeft  = 141, 
    CblasRight = 142
};

/*
 * ===========================================================================
 * Prototypes for level 1 BLAS functions (complex are recast as routines)
 * ===========================================================================
 */
float  cblas_sdsdot(const cblas_int_t N, const float alpha, const float *X,
                    const cblas_int_t incX, const float *Y, const cblas_int_t incY);
double cblas_dsdot(const cblas_int_t N, const float *X, const cblas_int_t incX, const float *Y,
                   const cblas_int_t incY);
float  cblas_sdot(const cblas_int_t N, const float  *X, const cblas_int_t incX,
                  const float  *Y, const cblas_int_t incY);
double cblas_ddot(const cblas_int_t N, const double *X, const cblas_int_t incX,
                  const double *Y, const cblas_int_t incY);
/*
 * Functions having prefixes Z and C only
 */
void   cblas_cdotu_sub(const cblas_int_t N, const void *X, const cblas_int_t incX,
                       const void *Y, const cblas_int_t incY, void *dotu);
void   cblas_cdotc_sub(const cblas_int_t N, const void *X, const cblas_int_t incX,
                       const void *Y, const cblas_int_t incY, void *dotc);

void   cblas_zdotu_sub(const cblas_int_t N, const void *X, const cblas_int_t incX,
                       const void *Y, const cblas_int_t incY, void *dotu);
void   cblas_zdotc_sub(const cblas_int_t N, const void *X, const cblas_int_t incX,
                       const void *Y, const cblas_int_t incY, void *dotc);


/*
 * Functions having prefixes S D SC DZ
 */
float  cblas_snrm2(const cblas_int_t N, const float *X, const cblas_int_t incX);
float  cblas_sasum(const cblas_int_t N, const float *X, const cblas_int_t incX);

double cblas_dnrm2(const cblas_int_t N, const double *X, const cblas_int_t incX);
double cblas_dasum(const cblas_int_t N, const double *X, const cblas_int_t incX);

float  cblas_scnrm2(const cblas_int_t N, const void *X, const cblas_int_t incX);
float  cblas_scasum(const cblas_int_t N, const void *X, const cblas_int_t incX);

double cblas_dznrm2(const cblas_int_t N, const void *X, const cblas_int_t incX);
double cblas_dzasum(const cblas_int_t N, const void *X, const cblas_int_t incX);


/*
 * Functions having standard 4 prefixes (S D C Z)
 */
cblas_int_t cblas_isamax(const cblas_int_t N, const float  *X, const cblas_int_t incX);
cblas_int_t cblas_idamax(const cblas_int_t N, const double *X, const cblas_int_t incX);
cblas_int_t cblas_icamax(const cblas_int_t N, const void   *X, const cblas_int_t incX);
cblas_int_t cblas_izamax(const cblas_int_t N, const void   *X, const cblas_int_t incX);

/*
 * ===========================================================================
 * Prototypes for level 1 BLAS routines
 * ===========================================================================
 */

/*
 * Routines with standard 4 prefixes (s, d, c, z)
 */
void cblas_sswap(const cblas_int_t N, float *X, const cblas_int_t incX,
                 float *Y, const cblas_int_t incY);
void cblas_scopy(const cblas_int_t N, const float *X, const cblas_int_t incX,
                 float *Y, const cblas_int_t incY);
void cblas_saxpy(const cblas_int_t N, const float alpha, const float *X,
                 const cblas_int_t incX, float *Y, const cblas_int_t incY);
void catlas_saxpby(const cblas_int_t N, const float alpha, const float *X,
                  const cblas_int_t incX, const float beta, float *Y, const cblas_int_t incY);
void catlas_sset
   (const cblas_int_t N, const float alpha, float *X, const cblas_int_t incX);

void cblas_dswap(const cblas_int_t N, double *X, const cblas_int_t incX,
                 double *Y, const cblas_int_t incY);
void cblas_dcopy(const cblas_int_t N, const double *X, const cblas_int_t incX,
                 double *Y, const cblas_int_t incY);
void cblas_daxpy(const cblas_int_t N, const double alpha, const double *X,
                 const cblas_int_t incX, double *Y, const cblas_int_t incY);
void catlas_daxpby(const cblas_int_t N, const double alpha, const double *X,
                  const cblas_int_t incX, const double beta, double *Y, const cblas_int_t incY);
void catlas_dset
   (const cblas_int_t N, const double alpha, double *X, const cblas_int_t incX);

void cblas_cswap(const cblas_int_t N, void *X, const cblas_int_t incX,
                 void *Y, const cblas_int_t incY);
void cblas_ccopy(const cblas_int_t N, const void *X, const cblas_int_t incX,
                 void *Y, const cblas_int_t incY);
void cblas_caxpy(const cblas_int_t N, const void *alpha, const void *X,
                 const cblas_int_t incX, void *Y, const cblas_int_t incY);
void catlas_caxpby(const cblas_int_t N, const void *alpha, const void *X,
                  const cblas_int_t incX, const void *beta, void *Y, const cblas_int_t incY);
void catlas_cset
   (const cblas_int_t N, const void *alpha, void *X, const cblas_int_t incX);

void cblas_zswap(const cblas_int_t N, void *X, const cblas_int_t incX,
                 void *Y, const cblas_int_t incY);
void cblas_zcopy(const cblas_int_t N, const void *X, const cblas_int_t incX,
                 void *Y, const cblas_int_t incY);
void cblas_zaxpy(const cblas_int_t N, const void *alpha, const void *X,
                 const cblas_int_t incX, void *Y, const cblas_int_t incY);
void catlas_zaxpby(const cblas_int_t N, const void *alpha, const void *X,
                  const cblas_int_t incX, const void *beta, void *Y, const cblas_int_t incY);
void catlas_zset
   (const cblas_int_t N, const void *alpha, void *X, const cblas_int_t incX);


/*
 * Routines with S and D prefix only
 */
void cblas_srotg(float *a, float *b, float *c, float *s);
void cblas_srotmg(float *d1, float *d2, float *b1, const float b2, float *P);
void cblas_srot(const cblas_int_t N, float *X, const cblas_int_t incX,
                float *Y, const cblas_int_t incY, const float c, const float s);
void cblas_srotm(const cblas_int_t N, float *X, const cblas_int_t incX,
                float *Y, const cblas_int_t incY, const float *P);

void cblas_drotg(double *a, double *b, double *c, double *s);
void cblas_drotmg(double *d1, double *d2, double *b1, const double b2, double *P);
void cblas_drot(const cblas_int_t N, double *X, const cblas_int_t incX,
                double *Y, const cblas_int_t incY, const double c, const double s);
void cblas_drotm(const cblas_int_t N, double *X, const cblas_int_t incX,
                double *Y, const cblas_int_t incY, const double *P);


/*
 * Routines with S D C Z CS and ZD prefixes
 */
void cblas_sscal(const cblas_int_t N, const float alpha, float *X, const cblas_int_t incX);
void cblas_dscal(const cblas_int_t N, const double alpha, double *X, const cblas_int_t incX);
void cblas_cscal(const cblas_int_t N, const void *alpha, void *X, const cblas_int_t incX);
void cblas_zscal(const cblas_int_t N, const void *alpha, void *X, const cblas_int_t incX);
void cblas_csscal(const cblas_int_t N, const float alpha, void *X, const cblas_int_t incX);
void cblas_zdscal(const cblas_int_t N, const double alpha, void *X, const cblas_int_t incX);

/*
 * Extra reference routines provided by ATLAS, but not mandated by the standard
 */
void cblas_crotg(void *a, void *b, void *c, void *s);
void cblas_zrotg(void *a, void *b, void *c, void *s);
void cblas_csrot(const cblas_int_t N, void *X, const cblas_int_t incX, void *Y, const cblas_int_t incY,
                 const float c, const float s);
void cblas_zdrot(const cblas_int_t N, void *X, const cblas_int_t incX, void *Y, const cblas_int_t incY,
                 const double c, const double s);

/*
 * ===========================================================================
 * Prototypes for level 2 BLAS
 * ===========================================================================
 */

/*
 * Routines with standard 4 prefixes (S, D, C, Z)
 */
void cblas_sgemv(const enum CBLAS_ORDER Order,
                 const enum CBLAS_TRANSPOSE TransA, const cblas_int_t M, const cblas_int_t N,
                 const float alpha, const float *A, const cblas_int_t lda,
                 const float *X, const cblas_int_t incX, const float beta,
                 float *Y, const cblas_int_t incY);
void cblas_sgbmv(const enum CBLAS_ORDER Order,
                 const enum CBLAS_TRANSPOSE TransA, const cblas_int_t M, const cblas_int_t N,
                 const cblas_int_t KL, const cblas_int_t KU, const float alpha,
                 const float *A, const cblas_int_t lda, const float *X,
                 const cblas_int_t incX, const float beta, float *Y, const cblas_int_t incY);
void cblas_strmv(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const cblas_int_t N, const float *A, const cblas_int_t lda,
                 float *X, const cblas_int_t incX);
void cblas_stbmv(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const cblas_int_t N, const cblas_int_t K, const float *A, const cblas_int_t lda,
                 float *X, const cblas_int_t incX);
void cblas_stpmv(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const cblas_int_t N, const float *Ap, float *X, const cblas_int_t incX);
void cblas_strsv(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const cblas_int_t N, const float *A, const cblas_int_t lda, float *X,
                 const cblas_int_t incX);
void cblas_stbsv(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const cblas_int_t N, const cblas_int_t K, const float *A, const cblas_int_t lda,
                 float *X, const cblas_int_t incX);
void cblas_stpsv(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const cblas_int_t N, const float *Ap, float *X, const cblas_int_t incX);

void cblas_dgemv(const enum CBLAS_ORDER Order,
                 const enum CBLAS_TRANSPOSE TransA, const cblas_int_t M, const cblas_int_t N,
                 const double alpha, const double *A, const cblas_int_t lda,
                 const double *X, const cblas_int_t incX, const double beta,
                 double *Y, const cblas_int_t incY);
void cblas_dgbmv(const enum CBLAS_ORDER Order,
                 const enum CBLAS_TRANSPOSE TransA, const cblas_int_t M, const cblas_int_t N,
                 const cblas_int_t KL, const cblas_int_t KU, const double alpha,
                 const double *A, const cblas_int_t lda, const double *X,
                 const cblas_int_t incX, const double beta, double *Y, const cblas_int_t incY);
void cblas_dtrmv(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const cblas_int_t N, const double *A, const cblas_int_t lda,
                 double *X, const cblas_int_t incX);
void cblas_dtbmv(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const cblas_int_t N, const cblas_int_t K, const double *A, const cblas_int_t lda,
                 double *X, const cblas_int_t incX);
void cblas_dtpmv(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const cblas_int_t N, const double *Ap, double *X, const cblas_int_t incX);
void cblas_dtrsv(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const cblas_int_t N, const double *A, const cblas_int_t lda, double *X,
                 const cblas_int_t incX);
void cblas_dtbsv(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const cblas_int_t N, const cblas_int_t K, const double *A, const cblas_int_t lda,
                 double *X, const cblas_int_t incX);
void cblas_dtpsv(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const cblas_int_t N, const double *Ap, double *X, const cblas_int_t incX);

void cblas_cgemv(const enum CBLAS_ORDER Order,
                 const enum CBLAS_TRANSPOSE TransA, const cblas_int_t M, const cblas_int_t N,
                 const void *alpha, const void *A, const cblas_int_t lda,
                 const void *X, const cblas_int_t incX, const void *beta,
                 void *Y, const cblas_int_t incY);
void cblas_cgbmv(const enum CBLAS_ORDER Order,
                 const enum CBLAS_TRANSPOSE TransA, const cblas_int_t M, const cblas_int_t N,
                 const cblas_int_t KL, const cblas_int_t KU, const void *alpha,
                 const void *A, const cblas_int_t lda, const void *X,
                 const cblas_int_t incX, const void *beta, void *Y, const cblas_int_t incY);
void cblas_ctrmv(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const cblas_int_t N, const void *A, const cblas_int_t lda,
                 void *X, const cblas_int_t incX);
void cblas_ctbmv(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const cblas_int_t N, const cblas_int_t K, const void *A, const cblas_int_t lda,
                 void *X, const cblas_int_t incX);
void cblas_ctpmv(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const cblas_int_t N, const void *Ap, void *X, const cblas_int_t incX);
void cblas_ctrsv(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const cblas_int_t N, const void *A, const cblas_int_t lda, void *X,
                 const cblas_int_t incX);
void cblas_ctbsv(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const cblas_int_t N, const cblas_int_t K, const void *A, const cblas_int_t lda,
                 void *X, const cblas_int_t incX);
void cblas_ctpsv(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const cblas_int_t N, const void *Ap, void *X, const cblas_int_t incX);

void cblas_zgemv(const enum CBLAS_ORDER Order,
                 const enum CBLAS_TRANSPOSE TransA, const cblas_int_t M, const cblas_int_t N,
                 const void *alpha, const void *A, const cblas_int_t lda,
                 const void *X, const cblas_int_t incX, const void *beta,
                 void *Y, const cblas_int_t incY);
void cblas_zgbmv(const enum CBLAS_ORDER Order,
                 const enum CBLAS_TRANSPOSE TransA, const cblas_int_t M, const cblas_int_t N,
                 const cblas_int_t KL, const cblas_int_t KU, const void *alpha,
                 const void *A, const cblas_int_t lda, const void *X,
                 const cblas_int_t incX, const void *beta, void *Y, const cblas_int_t incY);
void cblas_ztrmv(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const cblas_int_t N, const void *A, const cblas_int_t lda,
                 void *X, const cblas_int_t incX);
void cblas_ztbmv(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const cblas_int_t N, const cblas_int_t K, const void *A, const cblas_int_t lda,
                 void *X, const cblas_int_t incX);
void cblas_ztpmv(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const cblas_int_t N, const void *Ap, void *X, const cblas_int_t incX);
void cblas_ztrsv(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const cblas_int_t N, const void *A, const cblas_int_t lda, void *X,
                 const cblas_int_t incX);
void cblas_ztbsv(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const cblas_int_t N, const cblas_int_t K, const void *A, const cblas_int_t lda,
                 void *X, const cblas_int_t incX);
void cblas_ztpsv(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const cblas_int_t N, const void *Ap, void *X, const cblas_int_t incX);


/*
 * Routines with S and D prefixes only
 */
void cblas_ssymv(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                 const cblas_int_t N, const float alpha, const float *A,
                 const cblas_int_t lda, const float *X, const cblas_int_t incX,
                 const float beta, float *Y, const cblas_int_t incY);
void cblas_ssbmv(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                 const cblas_int_t N, const cblas_int_t K, const float alpha, const float *A,
                 const cblas_int_t lda, const float *X, const cblas_int_t incX,
                 const float beta, float *Y, const cblas_int_t incY);
void cblas_sspmv(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                 const cblas_int_t N, const float alpha, const float *Ap,
                 const float *X, const cblas_int_t incX,
                 const float beta, float *Y, const cblas_int_t incY);
void cblas_sger(const enum CBLAS_ORDER Order, const cblas_int_t M, const cblas_int_t N,
                const float alpha, const float *X, const cblas_int_t incX,
                const float *Y, const cblas_int_t incY, float *A, const cblas_int_t lda);
void cblas_ssyr(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                const cblas_int_t N, const float alpha, const float *X,
                const cblas_int_t incX, float *A, const cblas_int_t lda);
void cblas_sspr(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                const cblas_int_t N, const float alpha, const float *X,
                const cblas_int_t incX, float *Ap);
void cblas_ssyr2(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                const cblas_int_t N, const float alpha, const float *X,
                const cblas_int_t incX, const float *Y, const cblas_int_t incY, float *A,
                const cblas_int_t lda);
void cblas_sspr2(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                const cblas_int_t N, const float alpha, const float *X,
                const cblas_int_t incX, const float *Y, const cblas_int_t incY, float *A);

void cblas_dsymv(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                 const cblas_int_t N, const double alpha, const double *A,
                 const cblas_int_t lda, const double *X, const cblas_int_t incX,
                 const double beta, double *Y, const cblas_int_t incY);
void cblas_dsbmv(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                 const cblas_int_t N, const cblas_int_t K, const double alpha, const double *A,
                 const cblas_int_t lda, const double *X, const cblas_int_t incX,
                 const double beta, double *Y, const cblas_int_t incY);
void cblas_dspmv(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                 const cblas_int_t N, const double alpha, const double *Ap,
                 const double *X, const cblas_int_t incX,
                 const double beta, double *Y, const cblas_int_t incY);
void cblas_dger(const enum CBLAS_ORDER Order, const cblas_int_t M, const cblas_int_t N,
                const double alpha, const double *X, const cblas_int_t incX,
                const double *Y, const cblas_int_t incY, double *A, const cblas_int_t lda);
void cblas_dsyr(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                const cblas_int_t N, const double alpha, const double *X,
                const cblas_int_t incX, double *A, const cblas_int_t lda);
void cblas_dspr(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                const cblas_int_t N, const double alpha, const double *X,
                const cblas_int_t incX, double *Ap);
void cblas_dsyr2(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                const cblas_int_t N, const double alpha, const double *X,
                const cblas_int_t incX, const double *Y, const cblas_int_t incY, double *A,
                const cblas_int_t lda);
void cblas_dspr2(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                const cblas_int_t N, const double alpha, const double *X,
                const cblas_int_t incX, const double *Y, const cblas_int_t incY, double *A);


/*
 * Routines with C and Z prefixes only
 */
void cblas_chemv(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                 const cblas_int_t N, const void *alpha, const void *A,
                 const cblas_int_t lda, const void *X, const cblas_int_t incX,
                 const void *beta, void *Y, const cblas_int_t incY);
void cblas_chbmv(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                 const cblas_int_t N, const cblas_int_t K, const void *alpha, const void *A,
                 const cblas_int_t lda, const void *X, const cblas_int_t incX,
                 const void *beta, void *Y, const cblas_int_t incY);
void cblas_chpmv(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                 const cblas_int_t N, const void *alpha, const void *Ap,
                 const void *X, const cblas_int_t incX,
                 const void *beta, void *Y, const cblas_int_t incY);
void cblas_cgeru(const enum CBLAS_ORDER Order, const cblas_int_t M, const cblas_int_t N,
                 const void *alpha, const void *X, const cblas_int_t incX,
                 const void *Y, const cblas_int_t incY, void *A, const cblas_int_t lda);
void cblas_cgerc(const enum CBLAS_ORDER Order, const cblas_int_t M, const cblas_int_t N,
                 const void *alpha, const void *X, const cblas_int_t incX,
                 const void *Y, const cblas_int_t incY, void *A, const cblas_int_t lda);
void cblas_cher(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                const cblas_int_t N, const float alpha, const void *X, const cblas_int_t incX,
                void *A, const cblas_int_t lda);
void cblas_chpr(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                const cblas_int_t N, const float alpha, const void *X,
                const cblas_int_t incX, void *A);
void cblas_cher2(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo, const cblas_int_t N,
                const void *alpha, const void *X, const cblas_int_t incX,
                const void *Y, const cblas_int_t incY, void *A, const cblas_int_t lda);
void cblas_chpr2(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo, const cblas_int_t N,
                const void *alpha, const void *X, const cblas_int_t incX,
                const void *Y, const cblas_int_t incY, void *Ap);

void cblas_zhemv(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                 const cblas_int_t N, const void *alpha, const void *A,
                 const cblas_int_t lda, const void *X, const cblas_int_t incX,
                 const void *beta, void *Y, const cblas_int_t incY);
void cblas_zhbmv(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                 const cblas_int_t N, const cblas_int_t K, const void *alpha, const void *A,
                 const cblas_int_t lda, const void *X, const cblas_int_t incX,
                 const void *beta, void *Y, const cblas_int_t incY);
void cblas_zhpmv(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                 const cblas_int_t N, const void *alpha, const void *Ap,
                 const void *X, const cblas_int_t incX,
                 const void *beta, void *Y, const cblas_int_t incY);
void cblas_zgeru(const enum CBLAS_ORDER Order, const cblas_int_t M, const cblas_int_t N,
                 const void *alpha, const void *X, const cblas_int_t incX,
                 const void *Y, const cblas_int_t incY, void *A, const cblas_int_t lda);
void cblas_zgerc(const enum CBLAS_ORDER Order, const cblas_int_t M, const cblas_int_t N,
                 const void *alpha, const void *X, const cblas_int_t incX,
                 const void *Y, const cblas_int_t incY, void *A, const cblas_int_t lda);
void cblas_zher(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                const cblas_int_t N, const double alpha, const void *X, const cblas_int_t incX,
                void *A, const cblas_int_t lda);
void cblas_zhpr(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                const cblas_int_t N, const double alpha, const void *X,
                const cblas_int_t incX, void *A);
void cblas_zher2(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo, const cblas_int_t N,
                const void *alpha, const void *X, const cblas_int_t incX,
                const void *Y, const cblas_int_t incY, void *A, const cblas_int_t lda);
void cblas_zhpr2(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo, const cblas_int_t N,
                const void *alpha, const void *X, const cblas_int_t incX,
                const void *Y, const cblas_int_t incY, void *Ap);

/*
 * ===========================================================================
 * Prototypes for level 3 BLAS
 * ===========================================================================
 */

/*
 * Routines with standard 4 prefixes (S, D, C, Z)
 */
void cblas_sgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_TRANSPOSE TransB, const cblas_int_t M, const cblas_int_t N,
                 const cblas_int_t K, const float alpha, const float *A,
                 const cblas_int_t lda, const float *B, const cblas_int_t ldb,
                 const float beta, float *C, const cblas_int_t ldc);
void cblas_ssymm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo, const cblas_int_t M, const cblas_int_t N,
                 const float alpha, const float *A, const cblas_int_t lda,
                 const float *B, const cblas_int_t ldb, const float beta,
                 float *C, const cblas_int_t ldc);
void cblas_ssyrk(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE Trans, const cblas_int_t N, const cblas_int_t K,
                 const float alpha, const float *A, const cblas_int_t lda,
                 const float beta, float *C, const cblas_int_t ldc);
void cblas_ssyr2k(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                  const enum CBLAS_TRANSPOSE Trans, const cblas_int_t N, const cblas_int_t K,
                  const float alpha, const float *A, const cblas_int_t lda,
                  const float *B, const cblas_int_t ldb, const float beta,
                  float *C, const cblas_int_t ldc);
void cblas_strmm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_DIAG Diag, const cblas_int_t M, const cblas_int_t N,
                 const float alpha, const float *A, const cblas_int_t lda,
                 float *B, const cblas_int_t ldb);
void cblas_strsm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_DIAG Diag, const cblas_int_t M, const cblas_int_t N,
                 const float alpha, const float *A, const cblas_int_t lda,
                 float *B, const cblas_int_t ldb);

void cblas_dgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_TRANSPOSE TransB, const cblas_int_t M, const cblas_int_t N,
                 const cblas_int_t K, const double alpha, const double *A,
                 const cblas_int_t lda, const double *B, const cblas_int_t ldb,
                 const double beta, double *C, const cblas_int_t ldc);
void cblas_dsymm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo, const cblas_int_t M, const cblas_int_t N,
                 const double alpha, const double *A, const cblas_int_t lda,
                 const double *B, const cblas_int_t ldb, const double beta,
                 double *C, const cblas_int_t ldc);
void cblas_dsyrk(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE Trans, const cblas_int_t N, const cblas_int_t K,
                 const double alpha, const double *A, const cblas_int_t lda,
                 const double beta, double *C, const cblas_int_t ldc);
void cblas_dsyr2k(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                  const enum CBLAS_TRANSPOSE Trans, const cblas_int_t N, const cblas_int_t K,
                  const double alpha, const double *A, const cblas_int_t lda,
                  const double *B, const cblas_int_t ldb, const double beta,
                  double *C, const cblas_int_t ldc);
void cblas_dtrmm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_DIAG Diag, const cblas_int_t M, const cblas_int_t N,
                 const double alpha, const double *A, const cblas_int_t lda,
                 double *B, const cblas_int_t ldb);
void cblas_dtrsm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_DIAG Diag, const cblas_int_t M, const cblas_int_t N,
                 const double alpha, const double *A, const cblas_int_t lda,
                 double *B, const cblas_int_t ldb);

void cblas_cgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_TRANSPOSE TransB, const cblas_int_t M, const cblas_int_t N,
                 const cblas_int_t K, const void *alpha, const void *A,
                 const cblas_int_t lda, const void *B, const cblas_int_t ldb,
                 const void *beta, void *C, const cblas_int_t ldc);
void cblas_csymm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo, const cblas_int_t M, const cblas_int_t N,
                 const void *alpha, const void *A, const cblas_int_t lda,
                 const void *B, const cblas_int_t ldb, const void *beta,
                 void *C, const cblas_int_t ldc);
void cblas_csyrk(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE Trans, const cblas_int_t N, const cblas_int_t K,
                 const void *alpha, const void *A, const cblas_int_t lda,
                 const void *beta, void *C, const cblas_int_t ldc);
void cblas_csyr2k(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                  const enum CBLAS_TRANSPOSE Trans, const cblas_int_t N, const cblas_int_t K,
                  const void *alpha, const void *A, const cblas_int_t lda,
                  const void *B, const cblas_int_t ldb, const void *beta,
                  void *C, const cblas_int_t ldc);
void cblas_ctrmm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_DIAG Diag, const cblas_int_t M, const cblas_int_t N,
                 const void *alpha, const void *A, const cblas_int_t lda,
                 void *B, const cblas_int_t ldb);
void cblas_ctrsm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_DIAG Diag, const cblas_int_t M, const cblas_int_t N,
                 const void *alpha, const void *A, const cblas_int_t lda,
                 void *B, const cblas_int_t ldb);

void cblas_zgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_TRANSPOSE TransB, const cblas_int_t M, const cblas_int_t N,
                 const cblas_int_t K, const void *alpha, const void *A,
                 const cblas_int_t lda, const void *B, const cblas_int_t ldb,
                 const void *beta, void *C, const cblas_int_t ldc);
void cblas_zsymm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo, const cblas_int_t M, const cblas_int_t N,
                 const void *alpha, const void *A, const cblas_int_t lda,
                 const void *B, const cblas_int_t ldb, const void *beta,
                 void *C, const cblas_int_t ldc);
void cblas_zsyrk(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE Trans, const cblas_int_t N, const cblas_int_t K,
                 const void *alpha, const void *A, const cblas_int_t lda,
                 const void *beta, void *C, const cblas_int_t ldc);
void cblas_zsyr2k(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                  const enum CBLAS_TRANSPOSE Trans, const cblas_int_t N, const cblas_int_t K,
                  const void *alpha, const void *A, const cblas_int_t lda,
                  const void *B, const cblas_int_t ldb, const void *beta,
                  void *C, const cblas_int_t ldc);
void cblas_ztrmm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_DIAG Diag, const cblas_int_t M, const cblas_int_t N,
                 const void *alpha, const void *A, const cblas_int_t lda,
                 void *B, const cblas_int_t ldb);
void cblas_ztrsm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_DIAG Diag, const cblas_int_t M, const cblas_int_t N,
                 const void *alpha, const void *A, const cblas_int_t lda,
                 void *B, const cblas_int_t ldb);


/*
 * Routines with prefixes C and Z only
 */
void cblas_chemm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo, const cblas_int_t M, const cblas_int_t N,
                 const void *alpha, const void *A, const cblas_int_t lda,
                 const void *B, const cblas_int_t ldb, const void *beta,
                 void *C, const cblas_int_t ldc);
void cblas_cherk(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE Trans, const cblas_int_t N, const cblas_int_t K,
                 const float alpha, const void *A, const cblas_int_t lda,
                 const float beta, void *C, const cblas_int_t ldc);
void cblas_cher2k(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                  const enum CBLAS_TRANSPOSE Trans, const cblas_int_t N, const cblas_int_t K,
                  const void *alpha, const void *A, const cblas_int_t lda,
                  const void *B, const cblas_int_t ldb, const float beta,
                  void *C, const cblas_int_t ldc);
void cblas_zhemm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo, const cblas_int_t M, const cblas_int_t N,
                 const void *alpha, const void *A, const cblas_int_t lda,
                 const void *B, const cblas_int_t ldb, const void *beta,
                 void *C, const cblas_int_t ldc);
void cblas_zherk(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE Trans, const cblas_int_t N, const cblas_int_t K,
                 const double alpha, const void *A, const cblas_int_t lda,
                 const double beta, void *C, const cblas_int_t ldc);
void cblas_zher2k(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                  const enum CBLAS_TRANSPOSE Trans, const cblas_int_t N, const cblas_int_t K,
                  const void *alpha, const void *A, const cblas_int_t lda,
                  const void *B, const cblas_int_t ldb, const double beta,
                  void *C, const cblas_int_t ldc);

cblas_int_t cblas_errprn(cblas_int_t ierr, cblas_int_t info, char *form, ...);

#ifdef __cplusplus
}
#endif

#endif // _GLUE_CBLAS_H_INCLUDED_
