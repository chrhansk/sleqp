/** 
 * @file clapack.h
 * @brief      C LAPACK header file
 * @author     NETLIB Authors, Based on clapack.h from http://www.netlib.org/clapack/clapack.h
 * @date       2009 Jun 1
 */

#pragma once
#ifndef _GLUE_CLAPACK_H_INCLUDE_
#define _GLUE_CLAPACK_H_INCLUDE_

#ifdef __cplusplus
extern "C" {
#endif

#include "cblas.h"

typedef struct { float r, i; } complex;
typedef struct { double r, i; } doublecomplex;

/* Subroutine */ int cbdsqr_(const char *uplo, cblas_int_t *n, cblas_int_t *ncvt, cblas_int_t *
    nru, cblas_int_t *ncc, float *d__, float *e, complex *vt, cblas_int_t *ldvt, 
    complex *u, cblas_int_t *ldu, complex *c__, cblas_int_t *ldc, float *rwork, 
    cblas_int_t *info);

/* Subroutine */ int cgbbrd_(const char *vect, cblas_int_t *m, cblas_int_t *n, cblas_int_t *ncc,
     cblas_int_t *kl, cblas_int_t *ku, complex *ab, cblas_int_t *ldab, float *d__, 
    float *e, complex *q, cblas_int_t *ldq, complex *pt, cblas_int_t *ldpt, 
    complex *c__, cblas_int_t *ldc, complex *work, float *rwork, cblas_int_t *info);

/* Subroutine */ int cgbcon_(const char *norm, cblas_int_t *n, cblas_int_t *kl, cblas_int_t *ku,
     complex *ab, cblas_int_t *ldab, cblas_int_t *ipiv, float *anorm, float *rcond, 
    complex *work, float *rwork, cblas_int_t *info);

/* Subroutine */ int cgbequ_(cblas_int_t *m, cblas_int_t *n, cblas_int_t *kl, cblas_int_t *ku,
     complex *ab, cblas_int_t *ldab, float *r__, float *c__, float *rowcnd, float 
    *colcnd, float *amax, cblas_int_t *info);

/* Subroutine */ int cgbrfs_(const char *trans, cblas_int_t *n, cblas_int_t *kl, cblas_int_t *
    ku, cblas_int_t *nrhs, complex *ab, cblas_int_t *ldab, complex *afb, cblas_int_t *
    ldafb, cblas_int_t *ipiv, complex *b, cblas_int_t *ldb, complex *x, cblas_int_t *
    ldx, float *ferr, float *berr, complex *work, float *rwork, cblas_int_t *
    info);

/* Subroutine */ int cgbsv_(cblas_int_t *n, cblas_int_t *kl, cblas_int_t *ku, cblas_int_t *
    nrhs, complex *ab, cblas_int_t *ldab, cblas_int_t *ipiv, complex *b, cblas_int_t *
    ldb, cblas_int_t *info);

/* Subroutine */ int cgbsvx_(const char *fact, const char *trans, cblas_int_t *n, cblas_int_t *kl,
     cblas_int_t *ku, cblas_int_t *nrhs, complex *ab, cblas_int_t *ldab, complex *afb,
     cblas_int_t *ldafb, cblas_int_t *ipiv, const char *equed, float *r__, float *c__, 
    complex *b, cblas_int_t *ldb, complex *x, cblas_int_t *ldx, float *rcond, float 
    *ferr, float *berr, complex *work, float *rwork, cblas_int_t *info);

/* Subroutine */ int cgbtf2_(cblas_int_t *m, cblas_int_t *n, cblas_int_t *kl, cblas_int_t *ku,
     complex *ab, cblas_int_t *ldab, cblas_int_t *ipiv, cblas_int_t *info);

/* Subroutine */ int cgbtrf_(cblas_int_t *m, cblas_int_t *n, cblas_int_t *kl, cblas_int_t *ku,
     complex *ab, cblas_int_t *ldab, cblas_int_t *ipiv, cblas_int_t *info);

/* Subroutine */ int cgbtrs_(const char *trans, cblas_int_t *n, cblas_int_t *kl, cblas_int_t *
    ku, cblas_int_t *nrhs, complex *ab, cblas_int_t *ldab, cblas_int_t *ipiv, complex 
    *b, cblas_int_t *ldb, cblas_int_t *info);

/* Subroutine */ int cgebak_(const char *job, const char *side, cblas_int_t *n, cblas_int_t *ilo, 
    cblas_int_t *ihi, float *scale, cblas_int_t *m, complex *v, cblas_int_t *ldv, 
    cblas_int_t *info);

/* Subroutine */ int cgebal_(const char *job, cblas_int_t *n, complex *a, cblas_int_t *lda, 
    cblas_int_t *ilo, cblas_int_t *ihi, float *scale, cblas_int_t *info);

/* Subroutine */ int cgebd2_(cblas_int_t *m, cblas_int_t *n, complex *a, cblas_int_t *lda,
     float *d__, float *e, complex *tauq, complex *taup, complex *work, 
    cblas_int_t *info);

/* Subroutine */ int cgebrd_(cblas_int_t *m, cblas_int_t *n, complex *a, cblas_int_t *lda,
     float *d__, float *e, complex *tauq, complex *taup, complex *work, 
    cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int cgecon_(const char *norm, cblas_int_t *n, complex *a, cblas_int_t *lda,
     float *anorm, float *rcond, complex *work, float *rwork, cblas_int_t *info);

/* Subroutine */ int cgeequ_(cblas_int_t *m, cblas_int_t *n, complex *a, cblas_int_t *lda,
     float *r__, float *c__, float *rowcnd, float *colcnd, float *amax, 
    cblas_int_t *info);

/* Subroutine */ int cgees_(const char *jobvs, const char *sort, L_fp select, cblas_int_t *n, 
    complex *a, cblas_int_t *lda, cblas_int_t *sdim, complex *w, complex *vs, 
    cblas_int_t *ldvs, complex *work, cblas_int_t *lwork, float *rwork, logical *
    bwork, cblas_int_t *info);

/* Subroutine */ int cgeesx_(const char *jobvs, const char *sort, L_fp select, const char *
    sense, cblas_int_t *n, complex *a, cblas_int_t *lda, cblas_int_t *sdim, complex *
    w, complex *vs, cblas_int_t *ldvs, float *rconde, float *rcondv, complex *
    work, cblas_int_t *lwork, float *rwork, logical *bwork, cblas_int_t *info);

/* Subroutine */ int cgeev_(const char *jobvl, const char *jobvr, cblas_int_t *n, complex *a, 
    cblas_int_t *lda, complex *w, complex *vl, cblas_int_t *ldvl, complex *vr, 
    cblas_int_t *ldvr, complex *work, cblas_int_t *lwork, float *rwork, cblas_int_t *
    info);

/* Subroutine */ int cgeevx_(const char *balanc, const char *jobvl, const char *jobvr, const char *
    sense, cblas_int_t *n, complex *a, cblas_int_t *lda, complex *w, complex *vl, 
    cblas_int_t *ldvl, complex *vr, cblas_int_t *ldvr, cblas_int_t *ilo, cblas_int_t *ihi,
     float *scale, float *abnrm, float *rconde, float *rcondv, complex *work,
    cblas_int_t *lwork, float *rwork, cblas_int_t *info);

/* Subroutine */ int cgegs_(const char *jobvsl, const char *jobvsr, cblas_int_t *n, complex *
    a, cblas_int_t *lda, complex *b, cblas_int_t *ldb, complex *alpha, complex *
    beta, complex *vsl, cblas_int_t *ldvsl, complex *vsr, cblas_int_t *ldvsr,
    complex *work, cblas_int_t *lwork, float *rwork, cblas_int_t *info);

/* Subroutine */ int cgegv_(const char *jobvl, const char *jobvr, cblas_int_t *n, complex *a, 
    cblas_int_t *lda, complex *b, cblas_int_t *ldb, complex *alpha, complex *beta,
     complex *vl, cblas_int_t *ldvl, complex *vr, cblas_int_t *ldvr, complex *
    work, cblas_int_t *lwork, float *rwork, cblas_int_t *info);

/* Subroutine */ int cgehd2_(cblas_int_t *n, cblas_int_t *ilo, cblas_int_t *ihi, complex *
    a, cblas_int_t *lda, complex *tau, complex *work, cblas_int_t *info);

/* Subroutine */ int cgehrd_(cblas_int_t *n, cblas_int_t *ilo, cblas_int_t *ihi, complex *
    a, cblas_int_t *lda, complex *tau, complex *work, cblas_int_t *lwork, cblas_int_t 
    *info);

/* Subroutine */ int cgelq2_(cblas_int_t *m, cblas_int_t *n, complex *a, cblas_int_t *lda,
     complex *tau, complex *work, cblas_int_t *info);

/* Subroutine */ int cgelqf_(cblas_int_t *m, cblas_int_t *n, complex *a, cblas_int_t *lda,
     complex *tau, complex *work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int cgels_(const char *trans, cblas_int_t *m, cblas_int_t *n, cblas_int_t *
    nrhs, complex *a, cblas_int_t *lda, complex *b, cblas_int_t *ldb, complex *
    work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int cgelsd_(cblas_int_t *m, cblas_int_t *n, cblas_int_t *nrhs, complex *
    a, cblas_int_t *lda, complex *b, cblas_int_t *ldb, float *s, float *rcond, 
    cblas_int_t *rank, complex *work, cblas_int_t *lwork, float *rwork, cblas_int_t *
    iwork, cblas_int_t *info);

/* Subroutine */ int cgelss_(cblas_int_t *m, cblas_int_t *n, cblas_int_t *nrhs, complex *
    a, cblas_int_t *lda, complex *b, cblas_int_t *ldb, float *s, float *rcond, 
    cblas_int_t *rank, complex *work, cblas_int_t *lwork, float *rwork, cblas_int_t *
    info);

/* Subroutine */ int cgelsx_(cblas_int_t *m, cblas_int_t *n, cblas_int_t *nrhs, complex *
    a, cblas_int_t *lda, complex *b, cblas_int_t *ldb, cblas_int_t *jpvt, float *rcond,
     cblas_int_t *rank, complex *work, float *rwork, cblas_int_t *info);

/* Subroutine */ int cgelsy_(cblas_int_t *m, cblas_int_t *n, cblas_int_t *nrhs, complex *
    a, cblas_int_t *lda, complex *b, cblas_int_t *ldb, cblas_int_t *jpvt, float *rcond,
     cblas_int_t *rank, complex *work, cblas_int_t *lwork, float *rwork, cblas_int_t *
    info);

/* Subroutine */ int cgeql2_(cblas_int_t *m, cblas_int_t *n, complex *a, cblas_int_t *lda,
     complex *tau, complex *work, cblas_int_t *info);

/* Subroutine */ int cgeqlf_(cblas_int_t *m, cblas_int_t *n, complex *a, cblas_int_t *lda,
     complex *tau, complex *work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int cgeqp3_(cblas_int_t *m, cblas_int_t *n, complex *a, cblas_int_t *lda,
     cblas_int_t *jpvt, complex *tau, complex *work, cblas_int_t *lwork, float *
    rwork, cblas_int_t *info);

/* Subroutine */ int cgeqpf_(cblas_int_t *m, cblas_int_t *n, complex *a, cblas_int_t *lda,
     cblas_int_t *jpvt, complex *tau, complex *work, float *rwork, cblas_int_t *
    info);

/* Subroutine */ int cgeqr2_(cblas_int_t *m, cblas_int_t *n, complex *a, cblas_int_t *lda,
     complex *tau, complex *work, cblas_int_t *info);

/* Subroutine */ int cgeqrf_(cblas_int_t *m, cblas_int_t *n, complex *a, cblas_int_t *lda,
     complex *tau, complex *work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int cgerfs_(const char *trans, cblas_int_t *n, cblas_int_t *nrhs, complex *
    a, cblas_int_t *lda, complex *af, cblas_int_t *ldaf, cblas_int_t *ipiv, complex *
    b, cblas_int_t *ldb, complex *x, cblas_int_t *ldx, float *ferr, float *berr, 
    complex *work, float *rwork, cblas_int_t *info);

/* Subroutine */ int cgerq2_(cblas_int_t *m, cblas_int_t *n, complex *a, cblas_int_t *lda,
     complex *tau, complex *work, cblas_int_t *info);

/* Subroutine */ int cgerqf_(cblas_int_t *m, cblas_int_t *n, complex *a, cblas_int_t *lda,
     complex *tau, complex *work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int cgesc2_(cblas_int_t *n, complex *a, cblas_int_t *lda, complex *
    rhs, cblas_int_t *ipiv, cblas_int_t *jpiv, float *scale);

/* Subroutine */ int cgesdd_(const char *jobz, cblas_int_t *m, cblas_int_t *n, complex *a, 
    cblas_int_t *lda, float *s, complex *u, cblas_int_t *ldu, complex *vt, cblas_int_t 
    *ldvt, complex *work, cblas_int_t *lwork, float *rwork, cblas_int_t *iwork, 
    cblas_int_t *info);

/* Subroutine */ int cgesv_(cblas_int_t *n, cblas_int_t *nrhs, complex *a, cblas_int_t *
    lda, cblas_int_t *ipiv, complex *b, cblas_int_t *ldb, cblas_int_t *info);

/* Subroutine */ int cgesvd_(const char *jobu, const char *jobvt, cblas_int_t *m, cblas_int_t *n, 
    complex *a, cblas_int_t *lda, float *s, complex *u, cblas_int_t *ldu, complex *
    vt, cblas_int_t *ldvt, complex *work, cblas_int_t *lwork, float *rwork, 
    cblas_int_t *info);

/* Subroutine */ int cgesvx_(const char *fact, const char *trans, cblas_int_t *n, cblas_int_t *
    nrhs, complex *a, cblas_int_t *lda, complex *af, cblas_int_t *ldaf, cblas_int_t *
    ipiv, const char *equed, float *r__, float *c__, complex *b, cblas_int_t *ldb, 
    complex *x, cblas_int_t *ldx, float *rcond, float *ferr, float *berr, 
    complex *work, float *rwork, cblas_int_t *info);

/* Subroutine */ int cgetc2_(cblas_int_t *n, complex *a, cblas_int_t *lda, cblas_int_t *
    ipiv, cblas_int_t *jpiv, cblas_int_t *info);

/* Subroutine */ int cgetf2_(cblas_int_t *m, cblas_int_t *n, complex *a, cblas_int_t *lda,
     cblas_int_t *ipiv, cblas_int_t *info);

/* Subroutine */ int cgetrf_(cblas_int_t *m, cblas_int_t *n, complex *a, cblas_int_t *lda,
     cblas_int_t *ipiv, cblas_int_t *info);

/* Subroutine */ int cgetri_(cblas_int_t *n, complex *a, cblas_int_t *lda, cblas_int_t *
    ipiv, complex *work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int cgetrs_(const char *trans, cblas_int_t *n, cblas_int_t *nrhs, complex *
    a, cblas_int_t *lda, cblas_int_t *ipiv, complex *b, cblas_int_t *ldb, cblas_int_t *
    info);

/* Subroutine */ int cggbak_(const char *job, const char *side, cblas_int_t *n, cblas_int_t *ilo, 
    cblas_int_t *ihi, float *lscale, float *rscale, cblas_int_t *m, complex *v, 
    cblas_int_t *ldv, cblas_int_t *info);

/* Subroutine */ int cggbal_(const char *job, cblas_int_t *n, complex *a, cblas_int_t *lda, 
    complex *b, cblas_int_t *ldb, cblas_int_t *ilo, cblas_int_t *ihi, float *lscale, 
    float *rscale, float *work, cblas_int_t *info);

/* Subroutine */ int cgges_(const char *jobvsl, const char *jobvsr, const char *sort, L_fp 
    selctg, cblas_int_t *n, complex *a, cblas_int_t *lda, complex *b, cblas_int_t *
    ldb, cblas_int_t *sdim, complex *alpha, complex *beta, complex *vsl, 
    cblas_int_t *ldvsl, complex *vsr, cblas_int_t *ldvsr, complex *work, cblas_int_t *
    lwork, float *rwork, logical *bwork, cblas_int_t *info);

/* Subroutine */ int cggesx_(const char *jobvsl, const char *jobvsr, const char *sort, L_fp 
    selctg, const char *sense, cblas_int_t *n, complex *a, cblas_int_t *lda, complex *b,
     cblas_int_t *ldb, cblas_int_t *sdim, complex *alpha, complex *beta, complex *
    vsl, cblas_int_t *ldvsl, complex *vsr, cblas_int_t *ldvsr, float *rconde, float 
    *rcondv, complex *work, cblas_int_t *lwork, float *rwork, cblas_int_t *iwork, 
    cblas_int_t *liwork, logical *bwork, cblas_int_t *info);

/* Subroutine */ int cggev_(const char *jobvl, const char *jobvr, cblas_int_t *n, complex *a, 
    cblas_int_t *lda, complex *b, cblas_int_t *ldb, complex *alpha, complex *beta,
     complex *vl, cblas_int_t *ldvl, complex *vr, cblas_int_t *ldvr, complex *
    work, cblas_int_t *lwork, float *rwork, cblas_int_t *info);

/* Subroutine */ int cggevx_(const char *balanc, const char *jobvl, const char *jobvr, const char *
    sense, cblas_int_t *n, complex *a, cblas_int_t *lda, complex *b, cblas_int_t *ldb,
     complex *alpha, complex *beta, complex *vl, cblas_int_t *ldvl, complex *
    vr, cblas_int_t *ldvr, cblas_int_t *ilo, cblas_int_t *ihi, float *lscale, float *
    rscale, float *abnrm, float *bbnrm, float *rconde, float *rcondv, complex 
    *work, cblas_int_t *lwork, float *rwork, cblas_int_t *iwork, logical *bwork, 
    cblas_int_t *info);

/* Subroutine */ int cggglm_(cblas_int_t *n, cblas_int_t *m, cblas_int_t *p, complex *a, 
    cblas_int_t *lda, complex *b, cblas_int_t *ldb, complex *d__, complex *x, 
    complex *y, complex *work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int cgghrd_(const char *compq, const char *compz, cblas_int_t *n, cblas_int_t *
    ilo, cblas_int_t *ihi, complex *a, cblas_int_t *lda, complex *b, cblas_int_t *ldb,
     complex *q, cblas_int_t *ldq, complex *z__, cblas_int_t *ldz, cblas_int_t *info);

/* Subroutine */ int cgglse_(cblas_int_t *m, cblas_int_t *n, cblas_int_t *p, complex *a, 
    cblas_int_t *lda, complex *b, cblas_int_t *ldb, complex *c__, complex *d__, 
    complex *x, complex *work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int cggqrf_(cblas_int_t *n, cblas_int_t *m, cblas_int_t *p, complex *a, 
    cblas_int_t *lda, complex *taua, complex *b, cblas_int_t *ldb, complex *taub, 
    complex *work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int cggrqf_(cblas_int_t *m, cblas_int_t *p, cblas_int_t *n, complex *a, 
    cblas_int_t *lda, complex *taua, complex *b, cblas_int_t *ldb, complex *taub, 
    complex *work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int cggsvd_(const char *jobu, const char *jobv, const char *jobq, cblas_int_t *m, 
    cblas_int_t *n, cblas_int_t *p, cblas_int_t *k, cblas_int_t *l, complex *a, cblas_int_t *
    lda, complex *b, cblas_int_t *ldb, float *alpha, float *beta, complex *u, 
    cblas_int_t *ldu, complex *v, cblas_int_t *ldv, complex *q, cblas_int_t *ldq, 
    complex *work, float *rwork, cblas_int_t *iwork, cblas_int_t *info);

/* Subroutine */ int cggsvp_(const char *jobu, const char *jobv, const char *jobq, cblas_int_t *m, 
    cblas_int_t *p, cblas_int_t *n, complex *a, cblas_int_t *lda, complex *b, cblas_int_t 
    *ldb, float *tola, float *tolb, cblas_int_t *k, cblas_int_t *l, complex *u, 
    cblas_int_t *ldu, complex *v, cblas_int_t *ldv, complex *q, cblas_int_t *ldq, 
    cblas_int_t *iwork, float *rwork, complex *tau, complex *work, cblas_int_t *
    info);

/* Subroutine */ int cgtcon_(const char *norm, cblas_int_t *n, complex *dl, complex *
    d__, complex *du, complex *du2, cblas_int_t *ipiv, float *anorm, float *
    rcond, complex *work, cblas_int_t *info);

/* Subroutine */ int cgtrfs_(const char *trans, cblas_int_t *n, cblas_int_t *nrhs, complex *
    dl, complex *d__, complex *du, complex *dlf, complex *df, complex *
    duf, complex *du2, cblas_int_t *ipiv, complex *b, cblas_int_t *ldb, complex *
    x, cblas_int_t *ldx, float *ferr, float *berr, complex *work, float *rwork, 
    cblas_int_t *info);

/* Subroutine */ int cgtsv_(cblas_int_t *n, cblas_int_t *nrhs, complex *dl, complex *
    d__, complex *du, complex *b, cblas_int_t *ldb, cblas_int_t *info);

/* Subroutine */ int cgtsvx_(const char *fact, const char *trans, cblas_int_t *n, cblas_int_t *
    nrhs, complex *dl, complex *d__, complex *du, complex *dlf, complex *
    df, complex *duf, complex *du2, cblas_int_t *ipiv, complex *b, cblas_int_t *
    ldb, complex *x, cblas_int_t *ldx, float *rcond, float *ferr, float *berr, 
    complex *work, float *rwork, cblas_int_t *info);

/* Subroutine */ int cgttrf_(cblas_int_t *n, complex *dl, complex *d__, complex *
    du, complex *du2, cblas_int_t *ipiv, cblas_int_t *info);

/* Subroutine */ int cgttrs_(const char *trans, cblas_int_t *n, cblas_int_t *nrhs, complex *
    dl, complex *d__, complex *du, complex *du2, cblas_int_t *ipiv, complex *
    b, cblas_int_t *ldb, cblas_int_t *info);

/* Subroutine */ int cgtts2_(cblas_int_t *itrans, cblas_int_t *n, cblas_int_t *nrhs, 
    complex *dl, complex *d__, complex *du, complex *du2, cblas_int_t *ipiv, 
    complex *b, cblas_int_t *ldb);

/* Subroutine */ int chbev_(const char *jobz, const char *uplo, cblas_int_t *n, cblas_int_t *kd, 
    complex *ab, cblas_int_t *ldab, float *w, complex *z__, cblas_int_t *ldz, 
    complex *work, float *rwork, cblas_int_t *info);

/* Subroutine */ int chbevd_(const char *jobz, const char *uplo, cblas_int_t *n, cblas_int_t *kd, 
    complex *ab, cblas_int_t *ldab, float *w, complex *z__, cblas_int_t *ldz, 
    complex *work, cblas_int_t *lwork, float *rwork, cblas_int_t *lrwork, cblas_int_t *
    iwork, cblas_int_t *liwork, cblas_int_t *info);

/* Subroutine */ int chbevx_(const char *jobz, const char *range, const char *uplo, cblas_int_t *n, 
    cblas_int_t *kd, complex *ab, cblas_int_t *ldab, complex *q, cblas_int_t *ldq, 
    float *vl, float *vu, cblas_int_t *il, cblas_int_t *iu, float *abstol, cblas_int_t *
    m, float *w, complex *z__, cblas_int_t *ldz, complex *work, float *rwork, 
    cblas_int_t *iwork, cblas_int_t *ifail, cblas_int_t *info);

/* Subroutine */ int chbgst_(const char *vect, const char *uplo, cblas_int_t *n, cblas_int_t *ka, 
    cblas_int_t *kb, complex *ab, cblas_int_t *ldab, complex *bb, cblas_int_t *ldbb, 
    complex *x, cblas_int_t *ldx, complex *work, float *rwork, cblas_int_t *info);

/* Subroutine */ int chbgv_(const char *jobz, const char *uplo, cblas_int_t *n, cblas_int_t *ka, 
    cblas_int_t *kb, complex *ab, cblas_int_t *ldab, complex *bb, cblas_int_t *ldbb, 
    float *w, complex *z__, cblas_int_t *ldz, complex *work, float *rwork, 
    cblas_int_t *info);

/* Subroutine */ int chbgvd_(const char *jobz, const char *uplo, cblas_int_t *n, cblas_int_t *ka, 
    cblas_int_t *kb, complex *ab, cblas_int_t *ldab, complex *bb, cblas_int_t *ldbb, 
    float *w, complex *z__, cblas_int_t *ldz, complex *work, cblas_int_t *lwork, 
    float *rwork, cblas_int_t *lrwork, cblas_int_t *iwork, cblas_int_t *liwork, 
    cblas_int_t *info);

/* Subroutine */ int chbgvx_(const char *jobz, const char *range, const char *uplo, cblas_int_t *n, 
    cblas_int_t *ka, cblas_int_t *kb, complex *ab, cblas_int_t *ldab, complex *bb, 
    cblas_int_t *ldbb, complex *q, cblas_int_t *ldq, float *vl, float *vu, cblas_int_t *
    il, cblas_int_t *iu, float *abstol, cblas_int_t *m, float *w, complex *z__, 
    cblas_int_t *ldz, complex *work, float *rwork, cblas_int_t *iwork, cblas_int_t *
    ifail, cblas_int_t *info);

/* Subroutine */ int chbtrd_(const char *vect, const char *uplo, cblas_int_t *n, cblas_int_t *kd, 
    complex *ab, cblas_int_t *ldab, float *d__, float *e, complex *q, cblas_int_t *
    ldq, complex *work, cblas_int_t *info);

/* Subroutine */ int checon_(const char *uplo, cblas_int_t *n, complex *a, cblas_int_t *lda,
     cblas_int_t *ipiv, float *anorm, float *rcond, complex *work, cblas_int_t *
    info);

/* Subroutine */ int cheev_(const char *jobz, const char *uplo, cblas_int_t *n, complex *a, 
    cblas_int_t *lda, float *w, complex *work, cblas_int_t *lwork, float *rwork, 
    cblas_int_t *info);

/* Subroutine */ int cheevd_(const char *jobz, const char *uplo, cblas_int_t *n, complex *a, 
    cblas_int_t *lda, float *w, complex *work, cblas_int_t *lwork, float *rwork, 
    cblas_int_t *lrwork, cblas_int_t *iwork, cblas_int_t *liwork, cblas_int_t *info);

/* Subroutine */ int cheevr_(const char *jobz, const char *range, const char *uplo, cblas_int_t *n, 
    complex *a, cblas_int_t *lda, float *vl, float *vu, cblas_int_t *il, cblas_int_t *
    iu, float *abstol, cblas_int_t *m, float *w, complex *z__, cblas_int_t *ldz, 
    cblas_int_t *isuppz, complex *work, cblas_int_t *lwork, float *rwork, cblas_int_t *
    lrwork, cblas_int_t *iwork, cblas_int_t *liwork, cblas_int_t *info);

/* Subroutine */ int cheevx_(const char *jobz, const char *range, const char *uplo, cblas_int_t *n, 
    complex *a, cblas_int_t *lda, float *vl, float *vu, cblas_int_t *il, cblas_int_t *
    iu, float *abstol, cblas_int_t *m, float *w, complex *z__, cblas_int_t *ldz, 
    complex *work, cblas_int_t *lwork, float *rwork, cblas_int_t *iwork, cblas_int_t *
    ifail, cblas_int_t *info);

/* Subroutine */ int chegs2_(cblas_int_t *itype, const char *uplo, cblas_int_t *n, complex *
    a, cblas_int_t *lda, complex *b, cblas_int_t *ldb, cblas_int_t *info);

/* Subroutine */ int chegst_(cblas_int_t *itype, const char *uplo, cblas_int_t *n, complex *
    a, cblas_int_t *lda, complex *b, cblas_int_t *ldb, cblas_int_t *info);

/* Subroutine */ int chegv_(cblas_int_t *itype, const char *jobz, const char *uplo, cblas_int_t *
    n, complex *a, cblas_int_t *lda, complex *b, cblas_int_t *ldb, float *w,
    complex *work, cblas_int_t *lwork, float *rwork, cblas_int_t *info);

/* Subroutine */ int chegvd_(cblas_int_t *itype, const char *jobz, const char *uplo, cblas_int_t *
    n, complex *a, cblas_int_t *lda, complex *b, cblas_int_t *ldb, float *w, 
    complex *work, cblas_int_t *lwork, float *rwork, cblas_int_t *lrwork, cblas_int_t *
    iwork, cblas_int_t *liwork, cblas_int_t *info);

/* Subroutine */ int chegvx_(cblas_int_t *itype, const char *jobz, const char *range, const char *
    uplo, cblas_int_t *n, complex *a, cblas_int_t *lda, complex *b, cblas_int_t *ldb, 
    float *vl, float *vu, cblas_int_t *il, cblas_int_t *iu, float *abstol, cblas_int_t *
    m, float *w, complex *z__, cblas_int_t *ldz, complex *work, cblas_int_t *lwork,
     float *rwork, cblas_int_t *iwork, cblas_int_t *ifail, cblas_int_t *info);

/* Subroutine */ int cherfs_(const char *uplo, cblas_int_t *n, cblas_int_t *nrhs, complex *
    a, cblas_int_t *lda, complex *af, cblas_int_t *ldaf, cblas_int_t *ipiv, complex *
    b, cblas_int_t *ldb, complex *x, cblas_int_t *ldx, float *ferr, float *berr, 
    complex *work, float *rwork, cblas_int_t *info);

/* Subroutine */ int chesv_(const char *uplo, cblas_int_t *n, cblas_int_t *nrhs, complex *a,
     cblas_int_t *lda, cblas_int_t *ipiv, complex *b, cblas_int_t *ldb, complex *work,
     cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int chesvx_(const char *fact, const char *uplo, cblas_int_t *n, cblas_int_t *
    nrhs, complex *a, cblas_int_t *lda, complex *af, cblas_int_t *ldaf, cblas_int_t *
    ipiv, complex *b, cblas_int_t *ldb, complex *x, cblas_int_t *ldx, float *rcond,
     float *ferr, float *berr, complex *work, cblas_int_t *lwork, float *rwork, 
    cblas_int_t *info);

/* Subroutine */ int chetd2_(const char *uplo, cblas_int_t *n, complex *a, cblas_int_t *lda,
     float *d__, float *e, complex *tau, cblas_int_t *info);

/* Subroutine */ int chetf2_(const char *uplo, cblas_int_t *n, complex *a, cblas_int_t *lda,
     cblas_int_t *ipiv, cblas_int_t *info);

/* Subroutine */ int chetrd_(const char *uplo, cblas_int_t *n, complex *a, cblas_int_t *lda,
     float *d__, float *e, complex *tau, complex *work, cblas_int_t *lwork, 
    cblas_int_t *info);

/* Subroutine */ int chetrf_(const char *uplo, cblas_int_t *n, complex *a, cblas_int_t *lda,
     cblas_int_t *ipiv, complex *work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int chetri_(const char *uplo, cblas_int_t *n, complex *a, cblas_int_t *lda,
     cblas_int_t *ipiv, complex *work, cblas_int_t *info);

/* Subroutine */ int chetrs_(const char *uplo, cblas_int_t *n, cblas_int_t *nrhs, complex *
    a, cblas_int_t *lda, cblas_int_t *ipiv, complex *b, cblas_int_t *ldb, cblas_int_t *
    info);

/* Subroutine */ int chgeqz_(const char *job, const char *compq, const char *compz, cblas_int_t *n, 
    cblas_int_t *ilo, cblas_int_t *ihi, complex *h__, cblas_int_t *ldh, complex *t, 
    cblas_int_t *ldt, complex *alpha, complex *beta, complex *q, cblas_int_t *ldq,
     complex *z__, cblas_int_t *ldz, complex *work, cblas_int_t *lwork, float *
    rwork, cblas_int_t *info);

/* Subroutine */ int chpcon_(const char *uplo, cblas_int_t *n, complex *ap, cblas_int_t *
    ipiv, float *anorm, float *rcond, complex *work, cblas_int_t *info);

/* Subroutine */ int chpev_(const char *jobz, const char *uplo, cblas_int_t *n, complex *ap, 
    float *w, complex *z__, cblas_int_t *ldz, complex *work, float *rwork, 
    cblas_int_t *info);

/* Subroutine */ int chpevd_(const char *jobz, const char *uplo, cblas_int_t *n, complex *ap, 
    float *w, complex *z__, cblas_int_t *ldz, complex *work, cblas_int_t *lwork, 
    float *rwork, cblas_int_t *lrwork, cblas_int_t *iwork, cblas_int_t *liwork, 
    cblas_int_t *info);

/* Subroutine */ int chpevx_(const char *jobz, const char *range, const char *uplo, cblas_int_t *n, 
    complex *ap, float *vl, float *vu, cblas_int_t *il, cblas_int_t *iu, float *
    abstol, cblas_int_t *m, float *w, complex *z__, cblas_int_t *ldz, complex *
    work, float *rwork, cblas_int_t *iwork, cblas_int_t *ifail, cblas_int_t *info);

/* Subroutine */ int chpgst_(cblas_int_t *itype, const char *uplo, cblas_int_t *n, complex *
    ap, complex *bp, cblas_int_t *info);

/* Subroutine */ int chpgv_(cblas_int_t *itype, const char *jobz, const char *uplo, cblas_int_t *
    n, complex *ap, complex *bp, float *w, complex *z__, cblas_int_t *ldz, 
    complex *work, float *rwork, cblas_int_t *info);

/* Subroutine */ int chpgvd_(cblas_int_t *itype, const char *jobz, const char *uplo, cblas_int_t *
    n, complex *ap, complex *bp, float *w, complex *z__, cblas_int_t *ldz, 
    complex *work, cblas_int_t *lwork, float *rwork, cblas_int_t *lrwork, cblas_int_t *
    iwork, cblas_int_t *liwork, cblas_int_t *info);

/* Subroutine */ int chpgvx_(cblas_int_t *itype, const char *jobz, const char *range, const char *
    uplo, cblas_int_t *n, complex *ap, complex *bp, float *vl, float *vu, 
    cblas_int_t *il, cblas_int_t *iu, float *abstol, cblas_int_t *m, float *w, complex *
    z__, cblas_int_t *ldz, complex *work, float *rwork, cblas_int_t *iwork, 
    cblas_int_t *ifail, cblas_int_t *info);

/* Subroutine */ int chprfs_(const char *uplo, cblas_int_t *n, cblas_int_t *nrhs, complex *
    ap, complex *afp, cblas_int_t *ipiv, complex *b, cblas_int_t *ldb, complex *x,
     cblas_int_t *ldx, float *ferr, float *berr, complex *work, float *rwork, 
    cblas_int_t *info);

/* Subroutine */ int chpsv_(const char *uplo, cblas_int_t *n, cblas_int_t *nrhs, complex *
    ap, cblas_int_t *ipiv, complex *b, cblas_int_t *ldb, cblas_int_t *info);

/* Subroutine */ int chpsvx_(const char *fact, const char *uplo, cblas_int_t *n, cblas_int_t *
    nrhs, complex *ap, complex *afp, cblas_int_t *ipiv, complex *b, cblas_int_t *
    ldb, complex *x, cblas_int_t *ldx, float *rcond, float *ferr, float *berr, 
    complex *work, float *rwork, cblas_int_t *info);

/* Subroutine */ int chptrd_(const char *uplo, cblas_int_t *n, complex *ap, float *d__, 
    float *e, complex *tau, cblas_int_t *info);

/* Subroutine */ int chptrf_(const char *uplo, cblas_int_t *n, complex *ap, cblas_int_t *
    ipiv, cblas_int_t *info);

/* Subroutine */ int chptri_(const char *uplo, cblas_int_t *n, complex *ap, cblas_int_t *
    ipiv, complex *work, cblas_int_t *info);

/* Subroutine */ int chptrs_(const char *uplo, cblas_int_t *n, cblas_int_t *nrhs, complex *
    ap, cblas_int_t *ipiv, complex *b, cblas_int_t *ldb, cblas_int_t *info);

/* Subroutine */ int chsein_(const char *side, const char *eigsrc, const char *initv, logical *
    select, cblas_int_t *n, complex *h__, cblas_int_t *ldh, complex *w, complex *
    vl, cblas_int_t *ldvl, complex *vr, cblas_int_t *ldvr, cblas_int_t *mm, cblas_int_t *
    m, complex *work, float *rwork, cblas_int_t *ifaill, cblas_int_t *ifailr, 
    cblas_int_t *info);

/* Subroutine */ int chseqr_(const char *job, const char *compz, cblas_int_t *n, cblas_int_t *ilo,
    cblas_int_t *ihi, complex *h__, cblas_int_t *ldh, complex *w, complex *z__,
    cblas_int_t *ldz, complex *work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int clabrd_(cblas_int_t *m, cblas_int_t *n, cblas_int_t *nb, complex *a, 
    cblas_int_t *lda, float *d__, float *e, complex *tauq, complex *taup, 
    complex *x, cblas_int_t *ldx, complex *y, cblas_int_t *ldy);

/* Subroutine */ int clacgv_(cblas_int_t *n, complex *x, cblas_int_t *incx);

/* Subroutine */ int clacn2_(cblas_int_t *n, complex *v, complex *x, float *est, 
    cblas_int_t *kase, cblas_int_t *isave);

/* Subroutine */ int clacon_(cblas_int_t *n, complex *v, complex *x, float *est, 
    cblas_int_t *kase);

/* Subroutine */ int clacp2_(const char *uplo, cblas_int_t *m, cblas_int_t *n, float *a, 
    cblas_int_t *lda, complex *b, cblas_int_t *ldb);

/* Subroutine */ int clacpy_(const char *uplo, cblas_int_t *m, cblas_int_t *n, complex *a, 
    cblas_int_t *lda, complex *b, cblas_int_t *ldb);

/* Subroutine */ int clacrm_(cblas_int_t *m, cblas_int_t *n, complex *a, cblas_int_t *lda,
     float *b, cblas_int_t *ldb, complex *c__, cblas_int_t *ldc, float *rwork);

/* Subroutine */ int clacrt_(cblas_int_t *n, complex *cx, cblas_int_t *incx, complex *
    cy, cblas_int_t *incy, complex *c__, complex *s);

/* Complex */ VOID cladiv_(complex * ret_val, complex *x, complex *y);

/* Subroutine */ int claed0_(cblas_int_t *qsiz, cblas_int_t *n, float *d__, float *e, 
    complex *q, cblas_int_t *ldq, complex *qstore, cblas_int_t *ldqs, float *rwork,
     cblas_int_t *iwork, cblas_int_t *info);

/* Subroutine */ int claed7_(cblas_int_t *n, cblas_int_t *cutpnt, cblas_int_t *qsiz, 
    cblas_int_t *tlvls, cblas_int_t *curlvl, cblas_int_t *curpbm, float *d__, complex *
    q, cblas_int_t *ldq, float *rho, cblas_int_t *indxq, float *qstore, cblas_int_t *
    qptr, cblas_int_t *prmptr, cblas_int_t *perm, cblas_int_t *givptr, cblas_int_t *
    givcol, float *givnum, complex *work, float *rwork, cblas_int_t *iwork, 
    cblas_int_t *info);

/* Subroutine */ int claed8_(cblas_int_t *k, cblas_int_t *n, cblas_int_t *qsiz, complex *
    q, cblas_int_t *ldq, float *d__, float *rho, cblas_int_t *cutpnt, float *z__, 
    float *dlamda, complex *q2, cblas_int_t *ldq2, float *w, cblas_int_t *indxp, 
    cblas_int_t *indx, cblas_int_t *indxq, cblas_int_t *perm, cblas_int_t *givptr, 
    cblas_int_t *givcol, float *givnum, cblas_int_t *info);

/* Subroutine */ int claein_(logical *rightv, logical *noinit, cblas_int_t *n, 
    complex *h__, cblas_int_t *ldh, complex *w, complex *v, complex *b, 
    cblas_int_t *ldb, float *rwork, float *eps3, float *smlnum, cblas_int_t *info);

/* Subroutine */ int claesy_(complex *a, complex *b, complex *c__, complex *
    rt1, complex *rt2, complex *evscal, complex *cs1, complex *sn1);

/* Subroutine */ int claev2_(complex *a, complex *b, complex *c__, float *rt1, 
    float *rt2, float *cs1, complex *sn1);

/* Subroutine */ int clag2z_(cblas_int_t *m, cblas_int_t *n, complex *sa, cblas_int_t *
    ldsa, doublecomplex *a, cblas_int_t *lda, cblas_int_t *info);

/* Subroutine */ int clags2_(logical *upper, float *a1, complex *a2, float *a3, 
    float *b1, complex *b2, float *b3, float *csu, complex *snu, float *csv, 
    complex *snv, float *csq, complex *snq);

/* Subroutine */ int clagtm_(const char *trans, cblas_int_t *n, cblas_int_t *nrhs, float *
    alpha, complex *dl, complex *d__, complex *du, complex *x, cblas_int_t *
    ldx, float *beta, complex *b, cblas_int_t *ldb);

/* Subroutine */ int clahef_(const char *uplo, cblas_int_t *n, cblas_int_t *nb, cblas_int_t *kb,
     complex *a, cblas_int_t *lda, cblas_int_t *ipiv, complex *w, cblas_int_t *ldw, 
    cblas_int_t *info);

/* Subroutine */ int clahqr_(logical *wantt, logical *wantz, cblas_int_t *n, 
    cblas_int_t *ilo, cblas_int_t *ihi, complex *h__, cblas_int_t *ldh, complex *w, 
    cblas_int_t *iloz, cblas_int_t *ihiz, complex *z__, cblas_int_t *ldz, cblas_int_t *
    info);

/* Subroutine */ int clahr2_(cblas_int_t *n, cblas_int_t *k, cblas_int_t *nb, complex *a, 
    cblas_int_t *lda, complex *tau, complex *t, cblas_int_t *ldt, complex *y, 
    cblas_int_t *ldy);

/* Subroutine */ int clahrd_(cblas_int_t *n, cblas_int_t *k, cblas_int_t *nb, complex *a, 
    cblas_int_t *lda, complex *tau, complex *t, cblas_int_t *ldt, complex *y, 
    cblas_int_t *ldy);

/* Subroutine */ int claic1_(cblas_int_t *job, cblas_int_t *j, complex *x, float *sest,
     complex *w, complex *gamma, float *sestpr, complex *s, complex *c__);

/* Subroutine */ int clals0_(cblas_int_t *icompq, cblas_int_t *nl, cblas_int_t *nr, 
    cblas_int_t *sqre, cblas_int_t *nrhs, complex *b, cblas_int_t *ldb, complex *bx, 
    cblas_int_t *ldbx, cblas_int_t *perm, cblas_int_t *givptr, cblas_int_t *givcol, 
    cblas_int_t *ldgcol, float *givnum, cblas_int_t *ldgnum, float *poles, float *
    difl, float *difr, float *z__, cblas_int_t *k, float *c__, float *s, float *
    rwork, cblas_int_t *info);

/* Subroutine */ int clalsa_(cblas_int_t *icompq, cblas_int_t *smlsiz, cblas_int_t *n, 
    cblas_int_t *nrhs, complex *b, cblas_int_t *ldb, complex *bx, cblas_int_t *ldbx, 
    float *u, cblas_int_t *ldu, float *vt, cblas_int_t *k, float *difl, float *difr, 
    float *z__, float *poles, cblas_int_t *givptr, cblas_int_t *givcol, cblas_int_t *
    ldgcol, cblas_int_t *perm, float *givnum, float *c__, float *s, float *rwork, 
    cblas_int_t *iwork, cblas_int_t *info);

/* Subroutine */ int clalsd_(const char *uplo, cblas_int_t *smlsiz, cblas_int_t *n, cblas_int_t 
    *nrhs, float *d__, float *e, complex *b, cblas_int_t *ldb, float *rcond, 
    cblas_int_t *rank, complex *work, float *rwork, cblas_int_t *iwork, cblas_int_t *
    info);

/* Subroutine */ int clapll_(cblas_int_t *n, complex *x, cblas_int_t *incx, complex *
    y, cblas_int_t *incy, float *ssmin);

/* Subroutine */ int clapmt_(logical *forwrd, cblas_int_t *m, cblas_int_t *n, complex 
    *x, cblas_int_t *ldx, cblas_int_t *k);

/* Subroutine */ int claqgb_(cblas_int_t *m, cblas_int_t *n, cblas_int_t *kl, cblas_int_t *ku,
     complex *ab, cblas_int_t *ldab, float *r__, float *c__, float *rowcnd, float 
    *colcnd, float *amax, const char *equed);

/* Subroutine */ int claqge_(cblas_int_t *m, cblas_int_t *n, complex *a, cblas_int_t *lda,
     float *r__, float *c__, float *rowcnd, float *colcnd, float *amax, const char *
    equed);

/* Subroutine */ int claqhb_(const char *uplo, cblas_int_t *n, cblas_int_t *kd, complex *ab,
     cblas_int_t *ldab, float *s, float *scond, float *amax, const char *equed      );

/* Subroutine */ int claqhe_(const char *uplo, cblas_int_t *n, complex *a, cblas_int_t *lda,
     float *s, float *scond, float *amax, const char *equed);

/* Subroutine */ int claqhp_(const char *uplo, cblas_int_t *n, complex *ap, float *s, 
    float *scond, float *amax, const char *equed    );

/* Subroutine */ int claqp2_(cblas_int_t *m, cblas_int_t *n, cblas_int_t *offset, complex 
    *a, cblas_int_t *lda, cblas_int_t *jpvt, complex *tau, float *vn1, float *vn2, 
    complex *work);

/* Subroutine */ int claqps_(cblas_int_t *m, cblas_int_t *n, cblas_int_t *offset, cblas_int_t 
    *nb, cblas_int_t *kb, complex *a, cblas_int_t *lda, cblas_int_t *jpvt, complex *
    tau, float *vn1, float *vn2, complex *auxv, complex *f, cblas_int_t *ldf);

/* Subroutine */ int claqr0_(logical *wantt, logical *wantz, cblas_int_t *n, 
    cblas_int_t *ilo, cblas_int_t *ihi, complex *h__, cblas_int_t *ldh, complex *w, 
    cblas_int_t *iloz, cblas_int_t *ihiz, complex *z__, cblas_int_t *ldz, complex *
    work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int claqr1_(cblas_int_t *n, complex *h__, cblas_int_t *ldh, complex *
    s1, complex *s2, complex *v);

/* Subroutine */ int claqr2_(logical *wantt, logical *wantz, cblas_int_t *n, 
    cblas_int_t *ktop, cblas_int_t *kbot, cblas_int_t *nw, complex *h__, cblas_int_t *ldh,
     cblas_int_t *iloz, cblas_int_t *ihiz, complex *z__, cblas_int_t *ldz, cblas_int_t *
    ns, cblas_int_t *nd, complex *sh, complex *v, cblas_int_t *ldv, cblas_int_t *nh, 
    complex *t, cblas_int_t *ldt, cblas_int_t *nv, complex *wv, cblas_int_t *ldwv, 
    complex *work, cblas_int_t *lwork);

/* Subroutine */ int claqr3_(logical *wantt, logical *wantz, cblas_int_t *n, 
    cblas_int_t *ktop, cblas_int_t *kbot, cblas_int_t *nw, complex *h__, cblas_int_t *ldh,
     cblas_int_t *iloz, cblas_int_t *ihiz, complex *z__, cblas_int_t *ldz, cblas_int_t *
    ns, cblas_int_t *nd, complex *sh, complex *v, cblas_int_t *ldv, cblas_int_t *nh, 
    complex *t, cblas_int_t *ldt, cblas_int_t *nv, complex *wv, cblas_int_t *ldwv, 
    complex *work, cblas_int_t *lwork);

/* Subroutine */ int claqr4_(logical *wantt, logical *wantz, cblas_int_t *n, 
    cblas_int_t *ilo, cblas_int_t *ihi, complex *h__, cblas_int_t *ldh, complex *w, 
    cblas_int_t *iloz, cblas_int_t *ihiz, complex *z__, cblas_int_t *ldz, complex *
    work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int claqr5_(logical *wantt, logical *wantz, cblas_int_t *kacc22, 
    cblas_int_t *n, cblas_int_t *ktop, cblas_int_t *kbot, cblas_int_t *nshfts, complex *s,
     complex *h__, cblas_int_t *ldh, cblas_int_t *iloz, cblas_int_t *ihiz, complex *
    z__, cblas_int_t *ldz, complex *v, cblas_int_t *ldv, complex *u, cblas_int_t *ldu,
     cblas_int_t *nv, complex *wv, cblas_int_t *ldwv, cblas_int_t *nh, complex *wh, 
    cblas_int_t *ldwh);

/* Subroutine */ int claqsb_(const char *uplo, cblas_int_t *n, cblas_int_t *kd, complex *ab,
     cblas_int_t *ldab, float *s, float *scond, float *amax, const char *equed      );

/* Subroutine */ int claqsp_(const char *uplo, cblas_int_t *n, complex *ap, float *s, 
    float *scond, float *amax, const char *equed    );

/* Subroutine */ int claqsy_(const char *uplo, cblas_int_t *n, complex *a, cblas_int_t *lda,
     float *s, float *scond, float *amax, const char *equed);

/* Subroutine */ int clar1v_(cblas_int_t *n, cblas_int_t *b1, cblas_int_t *bn, float *
    lambda, float *d__, float *l, float *ld, float *lld, float *pivmin, float *
    gaptol, complex *z__, logical *wantnc, cblas_int_t *negcnt, float *ztz, 
    float *mingma, cblas_int_t *r__, cblas_int_t *isuppz, float *nrminv, float *
    resid, float *rqcorr, float *work);

/* Subroutine */ int clar2v_(cblas_int_t *n, complex *x, complex *y, complex *z__,
     cblas_int_t *incx, float *c__, complex *s, cblas_int_t *incc);

/* Subroutine */ int clarcm_(cblas_int_t *m, cblas_int_t *n, float *a, cblas_int_t *lda, 
    complex *b, cblas_int_t *ldb, complex *c__, cblas_int_t *ldc, float *rwork);

/* Subroutine */ int clarf_(const char *side, cblas_int_t *m, cblas_int_t *n, complex *v, 
    cblas_int_t *incv, complex *tau, complex *c__, cblas_int_t *ldc, complex *
    work);

/* Subroutine */ int clarfb_(const char *side, const char *trans, const char *direct, const char *
    storev, cblas_int_t *m, cblas_int_t *n, cblas_int_t *k, complex *v, cblas_int_t *ldv, 
    complex *t, cblas_int_t *ldt, complex *c__, cblas_int_t *ldc, complex *work, 
    cblas_int_t *ldwork);

/* Subroutine */ int clarfg_(cblas_int_t *n, complex *alpha, complex *x, cblas_int_t *
    incx, complex *tau);

/* Subroutine */ int clarft_(const char *direct, const char *storev, cblas_int_t *n, cblas_int_t *
    k, complex *v, cblas_int_t *ldv, complex *tau, complex *t, cblas_int_t *ldt);

/* Subroutine */ int clarfx_(const char *side, cblas_int_t *m, cblas_int_t *n, complex *v, 
    complex *tau, complex *c__, cblas_int_t *ldc, complex *work     );

/* Subroutine */ int clargv_(cblas_int_t *n, complex *x, cblas_int_t *incx, complex *
    y, cblas_int_t *incy, float *c__, cblas_int_t *incc);

/* Subroutine */ int clarnv_(cblas_int_t *idist, cblas_int_t *iseed, cblas_int_t *n, 
    complex *x);

/* Subroutine */ int clarrv_(cblas_int_t *n, float *vl, float *vu, float *d__, float *
    l, float *pivmin, cblas_int_t *isplit, cblas_int_t *m, cblas_int_t *dol, cblas_int_t *
    dou, float *minrgp, float *rtol1, float *rtol2, float *w, float *werr, 
    float *wgap, cblas_int_t *iblock, cblas_int_t *indexw, float *gers, complex *
    z__, cblas_int_t *ldz, cblas_int_t *isuppz, float *work, cblas_int_t *iwork, 
    cblas_int_t *info);

/* Subroutine */ int clartg_(complex *f, complex *g, float *cs, complex *sn, 
    complex *r__);

/* Subroutine */ int clartv_(cblas_int_t *n, complex *x, cblas_int_t *incx, complex *
    y, cblas_int_t *incy, float *c__, complex *s, cblas_int_t *incc);

/* Subroutine */ int clarz_(const char *side, cblas_int_t *m, cblas_int_t *n, cblas_int_t *l, 
    complex *v, cblas_int_t *incv, complex *tau, complex *c__, cblas_int_t *ldc, 
    complex *work);

/* Subroutine */ int clarzb_(const char *side, const char *trans, const char *direct, const char *
    storev, cblas_int_t *m, cblas_int_t *n, cblas_int_t *k, cblas_int_t *l, complex *v, 
    cblas_int_t *ldv, complex *t, cblas_int_t *ldt, complex *c__, cblas_int_t *ldc, 
    complex *work, cblas_int_t *ldwork);

/* Subroutine */ int clarzt_(const char *direct, const char *storev, cblas_int_t *n, cblas_int_t *
    k, complex *v, cblas_int_t *ldv, complex *tau, complex *t, cblas_int_t *ldt);

/* Subroutine */ int clascl_(const char *type__, cblas_int_t *kl, cblas_int_t *ku, float *
    cfrom, float *cto, cblas_int_t *m, cblas_int_t *n, complex *a, cblas_int_t *lda, 
    cblas_int_t *info);

/* Subroutine */ int claset_(const char *uplo, cblas_int_t *m, cblas_int_t *n, complex *
    alpha, complex *beta, complex *a, cblas_int_t *lda);

/* Subroutine */ int clasr_(const char *side, const char *pivot, const char *direct, cblas_int_t *m,
     cblas_int_t *n, float *c__, float *s, complex *a, cblas_int_t *lda     );

/* Subroutine */ int classq_(cblas_int_t *n, complex *x, cblas_int_t *incx, float *
    scale, float *sumsq);

/* Subroutine */ int claswp_(cblas_int_t *n, complex *a, cblas_int_t *lda, cblas_int_t *
    k1, cblas_int_t *k2, cblas_int_t *ipiv, cblas_int_t *incx);

/* Subroutine */ int clasyf_(const char *uplo, cblas_int_t *n, cblas_int_t *nb, cblas_int_t *kb,
     complex *a, cblas_int_t *lda, cblas_int_t *ipiv, complex *w, cblas_int_t *ldw, 
    cblas_int_t *info);

/* Subroutine */ int clatbs_(const char *uplo, const char *trans, const char *diag, const char *
    normin, cblas_int_t *n, cblas_int_t *kd, complex *ab, cblas_int_t *ldab, complex *
    x, float *scale, float *cnorm, cblas_int_t *info);

/* Subroutine */ int clatdf_(cblas_int_t *ijob, cblas_int_t *n, complex *z__, cblas_int_t 
    *ldz, complex *rhs, float *rdsum, float *rdscal, cblas_int_t *ipiv, cblas_int_t 
    *jpiv);

/* Subroutine */ int clatps_(const char *uplo, const char *trans, const char *diag, const char *
    normin, cblas_int_t *n, complex *ap, complex *x, float *scale, float *cnorm,
     cblas_int_t *info);

/* Subroutine */ int clatrd_(const char *uplo, cblas_int_t *n, cblas_int_t *nb, complex *a, 
    cblas_int_t *lda, float *e, complex *tau, complex *w, cblas_int_t *ldw      );

/* Subroutine */ int clatrs_(const char *uplo, const char *trans, const char *diag, const char *
    normin, cblas_int_t *n, complex *a, cblas_int_t *lda, complex *x, float *scale,
     float *cnorm, cblas_int_t *info);

/* Subroutine */ int clatrz_(cblas_int_t *m, cblas_int_t *n, cblas_int_t *l, complex *a, 
    cblas_int_t *lda, complex *tau, complex *work);

/* Subroutine */ int clatzm_(const char *side, cblas_int_t *m, cblas_int_t *n, complex *v, 
    cblas_int_t *incv, complex *tau, complex *c1, complex *c2, cblas_int_t *ldc, 
    complex *work);

/* Subroutine */ int clauu2_(const char *uplo, cblas_int_t *n, complex *a, cblas_int_t *lda,
     cblas_int_t *info);

/* Subroutine */ int clauum_(const char *uplo, cblas_int_t *n, complex *a, cblas_int_t *lda,
     cblas_int_t *info);

/* Subroutine */ int cpbcon_(const char *uplo, cblas_int_t *n, cblas_int_t *kd, complex *ab,
     cblas_int_t *ldab, float *anorm, float *rcond, complex *work, float *rwork, 
    cblas_int_t *info);

/* Subroutine */ int cpbequ_(const char *uplo, cblas_int_t *n, cblas_int_t *kd, complex *ab,
     cblas_int_t *ldab, float *s, float *scond, float *amax, cblas_int_t *info);

/* Subroutine */ int cpbrfs_(const char *uplo, cblas_int_t *n, cblas_int_t *kd, cblas_int_t *
    nrhs, complex *ab, cblas_int_t *ldab, complex *afb, cblas_int_t *ldafb, 
    complex *b, cblas_int_t *ldb, complex *x, cblas_int_t *ldx, float *ferr, float *
    berr, complex *work, float *rwork, cblas_int_t *info);

/* Subroutine */ int cpbstf_(const char *uplo, cblas_int_t *n, cblas_int_t *kd, complex *ab,
     cblas_int_t *ldab, cblas_int_t *info);

/* Subroutine */ int cpbsv_(const char *uplo, cblas_int_t *n, cblas_int_t *kd, cblas_int_t *
    nrhs, complex *ab, cblas_int_t *ldab, complex *b, cblas_int_t *ldb, cblas_int_t *
    info);

/* Subroutine */ int cpbsvx_(const char *fact, const char *uplo, cblas_int_t *n, cblas_int_t *kd, 
    cblas_int_t *nrhs, complex *ab, cblas_int_t *ldab, complex *afb, cblas_int_t *
    ldafb, const char *equed, float *s, complex *b, cblas_int_t *ldb, complex *x, 
    cblas_int_t *ldx, float *rcond, float *ferr, float *berr, complex *work, 
    float *rwork, cblas_int_t *info);

/* Subroutine */ int cpbtf2_(const char *uplo, cblas_int_t *n, cblas_int_t *kd, complex *ab,
     cblas_int_t *ldab, cblas_int_t *info);

/* Subroutine */ int cpbtrf_(const char *uplo, cblas_int_t *n, cblas_int_t *kd, complex *ab,
     cblas_int_t *ldab, cblas_int_t *info);

/* Subroutine */ int cpbtrs_(const char *uplo, cblas_int_t *n, cblas_int_t *kd, cblas_int_t *
    nrhs, complex *ab, cblas_int_t *ldab, complex *b, cblas_int_t *ldb, cblas_int_t *
    info);

/* Subroutine */ int cpocon_(const char *uplo, cblas_int_t *n, complex *a, cblas_int_t *lda,
     float *anorm, float *rcond, complex *work, float *rwork, cblas_int_t *info);

/* Subroutine */ int cpoequ_(cblas_int_t *n, complex *a, cblas_int_t *lda, float *s, 
    float *scond, float *amax, cblas_int_t *info);

/* Subroutine */ int cporfs_(const char *uplo, cblas_int_t *n, cblas_int_t *nrhs, complex *
    a, cblas_int_t *lda, complex *af, cblas_int_t *ldaf, complex *b, cblas_int_t *ldb,
     complex *x, cblas_int_t *ldx, float *ferr, float *berr, complex *work, 
    float *rwork, cblas_int_t *info);

/* Subroutine */ int cposv_(const char *uplo, cblas_int_t *n, cblas_int_t *nrhs, complex *a,
     cblas_int_t *lda, complex *b, cblas_int_t *ldb, cblas_int_t *info);

/* Subroutine */ int cposvx_(const char *fact, const char *uplo, cblas_int_t *n, cblas_int_t *
    nrhs, complex *a, cblas_int_t *lda, complex *af, cblas_int_t *ldaf, const char *
    equed, float *s, complex *b, cblas_int_t *ldb, complex *x, cblas_int_t *ldx, 
    float *rcond, float *ferr, float *berr, complex *work, float *rwork, 
    cblas_int_t *info);

/* Subroutine */ int cpotf2_(const char *uplo, cblas_int_t *n, complex *a, cblas_int_t *lda,
     cblas_int_t *info);

/* Subroutine */ int cpotrf_(const char *uplo, cblas_int_t *n, complex *a, cblas_int_t *lda,
     cblas_int_t *info);

/* Subroutine */ int cpotri_(const char *uplo, cblas_int_t *n, complex *a, cblas_int_t *lda,
     cblas_int_t *info);

/* Subroutine */ int cpotrs_(const char *uplo, cblas_int_t *n, cblas_int_t *nrhs, complex *
    a, cblas_int_t *lda, complex *b, cblas_int_t *ldb, cblas_int_t *info);

/* Subroutine */ int cppcon_(const char *uplo, cblas_int_t *n, complex *ap, float *anorm,
     float *rcond, complex *work, float *rwork, cblas_int_t *info);

/* Subroutine */ int cppequ_(const char *uplo, cblas_int_t *n, complex *ap, float *s, 
    float *scond, float *amax, cblas_int_t *info);

/* Subroutine */ int cpprfs_(const char *uplo, cblas_int_t *n, cblas_int_t *nrhs, complex *
    ap, complex *afp, complex *b, cblas_int_t *ldb, complex *x, cblas_int_t *ldx, 
    float *ferr, float *berr, complex *work, float *rwork, cblas_int_t *info);

/* Subroutine */ int cppsv_(const char *uplo, cblas_int_t *n, cblas_int_t *nrhs, complex *
    ap, complex *b, cblas_int_t *ldb, cblas_int_t *info);

/* Subroutine */ int cppsvx_(const char *fact, const char *uplo, cblas_int_t *n, cblas_int_t *
    nrhs, complex *ap, complex *afp, const char *equed, float *s, complex *b, 
    cblas_int_t *ldb, complex *x, cblas_int_t *ldx, float *rcond, float *ferr, float 
    *berr, complex *work, float *rwork, cblas_int_t *info);

/* Subroutine */ int cpptrf_(const char *uplo, cblas_int_t *n, complex *ap, cblas_int_t *
    info);

/* Subroutine */ int cpptri_(const char *uplo, cblas_int_t *n, complex *ap, cblas_int_t *
    info);

/* Subroutine */ int cpptrs_(const char *uplo, cblas_int_t *n, cblas_int_t *nrhs, complex *
    ap, complex *b, cblas_int_t *ldb, cblas_int_t *info);

/* Subroutine */ int cptcon_(cblas_int_t *n, float *d__, complex *e, float *anorm, 
    float *rcond, float *rwork, cblas_int_t *info);

/* Subroutine */ int cpteqr_(const char *compz, cblas_int_t *n, float *d__, float *e, 
    complex *z__, cblas_int_t *ldz, float *work, cblas_int_t *info);

/* Subroutine */ int cptrfs_(const char *uplo, cblas_int_t *n, cblas_int_t *nrhs, float *d__,
     complex *e, float *df, complex *ef, complex *b, cblas_int_t *ldb, complex 
    *x, cblas_int_t *ldx, float *ferr, float *berr, complex *work, float *rwork, 
    cblas_int_t *info);

/* Subroutine */ int cptsv_(cblas_int_t *n, cblas_int_t *nrhs, float *d__, complex *e, 
    complex *b, cblas_int_t *ldb, cblas_int_t *info);

/* Subroutine */ int cptsvx_(const char *fact, cblas_int_t *n, cblas_int_t *nrhs, float *d__,
     complex *e, float *df, complex *ef, complex *b, cblas_int_t *ldb, complex 
    *x, cblas_int_t *ldx, float *rcond, float *ferr, float *berr, complex *work, 
    float *rwork, cblas_int_t *info);

/* Subroutine */ int cpttrf_(cblas_int_t *n, float *d__, complex *e, cblas_int_t *info);

/* Subroutine */ int cpttrs_(const char *uplo, cblas_int_t *n, cblas_int_t *nrhs, float *d__,
     complex *e, complex *b, cblas_int_t *ldb, cblas_int_t *info);

/* Subroutine */ int cptts2_(cblas_int_t *iuplo, cblas_int_t *n, cblas_int_t *nrhs, float *
    d__, complex *e, complex *b, cblas_int_t *ldb);

/* Subroutine */ int crot_(cblas_int_t *n, complex *cx, cblas_int_t *incx, complex *
    cy, cblas_int_t *incy, float *c__, complex *s);

/* Subroutine */ int cspcon_(const char *uplo, cblas_int_t *n, complex *ap, cblas_int_t *
    ipiv, float *anorm, float *rcond, complex *work, cblas_int_t *info);

/* Subroutine */ int cspmv_(const char *uplo, cblas_int_t *n, complex *alpha, complex *
    ap, complex *x, cblas_int_t *incx, complex *beta, complex *y, cblas_int_t *
    incy);

/* Subroutine */ int cspr_(const char *uplo, cblas_int_t *n, complex *alpha, complex *x,
     cblas_int_t *incx, complex *ap);

/* Subroutine */ int csprfs_(const char *uplo, cblas_int_t *n, cblas_int_t *nrhs, complex *
    ap, complex *afp, cblas_int_t *ipiv, complex *b, cblas_int_t *ldb, complex *x,
     cblas_int_t *ldx, float *ferr, float *berr, complex *work, float *rwork, 
    cblas_int_t *info);

/* Subroutine */ int cspsv_(const char *uplo, cblas_int_t *n, cblas_int_t *nrhs, complex *
    ap, cblas_int_t *ipiv, complex *b, cblas_int_t *ldb, cblas_int_t *info);

/* Subroutine */ int cspsvx_(const char *fact, const char *uplo, cblas_int_t *n, cblas_int_t *
    nrhs, complex *ap, complex *afp, cblas_int_t *ipiv, complex *b, cblas_int_t *
    ldb, complex *x, cblas_int_t *ldx, float *rcond, float *ferr, float *berr, 
    complex *work, float *rwork, cblas_int_t *info);

/* Subroutine */ int csptrf_(const char *uplo, cblas_int_t *n, complex *ap, cblas_int_t *
    ipiv, cblas_int_t *info);

/* Subroutine */ int csptri_(const char *uplo, cblas_int_t *n, complex *ap, cblas_int_t *
    ipiv, complex *work, cblas_int_t *info);

/* Subroutine */ int csptrs_(const char *uplo, cblas_int_t *n, cblas_int_t *nrhs, complex *
    ap, cblas_int_t *ipiv, complex *b, cblas_int_t *ldb, cblas_int_t *info);

/* Subroutine */ int csrscl_(cblas_int_t *n, float *sa, complex *sx, cblas_int_t *incx);

/* Subroutine */ int cstedc_(const char *compz, cblas_int_t *n, float *d__, float *e, 
    complex *z__, cblas_int_t *ldz, complex *work, cblas_int_t *lwork, float *
    rwork, cblas_int_t *lrwork, cblas_int_t *iwork, cblas_int_t *liwork, cblas_int_t *
    info);

/* Subroutine */ int cstegr_(const char *jobz, const char *range, cblas_int_t *n, float *d__, 
    float *e, float *vl, float *vu, cblas_int_t *il, cblas_int_t *iu, float *abstol, 
    cblas_int_t *m, float *w, complex *z__, cblas_int_t *ldz, cblas_int_t *isuppz, 
    float *work, cblas_int_t *lwork, cblas_int_t *iwork, cblas_int_t *liwork, cblas_int_t *
    info);

/* Subroutine */ int cstein_(cblas_int_t *n, float *d__, float *e, cblas_int_t *m, float 
    *w, cblas_int_t *iblock, cblas_int_t *isplit, complex *z__, cblas_int_t *ldz, 
    float *work, cblas_int_t *iwork, cblas_int_t *ifail, cblas_int_t *info);

/* Subroutine */ int cstemr_(const char *jobz, const char *range, cblas_int_t *n, float *d__, 
    float *e, float *vl, float *vu, cblas_int_t *il, cblas_int_t *iu, cblas_int_t *m, 
    float *w, complex *z__, cblas_int_t *ldz, cblas_int_t *nzc, cblas_int_t *isuppz, 
    logical *tryrac, float *work, cblas_int_t *lwork, cblas_int_t *iwork, cblas_int_t *
    liwork, cblas_int_t *info);

/* Subroutine */ int csteqr_(const char *compz, cblas_int_t *n, float *d__, float *e, 
    complex *z__, cblas_int_t *ldz, float *work, cblas_int_t *info);

/* Subroutine */ int csycon_(const char *uplo, cblas_int_t *n, complex *a, cblas_int_t *lda,
     cblas_int_t *ipiv, float *anorm, float *rcond, complex *work, cblas_int_t *
    info);

/* Subroutine */ int csymv_(const char *uplo, cblas_int_t *n, complex *alpha, complex *
    a, cblas_int_t *lda, complex *x, cblas_int_t *incx, complex *beta, complex *y,
     cblas_int_t *incy);

/* Subroutine */ int csyr_(const char *uplo, cblas_int_t *n, complex *alpha, complex *x,
     cblas_int_t *incx, complex *a, cblas_int_t *lda);

/* Subroutine */ int csyrfs_(const char *uplo, cblas_int_t *n, cblas_int_t *nrhs, complex *
    a, cblas_int_t *lda, complex *af, cblas_int_t *ldaf, cblas_int_t *ipiv, complex *
    b, cblas_int_t *ldb, complex *x, cblas_int_t *ldx, float *ferr, float *berr, 
    complex *work, float *rwork, cblas_int_t *info);

/* Subroutine */ int csysv_(const char *uplo, cblas_int_t *n, cblas_int_t *nrhs, complex *a,
     cblas_int_t *lda, cblas_int_t *ipiv, complex *b, cblas_int_t *ldb, complex *work,
     cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int csysvx_(const char *fact, const char *uplo, cblas_int_t *n, cblas_int_t *
    nrhs, complex *a, cblas_int_t *lda, complex *af, cblas_int_t *ldaf, cblas_int_t *
    ipiv, complex *b, cblas_int_t *ldb, complex *x, cblas_int_t *ldx, float *rcond,
     float *ferr, float *berr, complex *work, cblas_int_t *lwork, float *rwork, 
    cblas_int_t *info);

/* Subroutine */ int csytf2_(const char *uplo, cblas_int_t *n, complex *a, cblas_int_t *lda,
     cblas_int_t *ipiv, cblas_int_t *info);

/* Subroutine */ int csytrf_(const char *uplo, cblas_int_t *n, complex *a, cblas_int_t *lda,
     cblas_int_t *ipiv, complex *work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int csytri_(const char *uplo, cblas_int_t *n, complex *a, cblas_int_t *lda,
     cblas_int_t *ipiv, complex *work, cblas_int_t *info);

/* Subroutine */ int csytrs_(const char *uplo, cblas_int_t *n, cblas_int_t *nrhs, complex *
    a, cblas_int_t *lda, cblas_int_t *ipiv, complex *b, cblas_int_t *ldb, cblas_int_t *
    info);

/* Subroutine */ int ctbcon_(const char *norm, const char *uplo, const char *diag, cblas_int_t *n, 
    cblas_int_t *kd, complex *ab, cblas_int_t *ldab, float *rcond, complex *work, 
    float *rwork, cblas_int_t *info);

/* Subroutine */ int ctbrfs_(const char *uplo, const char *trans, const char *diag, cblas_int_t *n, 
    cblas_int_t *kd, cblas_int_t *nrhs, complex *ab, cblas_int_t *ldab, complex *b, 
    cblas_int_t *ldb, complex *x, cblas_int_t *ldx, float *ferr, float *berr, 
    complex *work, float *rwork, cblas_int_t *info);

/* Subroutine */ int ctbtrs_(const char *uplo, const char *trans, const char *diag, cblas_int_t *n, 
    cblas_int_t *kd, cblas_int_t *nrhs, complex *ab, cblas_int_t *ldab, complex *b, 
    cblas_int_t *ldb, cblas_int_t *info);

/* Subroutine */ int ctgevc_(const char *side, const char *howmny, logical *select, 
    cblas_int_t *n, complex *s, cblas_int_t *lds, complex *p, cblas_int_t *ldp, 
    complex *vl, cblas_int_t *ldvl, complex *vr, cblas_int_t *ldvr, cblas_int_t *mm, 
    cblas_int_t *m, complex *work, float *rwork, cblas_int_t *info);

/* Subroutine */ int ctgex2_(logical *wantq, logical *wantz, cblas_int_t *n, 
    complex *a, cblas_int_t *lda, complex *b, cblas_int_t *ldb, complex *q, 
    cblas_int_t *ldq, complex *z__, cblas_int_t *ldz, cblas_int_t *j1, cblas_int_t *info);

/* Subroutine */ int ctgexc_(logical *wantq, logical *wantz, cblas_int_t *n, 
    complex *a, cblas_int_t *lda, complex *b, cblas_int_t *ldb, complex *q, 
    cblas_int_t *ldq, complex *z__, cblas_int_t *ldz, cblas_int_t *ifst, cblas_int_t *
    ilst, cblas_int_t *info);

/* Subroutine */ int ctgsen_(cblas_int_t *ijob, logical *wantq, logical *wantz, 
    logical *select, cblas_int_t *n, complex *a, cblas_int_t *lda, complex *b, 
    cblas_int_t *ldb, complex *alpha, complex *beta, complex *q, cblas_int_t *ldq,
     complex *z__, cblas_int_t *ldz, cblas_int_t *m, float *pl, float *pr, float *
    dif, complex *work, cblas_int_t *lwork, cblas_int_t *iwork, cblas_int_t *liwork, 
    cblas_int_t *info);

/* Subroutine */ int ctgsja_(const char *jobu, const char *jobv, const char *jobq, cblas_int_t *m, 
    cblas_int_t *p, cblas_int_t *n, cblas_int_t *k, cblas_int_t *l, complex *a, cblas_int_t *
    lda, complex *b, cblas_int_t *ldb, float *tola, float *tolb, float *alpha, 
    float *beta, complex *u, cblas_int_t *ldu, complex *v, cblas_int_t *ldv, 
    complex *q, cblas_int_t *ldq, complex *work, cblas_int_t *ncycle, cblas_int_t *
    info);

/* Subroutine */ int ctgsna_(const char *job, const char *howmny, logical *select, 
    cblas_int_t *n, complex *a, cblas_int_t *lda, complex *b, cblas_int_t *ldb, 
    complex *vl, cblas_int_t *ldvl, complex *vr, cblas_int_t *ldvr, float *s, float 
    *dif, cblas_int_t *mm, cblas_int_t *m, complex *work, cblas_int_t *lwork, cblas_int_t 
    *iwork, cblas_int_t *info);

/* Subroutine */ int ctgsy2_(const char *trans, cblas_int_t *ijob, cblas_int_t *m, cblas_int_t *
    n, complex *a, cblas_int_t *lda, complex *b, cblas_int_t *ldb, complex *c__, 
    cblas_int_t *ldc, complex *d__, cblas_int_t *ldd, complex *e, cblas_int_t *lde, 
    complex *f, cblas_int_t *ldf, float *scale, float *rdsum, float *rdscal, 
    cblas_int_t *info);

/* Subroutine */ int ctgsyl_(const char *trans, cblas_int_t *ijob, cblas_int_t *m, cblas_int_t *
    n, complex *a, cblas_int_t *lda, complex *b, cblas_int_t *ldb, complex *c__, 
    cblas_int_t *ldc, complex *d__, cblas_int_t *ldd, complex *e, cblas_int_t *lde, 
    complex *f, cblas_int_t *ldf, float *scale, float *dif, complex *work, 
    cblas_int_t *lwork, cblas_int_t *iwork, cblas_int_t *info);

/* Subroutine */ int ctpcon_(const char *norm, const char *uplo, const char *diag, cblas_int_t *n, 
    complex *ap, float *rcond, complex *work, float *rwork, cblas_int_t *info);

/* Subroutine */ int ctprfs_(const char *uplo, const char *trans, const char *diag, cblas_int_t *n, 
    cblas_int_t *nrhs, complex *ap, complex *b, cblas_int_t *ldb, complex *x, 
    cblas_int_t *ldx, float *ferr, float *berr, complex *work, float *rwork, 
    cblas_int_t *info);

/* Subroutine */ int ctptri_(const char *uplo, const char *diag, cblas_int_t *n, complex *ap, 
    cblas_int_t *info);

/* Subroutine */ int ctptrs_(const char *uplo, const char *trans, const char *diag, cblas_int_t *n, 
    cblas_int_t *nrhs, complex *ap, complex *b, cblas_int_t *ldb, cblas_int_t *info);

/* Subroutine */ int ctrcon_(const char *norm, const char *uplo, const char *diag, cblas_int_t *n, 
    complex *a, cblas_int_t *lda, float *rcond, complex *work, float *rwork, 
    cblas_int_t *info);

/* Subroutine */ int ctrevc_(const char *side, const char *howmny, logical *select, 
    cblas_int_t *n, complex *t, cblas_int_t *ldt, complex *vl, cblas_int_t *ldvl, 
    complex *vr, cblas_int_t *ldvr, cblas_int_t *mm, cblas_int_t *m, complex *work, 
    float *rwork, cblas_int_t *info);

/* Subroutine */ int ctrexc_(const char *compq, cblas_int_t *n, complex *t, cblas_int_t *
    ldt, complex *q, cblas_int_t *ldq, cblas_int_t *ifst, cblas_int_t *ilst, cblas_int_t *
    info);

/* Subroutine */ int ctrrfs_(const char *uplo, const char *trans, const char *diag, cblas_int_t *n, 
    cblas_int_t *nrhs, complex *a, cblas_int_t *lda, complex *b, cblas_int_t *ldb, 
    complex *x, cblas_int_t *ldx, float *ferr, float *berr, complex *work, float 
    *rwork, cblas_int_t *info);

/* Subroutine */ int ctrsen_(const char *job, const char *compq, logical *select, cblas_int_t 
    *n, complex *t, cblas_int_t *ldt, complex *q, cblas_int_t *ldq, complex *w, 
    cblas_int_t *m, float *s, float *sep, complex *work, cblas_int_t *lwork, 
    cblas_int_t *info);

/* Subroutine */ int ctrsna_(const char *job, const char *howmny, logical *select, 
    cblas_int_t *n, complex *t, cblas_int_t *ldt, complex *vl, cblas_int_t *ldvl, 
    complex *vr, cblas_int_t *ldvr, float *s, float *sep, cblas_int_t *mm, cblas_int_t *
    m, complex *work, cblas_int_t *ldwork, float *rwork, cblas_int_t *info);

/* Subroutine */ int ctrsyl_(const char *trana, const char *tranb, cblas_int_t *isgn, cblas_int_t 
    *m, cblas_int_t *n, complex *a, cblas_int_t *lda, complex *b, cblas_int_t *ldb, 
    complex *c__, cblas_int_t *ldc, float *scale, cblas_int_t *info);

/* Subroutine */ int ctrti2_(const char *uplo, const char *diag, cblas_int_t *n, complex *a, 
    cblas_int_t *lda, cblas_int_t *info);

/* Subroutine */ int ctrtri_(const char *uplo, const char *diag, cblas_int_t *n, complex *a, 
    cblas_int_t *lda, cblas_int_t *info);

/* Subroutine */ int ctrtrs_(const char *uplo, const char *trans, const char *diag, cblas_int_t *n, 
    cblas_int_t *nrhs, complex *a, cblas_int_t *lda, complex *b, cblas_int_t *ldb, 
    cblas_int_t *info);

/* Subroutine */ int ctzrqf_(cblas_int_t *m, cblas_int_t *n, complex *a, cblas_int_t *lda,
     complex *tau, cblas_int_t *info);

/* Subroutine */ int ctzrzf_(cblas_int_t *m, cblas_int_t *n, complex *a, cblas_int_t *lda,
     complex *tau, complex *work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int cung2l_(cblas_int_t *m, cblas_int_t *n, cblas_int_t *k, complex *a, 
    cblas_int_t *lda, complex *tau, complex *work, cblas_int_t *info);

/* Subroutine */ int cung2r_(cblas_int_t *m, cblas_int_t *n, cblas_int_t *k, complex *a, 
    cblas_int_t *lda, complex *tau, complex *work, cblas_int_t *info);

/* Subroutine */ int cungbr_(const char *vect, cblas_int_t *m, cblas_int_t *n, cblas_int_t *k, 
    complex *a, cblas_int_t *lda, complex *tau, complex *work, cblas_int_t *lwork,
     cblas_int_t *info);

/* Subroutine */ int cunghr_(cblas_int_t *n, cblas_int_t *ilo, cblas_int_t *ihi, complex *
    a, cblas_int_t *lda, complex *tau, complex *work, cblas_int_t *lwork, cblas_int_t 
    *info);

/* Subroutine */ int cungl2_(cblas_int_t *m, cblas_int_t *n, cblas_int_t *k, complex *a, 
    cblas_int_t *lda, complex *tau, complex *work, cblas_int_t *info);

/* Subroutine */ int cunglq_(cblas_int_t *m, cblas_int_t *n, cblas_int_t *k, complex *a, 
    cblas_int_t *lda, complex *tau, complex *work, cblas_int_t *lwork, cblas_int_t *
    info);

/* Subroutine */ int cungql_(cblas_int_t *m, cblas_int_t *n, cblas_int_t *k, complex *a, 
    cblas_int_t *lda, complex *tau, complex *work, cblas_int_t *lwork, cblas_int_t *
    info);

/* Subroutine */ int cungqr_(cblas_int_t *m, cblas_int_t *n, cblas_int_t *k, complex *a, 
    cblas_int_t *lda, complex *tau, complex *work, cblas_int_t *lwork, cblas_int_t *
    info);

/* Subroutine */ int cungr2_(cblas_int_t *m, cblas_int_t *n, cblas_int_t *k, complex *a, 
    cblas_int_t *lda, complex *tau, complex *work, cblas_int_t *info);

/* Subroutine */ int cungrq_(cblas_int_t *m, cblas_int_t *n, cblas_int_t *k, complex *a, 
    cblas_int_t *lda, complex *tau, complex *work, cblas_int_t *lwork, cblas_int_t *
    info);

/* Subroutine */ int cungtr_(const char *uplo, cblas_int_t *n, complex *a, cblas_int_t *lda,
     complex *tau, complex *work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int cunm2l_(const char *side, const char *trans, cblas_int_t *m, cblas_int_t *n, 
    cblas_int_t *k, complex *a, cblas_int_t *lda, complex *tau, complex *c__, 
    cblas_int_t *ldc, complex *work, cblas_int_t *info);

/* Subroutine */ int cunm2r_(const char *side, const char *trans, cblas_int_t *m, cblas_int_t *n, 
    cblas_int_t *k, complex *a, cblas_int_t *lda, complex *tau, complex *c__, 
    cblas_int_t *ldc, complex *work, cblas_int_t *info);

/* Subroutine */ int cunmbr_(const char *vect, const char *side, const char *trans, cblas_int_t *m, 
    cblas_int_t *n, cblas_int_t *k, complex *a, cblas_int_t *lda, complex *tau, 
    complex *c__, cblas_int_t *ldc, complex *work, cblas_int_t *lwork, cblas_int_t *
    info);

/* Subroutine */ int cunmhr_(const char *side, const char *trans, cblas_int_t *m, cblas_int_t *n, 
    cblas_int_t *ilo, cblas_int_t *ihi, complex *a, cblas_int_t *lda, complex *tau, 
    complex *c__, cblas_int_t *ldc, complex *work, cblas_int_t *lwork, cblas_int_t *
    info);

/* Subroutine */ int cunml2_(const char *side, const char *trans, cblas_int_t *m, cblas_int_t *n, 
    cblas_int_t *k, complex *a, cblas_int_t *lda, complex *tau, complex *c__, 
    cblas_int_t *ldc, complex *work, cblas_int_t *info);

/* Subroutine */ int cunmlq_(const char *side, const char *trans, cblas_int_t *m, cblas_int_t *n, 
    cblas_int_t *k, complex *a, cblas_int_t *lda, complex *tau, complex *c__,
        cblas_int_t *ldc, complex *work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int cunmql_(const char *side, const char *trans, cblas_int_t *m, cblas_int_t *n, 
    cblas_int_t *k, complex *a, cblas_int_t *lda, complex *tau, complex *c__,
        cblas_int_t *ldc, complex *work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int cunmqr_(const char *side, const char *trans, cblas_int_t *m, cblas_int_t *n, 
    cblas_int_t *k, complex *a, cblas_int_t *lda, complex *tau, complex *c__,
        cblas_int_t *ldc, complex *work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int cunmr2_(const char *side, const char *trans, cblas_int_t *m, cblas_int_t *n, 
    cblas_int_t *k, complex *a, cblas_int_t *lda, complex *tau, complex *c__, 
    cblas_int_t *ldc, complex *work, cblas_int_t *info);

/* Subroutine */ int cunmr3_(const char *side, const char *trans, cblas_int_t *m, cblas_int_t *n, 
    cblas_int_t *k, cblas_int_t *l, complex *a, cblas_int_t *lda, complex *tau, 
    complex *c__, cblas_int_t *ldc, complex *work, cblas_int_t *info);

/* Subroutine */ int cunmrq_(const char *side, const char *trans, cblas_int_t *m, cblas_int_t *n, 
    cblas_int_t *k, complex *a, cblas_int_t *lda, complex *tau, complex *c__,   
    cblas_int_t *ldc, complex *work, cblas_int_t *lwork, cblas_int_t *iinfo);

/* Subroutine */ int cupgtr_(const char *uplo, cblas_int_t *n, complex *ap, complex *
    tau, complex *q, cblas_int_t *ldq, complex *work, cblas_int_t *info);

/* Subroutine */ int cupmtr_(const char *side, const char *uplo, const char *trans, cblas_int_t *m, 
    cblas_int_t *n, complex *ap, complex *tau, complex *c__, cblas_int_t *ldc, 
    complex *work, cblas_int_t *info);

/* Subroutine */ int dbdsdc_(const char *uplo, const char *compq, cblas_int_t *n, double *
    d__, double *e, double *u, cblas_int_t *ldu, double *vt, 
    cblas_int_t *ldvt, double *q, cblas_int_t *iq, double *work, cblas_int_t *
    iwork, cblas_int_t *info);

/* Subroutine */ int dbdsqr_(const char *uplo, cblas_int_t *n, cblas_int_t *ncvt, cblas_int_t *
    nru, cblas_int_t *ncc, double *d__, double *e, double *vt, 
    cblas_int_t *ldvt, double *u, cblas_int_t *ldu, double *c__, cblas_int_t *
    ldc, double *work, cblas_int_t *info);

/* Subroutine */ int ddisna_(const char *job, cblas_int_t *m, cblas_int_t *n, double *
    d__, double *sep, cblas_int_t *info);

/* Subroutine */ int dgbbrd_(const char *vect, cblas_int_t *m, cblas_int_t *n, cblas_int_t *ncc,
     cblas_int_t *kl, cblas_int_t *ku, double *ab, cblas_int_t *ldab, double *
    d__, double *e, double *q, cblas_int_t *ldq, double *pt, 
    cblas_int_t *ldpt, double *c__, cblas_int_t *ldc, double *work, 
    cblas_int_t *info);

/* Subroutine */ int dgbcon_(const char *norm, cblas_int_t *n, cblas_int_t *kl, cblas_int_t *ku,
     double *ab, cblas_int_t *ldab, cblas_int_t *ipiv, double *anorm, 
    double *rcond, double *work, cblas_int_t *iwork, cblas_int_t *info);

/* Subroutine */ int dgbequ_(cblas_int_t *m, cblas_int_t *n, cblas_int_t *kl, cblas_int_t *ku,
     double *ab, cblas_int_t *ldab, double *r__, double *c__, 
    double *rowcnd, double *colcnd, double *amax, cblas_int_t *
    info);

/* Subroutine */ int dgbrfs_(const char *trans, cblas_int_t *n, cblas_int_t *kl, cblas_int_t *
    ku, cblas_int_t *nrhs, double *ab, cblas_int_t *ldab, double *afb, 
    cblas_int_t *ldafb, cblas_int_t *ipiv, double *b, cblas_int_t *ldb, 
    double *x, cblas_int_t *ldx, double *ferr, double *berr, 
    double *work, cblas_int_t *iwork, cblas_int_t *info);

/* Subroutine */ int dgbsv_(cblas_int_t *n, cblas_int_t *kl, cblas_int_t *ku, cblas_int_t *
    nrhs, double *ab, cblas_int_t *ldab, cblas_int_t *ipiv, double *b, 
    cblas_int_t *ldb, cblas_int_t *info);

/* Subroutine */ int dgbsvx_(const char *fact, const char *trans, cblas_int_t *n, cblas_int_t *kl,
     cblas_int_t *ku, cblas_int_t *nrhs, double *ab, cblas_int_t *ldab, 
    double *afb, cblas_int_t *ldafb, cblas_int_t *ipiv, const char *equed, 
    double *r__, double *c__, double *b, cblas_int_t *ldb, 
    double *x, cblas_int_t *ldx, double *rcond, double *ferr, 
    double *berr, double *work, cblas_int_t *iwork, cblas_int_t *info);

/* Subroutine */ int dgbtf2_(cblas_int_t *m, cblas_int_t *n, cblas_int_t *kl, cblas_int_t *ku,
     double *ab, cblas_int_t *ldab, cblas_int_t *ipiv, cblas_int_t *info);

/* Subroutine */ int dgbtrf_(cblas_int_t *m, cblas_int_t *n, cblas_int_t *kl, cblas_int_t *ku,
     double *ab, cblas_int_t *ldab, cblas_int_t *ipiv, cblas_int_t *info);

/* Subroutine */ int dgbtrs_(const char *trans, cblas_int_t *n, cblas_int_t *kl, cblas_int_t *
    ku, cblas_int_t *nrhs, double *ab, cblas_int_t *ldab, cblas_int_t *ipiv, 
    double *b, cblas_int_t *ldb, cblas_int_t *info);

/* Subroutine */ int dgebak_(const char *job, const char *side, cblas_int_t *n, cblas_int_t *ilo, 
    cblas_int_t *ihi, double *scale, cblas_int_t *m, double *v, cblas_int_t *
    ldv, cblas_int_t *info);

/* Subroutine */ int dgebal_(const char *job, cblas_int_t *n, double *a, cblas_int_t *
    lda, cblas_int_t *ilo, cblas_int_t *ihi, double *scale, cblas_int_t *info);

/* Subroutine */ int dgebd2_(cblas_int_t *m, cblas_int_t *n, double *a, cblas_int_t *
    lda, double *d__, double *e, double *tauq, double *
    taup, double *work, cblas_int_t *info);

/* Subroutine */ int dgebrd_(cblas_int_t *m, cblas_int_t *n, double *a, cblas_int_t *
    lda, double *d__, double *e, double *tauq, double *
    taup, double *work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int dgecon_(const char *norm, cblas_int_t *n, double *a, cblas_int_t *
    lda, double *anorm, double *rcond, double *work, cblas_int_t *
    iwork, cblas_int_t *info);

/* Subroutine */ int dgeequ_(cblas_int_t *m, cblas_int_t *n, double *a, cblas_int_t *
    lda, double *r__, double *c__, double *rowcnd, double 
    *colcnd, double *amax, cblas_int_t *info);

/* Subroutine */ int dgees_(const char *jobvs, const char *sort, L_fp select, cblas_int_t *n, 
    double *a, cblas_int_t *lda, cblas_int_t *sdim, double *wr, 
    double *wi, double *vs, cblas_int_t *ldvs, double *work, 
    cblas_int_t *lwork, logical *bwork, cblas_int_t *info);

/* Subroutine */ int dgeesx_(const char *jobvs, const char *sort, L_fp select, const char *
    sense, cblas_int_t *n, double *a, cblas_int_t *lda, cblas_int_t *sdim, 
    double *wr, double *wi, double *vs, cblas_int_t *ldvs, 
    double *rconde, double *rcondv, double *work, cblas_int_t *
    lwork, cblas_int_t *iwork, cblas_int_t *liwork, logical *bwork, cblas_int_t *info);

/* Subroutine */ int dgeev_(const char *jobvl, const char *jobvr, cblas_int_t *n, double *
    a, cblas_int_t *lda, double *wr, double *wi, double *vl, 
    cblas_int_t *ldvl, double *vr, cblas_int_t *ldvr, double *work, 
    cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int dgeevx_(const char *balanc, const char *jobvl, const char *jobvr, const char *
    sense, cblas_int_t *n, double *a, cblas_int_t *lda, double *wr, 
    double *wi, double *vl, cblas_int_t *ldvl, double *vr, 
    cblas_int_t *ldvr, cblas_int_t *ilo, cblas_int_t *ihi, double *scale, 
    double *abnrm, double *rconde, double *rcondv, double   
    *work, cblas_int_t *lwork, cblas_int_t *iwork, cblas_int_t *info);

/* Subroutine */ int dgegs_(const char *jobvsl, const char *jobvsr, cblas_int_t *n, 
    double *a, cblas_int_t *lda, double *b, cblas_int_t *ldb, double *
    alphar, double *alphai, double *beta, double *vsl, 
    cblas_int_t *ldvsl, double *vsr, cblas_int_t *ldvsr, double *work, 
    cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int dgegv_(const char *jobvl, const char *jobvr, cblas_int_t *n, double *
    a, cblas_int_t *lda, double *b, cblas_int_t *ldb, double *alphar, 
    double *alphai, double *beta, double *vl, cblas_int_t *ldvl, 
    double *vr, cblas_int_t *ldvr, double *work, cblas_int_t *lwork, 
    cblas_int_t *info);

/* Subroutine */ int dgehd2_(cblas_int_t *n, cblas_int_t *ilo, cblas_int_t *ihi, 
    double *a, cblas_int_t *lda, double *tau, double *work, 
    cblas_int_t *info);

/* Subroutine */ int dgehrd_(cblas_int_t *n, cblas_int_t *ilo, cblas_int_t *ihi, 
    double *a, cblas_int_t *lda, double *tau, double *work, 
    cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int dgelq2_(cblas_int_t *m, cblas_int_t *n, double *a, cblas_int_t *
    lda, double *tau, double *work, cblas_int_t *info);

/* Subroutine */ int dgelqf_(cblas_int_t *m, cblas_int_t *n, double *a, cblas_int_t *
    lda, double *tau, double *work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int dgels_(const char *trans, cblas_int_t *m, cblas_int_t *n, cblas_int_t *
    nrhs, double *a, cblas_int_t *lda, double *b, cblas_int_t *ldb, 
    double *work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int dgelsd_(cblas_int_t *m, cblas_int_t *n, cblas_int_t *nrhs, 
    double *a, cblas_int_t *lda, double *b, cblas_int_t *ldb, double *
    s, double *rcond, cblas_int_t *rank, double *work, cblas_int_t *lwork,
     cblas_int_t *iwork, cblas_int_t *info);

/* Subroutine */ int dgelss_(cblas_int_t *m, cblas_int_t *n, cblas_int_t *nrhs, 
    double *a, cblas_int_t *lda, double *b, cblas_int_t *ldb, double *
    s, double *rcond, cblas_int_t *rank, double *work, cblas_int_t *lwork,
     cblas_int_t *info);

/* Subroutine */ int dgelsx_(cblas_int_t *m, cblas_int_t *n, cblas_int_t *nrhs, 
    double *a, cblas_int_t *lda, double *b, cblas_int_t *ldb, cblas_int_t *
    jpvt, double *rcond, cblas_int_t *rank, double *work, cblas_int_t *
    info);

/* Subroutine */ int dgelsy_(cblas_int_t *m, cblas_int_t *n, cblas_int_t *nrhs, 
    double *a, cblas_int_t *lda, double *b, cblas_int_t *ldb, cblas_int_t *
    jpvt, double *rcond, cblas_int_t *rank, double *work, cblas_int_t *
    lwork, cblas_int_t *info);

/* Subroutine */ int dgeql2_(cblas_int_t *m, cblas_int_t *n, double *a, cblas_int_t *
    lda, double *tau, double *work, cblas_int_t *info);

/* Subroutine */ int dgeqlf_(cblas_int_t *m, cblas_int_t *n, double *a, cblas_int_t *
    lda, double *tau, double *work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int dgeqp3_(cblas_int_t *m, cblas_int_t *n, double *a, cblas_int_t *
    lda, cblas_int_t *jpvt, double *tau, double *work, cblas_int_t *lwork,
     cblas_int_t *info);

/* Subroutine */ int dgeqpf_(cblas_int_t *m, cblas_int_t *n, double *a, cblas_int_t *
    lda, cblas_int_t *jpvt, double *tau, double *work, cblas_int_t *info);

/* Subroutine */ int dgeqr2_(cblas_int_t *m, cblas_int_t *n, double *a, cblas_int_t *
    lda, double *tau, double *work, cblas_int_t *info);

/* Subroutine */ int dgeqrf_(cblas_int_t *m, cblas_int_t *n, double *a, cblas_int_t *
    lda, double *tau, double *work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int dgerfs_(const char *trans, cblas_int_t *n, cblas_int_t *nrhs, 
    double *a, cblas_int_t *lda, double *af, cblas_int_t *ldaf, cblas_int_t *
    ipiv, double *b, cblas_int_t *ldb, double *x, cblas_int_t *ldx, 
    double *ferr, double *berr, double *work, cblas_int_t *iwork, 
    cblas_int_t *info);

/* Subroutine */ int dgerq2_(cblas_int_t *m, cblas_int_t *n, double *a, cblas_int_t *
    lda, double *tau, double *work, cblas_int_t *info);

/* Subroutine */ int dgerqf_(cblas_int_t *m, cblas_int_t *n, double *a, cblas_int_t *
    lda, double *tau, double *work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int dgesc2_(cblas_int_t *n, double *a, cblas_int_t *lda, 
    double *rhs, cblas_int_t *ipiv, cblas_int_t *jpiv, double *scale);

/* Subroutine */ int dgesdd_(const char *jobz, cblas_int_t *m, cblas_int_t *n, double *
    a, cblas_int_t *lda, double *s, double *u, cblas_int_t *ldu, 
    double *vt, cblas_int_t *ldvt, double *work, cblas_int_t *lwork, 
    cblas_int_t *iwork, cblas_int_t *info);

/* Subroutine */ int dgesv_(cblas_int_t *n, cblas_int_t *nrhs, double *a, cblas_int_t 
    *lda, cblas_int_t *ipiv, double *b, cblas_int_t *ldb, cblas_int_t *info);

/* Subroutine */ int dgesvd_(const char *jobu, const char *jobvt, cblas_int_t *m, cblas_int_t *n, 
    double *a, cblas_int_t *lda, double *s, double *u, cblas_int_t *
    ldu, double *vt, cblas_int_t *ldvt, double *work, cblas_int_t *lwork, 
    cblas_int_t *info);

/* Subroutine */ int dgesvx_(const char *fact, const char *trans, cblas_int_t *n, cblas_int_t *
    nrhs, double *a, cblas_int_t *lda, double *af, cblas_int_t *ldaf, 
    cblas_int_t *ipiv, const char *equed, double *r__, double *c__, 
    double *b, cblas_int_t *ldb, double *x, cblas_int_t *ldx, double *
    rcond, double *ferr, double *berr, double *work, cblas_int_t *
    iwork, cblas_int_t *info);

/* Subroutine */ int dgetc2_(cblas_int_t *n, double *a, cblas_int_t *lda, cblas_int_t 
    *ipiv, cblas_int_t *jpiv, cblas_int_t *info);

/* Subroutine */ int dgetf2_(cblas_int_t *m, cblas_int_t *n, double *a, cblas_int_t *
    lda, cblas_int_t *ipiv, cblas_int_t *info);

/* Subroutine */ int dgetrf_(cblas_int_t *m, cblas_int_t *n, double *a, cblas_int_t *
    lda, cblas_int_t *ipiv, cblas_int_t *info);

/* Subroutine */ int dgetri_(cblas_int_t *n, double *a, cblas_int_t *lda, cblas_int_t 
    *ipiv, double *work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int dgetrs_(const char *trans, cblas_int_t *n, cblas_int_t *nrhs, 
    double *a, cblas_int_t *lda, cblas_int_t *ipiv, double *b, cblas_int_t *
    ldb, cblas_int_t *info);

/* Subroutine */ int dggbak_(const char *job, const char *side, cblas_int_t *n, cblas_int_t *ilo, 
    cblas_int_t *ihi, double *lscale, double *rscale, cblas_int_t *m, 
    double *v, cblas_int_t *ldv, cblas_int_t *info);

/* Subroutine */ int dggbal_(const char *job, cblas_int_t *n, double *a, cblas_int_t *
    lda, double *b, cblas_int_t *ldb, cblas_int_t *ilo, cblas_int_t *ihi, 
    double *lscale, double *rscale, double *work, cblas_int_t *
    info);

/* Subroutine */ int dgges_(const char *jobvsl, const char *jobvsr, const char *sort, L_fp 
    selctg, cblas_int_t *n, double *a, cblas_int_t *lda, double *b, 
    cblas_int_t *ldb, cblas_int_t *sdim, double *alphar, double *alphai, 
    double *beta, double *vsl, cblas_int_t *ldvsl, double *vsr, 
    cblas_int_t *ldvsr, double *work, cblas_int_t *lwork, logical *bwork, 
    cblas_int_t *info);

/* Subroutine */ int dggesx_(const char *jobvsl, const char *jobvsr, const char *sort, L_fp 
    selctg, const char *sense, cblas_int_t *n, double *a, cblas_int_t *lda, 
    double *b, cblas_int_t *ldb, cblas_int_t *sdim, double *alphar, 
    double *alphai, double *beta, double *vsl, cblas_int_t *ldvsl,
     double *vsr, cblas_int_t *ldvsr, double *rconde, double *
    rcondv, double *work, cblas_int_t *lwork, cblas_int_t *iwork, cblas_int_t *     
    liwork, logical *bwork, cblas_int_t *info);

/* Subroutine */ int dggev_(const char *jobvl, const char *jobvr, cblas_int_t *n, double *
    a, cblas_int_t *lda, double *b, cblas_int_t *ldb, double *alphar, 
    double *alphai, double *beta, double *vl, cblas_int_t *ldvl, 
    double *vr, cblas_int_t *ldvr, double *work, cblas_int_t *lwork, 
    cblas_int_t *info);

/* Subroutine */ int dggevx_(const char *balanc, const char *jobvl, const char *jobvr, const char *
    sense, cblas_int_t *n, double *a, cblas_int_t *lda, double *b, 
    cblas_int_t *ldb, double *alphar, double *alphai, double *
    beta, double *vl, cblas_int_t *ldvl, double *vr, cblas_int_t *ldvr, 
    cblas_int_t *ilo, cblas_int_t *ihi, double *lscale, double *rscale, 
    double *abnrm, double *bbnrm, double *rconde, double *
    rcondv, double *work, cblas_int_t *lwork, cblas_int_t *iwork, logical *     
    bwork, cblas_int_t *info);

/* Subroutine */ int dggglm_(cblas_int_t *n, cblas_int_t *m, cblas_int_t *p, double *
    a, cblas_int_t *lda, double *b, cblas_int_t *ldb, double *d__, 
    double *x, double *y, double *work, cblas_int_t *lwork, 
    cblas_int_t *info);

/* Subroutine */ int dgghrd_(const char *compq, const char *compz, cblas_int_t *n, cblas_int_t *
    ilo, cblas_int_t *ihi, double *a, cblas_int_t *lda, double *b, 
    cblas_int_t *ldb, double *q, cblas_int_t *ldq, double *z__, cblas_int_t *
    ldz, cblas_int_t *info);

/* Subroutine */ int dgglse_(cblas_int_t *m, cblas_int_t *n, cblas_int_t *p, double *
    a, cblas_int_t *lda, double *b, cblas_int_t *ldb, double *c__, 
    double *d__, double *x, double *work, cblas_int_t *lwork, 
    cblas_int_t *info);

/* Subroutine */ int dggqrf_(cblas_int_t *n, cblas_int_t *m, cblas_int_t *p, double *
    a, cblas_int_t *lda, double *taua, double *b, cblas_int_t *ldb, 
    double *taub, double *work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int dggrqf_(cblas_int_t *m, cblas_int_t *p, cblas_int_t *n, double *
    a, cblas_int_t *lda, double *taua, double *b, cblas_int_t *ldb, 
    double *taub, double *work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int dggsvd_(const char *jobu, const char *jobv, const char *jobq, cblas_int_t *m, 
    cblas_int_t *n, cblas_int_t *p, cblas_int_t *k, cblas_int_t *l, double *a, 
    cblas_int_t *lda, double *b, cblas_int_t *ldb, double *alpha, 
    double *beta, double *u, cblas_int_t *ldu, double *v, cblas_int_t 
    *ldv, double *q, cblas_int_t *ldq, double *work, cblas_int_t *iwork, 
    cblas_int_t *info);

/* Subroutine */ int dggsvp_(const char *jobu, const char *jobv, const char *jobq, cblas_int_t *m, 
    cblas_int_t *p, cblas_int_t *n, double *a, cblas_int_t *lda, double *b, 
    cblas_int_t *ldb, double *tola, double *tolb, cblas_int_t *k, cblas_int_t 
    *l, double *u, cblas_int_t *ldu, double *v, cblas_int_t *ldv, 
    double *q, cblas_int_t *ldq, cblas_int_t *iwork, double *tau, 
    double *work, cblas_int_t *info);

/* Subroutine */ int dgtcon_(const char *norm, cblas_int_t *n, double *dl, 
    double *d__, double *du, double *du2, cblas_int_t *ipiv, 
    double *anorm, double *rcond, double *work, cblas_int_t *
    iwork, cblas_int_t *info);

/* Subroutine */ int dgtrfs_(const char *trans, cblas_int_t *n, cblas_int_t *nrhs, 
    double *dl, double *d__, double *du, double *dlf, 
    double *df, double *duf, double *du2, cblas_int_t *ipiv, 
    double *b, cblas_int_t *ldb, double *x, cblas_int_t *ldx, double *
    ferr, double *berr, double *work, cblas_int_t *iwork, cblas_int_t *
    info);

/* Subroutine */ int dgtsv_(cblas_int_t *n, cblas_int_t *nrhs, double *dl, 
    double *d__, double *du, double *b, cblas_int_t *ldb, cblas_int_t 
    *info);

/* Subroutine */ int dgtsvx_(const char *fact, const char *trans, cblas_int_t *n, cblas_int_t *
    nrhs, double *dl, double *d__, double *du, double *
    dlf, double *df, double *duf, double *du2, cblas_int_t *ipiv, 
    double *b, cblas_int_t *ldb, double *x, cblas_int_t *ldx, double *
    rcond, double *ferr, double *berr, double *work, cblas_int_t *
    iwork, cblas_int_t *info);

/* Subroutine */ int dgttrf_(cblas_int_t *n, double *dl, double *d__, 
    double *du, double *du2, cblas_int_t *ipiv, cblas_int_t *info);

/* Subroutine */ int dgttrs_(const char *trans, cblas_int_t *n, cblas_int_t *nrhs, 
    double *dl, double *d__, double *du, double *du2,   
    cblas_int_t *ipiv, double *b, cblas_int_t *ldb, cblas_int_t *info);

/* Subroutine */ int dgtts2_(cblas_int_t *itrans, cblas_int_t *n, cblas_int_t *nrhs, 
    double *dl, double *d__, double *du, double *du2, 
    cblas_int_t *ipiv, double *b, cblas_int_t *ldb);

/* Subroutine */ int dhgeqz_(const char *job, const char *compq, const char *compz, cblas_int_t *n, 
    cblas_int_t *ilo, cblas_int_t *ihi, double *h__, cblas_int_t *ldh, double 
    *t, cblas_int_t *ldt, double *alphar, double *alphai, double *
    beta, double *q, cblas_int_t *ldq, double *z__, cblas_int_t *ldz, 
    double *work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int dhsein_(const char *side, const char *eigsrc, const char *initv, logical *
    select, cblas_int_t *n, double *h__, cblas_int_t *ldh, double *wr, 
    double *wi, double *vl, cblas_int_t *ldvl, double *vr, 
    cblas_int_t *ldvr, cblas_int_t *mm, cblas_int_t *m, double *work, cblas_int_t *
    ifaill, cblas_int_t *ifailr, cblas_int_t *info);

/* Subroutine */ int dhseqr_(const char *job, const char *compz, cblas_int_t *n, cblas_int_t *ilo,
     cblas_int_t *ihi, double *h__, cblas_int_t *ldh, double *wr, 
    double *wi, double *z__, cblas_int_t *ldz, double *work, 
    cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int dlabad_(double *small, double *large);

/* Subroutine */ int dlabrd_(cblas_int_t *m, cblas_int_t *n, cblas_int_t *nb, double *
    a, cblas_int_t *lda, double *d__, double *e, double *tauq, 
    double *taup, double *x, cblas_int_t *ldx, double *y, cblas_int_t 
    *ldy);

/* Subroutine */ int dlacn2_(cblas_int_t *n, double *v, double *x, 
    cblas_int_t *isgn, double *est, cblas_int_t *kase, cblas_int_t *isave);

/* Subroutine */ int dlacon_(cblas_int_t *n, double *v, double *x, 
    cblas_int_t *isgn, double *est, cblas_int_t *kase);

/* Subroutine */ int dlacpy_(const char *uplo, cblas_int_t *m, cblas_int_t *n, double *
    a, cblas_int_t *lda, double *b, cblas_int_t *ldb);

/* Subroutine */ int dladiv_(double *a, double *b, double *c__, 
    double *d__, double *p, double *q);

/* Subroutine */ int dlae2_(double *a, double *b, double *c__, 
    double *rt1, double *rt2);

/* Subroutine */ int dlaebz_(cblas_int_t *ijob, cblas_int_t *nitmax, cblas_int_t *n, 
    cblas_int_t *mmax, cblas_int_t *minp, cblas_int_t *nbmin, double *abstol, 
    double *reltol, double *pivmin, double *d__, double *
    e, double *e2, cblas_int_t *nval, double *ab, double *c__, 
    cblas_int_t *mout, cblas_int_t *nab, double *work, cblas_int_t *iwork, 
    cblas_int_t *info);

/* Subroutine */ int dlaed0_(cblas_int_t *icompq, cblas_int_t *qsiz, cblas_int_t *n, 
    double *d__, double *e, double *q, cblas_int_t *ldq, 
    double *qstore, cblas_int_t *ldqs, double *work, cblas_int_t *iwork, 
    cblas_int_t *info);

/* Subroutine */ int dlaed1_(cblas_int_t *n, double *d__, double *q, 
    cblas_int_t *ldq, cblas_int_t *indxq, double *rho, cblas_int_t *cutpnt, 
    double *work, cblas_int_t *iwork, cblas_int_t *info);

/* Subroutine */ int dlaed2_(cblas_int_t *k, cblas_int_t *n, cblas_int_t *n1, double *
    d__, double *q, cblas_int_t *ldq, cblas_int_t *indxq, double *rho, 
    double *z__, double *dlamda, double *w, double *q2, 
    cblas_int_t *indx, cblas_int_t *indxc, cblas_int_t *indxp, cblas_int_t *coltyp, 
    cblas_int_t *info);

/* Subroutine */ int dlaed3_(cblas_int_t *k, cblas_int_t *n, cblas_int_t *n1, double *
    d__, double *q, cblas_int_t *ldq, double *rho, double *dlamda,
     double *q2, cblas_int_t *indx, cblas_int_t *ctot, double *w, 
    double *s, cblas_int_t *info);

/* Subroutine */ int dlaed4_(cblas_int_t *n, cblas_int_t *i__, double *d__, 
    double *z__, double *delta, double *rho, double *dlam,
     cblas_int_t *info);

/* Subroutine */ int dlaed5_(cblas_int_t *i__, double *d__, double *z__, 
    double *delta, double *rho, double *dlam);

/* Subroutine */ int dlaed6_(cblas_int_t *kniter, logical *orgati, double *
    rho, double *d__, double *z__, double *finit, double *
    tau, cblas_int_t *info);

/* Subroutine */ int dlaed7_(cblas_int_t *icompq, cblas_int_t *n, cblas_int_t *qsiz, 
    cblas_int_t *tlvls, cblas_int_t *curlvl, cblas_int_t *curpbm, double *d__, 
    double *q, cblas_int_t *ldq, cblas_int_t *indxq, double *rho, cblas_int_t 
    *cutpnt, double *qstore, cblas_int_t *qptr, cblas_int_t *prmptr, cblas_int_t *
    perm, cblas_int_t *givptr, cblas_int_t *givcol, double *givnum, 
    double *work, cblas_int_t *iwork, cblas_int_t *info);

/* Subroutine */ int dlaed8_(cblas_int_t *icompq, cblas_int_t *k, cblas_int_t *n, cblas_int_t 
    *qsiz, double *d__, double *q, cblas_int_t *ldq, cblas_int_t *indxq, 
    double *rho, cblas_int_t *cutpnt, double *z__, double *dlamda,
     double *q2, cblas_int_t *ldq2, double *w, cblas_int_t *perm, cblas_int_t 
    *givptr, cblas_int_t *givcol, double *givnum, cblas_int_t *indxp, cblas_int_t 
    *indx, cblas_int_t *info);

/* Subroutine */ int dlaed9_(cblas_int_t *k, cblas_int_t *kstart, cblas_int_t *kstop, 
    cblas_int_t *n, double *d__, double *q, cblas_int_t *ldq, double *
    rho, double *dlamda, double *w, double *s, cblas_int_t *lds, 
    cblas_int_t *info);

/* Subroutine */ int dlaeda_(cblas_int_t *n, cblas_int_t *tlvls, cblas_int_t *curlvl, 
    cblas_int_t *curpbm, cblas_int_t *prmptr, cblas_int_t *perm, cblas_int_t *givptr, 
    cblas_int_t *givcol, double *givnum, double *q, cblas_int_t *qptr, 
    double *z__, double *ztemp, cblas_int_t *info);

/* Subroutine */ int dlaein_(logical *rightv, logical *noinit, cblas_int_t *n, 
    double *h__, cblas_int_t *ldh, double *wr, double *wi, 
    double *vr, double *vi, double *b, cblas_int_t *ldb, 
    double *work, double *eps3, double *smlnum, double *
    bignum, cblas_int_t *info);

/* Subroutine */ int dlaev2_(double *a, double *b, double *c__, 
    double *rt1, double *rt2, double *cs1, double *sn1);

/* Subroutine */ int dlaexc_(logical *wantq, cblas_int_t *n, double *t, 
    cblas_int_t *ldt, double *q, cblas_int_t *ldq, cblas_int_t *j1, cblas_int_t *n1, 
    cblas_int_t *n2, double *work, cblas_int_t *info);

/* Subroutine */ int dlag2_(double *a, cblas_int_t *lda, double *b, 
    cblas_int_t *ldb, double *safmin, double *scale1, double *
    scale2, double *wr1, double *wr2, double *wi);

/* Subroutine */ int dlag2s_(cblas_int_t *m, cblas_int_t *n, double *a, cblas_int_t *
    lda, float *sa, cblas_int_t *ldsa, cblas_int_t *info);

/* Subroutine */ int dlags2_(logical *upper, double *a1, double *a2, 
    double *a3, double *b1, double *b2, double *b3, 
    double *csu, double *snu, double *csv, double *snv, 
    double *csq, double *snq);

/* Subroutine */ int dlagtf_(cblas_int_t *n, double *a, double *lambda, 
    double *b, double *c__, double *tol, double *d__, 
    cblas_int_t *in, cblas_int_t *info);

/* Subroutine */ int dlagtm_(const char *trans, cblas_int_t *n, cblas_int_t *nrhs, 
    double *alpha, double *dl, double *d__, double *du, 
    double *x, cblas_int_t *ldx, double *beta, double *b, cblas_int_t 
    *ldb);

/* Subroutine */ int dlagts_(cblas_int_t *job, cblas_int_t *n, double *a, 
    double *b, double *c__, double *d__, cblas_int_t *in, 
    double *y, double *tol, cblas_int_t *info);

/* Subroutine */ int dlagv2_(double *a, cblas_int_t *lda, double *b, 
    cblas_int_t *ldb, double *alphar, double *alphai, double *
    beta, double *csl, double *snl, double *csr, double *
    snr);

/* Subroutine */ int dlahqr_(logical *wantt, logical *wantz, cblas_int_t *n, 
    cblas_int_t *ilo, cblas_int_t *ihi, double *h__, cblas_int_t *ldh, double 
    *wr, double *wi, cblas_int_t *iloz, cblas_int_t *ihiz, double *z__, 
    cblas_int_t *ldz, cblas_int_t *info);

/* Subroutine */ int dlahr2_(cblas_int_t *n, cblas_int_t *k, cblas_int_t *nb, double *
    a, cblas_int_t *lda, double *tau, double *t, cblas_int_t *ldt, 
    double *y, cblas_int_t *ldy);

/* Subroutine */ int dlahrd_(cblas_int_t *n, cblas_int_t *k, cblas_int_t *nb, double *
    a, cblas_int_t *lda, double *tau, double *t, cblas_int_t *ldt, 
    double *y, cblas_int_t *ldy);

/* Subroutine */ int dlaic1_(cblas_int_t *job, cblas_int_t *j, double *x, 
    double *sest, double *w, double *gamma, double *
    sestpr, double *s, double *c__);

/* Subroutine */ int dlaln2_(logical *ltrans, cblas_int_t *na, cblas_int_t *nw, 
    double *smin, double *ca, double *a, cblas_int_t *lda, 
    double *d1, double *d2, double *b, cblas_int_t *ldb, 
    double *wr, double *wi, double *x, cblas_int_t *ldx, 
    double *scale, double *xnorm, cblas_int_t *info);

/* Subroutine */ int dlals0_(cblas_int_t *icompq, cblas_int_t *nl, cblas_int_t *nr, 
    cblas_int_t *sqre, cblas_int_t *nrhs, double *b, cblas_int_t *ldb, double 
    *bx, cblas_int_t *ldbx, cblas_int_t *perm, cblas_int_t *givptr, cblas_int_t *givcol, 
    cblas_int_t *ldgcol, double *givnum, cblas_int_t *ldgnum, double *
    poles, double *difl, double *difr, double *z__, cblas_int_t *
    k, double *c__, double *s, double *work, cblas_int_t *info);

/* Subroutine */ int dlalsa_(cblas_int_t *icompq, cblas_int_t *smlsiz, cblas_int_t *n, 
    cblas_int_t *nrhs, double *b, cblas_int_t *ldb, double *bx, cblas_int_t *
    ldbx, double *u, cblas_int_t *ldu, double *vt, cblas_int_t *k, 
    double *difl, double *difr, double *z__, double *
    poles, cblas_int_t *givptr, cblas_int_t *givcol, cblas_int_t *ldgcol, cblas_int_t *
    perm, double *givnum, double *c__, double *s, double *
    work, cblas_int_t *iwork, cblas_int_t *info);

/* Subroutine */ int dlalsd_(const char *uplo, cblas_int_t *smlsiz, cblas_int_t *n, cblas_int_t 
    *nrhs, double *d__, double *e, double *b, cblas_int_t *ldb, 
    double *rcond, cblas_int_t *rank, double *work, cblas_int_t *iwork, 
    cblas_int_t *info);

/* Subroutine */ int dlamrg_(cblas_int_t *n1, cblas_int_t *n2, double *a, cblas_int_t 
    *dtrd1, cblas_int_t *dtrd2, cblas_int_t *index);

/* Subroutine */ int dlanv2_(double *a, double *b, double *c__, 
    double *d__, double *rt1r, double *rt1i, double *rt2r,
     double *rt2i, double *cs, double *sn);

/* Subroutine */ int dlapll_(cblas_int_t *n, double *x, cblas_int_t *incx, 
    double *y, cblas_int_t *incy, double *ssmin);

/* Subroutine */ int dlapmt_(logical *forwrd, cblas_int_t *m, cblas_int_t *n, 
    double *x, cblas_int_t *ldx, cblas_int_t *k);

/* Subroutine */ int dlaqgb_(cblas_int_t *m, cblas_int_t *n, cblas_int_t *kl, cblas_int_t *ku,
     double *ab, cblas_int_t *ldab, double *r__, double *c__, 
    double *rowcnd, double *colcnd, double *amax, const char *equed);

/* Subroutine */ int dlaqge_(cblas_int_t *m, cblas_int_t *n, double *a, cblas_int_t *
    lda, double *r__, double *c__, double *rowcnd, double 
    *colcnd, double *amax, const char *equed);

/* Subroutine */ int dlaqp2_(cblas_int_t *m, cblas_int_t *n, cblas_int_t *offset, 
    double *a, cblas_int_t *lda, cblas_int_t *jpvt, double *tau, 
    double *vn1, double *vn2, double *work);

/* Subroutine */ int dlaqps_(cblas_int_t *m, cblas_int_t *n, cblas_int_t *offset, cblas_int_t 
    *nb, cblas_int_t *kb, double *a, cblas_int_t *lda, cblas_int_t *jpvt, 
    double *tau, double *vn1, double *vn2, double *auxv, 
    double *f, cblas_int_t *ldf);

/* Subroutine */ int dlaqr0_(logical *wantt, logical *wantz, cblas_int_t *n, 
    cblas_int_t *ilo, cblas_int_t *ihi, double *h__, cblas_int_t *ldh, double 
    *wr, double *wi, cblas_int_t *iloz, cblas_int_t *ihiz, double *z__, 
    cblas_int_t *ldz, double *work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int dlaqr1_(cblas_int_t *n, double *h__, cblas_int_t *ldh, 
    double *sr1, double *si1, double *sr2, double *si2, 
    double *v);

/* Subroutine */ int dlaqr2_(logical *wantt, logical *wantz, cblas_int_t *n, 
    cblas_int_t *ktop, cblas_int_t *kbot, cblas_int_t *nw, double *h__, cblas_int_t *
    ldh, cblas_int_t *iloz, cblas_int_t *ihiz, double *z__, cblas_int_t *ldz, 
    cblas_int_t *ns, cblas_int_t *nd, double *sr, double *si, double *
    v, cblas_int_t *ldv, cblas_int_t *nh, double *t, cblas_int_t *ldt, cblas_int_t *
    nv, double *wv, cblas_int_t *ldwv, double *work, cblas_int_t *lwork);

/* Subroutine */ int dlaqr3_(logical *wantt, logical *wantz, cblas_int_t *n, 
    cblas_int_t *ktop, cblas_int_t *kbot, cblas_int_t *nw, double *h__, cblas_int_t *
    ldh, cblas_int_t *iloz, cblas_int_t *ihiz, double *z__, cblas_int_t *ldz, 
    cblas_int_t *ns, cblas_int_t *nd, double *sr, double *si, double *
    v, cblas_int_t *ldv, cblas_int_t *nh, double *t, cblas_int_t *ldt, cblas_int_t *
    nv, double *wv, cblas_int_t *ldwv, double *work, cblas_int_t *lwork);

/* Subroutine */ int dlaqr4_(logical *wantt, logical *wantz, cblas_int_t *n, 
    cblas_int_t *ilo, cblas_int_t *ihi, double *h__, cblas_int_t *ldh, double 
    *wr, double *wi, cblas_int_t *iloz, cblas_int_t *ihiz, double *z__, 
    cblas_int_t *ldz, double *work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int dlaqr5_(logical *wantt, logical *wantz, cblas_int_t *kacc22, 
    cblas_int_t *n, cblas_int_t *ktop, cblas_int_t *kbot, cblas_int_t *nshfts, double 
    *sr, double *si, double *h__, cblas_int_t *ldh, cblas_int_t *iloz, 
    cblas_int_t *ihiz, double *z__, cblas_int_t *ldz, double *v, cblas_int_t *
    ldv, double *u, cblas_int_t *ldu, cblas_int_t *nv, double *wv, 
    cblas_int_t *ldwv, cblas_int_t *nh, double *wh, cblas_int_t *ldwh);

/* Subroutine */ int dlaqsb_(const char *uplo, cblas_int_t *n, cblas_int_t *kd, double *
    ab, cblas_int_t *ldab, double *s, double *scond, double *amax,
     const char *equed);

/* Subroutine */ int dlaqsp_(const char *uplo, cblas_int_t *n, double *ap, 
    double *s, double *scond, double *amax, const char *equed);

/* Subroutine */ int dlaqsy_(const char *uplo, cblas_int_t *n, double *a, cblas_int_t *
    lda, double *s, double *scond, double *amax, const char *equed);

/* Subroutine */ int dlaqtr_(logical *ltran, logical *lfloat, cblas_int_t *n, 
    double *t, cblas_int_t *ldt, double *b, double *w, double 
    *scale, double *x, double *work, cblas_int_t *info);

/* Subroutine */ int dlar1v_(cblas_int_t *n, cblas_int_t *b1, cblas_int_t *bn, double 
    *lambda, double *d__, double *l, double *ld, double *
    lld, double *pivmin, double *gaptol, double *z__, logical 
    *wantnc, cblas_int_t *negcnt, double *ztz, double *mingma, 
    cblas_int_t *r__, cblas_int_t *isuppz, double *nrminv, double *resid, 
    double *rqcorr, double *work);

/* Subroutine */ int dlar2v_(cblas_int_t *n, double *x, double *y, 
    double *z__, cblas_int_t *incx, double *c__, double *s, 
    cblas_int_t *incc);

/* Subroutine */ int dlarf_(const char *side, cblas_int_t *m, cblas_int_t *n, double *v,
     cblas_int_t *incv, double *tau, double *c__, cblas_int_t *ldc, 
    double *work);

/* Subroutine */ int dlarfb_(const char *side, const char *trans, const char *direct, const char *
    storev, cblas_int_t *m, cblas_int_t *n, cblas_int_t *k, double *v, cblas_int_t *
    ldv, double *t, cblas_int_t *ldt, double *c__, cblas_int_t *ldc, 
    double *work, cblas_int_t *ldwork);

/* Subroutine */ int dlarfg_(cblas_int_t *n, double *alpha, double *x, 
    cblas_int_t *incx, double *tau);

/* Subroutine */ int dlarft_(const char *direct, const char *storev, cblas_int_t *n, cblas_int_t *
    k, double *v, cblas_int_t *ldv, double *tau, double *t, 
    cblas_int_t *ldt);

/* Subroutine */ int dlarfx_(const char *side, cblas_int_t *m, cblas_int_t *n, double *
    v, double *tau, double *c__, cblas_int_t *ldc, double *work);

/* Subroutine */ int dlargv_(cblas_int_t *n, double *x, cblas_int_t *incx, 
    double *y, cblas_int_t *incy, double *c__, cblas_int_t *incc);

/* Subroutine */ int dlarnv_(cblas_int_t *idist, cblas_int_t *iseed, cblas_int_t *n, 
    double *x);

/* Subroutine */ int dlarra_(cblas_int_t *n, double *d__, double *e, 
    double *e2, double *spltol, double *tnrm, cblas_int_t *nsplit,
     cblas_int_t *isplit, cblas_int_t *info);

/* Subroutine */ int dlarrb_(cblas_int_t *n, double *d__, double *lld, 
    cblas_int_t *ifirst, cblas_int_t *ilast, double *rtol1, double *rtol2,
     cblas_int_t *offset, double *w, double *wgap, double *werr, 
    double *work, cblas_int_t *iwork, double *pivmin, double *
    spdiam, cblas_int_t *twist, cblas_int_t *info);

/* Subroutine */ int dlarrc_(const char *jobt, cblas_int_t *n, double *vl, 
    double *vu, double *d__, double *e, double *pivmin, 
    cblas_int_t *eigcnt, cblas_int_t *lcnt, cblas_int_t *rcnt, cblas_int_t *info);

/* Subroutine */ int dlarrd_(const char *range, const char *order, cblas_int_t *n, double 
    *vl, double *vu, cblas_int_t *il, cblas_int_t *iu, double *gers, 
    double *reltol, double *d__, double *e, double *e2, 
    double *pivmin, cblas_int_t *nsplit, cblas_int_t *isplit, cblas_int_t *m, 
    double *w, double *werr, double *wl, double *wu, 
    cblas_int_t *iblock, cblas_int_t *indexw, double *work, cblas_int_t *iwork, 
    cblas_int_t *info);

/* Subroutine */ int dlarre_(const char *range, cblas_int_t *n, double *vl, 
    double *vu, cblas_int_t *il, cblas_int_t *iu, double *d__, double 
    *e, double *e2, double *rtol1, double *rtol2, double *
    spltol, cblas_int_t *nsplit, cblas_int_t *isplit, cblas_int_t *m, double *w, 
    double *werr, double *wgap, cblas_int_t *iblock, cblas_int_t *indexw, 
    double *gers, double *pivmin, double *work, cblas_int_t *
    iwork, cblas_int_t *info);

/* Subroutine */ int dlarrf_(cblas_int_t *n, double *d__, double *l, 
    double *ld, cblas_int_t *clstrt, cblas_int_t *clend, double *w, 
    double *wgap, double *werr, double *spdiam, double *
    clgapl, double *clgapr, double *pivmin, double *sigma, 
    double *dplus, double *lplus, double *work, cblas_int_t *info);

/* Subroutine */ int dlarrj_(cblas_int_t *n, double *d__, double *e2, 
    cblas_int_t *ifirst, cblas_int_t *ilast, double *rtol, cblas_int_t *offset, 
    double *w, double *werr, double *work, cblas_int_t *iwork, 
    double *pivmin, double *spdiam, cblas_int_t *info);

/* Subroutine */ int dlarrk_(cblas_int_t *n, cblas_int_t *iw, double *gl, 
    double *gu, double *d__, double *e2, double *pivmin, 
    double *reltol, double *w, double *werr, cblas_int_t *info);

/* Subroutine */ int dlarrr_(cblas_int_t *n, double *d__, double *e, 
    cblas_int_t *info);

/* Subroutine */ int dlarrv_(cblas_int_t *n, double *vl, double *vu, 
    double *d__, double *l, double *pivmin, cblas_int_t *isplit, 
    cblas_int_t *m, cblas_int_t *dol, cblas_int_t *dou, double *minrgp, 
    double *rtol1, double *rtol2, double *w, double *werr,
     double *wgap, cblas_int_t *iblock, cblas_int_t *indexw, double *gers,
     double *z__, cblas_int_t *ldz, cblas_int_t *isuppz, double *work, 
    cblas_int_t *iwork, cblas_int_t *info);

/* Subroutine */ int dlartg_(double *f, double *g, double *cs, 
    double *sn, double *r__);

/* Subroutine */ int dlartv_(cblas_int_t *n, double *x, cblas_int_t *incx, 
    double *y, cblas_int_t *incy, double *c__, double *s, cblas_int_t 
    *incc);

/* Subroutine */ int dlaruv_(cblas_int_t *iseed, cblas_int_t *n, double *x);

/* Subroutine */ int dlarz_(const char *side, cblas_int_t *m, cblas_int_t *n, cblas_int_t *l, 
    double *v, cblas_int_t *incv, double *tau, double *c__, 
    cblas_int_t *ldc, double *work);

/* Subroutine */ int dlarzb_(const char *side, const char *trans, const char *direct, const char *
    storev, cblas_int_t *m, cblas_int_t *n, cblas_int_t *k, cblas_int_t *l, double *v,
     cblas_int_t *ldv, double *t, cblas_int_t *ldt, double *c__, cblas_int_t *
    ldc, double *work, cblas_int_t *ldwork      );

/* Subroutine */ int dlarzt_(const char *direct, const char *storev, cblas_int_t *n, cblas_int_t *
    k, double *v, cblas_int_t *ldv, double *tau, double *t, 
    cblas_int_t *ldt);

/* Subroutine */ int dlas2_(double *f, double *g, double *h__, 
    double *ssmin, double *ssmax);

/* Subroutine */ int dlascl_(const char *type__, cblas_int_t *kl, cblas_int_t *ku, 
    double *cfrom, double *cto, cblas_int_t *m, cblas_int_t *n, 
    double *a, cblas_int_t *lda, cblas_int_t *info);

/* Subroutine */ int dlasd0_(cblas_int_t *n, cblas_int_t *sqre, double *d__, 
    double *e, double *u, cblas_int_t *ldu, double *vt, cblas_int_t *
    ldvt, cblas_int_t *smlsiz, cblas_int_t *iwork, double *work, cblas_int_t *
    info);

/* Subroutine */ int dlasd1_(cblas_int_t *nl, cblas_int_t *nr, cblas_int_t *sqre, 
    double *d__, double *alpha, double *beta, double *u, 
    cblas_int_t *ldu, double *vt, cblas_int_t *ldvt, cblas_int_t *idxq, cblas_int_t *
    iwork, double *work, cblas_int_t *info);

/* Subroutine */ int dlasd2_(cblas_int_t *nl, cblas_int_t *nr, cblas_int_t *sqre, cblas_int_t 
    *k, double *d__, double *z__, double *alpha, double *
    beta, double *u, cblas_int_t *ldu, double *vt, cblas_int_t *ldvt, 
    double *dsigma, double *u2, cblas_int_t *ldu2, double *vt2, 
    cblas_int_t *ldvt2, cblas_int_t *idxp, cblas_int_t *idx, cblas_int_t *idxc, cblas_int_t *
    idxq, cblas_int_t *coltyp, cblas_int_t *info);

/* Subroutine */ int dlasd3_(cblas_int_t *nl, cblas_int_t *nr, cblas_int_t *sqre, cblas_int_t 
    *k, double *d__, double *q, cblas_int_t *ldq, double *dsigma, 
    double *u, cblas_int_t *ldu, double *u2, cblas_int_t *ldu2, 
    double *vt, cblas_int_t *ldvt, double *vt2, cblas_int_t *ldvt2, 
    cblas_int_t *idxc, cblas_int_t *ctot, double *z__, cblas_int_t *info);

/* Subroutine */ int dlasd4_(cblas_int_t *n, cblas_int_t *i__, double *d__, 
    double *z__, double *delta, double *rho, double *
    sigma, double *work, cblas_int_t *info);

/* Subroutine */ int dlasd5_(cblas_int_t *i__, double *d__, double *z__, 
    double *delta, double *rho, double *dsigma, double *
    work);

/* Subroutine */ int dlasd6_(cblas_int_t *icompq, cblas_int_t *nl, cblas_int_t *nr, 
    cblas_int_t *sqre, double *d__, double *vf, double *vl, 
    double *alpha, double *beta, cblas_int_t *idxq, cblas_int_t *perm, 
    cblas_int_t *givptr, cblas_int_t *givcol, cblas_int_t *ldgcol, double *givnum,
     cblas_int_t *ldgnum, double *poles, double *difl, double *
    difr, double *z__, cblas_int_t *k, double *c__, double *s, 
    double *work, cblas_int_t *iwork, cblas_int_t *info);

/* Subroutine */ int dlasd7_(cblas_int_t *icompq, cblas_int_t *nl, cblas_int_t *nr, 
    cblas_int_t *sqre, cblas_int_t *k, double *d__, double *z__, 
    double *zw, double *vf, double *vfw, double *vl, 
    double *vlw, double *alpha, double *beta, double *
    dsigma, cblas_int_t *idx, cblas_int_t *idxp, cblas_int_t *idxq, cblas_int_t *perm, 
    cblas_int_t *givptr, cblas_int_t *givcol, cblas_int_t *ldgcol, double *givnum,
     cblas_int_t *ldgnum, double *c__, double *s, cblas_int_t *info);

/* Subroutine */ int dlasd8_(cblas_int_t *icompq, cblas_int_t *k, double *d__, 
    double *z__, double *vf, double *vl, double *difl, 
    double *difr, cblas_int_t *lddifr, double *dsigma, double *
    work, cblas_int_t *info);

/* Subroutine */ int dlasda_(cblas_int_t *icompq, cblas_int_t *smlsiz, cblas_int_t *n, 
    cblas_int_t *sqre, double *d__, double *e, double *u, cblas_int_t 
    *ldu, double *vt, cblas_int_t *k, double *difl, double *difr, 
    double *z__, double *poles, cblas_int_t *givptr, cblas_int_t *givcol, 
    cblas_int_t *ldgcol, cblas_int_t *perm, double *givnum, double *c__, 
    double *s, double *work, cblas_int_t *iwork, cblas_int_t *info);

/* Subroutine */ int dlasdq_(const char *uplo, cblas_int_t *sqre, cblas_int_t *n, cblas_int_t *
    ncvt, cblas_int_t *nru, cblas_int_t *ncc, double *d__, double *e, 
    double *vt, cblas_int_t *ldvt, double *u, cblas_int_t *ldu, 
    double *c__, cblas_int_t *ldc, double *work, cblas_int_t *info);

/* Subroutine */ int dlasdt_(cblas_int_t *n, cblas_int_t *lvl, cblas_int_t *nd, cblas_int_t *
    inode, cblas_int_t *ndiml, cblas_int_t *ndimr, cblas_int_t *msub);

/* Subroutine */ int dlaset_(const char *uplo, cblas_int_t *m, cblas_int_t *n, double *
    alpha, double *beta, double *a, cblas_int_t *lda);

/* Subroutine */ int dlasq1_(cblas_int_t *n, double *d__, double *e, 
    double *work, cblas_int_t *info);

/* Subroutine */ int dlasq2_(cblas_int_t *n, double *z__, cblas_int_t *info);

/* Subroutine */ int dlasq3_(cblas_int_t *i0, cblas_int_t *n0, double *z__, 
    cblas_int_t *pp, double *dmin__, double *sigma, double *desig,
     double *qmax, cblas_int_t *nfail, cblas_int_t *iter, cblas_int_t *ndiv, 
    logical *ieee);

/* Subroutine */ int dlasq4_(cblas_int_t *i0, cblas_int_t *n0, double *z__, 
    cblas_int_t *pp, cblas_int_t *n0in, double *dmin__, double *dmin1, 
    double *dmin2, double *dn, double *dn1, double *dn2, 
    double *tau, cblas_int_t *ttype);

/* Subroutine */ int dlasq5_(cblas_int_t *i0, cblas_int_t *n0, double *z__, 
    cblas_int_t *pp, double *tau, double *dmin__, double *dmin1, 
    double *dmin2, double *dn, double *dnm1, double *dnm2,
     logical *ieee);

/* Subroutine */ int dlasq6_(cblas_int_t *i0, cblas_int_t *n0, double *z__, 
    cblas_int_t *pp, double *dmin__, double *dmin1, double *dmin2,
     double *dn, double *dnm1, double *dnm2);

/* Subroutine */ int dlasr_(const char *side, const char *pivot, const char *direct, cblas_int_t *m,
     cblas_int_t *n, double *c__, double *s, double *a, cblas_int_t *
    lda);

/* Subroutine */ int dlasrt_(const char *id, cblas_int_t *n, double *d__, cblas_int_t *
    info);

/* Subroutine */ int dlassq_(cblas_int_t *n, double *x, cblas_int_t *incx, 
    double *scale, double *sumsq);

/* Subroutine */ int dlasv2_(double *f, double *g, double *h__, 
    double *ssmin, double *ssmax, double *snr, double *
    csr, double *snl, double *csl);

/* Subroutine */ int dlaswp_(cblas_int_t *n, double *a, cblas_int_t *lda, cblas_int_t 
    *k1, cblas_int_t *k2, cblas_int_t *ipiv, cblas_int_t *incx);

/* Subroutine */ int dlasy2_(logical *ltranl, logical *ltranr, cblas_int_t *isgn, 
    cblas_int_t *n1, cblas_int_t *n2, double *tl, cblas_int_t *ldtl, double *
    tr, cblas_int_t *ldtr, double *b, cblas_int_t *ldb, double *scale, 
    double *x, cblas_int_t *ldx, double *xnorm, cblas_int_t *info);

/* Subroutine */ int dlasyf_(const char *uplo, cblas_int_t *n, cblas_int_t *nb, cblas_int_t *kb,
     double *a, cblas_int_t *lda, cblas_int_t *ipiv, double *w, cblas_int_t *
    ldw, cblas_int_t *info);

/* Subroutine */ int dlatbs_(const char *uplo, const char *trans, const char *diag, const char *
    normin, cblas_int_t *n, cblas_int_t *kd, double *ab, cblas_int_t *ldab, 
    double *x, double *scale, double *cnorm, cblas_int_t *info);

/* Subroutine */ int dlatdf_(cblas_int_t *ijob, cblas_int_t *n, double *z__, 
    cblas_int_t *ldz, double *rhs, double *rdsum, double *rdscal, 
    cblas_int_t *ipiv, cblas_int_t *jpiv);

/* Subroutine */ int dlatps_(const char *uplo, const char *trans, const char *diag, const char *
    normin, cblas_int_t *n, double *ap, double *x, double *scale, 
    double *cnorm, cblas_int_t *info);

/* Subroutine */ int dlatrd_(const char *uplo, cblas_int_t *n, cblas_int_t *nb, double *
    a, cblas_int_t *lda, double *e, double *tau, double *w, 
    cblas_int_t *ldw);

/* Subroutine */ int dlatrs_(const char *uplo, const char *trans, const char *diag, const char *
    normin, cblas_int_t *n, double *a, cblas_int_t *lda, double *x, 
    double *scale, double *cnorm, cblas_int_t *info);

/* Subroutine */ int dlatrz_(cblas_int_t *m, cblas_int_t *n, cblas_int_t *l, double *
    a, cblas_int_t *lda, double *tau, double *work);

/* Subroutine */ int dlatzm_(const char *side, cblas_int_t *m, cblas_int_t *n, double *
    v, cblas_int_t *incv, double *tau, double *c1, double *c2, 
    cblas_int_t *ldc, double *work);

/* Subroutine */ int dlauu2_(const char *uplo, cblas_int_t *n, double *a, cblas_int_t *
    lda, cblas_int_t *info);

/* Subroutine */ int dlauum_(const char *uplo, cblas_int_t *n, double *a, cblas_int_t *
    lda, cblas_int_t *info);

/* Subroutine */ int dlazq3_(cblas_int_t *i0, cblas_int_t *n0, double *z__, 
    cblas_int_t *pp, double *dmin__, double *sigma, double *desig,
     double *qmax, cblas_int_t *nfail, cblas_int_t *iter, cblas_int_t *ndiv, 
    logical *ieee, cblas_int_t *ttype, double *dmin1, double *dmin2, 
    double *dn, double *dn1, double *dn2, double *tau);

/* Subroutine */ int dlazq4_(cblas_int_t *i0, cblas_int_t *n0, double *z__, 
    cblas_int_t *pp, cblas_int_t *n0in, double *dmin__, double *dmin1, 
    double *dmin2, double *dn, double *dn1, double *dn2, 
    double *tau, cblas_int_t *ttype, double *g);

/* Subroutine */ int dopgtr_(const char *uplo, cblas_int_t *n, double *ap, 
    double *tau, double *q, cblas_int_t *ldq, double *work, 
    cblas_int_t *info);

/* Subroutine */ int dopmtr_(const char *side, const char *uplo, const char *trans, cblas_int_t *m, 
    cblas_int_t *n, double *ap, double *tau, double *c__, cblas_int_t 
    *ldc, double *work, cblas_int_t *info);

/* Subroutine */ int dorg2l_(cblas_int_t *m, cblas_int_t *n, cblas_int_t *k, double *
    a, cblas_int_t *lda, double *tau, double *work, cblas_int_t *info);

/* Subroutine */ int dorg2r_(cblas_int_t *m, cblas_int_t *n, cblas_int_t *k, double *
    a, cblas_int_t *lda, double *tau, double *work, cblas_int_t *info);

/* Subroutine */ int dorgbr_(const char *vect, cblas_int_t *m, cblas_int_t *n, cblas_int_t *k, 
    double *a, cblas_int_t *lda, double *tau, double *work, 
    cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int dorghr_(cblas_int_t *n, cblas_int_t *ilo, cblas_int_t *ihi, 
    double *a, cblas_int_t *lda, double *tau, double *work, 
    cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int dorgl2_(cblas_int_t *m, cblas_int_t *n, cblas_int_t *k, double *
    a, cblas_int_t *lda, double *tau, double *work, cblas_int_t *info);

/* Subroutine */ int dorglq_(cblas_int_t *m, cblas_int_t *n, cblas_int_t *k, double *
    a, cblas_int_t *lda, double *tau, double *work, cblas_int_t *lwork, 
    cblas_int_t *info);

/* Subroutine */ int dorgql_(cblas_int_t *m, cblas_int_t *n, cblas_int_t *k, double *
    a, cblas_int_t *lda, double *tau, double *work, cblas_int_t *lwork, 
    cblas_int_t *info);

/* Subroutine */ int dorgqr_(cblas_int_t *m, cblas_int_t *n, cblas_int_t *k, double *
    a, cblas_int_t *lda, double *tau, double *work, cblas_int_t *lwork, 
    cblas_int_t *info);

/* Subroutine */ int dorgr2_(cblas_int_t *m, cblas_int_t *n, cblas_int_t *k, double *
    a, cblas_int_t *lda, double *tau, double *work, cblas_int_t *info);

/* Subroutine */ int dorgrq_(cblas_int_t *m, cblas_int_t *n, cblas_int_t *k, double *
    a, cblas_int_t *lda, double *tau, double *work, cblas_int_t *lwork, 
    cblas_int_t *info);

/* Subroutine */ int dorgtr_(const char *uplo, cblas_int_t *n, double *a, cblas_int_t *
    lda, double *tau, double *work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int dorm2l_(const char *side, const char *trans, cblas_int_t *m, cblas_int_t *n, 
    cblas_int_t *k, double *a, cblas_int_t *lda, double *tau, double *
    c__, cblas_int_t *ldc, double *work, cblas_int_t *info);

/* Subroutine */ int dorm2r_(const char *side, const char *trans, cblas_int_t *m, cblas_int_t *n, 
    cblas_int_t *k, double *a, cblas_int_t *lda, double *tau, double *
    c__, cblas_int_t *ldc, double *work, cblas_int_t *info);

/* Subroutine */ int dormbr_(const char *vect, const char *side, const char *trans, cblas_int_t *m, 
    cblas_int_t *n, cblas_int_t *k, double *a, cblas_int_t *lda, double *tau, 
    double *c__, cblas_int_t *ldc, double *work, cblas_int_t *lwork, 
    cblas_int_t *info);

/* Subroutine */ int dormhr_(const char *side, const char *trans, cblas_int_t *m, cblas_int_t *n, 
    cblas_int_t *ilo, cblas_int_t *ihi, double *a, cblas_int_t *lda, double *
    tau, double *c__, cblas_int_t *ldc, double *work, cblas_int_t *lwork, 
    cblas_int_t *info);

/* Subroutine */ int dorml2_(const char *side, const char *trans, cblas_int_t *m, cblas_int_t *n, 
    cblas_int_t *k, double *a, cblas_int_t *lda, double *tau, double *
    c__, cblas_int_t *ldc, double *work, cblas_int_t *info);

/* Subroutine */ int dormlq_(const char *side, const char *trans, cblas_int_t *m, cblas_int_t *n, 
    cblas_int_t *k, double *a, cblas_int_t *lda, double *tau, double *
    c__, cblas_int_t *ldc, double *work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int dormql_(const char *side, const char *trans, cblas_int_t *m, cblas_int_t *n, 
    cblas_int_t *k, double *a, cblas_int_t *lda, double *tau, double *
    c__, cblas_int_t *ldc, double *work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int dormqr_(const char *side, const char *trans, cblas_int_t *m, cblas_int_t *n, 
    cblas_int_t *k, double *a, cblas_int_t *lda, double *tau, double *
    c__, cblas_int_t *ldc, double *work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int dormr2_(const char *side, const char *trans, cblas_int_t *m, cblas_int_t *n, 
    cblas_int_t *k, double *a, cblas_int_t *lda, double *tau, double *
    c__, cblas_int_t *ldc, double *work, cblas_int_t *info);

/* Subroutine */ int dormr3_(const char *side, const char *trans, cblas_int_t *m, cblas_int_t *n, 
    cblas_int_t *k, cblas_int_t *l, double *a, cblas_int_t *lda, double *tau, 
    double *c__, cblas_int_t *ldc, double *work, cblas_int_t *info);

/* Subroutine */ int dormrq_(const char *side, const char *trans, cblas_int_t *m, cblas_int_t *n, 
    cblas_int_t *k, double *a, cblas_int_t *lda, double *tau, double *
    c__, cblas_int_t *ldc, double *work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int dormrz_(const char *side, const char *trans, cblas_int_t *m, cblas_int_t *n, 
    cblas_int_t *k, cblas_int_t *l, double *a, cblas_int_t *lda, double *tau, 
    double *c__, cblas_int_t *ldc, double *work, cblas_int_t *lwork, 
    cblas_int_t *info);

/* Subroutine */ int dormtr_(const char *side, const char *uplo, const char *trans, cblas_int_t *m, 
    cblas_int_t *n, double *a, cblas_int_t *lda, double *tau, double *
    c__, cblas_int_t *ldc, double *work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int dpbcon_(const char *uplo, cblas_int_t *n, cblas_int_t *kd, double *
    ab, cblas_int_t *ldab, double *anorm, double *rcond, double *
    work, cblas_int_t *iwork, cblas_int_t *info);

/* Subroutine */ int dpbequ_(const char *uplo, cblas_int_t *n, cblas_int_t *kd, double *
    ab, cblas_int_t *ldab, double *s, double *scond, double *amax,
     cblas_int_t *info);

/* Subroutine */ int dpbrfs_(const char *uplo, cblas_int_t *n, cblas_int_t *kd, cblas_int_t *
    nrhs, double *ab, cblas_int_t *ldab, double *afb, cblas_int_t *ldafb, 
    double *b, cblas_int_t *ldb, double *x, cblas_int_t *ldx, double *
    ferr, double *berr, double *work, cblas_int_t *iwork, cblas_int_t *
    info);

/* Subroutine */ int dpbstf_(const char *uplo, cblas_int_t *n, cblas_int_t *kd, double *
    ab, cblas_int_t *ldab, cblas_int_t *info);

/* Subroutine */ int dpbsv_(const char *uplo, cblas_int_t *n, cblas_int_t *kd, cblas_int_t *
    nrhs, double *ab, cblas_int_t *ldab, double *b, cblas_int_t *ldb, 
    cblas_int_t *info);

/* Subroutine */ int dpbsvx_(const char *fact, const char *uplo, cblas_int_t *n, cblas_int_t *kd, 
    cblas_int_t *nrhs, double *ab, cblas_int_t *ldab, double *afb, 
    cblas_int_t *ldafb, const char *equed, double *s, double *b, cblas_int_t *
    ldb, double *x, cblas_int_t *ldx, double *rcond, double *ferr,
     double *berr, double *work, cblas_int_t *iwork, cblas_int_t *info);

/* Subroutine */ int dpbtf2_(const char *uplo, cblas_int_t *n, cblas_int_t *kd, double *
    ab, cblas_int_t *ldab, cblas_int_t *info);

/* Subroutine */ int dpbtrf_(const char *uplo, cblas_int_t *n, cblas_int_t *kd, double *
    ab, cblas_int_t *ldab, cblas_int_t *info);

/* Subroutine */ int dpbtrs_(const char *uplo, cblas_int_t *n, cblas_int_t *kd, cblas_int_t *
    nrhs, double *ab, cblas_int_t *ldab, double *b, cblas_int_t *ldb, 
    cblas_int_t *info);

/* Subroutine */ int dpocon_(const char *uplo, cblas_int_t *n, double *a, cblas_int_t *
    lda, double *anorm, double *rcond, double *work, cblas_int_t *
    iwork, cblas_int_t *info);

/* Subroutine */ int dpoequ_(cblas_int_t *n, double *a, cblas_int_t *lda, 
    double *s, double *scond, double *amax, cblas_int_t *info);

/* Subroutine */ int dporfs_(const char *uplo, cblas_int_t *n, cblas_int_t *nrhs, 
    double *a, cblas_int_t *lda, double *af, cblas_int_t *ldaf, 
    double *b, cblas_int_t *ldb, double *x, cblas_int_t *ldx, double *
    ferr, double *berr, double *work, cblas_int_t *iwork, cblas_int_t *
    info);

/* Subroutine */ int dposv_(const char *uplo, cblas_int_t *n, cblas_int_t *nrhs, double 
    *a, cblas_int_t *lda, double *b, cblas_int_t *ldb, cblas_int_t *info);

/* Subroutine */ int dposvx_(const char *fact, const char *uplo, cblas_int_t *n, cblas_int_t *
    nrhs, double *a, cblas_int_t *lda, double *af, cblas_int_t *ldaf, 
    const char *equed, double *s, double *b, cblas_int_t *ldb, double *
    x, cblas_int_t *ldx, double *rcond, double *ferr, double *
    berr, double *work, cblas_int_t *iwork, cblas_int_t *info);

/* Subroutine */ int dpotf2_(const char *uplo, cblas_int_t *n, double *a, cblas_int_t *
    lda, cblas_int_t *info);

/* Subroutine */ int dpotrf_(const char *uplo, cblas_int_t *n, double *a, cblas_int_t *
    lda, cblas_int_t *info);

/* Subroutine */ int dpotri_(const char *uplo, cblas_int_t *n, double *a, cblas_int_t *
    lda, cblas_int_t *info);

/* Subroutine */ int dpotrs_(const char *uplo, cblas_int_t *n, cblas_int_t *nrhs, 
    double *a, cblas_int_t *lda, double *b, cblas_int_t *ldb, cblas_int_t *
    info);

/* Subroutine */ int dppcon_(const char *uplo, cblas_int_t *n, double *ap, 
    double *anorm, double *rcond, double *work, cblas_int_t *
    iwork, cblas_int_t *info);

/* Subroutine */ int dppequ_(const char *uplo, cblas_int_t *n, double *ap, 
    double *s, double *scond, double *amax, cblas_int_t *info);

/* Subroutine */ int dpprfs_(const char *uplo, cblas_int_t *n, cblas_int_t *nrhs, 
    double *ap, double *afp, double *b, cblas_int_t *ldb, 
    double *x, cblas_int_t *ldx, double *ferr, double *berr, 
    double *work, cblas_int_t *iwork, cblas_int_t *info);

/* Subroutine */ int dppsv_(const char *uplo, cblas_int_t *n, cblas_int_t *nrhs, double 
    *ap, double *b, cblas_int_t *ldb, cblas_int_t *info);

/* Subroutine */ int dppsvx_(const char *fact, const char *uplo, cblas_int_t *n, cblas_int_t *
    nrhs, double *ap, double *afp, const char *equed, double *s, 
    double *b, cblas_int_t *ldb, double *x, cblas_int_t *ldx, double *
    rcond, double *ferr, double *berr, double *work, cblas_int_t *
    iwork, cblas_int_t *info);

/* Subroutine */ int dpptrf_(const char *uplo, cblas_int_t *n, double *ap, cblas_int_t *
    info);

/* Subroutine */ int dpptri_(const char *uplo, cblas_int_t *n, double *ap, cblas_int_t *
    info);

/* Subroutine */ int dpptrs_(const char *uplo, cblas_int_t *n, cblas_int_t *nrhs, 
    double *ap, double *b, cblas_int_t *ldb, cblas_int_t *info);

/* Subroutine */ int dptcon_(cblas_int_t *n, double *d__, double *e, 
    double *anorm, double *rcond, double *work, cblas_int_t *info);

/* Subroutine */ int dpteqr_(const char *compz, cblas_int_t *n, double *d__, 
    double *e, double *z__, cblas_int_t *ldz, double *work, 
    cblas_int_t *info);

/* Subroutine */ int dptrfs_(cblas_int_t *n, cblas_int_t *nrhs, double *d__, 
    double *e, double *df, double *ef, double *b, cblas_int_t 
    *ldb, double *x, cblas_int_t *ldx, double *ferr, double *berr,
     double *work, cblas_int_t *info);

/* Subroutine */ int dptsv_(cblas_int_t *n, cblas_int_t *nrhs, double *d__, 
    double *e, double *b, cblas_int_t *ldb, cblas_int_t *info);

/* Subroutine */ int dptsvx_(const char *fact, cblas_int_t *n, cblas_int_t *nrhs, 
    double *d__, double *e, double *df, double *ef, 
    double *b, cblas_int_t *ldb, double *x, cblas_int_t *ldx, double *
    rcond, double *ferr, double *berr, double *work, cblas_int_t *
    info);

/* Subroutine */ int dpttrf_(cblas_int_t *n, double *d__, double *e, 
    cblas_int_t *info);

/* Subroutine */ int dpttrs_(cblas_int_t *n, cblas_int_t *nrhs, double *d__, 
    double *e, double *b, cblas_int_t *ldb, cblas_int_t *info);

/* Subroutine */ int dptts2_(cblas_int_t *n, cblas_int_t *nrhs, double *d__, 
    double *e, double *b, cblas_int_t *ldb);

/* Subroutine */ int drscl_(cblas_int_t *n, double *sa, double *sx, 
    cblas_int_t *incx);

/* Subroutine */ int dsbev_(const char *jobz, const char *uplo, cblas_int_t *n, cblas_int_t *kd, 
    double *ab, cblas_int_t *ldab, double *w, double *z__, 
    cblas_int_t *ldz, double *work, cblas_int_t *info);

/* Subroutine */ int dsbevd_(const char *jobz, const char *uplo, cblas_int_t *n, cblas_int_t *kd, 
    double *ab, cblas_int_t *ldab, double *w, double *z__, 
    cblas_int_t *ldz, double *work, cblas_int_t *lwork, cblas_int_t *iwork, 
    cblas_int_t *liwork, cblas_int_t *info);

/* Subroutine */ int dsbevx_(const char *jobz, const char *range, const char *uplo, cblas_int_t *n, 
    cblas_int_t *kd, double *ab, cblas_int_t *ldab, double *q, cblas_int_t *
    ldq, double *vl, double *vu, cblas_int_t *il, cblas_int_t *iu, 
    double *abstol, cblas_int_t *m, double *w, double *z__, 
    cblas_int_t *ldz, double *work, cblas_int_t *iwork, cblas_int_t *ifail, 
    cblas_int_t *info);

/* Subroutine */ int dsbgst_(const char *vect, const char *uplo, cblas_int_t *n, cblas_int_t *ka, 
    cblas_int_t *kb, double *ab, cblas_int_t *ldab, double *bb, cblas_int_t *
    ldbb, double *x, cblas_int_t *ldx, double *work, cblas_int_t *info);

/* Subroutine */ int dsbgv_(const char *jobz, const char *uplo, cblas_int_t *n, cblas_int_t *ka, 
    cblas_int_t *kb, double *ab, cblas_int_t *ldab, double *bb, cblas_int_t *
    ldbb, double *w, double *z__, cblas_int_t *ldz, double *work, 
    cblas_int_t *info);

/* Subroutine */ int dsbgvd_(const char *jobz, const char *uplo, cblas_int_t *n, cblas_int_t *ka, 
    cblas_int_t *kb, double *ab, cblas_int_t *ldab, double *bb, cblas_int_t *
    ldbb, double *w, double *z__, cblas_int_t *ldz, double *work, 
    cblas_int_t *lwork, cblas_int_t *iwork, cblas_int_t *liwork, cblas_int_t *info);

/* Subroutine */ int dsbgvx_(const char *jobz, const char *range, const char *uplo, cblas_int_t *n, 
    cblas_int_t *ka, cblas_int_t *kb, double *ab, cblas_int_t *ldab, double *
    bb, cblas_int_t *ldbb, double *q, cblas_int_t *ldq, double *vl, 
    double *vu, cblas_int_t *il, cblas_int_t *iu, double *abstol, cblas_int_t 
    *m, double *w, double *z__, cblas_int_t *ldz, double *work, 
    cblas_int_t *iwork, cblas_int_t *ifail, cblas_int_t *info);

/* Subroutine */ int dsbtrd_(const char *vect, const char *uplo, cblas_int_t *n, cblas_int_t *kd, 
    double *ab, cblas_int_t *ldab, double *d__, double *e, 
    double *q, cblas_int_t *ldq, double *work, cblas_int_t *info);

/* Subroutine */ int dsgesv_(cblas_int_t *n, cblas_int_t *nrhs, double *a, 
    cblas_int_t *lda, cblas_int_t *ipiv, double *b, cblas_int_t *ldb, double *
    x, cblas_int_t *ldx, double *work, float *swork, cblas_int_t *iter, 
    cblas_int_t *info);

/* Subroutine */ int dspcon_(const char *uplo, cblas_int_t *n, double *ap, cblas_int_t *
    ipiv, double *anorm, double *rcond, double *work, cblas_int_t 
    *iwork, cblas_int_t *info);

/* Subroutine */ int dspev_(const char *jobz, const char *uplo, cblas_int_t *n, double *
    ap, double *w, double *z__, cblas_int_t *ldz, double *work, 
    cblas_int_t *info);

/* Subroutine */ int dspevd_(const char *jobz, const char *uplo, cblas_int_t *n, double *
    ap, double *w, double *z__, cblas_int_t *ldz, double *work, 
    cblas_int_t *lwork, cblas_int_t *iwork, cblas_int_t *liwork, cblas_int_t *info);

/* Subroutine */ int dspevx_(const char *jobz, const char *range, const char *uplo, cblas_int_t *n, 
    double *ap, double *vl, double *vu, cblas_int_t *il, cblas_int_t *
    iu, double *abstol, cblas_int_t *m, double *w, double *z__, 
    cblas_int_t *ldz, double *work, cblas_int_t *iwork, cblas_int_t *ifail, 
    cblas_int_t *info);

/* Subroutine */ int dspgst_(cblas_int_t *itype, const char *uplo, cblas_int_t *n, 
    double *ap, double *bp, cblas_int_t *info);

/* Subroutine */ int dspgv_(cblas_int_t *itype, const char *jobz, const char *uplo, cblas_int_t *
    n, double *ap, double *bp, double *w, double *z__, 
    cblas_int_t *ldz, double *work, cblas_int_t *info);

/* Subroutine */ int dspgvd_(cblas_int_t *itype, const char *jobz, const char *uplo, cblas_int_t *
    n, double *ap, double *bp, double *w, double *z__, 
    cblas_int_t *ldz, double *work, cblas_int_t *lwork, cblas_int_t *iwork, 
    cblas_int_t *liwork, cblas_int_t *info);

/* Subroutine */ int dspgvx_(cblas_int_t *itype, const char *jobz, const char *range, const char *
    uplo, cblas_int_t *n, double *ap, double *bp, double *vl, 
    double *vu, cblas_int_t *il, cblas_int_t *iu, double *abstol, cblas_int_t 
    *m, double *w, double *z__, cblas_int_t *ldz, double *work, 
    cblas_int_t *iwork, cblas_int_t *ifail, cblas_int_t *info);

/* Subroutine */ int dsprfs_(const char *uplo, cblas_int_t *n, cblas_int_t *nrhs, 
    double *ap, double *afp, cblas_int_t *ipiv, double *b, 
    cblas_int_t *ldb, double *x, cblas_int_t *ldx, double *ferr, 
    double *berr, double *work, cblas_int_t *iwork, cblas_int_t *info);

/* Subroutine */ int dspsv_(const char *uplo, cblas_int_t *n, cblas_int_t *nrhs, double 
    *ap, cblas_int_t *ipiv, double *b, cblas_int_t *ldb, cblas_int_t *info);

/* Subroutine */ int dspsvx_(const char *fact, const char *uplo, cblas_int_t *n, cblas_int_t *
    nrhs, double *ap, double *afp, cblas_int_t *ipiv, double *b, 
    cblas_int_t *ldb, double *x, cblas_int_t *ldx, double *rcond, 
    double *ferr, double *berr, double *work, cblas_int_t *iwork, 
    cblas_int_t *info);

/* Subroutine */ int dsptrd_(const char *uplo, cblas_int_t *n, double *ap, 
    double *d__, double *e, double *tau, cblas_int_t *info);

/* Subroutine */ int dsptrf_(const char *uplo, cblas_int_t *n, double *ap, cblas_int_t *
    ipiv, cblas_int_t *info);

/* Subroutine */ int dsptri_(const char *uplo, cblas_int_t *n, double *ap, cblas_int_t *
    ipiv, double *work, cblas_int_t *info);

/* Subroutine */ int dsptrs_(const char *uplo, cblas_int_t *n, cblas_int_t *nrhs, 
    double *ap, cblas_int_t *ipiv, double *b, cblas_int_t *ldb, cblas_int_t *
    info);

/* Subroutine */ int dstebz_(const char *range, const char *order, cblas_int_t *n, double 
    *vl, double *vu, cblas_int_t *il, cblas_int_t *iu, double *abstol, 
    double *d__, double *e, cblas_int_t *m, cblas_int_t *nsplit, 
    double *w, cblas_int_t *iblock, cblas_int_t *isplit, double *work, 
    cblas_int_t *iwork, cblas_int_t *info);

/* Subroutine */ int dstedc_(const char *compz, cblas_int_t *n, double *d__, 
    double *e, double *z__, cblas_int_t *ldz, double *work, 
    cblas_int_t *lwork, cblas_int_t *iwork, cblas_int_t *liwork, cblas_int_t *info);

/* Subroutine */ int dstegr_(const char *jobz, const char *range, cblas_int_t *n, double *
    d__, double *e, double *vl, double *vu, cblas_int_t *il, 
    cblas_int_t *iu, double *abstol, cblas_int_t *m, double *w, 
    double *z__, cblas_int_t *ldz, cblas_int_t *isuppz, double *work, 
    cblas_int_t *lwork, cblas_int_t *iwork, cblas_int_t *liwork, cblas_int_t *info);

/* Subroutine */ int dstein_(cblas_int_t *n, double *d__, double *e, 
    cblas_int_t *m, double *w, cblas_int_t *iblock, cblas_int_t *isplit, 
    double *z__, cblas_int_t *ldz, double *work, cblas_int_t *iwork, 
    cblas_int_t *ifail, cblas_int_t *info);

/* Subroutine */ int dstemr_(const char *jobz, const char *range, cblas_int_t *n, double *
    d__, double *e, double *vl, double *vu, cblas_int_t *il, 
    cblas_int_t *iu, cblas_int_t *m, double *w, double *z__, cblas_int_t *ldz,
     cblas_int_t *nzc, cblas_int_t *isuppz, logical *tryrac, double *work, 
    cblas_int_t *lwork, cblas_int_t *iwork, cblas_int_t *liwork, cblas_int_t *info);

/* Subroutine */ int dsteqr_(const char *compz, cblas_int_t *n, double *d__, 
    double *e, double *z__, cblas_int_t *ldz, double *work, 
    cblas_int_t *info);

/* Subroutine */ int dsterf_(cblas_int_t *n, double *d__, double *e, 
    cblas_int_t *info);

/* Subroutine */ int dstev_(const char *jobz, cblas_int_t *n, double *d__, 
    double *e, double *z__, cblas_int_t *ldz, double *work, 
    cblas_int_t *info);

/* Subroutine */ int dstevd_(const char *jobz, cblas_int_t *n, double *d__, 
    double *e, double *z__, cblas_int_t *ldz, double *work, 
    cblas_int_t *lwork, cblas_int_t *iwork, cblas_int_t *liwork, cblas_int_t *info);

/* Subroutine */ int dstevr_(const char *jobz, const char *range, cblas_int_t *n, double *
    d__, double *e, double *vl, double *vu, cblas_int_t *il, 
    cblas_int_t *iu, double *abstol, cblas_int_t *m, double *w, 
    double *z__, cblas_int_t *ldz, cblas_int_t *isuppz, double *work, 
    cblas_int_t *lwork, cblas_int_t *iwork, cblas_int_t *liwork, cblas_int_t *info);

/* Subroutine */ int dstevx_(const char *jobz, const char *range, cblas_int_t *n, double *
    d__, double *e, double *vl, double *vu, cblas_int_t *il, 
    cblas_int_t *iu, double *abstol, cblas_int_t *m, double *w, 
    double *z__, cblas_int_t *ldz, double *work, cblas_int_t *iwork, 
    cblas_int_t *ifail, cblas_int_t *info);

/* Subroutine */ int dsycon_(const char *uplo, cblas_int_t *n, double *a, cblas_int_t *
    lda, cblas_int_t *ipiv, double *anorm, double *rcond, double *
    work, cblas_int_t *iwork, cblas_int_t *info);

/* Subroutine */ int dsyev_(const char *jobz, const char *uplo, cblas_int_t *n, double *a,
     cblas_int_t *lda, double *w, double *work, cblas_int_t *lwork, 
    cblas_int_t *info);

/* Subroutine */ int dsyevd_(const char *jobz, const char *uplo, cblas_int_t *n, double *
    a, cblas_int_t *lda, double *w, double *work, cblas_int_t *lwork, 
    cblas_int_t *iwork, cblas_int_t *liwork, cblas_int_t *info);

/* Subroutine */ int dsyevr_(const char *jobz, const char *range, const char *uplo, cblas_int_t *n, 
    double *a, cblas_int_t *lda, double *vl, double *vu, cblas_int_t *
    il, cblas_int_t *iu, double *abstol, cblas_int_t *m, double *w, 
    double *z__, cblas_int_t *ldz, cblas_int_t *isuppz, double *work, 
    cblas_int_t *lwork, cblas_int_t *iwork, cblas_int_t *liwork, cblas_int_t *info);

/* Subroutine */ int dsyevx_(const char *jobz, const char *range, const char *uplo, cblas_int_t *n, 
    double *a, cblas_int_t *lda, double *vl, double *vu, cblas_int_t *
    il, cblas_int_t *iu, double *abstol, cblas_int_t *m, double *w, 
    double *z__, cblas_int_t *ldz, double *work, cblas_int_t *lwork, 
    cblas_int_t *iwork, cblas_int_t *ifail, cblas_int_t *info);

/* Subroutine */ int dsygs2_(cblas_int_t *itype, const char *uplo, cblas_int_t *n, 
    double *a, cblas_int_t *lda, double *b, cblas_int_t *ldb, cblas_int_t *
    info);

/* Subroutine */ int dsygst_(cblas_int_t *itype, const char *uplo, cblas_int_t *n, 
    double *a, cblas_int_t *lda, double *b, cblas_int_t *ldb, cblas_int_t *
    info);

/* Subroutine */ int dsygv_(cblas_int_t *itype, const char *jobz, const char *uplo, cblas_int_t *
    n, double *a, cblas_int_t *lda, double *b, cblas_int_t *ldb, 
    double *w, double *work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int dsygvd_(cblas_int_t *itype, const char *jobz, const char *uplo, cblas_int_t *
    n, double *a, cblas_int_t *lda, double *b, cblas_int_t *ldb, 
    double *w, double *work, cblas_int_t *lwork, cblas_int_t *iwork, 
    cblas_int_t *liwork, cblas_int_t *info);

/* Subroutine */ int dsygvx_(cblas_int_t *itype, const char *jobz, const char *range, const char *
    uplo, cblas_int_t *n, double *a, cblas_int_t *lda, double *b, cblas_int_t 
    *ldb, double *vl, double *vu, cblas_int_t *il, cblas_int_t *iu, 
    double *abstol, cblas_int_t *m, double *w, double *z__, 
    cblas_int_t *ldz, double *work, cblas_int_t *lwork, cblas_int_t *iwork, 
    cblas_int_t *ifail, cblas_int_t *info);

/* Subroutine */ int dsyrfs_(const char *uplo, cblas_int_t *n, cblas_int_t *nrhs, 
    double *a, cblas_int_t *lda, double *af, cblas_int_t *ldaf, cblas_int_t *
    ipiv, double *b, cblas_int_t *ldb, double *x, cblas_int_t *ldx, 
    double *ferr, double *berr, double *work, cblas_int_t *iwork, 
    cblas_int_t *info);

/* Subroutine */ int dsysv_(const char *uplo, cblas_int_t *n, cblas_int_t *nrhs, double 
    *a, cblas_int_t *lda, cblas_int_t *ipiv, double *b, cblas_int_t *ldb, 
    double *work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int dsysvx_(const char *fact, const char *uplo, cblas_int_t *n, cblas_int_t *
    nrhs, double *a, cblas_int_t *lda, double *af, cblas_int_t *ldaf, 
    cblas_int_t *ipiv, double *b, cblas_int_t *ldb, double *x, cblas_int_t *
    ldx, double *rcond, double *ferr, double *berr, 
    double *work, cblas_int_t *lwork, cblas_int_t *iwork, cblas_int_t *info);

/* Subroutine */ int dsytd2_(const char *uplo, cblas_int_t *n, double *a, cblas_int_t *
    lda, double *d__, double *e, double *tau, cblas_int_t *info);

/* Subroutine */ int dsytf2_(const char *uplo, cblas_int_t *n, double *a, cblas_int_t *
    lda, cblas_int_t *ipiv, cblas_int_t *info);

/* Subroutine */ int dsytrd_(const char *uplo, cblas_int_t *n, double *a, cblas_int_t *
    lda, double *d__, double *e, double *tau, double *
    work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int dsytrf_(const char *uplo, cblas_int_t *n, double *a, cblas_int_t *
    lda, cblas_int_t *ipiv, double *work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int dsytri_(const char *uplo, cblas_int_t *n, double *a, cblas_int_t *
    lda, cblas_int_t *ipiv, double *work, cblas_int_t *info);

/* Subroutine */ int dsytrs_(const char *uplo, cblas_int_t *n, cblas_int_t *nrhs, 
    double *a, cblas_int_t *lda, cblas_int_t *ipiv, double *b, cblas_int_t *
    ldb, cblas_int_t *info);

/* Subroutine */ int dtbcon_(const char *norm, const char *uplo, const char *diag, cblas_int_t *n, 
    cblas_int_t *kd, double *ab, cblas_int_t *ldab, double *rcond, 
    double *work, cblas_int_t *iwork, cblas_int_t *info);

/* Subroutine */ int dtbrfs_(const char *uplo, const char *trans, const char *diag, cblas_int_t *n, 
    cblas_int_t *kd, cblas_int_t *nrhs, double *ab, cblas_int_t *ldab, double 
    *b, cblas_int_t *ldb, double *x, cblas_int_t *ldx, double *ferr, 
    double *berr, double *work, cblas_int_t *iwork, cblas_int_t *info);

/* Subroutine */ int dtbtrs_(const char *uplo, const char *trans, const char *diag, cblas_int_t *n, 
    cblas_int_t *kd, cblas_int_t *nrhs, double *ab, cblas_int_t *ldab, double 
    *b, cblas_int_t *ldb, cblas_int_t *info);

/* Subroutine */ int dtgevc_(const char *side, const char *howmny, logical *select, 
    cblas_int_t *n, double *s, cblas_int_t *lds, double *p, cblas_int_t *ldp, 
    double *vl, cblas_int_t *ldvl, double *vr, cblas_int_t *ldvr, cblas_int_t 
    *mm, cblas_int_t *m, double *work, cblas_int_t *info);

/* Subroutine */ int dtgex2_(logical *wantq, logical *wantz, cblas_int_t *n, 
    double *a, cblas_int_t *lda, double *b, cblas_int_t *ldb, double *
    q, cblas_int_t *ldq, double *z__, cblas_int_t *ldz, cblas_int_t *j1, cblas_int_t *
    n1, cblas_int_t *n2, double *work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int dtgexc_(logical *wantq, logical *wantz, cblas_int_t *n, 
    double *a, cblas_int_t *lda, double *b, cblas_int_t *ldb, double *
    q, cblas_int_t *ldq, double *z__, cblas_int_t *ldz, cblas_int_t *ifst, 
    cblas_int_t *ilst, double *work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int dtgsen_(cblas_int_t *ijob, logical *wantq, logical *wantz, 
    logical *select, cblas_int_t *n, double *a, cblas_int_t *lda, double *
    b, cblas_int_t *ldb, double *alphar, double *alphai, double *
    beta, double *q, cblas_int_t *ldq, double *z__, cblas_int_t *ldz, 
    cblas_int_t *m, double *pl, double *pr, double *dif, 
    double *work, cblas_int_t *lwork, cblas_int_t *iwork, cblas_int_t *liwork, 
    cblas_int_t *info);

/* Subroutine */ int dtgsja_(const char *jobu, const char *jobv, const char *jobq, cblas_int_t *m, 
    cblas_int_t *p, cblas_int_t *n, cblas_int_t *k, cblas_int_t *l, double *a, 
    cblas_int_t *lda, double *b, cblas_int_t *ldb, double *tola, 
    double *tolb, double *alpha, double *beta, double *u, 
    cblas_int_t *ldu, double *v, cblas_int_t *ldv, double *q, cblas_int_t *
    ldq, double *work, cblas_int_t *ncycle, cblas_int_t *info);

/* Subroutine */ int dtgsna_(const char *job, const char *howmny, logical *select, 
    cblas_int_t *n, double *a, cblas_int_t *lda, double *b, cblas_int_t *ldb, 
    double *vl, cblas_int_t *ldvl, double *vr, cblas_int_t *ldvr, 
    double *s, double *dif, cblas_int_t *mm, cblas_int_t *m, double *
    work, cblas_int_t *lwork, cblas_int_t *iwork, cblas_int_t *info);

/* Subroutine */ int dtgsy2_(const char *trans, cblas_int_t *ijob, cblas_int_t *m, cblas_int_t *
    n, double *a, cblas_int_t *lda, double *b, cblas_int_t *ldb, 
    double *c__, cblas_int_t *ldc, double *d__, cblas_int_t *ldd, 
    double *e, cblas_int_t *lde, double *f, cblas_int_t *ldf, double *
    scale, double *rdsum, double *rdscal, cblas_int_t *iwork, cblas_int_t 
    *pq, cblas_int_t *info);

/* Subroutine */ int dtgsyl_(const char *trans, cblas_int_t *ijob, cblas_int_t *m, cblas_int_t *
    n, double *a, cblas_int_t *lda, double *b, cblas_int_t *ldb, 
    double *c__, cblas_int_t *ldc, double *d__, cblas_int_t *ldd, 
    double *e, cblas_int_t *lde, double *f, cblas_int_t *ldf, double *
    scale, double *dif, double *work, cblas_int_t *lwork, cblas_int_t *
    iwork, cblas_int_t *info);

/* Subroutine */ int dtpcon_(const char *norm, const char *uplo, const char *diag, cblas_int_t *n, 
    double *ap, double *rcond, double *work, cblas_int_t *iwork, 
    cblas_int_t *info);

/* Subroutine */ int dtprfs_(const char *uplo, const char *trans, const char *diag, cblas_int_t *n, 
    cblas_int_t *nrhs, double *ap, double *b, cblas_int_t *ldb, 
    double *x, cblas_int_t *ldx, double *ferr, double *berr, 
    double *work, cblas_int_t *iwork, cblas_int_t *info);

/* Subroutine */ int dtptri_(const char *uplo, const char *diag, cblas_int_t *n, double *
    ap, cblas_int_t *info);

/* Subroutine */ int dtptrs_(const char *uplo, const char *trans, const char *diag, cblas_int_t *n, 
    cblas_int_t *nrhs, double *ap, double *b, cblas_int_t *ldb, cblas_int_t *
    info);

/* Subroutine */ int dtrcon_(const char *norm, const char *uplo, const char *diag, cblas_int_t *n, 
    double *a, cblas_int_t *lda, double *rcond, double *work, 
    cblas_int_t *iwork, cblas_int_t *info);

/* Subroutine */ int dtrevc_(const char *side, const char *howmny, logical *select, 
    cblas_int_t *n, double *t, cblas_int_t *ldt, double *vl, cblas_int_t *
    ldvl, double *vr, cblas_int_t *ldvr, cblas_int_t *mm, cblas_int_t *m, 
    double *work, cblas_int_t *info);

/* Subroutine */ int dtrexc_(const char *compq, cblas_int_t *n, double *t, cblas_int_t *
    ldt, double *q, cblas_int_t *ldq, cblas_int_t *ifst, cblas_int_t *ilst, 
    double *work, cblas_int_t *info);

/* Subroutine */ int dtrrfs_(const char *uplo, const char *trans, const char *diag, cblas_int_t *n, 
    cblas_int_t *nrhs, double *a, cblas_int_t *lda, double *b, cblas_int_t *
    ldb, double *x, cblas_int_t *ldx, double *ferr, double *berr, 
    double *work, cblas_int_t *iwork, cblas_int_t *info);

/* Subroutine */ int dtrsen_(const char *job, const char *compq, logical *select, cblas_int_t 
    *n, double *t, cblas_int_t *ldt, double *q, cblas_int_t *ldq, 
    double *wr, double *wi, cblas_int_t *m, double *s, double 
    *sep, double *work, cblas_int_t *lwork, cblas_int_t *iwork, cblas_int_t *
    liwork, cblas_int_t *info);

/* Subroutine */ int dtrsna_(const char *job, const char *howmny, logical *select, 
    cblas_int_t *n, double *t, cblas_int_t *ldt, double *vl, cblas_int_t *
    ldvl, double *vr, cblas_int_t *ldvr, double *s, double *sep, 
    cblas_int_t *mm, cblas_int_t *m, double *work, cblas_int_t *ldwork, cblas_int_t *
    iwork, cblas_int_t *info);

/* Subroutine */ int dtrsyl_(const char *trana, const char *tranb, cblas_int_t *isgn, cblas_int_t 
    *m, cblas_int_t *n, double *a, cblas_int_t *lda, double *b, cblas_int_t *
    ldb, double *c__, cblas_int_t *ldc, double *scale, cblas_int_t *info);

/* Subroutine */ int dtrti2_(const char *uplo, const char *diag, cblas_int_t *n, double *
    a, cblas_int_t *lda, cblas_int_t *info);

/* Subroutine */ int dtrtri_(const char *uplo, const char *diag, cblas_int_t *n, double *
    a, cblas_int_t *lda, cblas_int_t *info);

/* Subroutine */ int dtrtrs_(const char *uplo, const char *trans, const char *diag, cblas_int_t *n, 
    cblas_int_t *nrhs, double *a, cblas_int_t *lda, double *b, cblas_int_t *
    ldb, cblas_int_t *info);

/* Subroutine */ int dtzrqf_(cblas_int_t *m, cblas_int_t *n, double *a, cblas_int_t *
    lda, double *tau, cblas_int_t *info);

/* Subroutine */ int dtzrzf_(cblas_int_t *m, cblas_int_t *n, double *a, cblas_int_t *
    lda, double *tau, double *work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int ilaver_(cblas_int_t *vers_major__, cblas_int_t *vers_minor__, 
    cblas_int_t *vers_patch__);

/* Subroutine */ int dgesv_(cblas_int_t *n, cblas_int_t *nrhs, double *a, cblas_int_t 
    *lda, cblas_int_t *ipiv, double *b, cblas_int_t *ldb, cblas_int_t *info);

/* Subroutine */ int sbdsdc_(const char *uplo, const char *compq, cblas_int_t *n, float *d__, 
    float *e, float *u, cblas_int_t *ldu, float *vt, cblas_int_t *ldvt, float *q,   
    cblas_int_t *iq, float *work, cblas_int_t *iwork, cblas_int_t *info);

/* Subroutine */ int sbdsqr_(const char *uplo, cblas_int_t *n, cblas_int_t *ncvt, cblas_int_t *
    nru, cblas_int_t *ncc, float *d__, float *e, float *vt, cblas_int_t *ldvt, float *
    u, cblas_int_t *ldu, float *c__, cblas_int_t *ldc, float *work, cblas_int_t *info);

/* Subroutine */ int sdisna_(const char *job, cblas_int_t *m, cblas_int_t *n, float *d__, 
    float *sep, cblas_int_t *info);

/* Subroutine */ int sgbbrd_(const char *vect, cblas_int_t *m, cblas_int_t *n, cblas_int_t *ncc,
     cblas_int_t *kl, cblas_int_t *ku, float *ab, cblas_int_t *ldab, float *d__, float *
    e, float *q, cblas_int_t *ldq, float *pt, cblas_int_t *ldpt, float *c__, cblas_int_t 
    *ldc, float *work, cblas_int_t *info);

/* Subroutine */ int sgbcon_(const char *norm, cblas_int_t *n, cblas_int_t *kl, cblas_int_t *ku,
     float *ab, cblas_int_t *ldab, cblas_int_t *ipiv, float *anorm, float *rcond, 
    float *work, cblas_int_t *iwork, cblas_int_t *info);

/* Subroutine */ int sgbequ_(cblas_int_t *m, cblas_int_t *n, cblas_int_t *kl, cblas_int_t *ku,
     float *ab, cblas_int_t *ldab, float *r__, float *c__, float *rowcnd, float *
    colcnd, float *amax, cblas_int_t *info);

/* Subroutine */ int sgbrfs_(const char *trans, cblas_int_t *n, cblas_int_t *kl, cblas_int_t *
    ku, cblas_int_t *nrhs, float *ab, cblas_int_t *ldab, float *afb, cblas_int_t *ldafb,
     cblas_int_t *ipiv, float *b, cblas_int_t *ldb, float *x, cblas_int_t *ldx, float *
    ferr, float *berr, float *work, cblas_int_t *iwork, cblas_int_t *info);

/* Subroutine */ int sgbsv_(cblas_int_t *n, cblas_int_t *kl, cblas_int_t *ku, cblas_int_t *
    nrhs, float *ab, cblas_int_t *ldab, cblas_int_t *ipiv, float *b, cblas_int_t *ldb, 
    cblas_int_t *info);

/* Subroutine */ int sgbsvx_(const char *fact, const char *trans, cblas_int_t *n, cblas_int_t *kl,
     cblas_int_t *ku, cblas_int_t *nrhs, float *ab, cblas_int_t *ldab, float *afb, 
    cblas_int_t *ldafb, cblas_int_t *ipiv, const char *equed, float *r__, float *c__, 
    float *b, cblas_int_t *ldb, float *x, cblas_int_t *ldx, float *rcond, float *ferr,
     float *berr, float *work, cblas_int_t *iwork, cblas_int_t *info);

/* Subroutine */ int sgbtf2_(cblas_int_t *m, cblas_int_t *n, cblas_int_t *kl, cblas_int_t *ku,
     float *ab, cblas_int_t *ldab, cblas_int_t *ipiv, cblas_int_t *info);

/* Subroutine */ int sgbtrf_(cblas_int_t *m, cblas_int_t *n, cblas_int_t *kl, cblas_int_t *ku,
     float *ab, cblas_int_t *ldab, cblas_int_t *ipiv, cblas_int_t *info);

/* Subroutine */ int sgbtrs_(const char *trans, cblas_int_t *n, cblas_int_t *kl, cblas_int_t *
    ku, cblas_int_t *nrhs, float *ab, cblas_int_t *ldab, cblas_int_t *ipiv, float *b, 
    cblas_int_t *ldb, cblas_int_t *info);

/* Subroutine */ int sgebak_(const char *job, const char *side, cblas_int_t *n, cblas_int_t *ilo, 
    cblas_int_t *ihi, float *scale, cblas_int_t *m, float *v, cblas_int_t *ldv, cblas_int_t 
    *info);

/* Subroutine */ int sgebal_(const char *job, cblas_int_t *n, float *a, cblas_int_t *lda, 
    cblas_int_t *ilo, cblas_int_t *ihi, float *scale, cblas_int_t *info);

/* Subroutine */ int sgebd2_(cblas_int_t *m, cblas_int_t *n, float *a, cblas_int_t *lda, 
    float *d__, float *e, float *tauq, float *taup, float *work, cblas_int_t *info);

/* Subroutine */ int sgebrd_(cblas_int_t *m, cblas_int_t *n, float *a, cblas_int_t *lda, 
    float *d__, float *e, float *tauq, float *taup, float *work, cblas_int_t *
    lwork, cblas_int_t *info);

/* Subroutine */ int sgecon_(const char *norm, cblas_int_t *n, float *a, cblas_int_t *lda, 
    float *anorm, float *rcond, float *work, cblas_int_t *iwork, cblas_int_t *info);

/* Subroutine */ int sgeequ_(cblas_int_t *m, cblas_int_t *n, float *a, cblas_int_t *lda, 
    float *r__, float *c__, float *rowcnd, float *colcnd, float *amax, cblas_int_t 
    *info);

/* Subroutine */ int sgees_(const char *jobvs, const char *sort, L_fp select, cblas_int_t *n, 
    float *a, cblas_int_t *lda, cblas_int_t *sdim, float *wr, float *wi, float *vs, 
    cblas_int_t *ldvs, float *work, cblas_int_t *lwork, logical *bwork, cblas_int_t *
    info);

/* Subroutine */ int sgeesx_(const char *jobvs, const char *sort, L_fp select, const char *
    sense, cblas_int_t *n, float *a, cblas_int_t *lda, cblas_int_t *sdim, float *wr, 
    float *wi, float *vs, cblas_int_t *ldvs, float *rconde, float *rcondv, float *
    work, cblas_int_t *lwork, cblas_int_t *iwork, cblas_int_t *liwork, logical *bwork,
     cblas_int_t *info);

/* Subroutine */ int sgeev_(const char *jobvl, const char *jobvr, cblas_int_t *n, float *a, 
    cblas_int_t *lda, float *wr, float *wi, float *vl, cblas_int_t *ldvl, float *vr,    
    cblas_int_t *ldvr, float *work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int sgeevx_(const char *balanc, const char *jobvl, const char *jobvr, const char *
    sense, cblas_int_t *n, float *a, cblas_int_t *lda, float *wr, float *wi, float *
    vl, cblas_int_t *ldvl, float *vr, cblas_int_t *ldvr, cblas_int_t *ilo, cblas_int_t *
    ihi, float *scale, float *abnrm, float *rconde, float *rcondv, float *work,
     cblas_int_t *lwork, cblas_int_t *iwork, cblas_int_t *info);

/* Subroutine */ int sgegs_(const char *jobvsl, const char *jobvsr, cblas_int_t *n, float *a, 
    cblas_int_t *lda, float *b, cblas_int_t *ldb, float *alphar, float *alphai, float 
    *beta, float *vsl, cblas_int_t *ldvsl, float *vsr, cblas_int_t *ldvsr, float *  
    work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int sgegv_(const char *jobvl, const char *jobvr, cblas_int_t *n, float *a, 
    cblas_int_t *lda, float *b, cblas_int_t *ldb, float *alphar, float *alphai, float 
    *beta, float *vl, cblas_int_t *ldvl, float *vr, cblas_int_t *ldvr, float *work, 
    cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int sgehd2_(cblas_int_t *n, cblas_int_t *ilo, cblas_int_t *ihi, float *a, 
    cblas_int_t *lda, float *tau, float *work, cblas_int_t *info);

/* Subroutine */ int sgehrd_(cblas_int_t *n, cblas_int_t *ilo, cblas_int_t *ihi, float *a, 
    cblas_int_t *lda, float *tau, float *work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int sgelq2_(cblas_int_t *m, cblas_int_t *n, float *a, cblas_int_t *lda, 
    float *tau, float *work, cblas_int_t *info);

/* Subroutine */ int sgelqf_(cblas_int_t *m, cblas_int_t *n, float *a, cblas_int_t *lda, 
    float *tau, float *work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int sgels_(const char *trans, cblas_int_t *m, cblas_int_t *n, cblas_int_t *
    nrhs, float *a, cblas_int_t *lda, float *b, cblas_int_t *ldb, float *work, 
    cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int sgelsd_(cblas_int_t *m, cblas_int_t *n, cblas_int_t *nrhs, float *a, 
    cblas_int_t *lda, float *b, cblas_int_t *ldb, float *s, float *rcond, cblas_int_t *
    rank, float *work, cblas_int_t *lwork, cblas_int_t *iwork, cblas_int_t *info);

/* Subroutine */ int sgelss_(cblas_int_t *m, cblas_int_t *n, cblas_int_t *nrhs, float *a, 
    cblas_int_t *lda, float *b, cblas_int_t *ldb, float *s, float *rcond, cblas_int_t *
    rank, float *work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int sgelsx_(cblas_int_t *m, cblas_int_t *n, cblas_int_t *nrhs, float *a, 
    cblas_int_t *lda, float *b, cblas_int_t *ldb, cblas_int_t *jpvt, float *rcond, 
    cblas_int_t *rank, float *work, cblas_int_t *info);

/* Subroutine */ int sgelsy_(cblas_int_t *m, cblas_int_t *n, cblas_int_t *nrhs, float *a, 
    cblas_int_t *lda, float *b, cblas_int_t *ldb, cblas_int_t *jpvt, float *rcond, 
    cblas_int_t *rank, float *work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int sgeql2_(cblas_int_t *m, cblas_int_t *n, float *a, cblas_int_t *lda, 
    float *tau, float *work, cblas_int_t *info);

/* Subroutine */ int sgeqlf_(cblas_int_t *m, cblas_int_t *n, float *a, cblas_int_t *lda, 
    float *tau, float *work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int sgeqp3_(cblas_int_t *m, cblas_int_t *n, float *a, cblas_int_t *lda, 
    cblas_int_t *jpvt, float *tau, float *work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int sgeqpf_(cblas_int_t *m, cblas_int_t *n, float *a, cblas_int_t *lda, 
    cblas_int_t *jpvt, float *tau, float *work, cblas_int_t *info);

/* Subroutine */ int sgeqr2_(cblas_int_t *m, cblas_int_t *n, float *a, cblas_int_t *lda, 
    float *tau, float *work, cblas_int_t *info);

/* Subroutine */ int sgeqrf_(cblas_int_t *m, cblas_int_t *n, float *a, cblas_int_t *lda, 
    float *tau, float *work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int sgerfs_(const char *trans, cblas_int_t *n, cblas_int_t *nrhs, float *a, 
    cblas_int_t *lda, float *af, cblas_int_t *ldaf, cblas_int_t *ipiv, float *b, 
    cblas_int_t *ldb, float *x, cblas_int_t *ldx, float *ferr, float *berr, float *
    work, cblas_int_t *iwork, cblas_int_t *info);

/* Subroutine */ int sgerq2_(cblas_int_t *m, cblas_int_t *n, float *a, cblas_int_t *lda, 
    float *tau, float *work, cblas_int_t *info);

/* Subroutine */ int sgerqf_(cblas_int_t *m, cblas_int_t *n, float *a, cblas_int_t *lda, 
    float *tau, float *work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int sgesc2_(cblas_int_t *n, float *a, cblas_int_t *lda, float *rhs, 
    cblas_int_t *ipiv, cblas_int_t *jpiv, float *scale);

/* Subroutine */ int sgesdd_(const char *jobz, cblas_int_t *m, cblas_int_t *n, float *a, 
    cblas_int_t *lda, float *s, float *u, cblas_int_t *ldu, float *vt, cblas_int_t *ldvt,    
    float *work, cblas_int_t *lwork, cblas_int_t *iwork, cblas_int_t *info);

/* Subroutine */ int sgesv_(cblas_int_t *n, cblas_int_t *nrhs, float *a, cblas_int_t *lda, 
    cblas_int_t *ipiv, float *b, cblas_int_t *ldb, cblas_int_t *info);

/* Subroutine */ int sgesvd_(const char *jobu, const char *jobvt, cblas_int_t *m, cblas_int_t *n, 
    float *a, cblas_int_t *lda, float *s, float *u, cblas_int_t *ldu, float *vt,    
    cblas_int_t *ldvt, float *work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int sgesvx_(const char *fact, const char *trans, cblas_int_t *n, cblas_int_t *
    nrhs, float *a, cblas_int_t *lda, float *af, cblas_int_t *ldaf, cblas_int_t *ipiv, 
    const char *equed, float *r__, float *c__, float *b, cblas_int_t *ldb, float *x, 
    cblas_int_t *ldx, float *rcond, float *ferr, float *berr, float *work, 
    cblas_int_t *iwork, cblas_int_t *info);

/* Subroutine */ int sgetc2_(cblas_int_t *n, float *a, cblas_int_t *lda, cblas_int_t *ipiv,
     cblas_int_t *jpiv, cblas_int_t *info);

/* Subroutine */ int sgetf2_(cblas_int_t *m, cblas_int_t *n, float *a, cblas_int_t *lda, 
    cblas_int_t *ipiv, cblas_int_t *info);

/* Subroutine */ int sgetrf_(cblas_int_t *m, cblas_int_t *n, float *a, cblas_int_t *lda, 
    cblas_int_t *ipiv, cblas_int_t *info);

/* Subroutine */ int sgetri_(cblas_int_t *n, float *a, cblas_int_t *lda, cblas_int_t *ipiv,
     float *work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int sgetrs_(const char *trans, cblas_int_t *n, cblas_int_t *nrhs, float *a, 
    cblas_int_t *lda, cblas_int_t *ipiv, float *b, cblas_int_t *ldb, cblas_int_t *info);

/* Subroutine */ int sggbak_(const char *job, const char *side, cblas_int_t *n, cblas_int_t *ilo, 
    cblas_int_t *ihi, float *lscale, float *rscale, cblas_int_t *m, float *v, 
    cblas_int_t *ldv, cblas_int_t *info);

/* Subroutine */ int sggbal_(const char *job, cblas_int_t *n, float *a, cblas_int_t *lda, 
    float *b, cblas_int_t *ldb, cblas_int_t *ilo, cblas_int_t *ihi, float *lscale, float 
    *rscale, float *work, cblas_int_t *info);

/* Subroutine */ int sgges_(const char *jobvsl, const char *jobvsr, const char *sort, L_fp 
    selctg, cblas_int_t *n, float *a, cblas_int_t *lda, float *b, cblas_int_t *ldb, 
    cblas_int_t *sdim, float *alphar, float *alphai, float *beta, float *vsl, 
    cblas_int_t *ldvsl, float *vsr, cblas_int_t *ldvsr, float *work, cblas_int_t *lwork,
     logical *bwork, cblas_int_t *info);

/* Subroutine */ int sggesx_(const char *jobvsl, const char *jobvsr, const char *sort, L_fp 
    selctg, const char *sense, cblas_int_t *n, float *a, cblas_int_t *lda, float *b, 
    cblas_int_t *ldb, cblas_int_t *sdim, float *alphar, float *alphai, float *beta, 
    float *vsl, cblas_int_t *ldvsl, float *vsr, cblas_int_t *ldvsr, float *rconde, 
    float *rcondv, float *work, cblas_int_t *lwork, cblas_int_t *iwork, cblas_int_t *   
    liwork, logical *bwork, cblas_int_t *info);

/* Subroutine */ int sggev_(const char *jobvl, const char *jobvr, cblas_int_t *n, float *a, 
    cblas_int_t *lda, float *b, cblas_int_t *ldb, float *alphar, float *alphai, float 
    *beta, float *vl, cblas_int_t *ldvl, float *vr, cblas_int_t *ldvr, float *work, 
    cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int sggevx_(const char *balanc, const char *jobvl, const char *jobvr, const char *
    sense, cblas_int_t *n, float *a, cblas_int_t *lda, float *b, cblas_int_t *ldb, float 
    *alphar, float *alphai, float *beta, float *vl, cblas_int_t *ldvl, float *vr, 
    cblas_int_t *ldvr, cblas_int_t *ilo, cblas_int_t *ihi, float *lscale, float *rscale,
     float *abnrm, float *bbnrm, float *rconde, float *rcondv, float *work,     
    cblas_int_t *lwork, cblas_int_t *iwork, logical *bwork, cblas_int_t *info);

/* Subroutine */ int sggglm_(cblas_int_t *n, cblas_int_t *m, cblas_int_t *p, float *a, 
    cblas_int_t *lda, float *b, cblas_int_t *ldb, float *d__, float *x, float *y, 
    float *work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int sgghrd_(const char *compq, const char *compz, cblas_int_t *n, cblas_int_t *
    ilo, cblas_int_t *ihi, float *a, cblas_int_t *lda, float *b, cblas_int_t *ldb, float 
    *q, cblas_int_t *ldq, float *z__, cblas_int_t *ldz, cblas_int_t *info);

/* Subroutine */ int sgglse_(cblas_int_t *m, cblas_int_t *n, cblas_int_t *p, float *a, 
    cblas_int_t *lda, float *b, cblas_int_t *ldb, float *c__, float *d__, float *x, 
    float *work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int sggqrf_(cblas_int_t *n, cblas_int_t *m, cblas_int_t *p, float *a, 
    cblas_int_t *lda, float *taua, float *b, cblas_int_t *ldb, float *taub, float *
    work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int sggrqf_(cblas_int_t *m, cblas_int_t *p, cblas_int_t *n, float *a, 
    cblas_int_t *lda, float *taua, float *b, cblas_int_t *ldb, float *taub, float *
    work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int sggsvd_(const char *jobu, const char *jobv, const char *jobq, cblas_int_t *m, 
    cblas_int_t *n, cblas_int_t *p, cblas_int_t *k, cblas_int_t *l, float *a, cblas_int_t *lda,
     float *b, cblas_int_t *ldb, float *alpha, float *beta, float *u, cblas_int_t *
    ldu, float *v, cblas_int_t *ldv, float *q, cblas_int_t *ldq, float *work, 
    cblas_int_t *iwork, cblas_int_t *info);

/* Subroutine */ int sggsvp_(const char *jobu, const char *jobv, const char *jobq, cblas_int_t *m, 
    cblas_int_t *p, cblas_int_t *n, float *a, cblas_int_t *lda, float *b, cblas_int_t *ldb, 
    float *tola, float *tolb, cblas_int_t *k, cblas_int_t *l, float *u, cblas_int_t *ldu,
     float *v, cblas_int_t *ldv, float *q, cblas_int_t *ldq, cblas_int_t *iwork, float *
    tau, float *work, cblas_int_t *info);

/* Subroutine */ int sgtcon_(const char *norm, cblas_int_t *n, float *dl, float *d__, 
    float *du, float *du2, cblas_int_t *ipiv, float *anorm, float *rcond, float *
    work, cblas_int_t *iwork, cblas_int_t *info);

/* Subroutine */ int sgtrfs_(const char *trans, cblas_int_t *n, cblas_int_t *nrhs, float *dl,
     float *d__, float *du, float *dlf, float *df, float *duf, float *du2, 
    cblas_int_t *ipiv, float *b, cblas_int_t *ldb, float *x, cblas_int_t *ldx, float *
    ferr, float *berr, float *work, cblas_int_t *iwork, cblas_int_t *info);

/* Subroutine */ int sgtsv_(cblas_int_t *n, cblas_int_t *nrhs, float *dl, float *d__, 
    float *du, float *b, cblas_int_t *ldb, cblas_int_t *info);

/* Subroutine */ int sgtsvx_(const char *fact, const char *trans, cblas_int_t *n, cblas_int_t *
    nrhs, float *dl, float *d__, float *du, float *dlf, float *df, float *duf, 
    float *du2, cblas_int_t *ipiv, float *b, cblas_int_t *ldb, float *x, cblas_int_t *
    ldx, float *rcond, float *ferr, float *berr, float *work, cblas_int_t *iwork, 
    cblas_int_t *info);

/* Subroutine */ int sgttrf_(cblas_int_t *n, float *dl, float *d__, float *du, float *
    du2, cblas_int_t *ipiv, cblas_int_t *info);

/* Subroutine */ int sgttrs_(const char *trans, cblas_int_t *n, cblas_int_t *nrhs, float *dl,
     float *d__, float *du, float *du2, cblas_int_t *ipiv, float *b, cblas_int_t *ldb,
     cblas_int_t *info);

/* Subroutine */ int sgtts2_(cblas_int_t *itrans, cblas_int_t *n, cblas_int_t *nrhs, float 
    *dl, float *d__, float *du, float *du2, cblas_int_t *ipiv, float *b, cblas_int_t *
    ldb);

/* Subroutine */ int shgeqz_(const char *job, const char *compq, const char *compz, cblas_int_t *n, 
    cblas_int_t *ilo, cblas_int_t *ihi, float *h__, cblas_int_t *ldh, float *t, cblas_int_t 
    *ldt, float *alphar, float *alphai, float *beta, float *q, cblas_int_t *ldq, 
    float *z__, cblas_int_t *ldz, float *work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int shsein_(const char *side, const char *eigsrc, const char *initv, logical *
    select, cblas_int_t *n, float *h__, cblas_int_t *ldh, float *wr, float *wi, float 
    *vl, cblas_int_t *ldvl, float *vr, cblas_int_t *ldvr, cblas_int_t *mm, cblas_int_t *m, 
    float *work, cblas_int_t *ifaill, cblas_int_t *ifailr, cblas_int_t *info);

/* Subroutine */ int shseqr_(const char *job, const char *compz, cblas_int_t *n, cblas_int_t *ilo,
     cblas_int_t *ihi, float *h__, cblas_int_t *ldh, float *wr, float *wi, float *z__,   
    cblas_int_t *ldz, float *work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int slabad_(float *small, float *large);

/* Subroutine */ int slabrd_(cblas_int_t *m, cblas_int_t *n, cblas_int_t *nb, float *a, 
    cblas_int_t *lda, float *d__, float *e, float *tauq, float *taup, float *x, 
    cblas_int_t *ldx, float *y, cblas_int_t *ldy);

/* Subroutine */ int slacn2_(cblas_int_t *n, float *v, float *x, cblas_int_t *isgn, 
    float *est, cblas_int_t *kase, cblas_int_t *isave);

/* Subroutine */ int slacon_(cblas_int_t *n, float *v, float *x, cblas_int_t *isgn, 
    float *est, cblas_int_t *kase);

/* Subroutine */ int slacpy_(const char *uplo, cblas_int_t *m, cblas_int_t *n, float *a, 
    cblas_int_t *lda, float *b, cblas_int_t *ldb);

/* Subroutine */ int sladiv_(float *a, float *b, float *c__, float *d__, float *p, 
    float *q);

/* Subroutine */ int slae2_(float *a, float *b, float *c__, float *rt1, float *rt2);

/* Subroutine */ int slaebz_(cblas_int_t *ijob, cblas_int_t *nitmax, cblas_int_t *n, 
    cblas_int_t *mmax, cblas_int_t *minp, cblas_int_t *nbmin, float *abstol, float *
    reltol, float *pivmin, float *d__, float *e, float *e2, cblas_int_t *nval, 
    float *ab, float *c__, cblas_int_t *mout, cblas_int_t *nab, float *work, cblas_int_t 
    *iwork, cblas_int_t *info);

/* Subroutine */ int slaed0_(cblas_int_t *icompq, cblas_int_t *qsiz, cblas_int_t *n, float 
    *d__, float *e, float *q, cblas_int_t *ldq, float *qstore, cblas_int_t *ldqs, 
    float *work, cblas_int_t *iwork, cblas_int_t *info);

/* Subroutine */ int slaed1_(cblas_int_t *n, float *d__, float *q, cblas_int_t *ldq, 
    cblas_int_t *indxq, float *rho, cblas_int_t *cutpnt, float *work, cblas_int_t *
    iwork, cblas_int_t *info);

/* Subroutine */ int slaed2_(cblas_int_t *k, cblas_int_t *n, cblas_int_t *n1, float *d__, 
    float *q, cblas_int_t *ldq, cblas_int_t *indxq, float *rho, float *z__, float *
    dlamda, float *w, float *q2, cblas_int_t *indx, cblas_int_t *indxc, cblas_int_t *
    indxp, cblas_int_t *coltyp, cblas_int_t *info);

/* Subroutine */ int slaed3_(cblas_int_t *k, cblas_int_t *n, cblas_int_t *n1, float *d__, 
    float *q, cblas_int_t *ldq, float *rho, float *dlamda, float *q2, cblas_int_t *
    indx, cblas_int_t *ctot, float *w, float *s, cblas_int_t *info);

/* Subroutine */ int slaed4_(cblas_int_t *n, cblas_int_t *i__, float *d__, float *z__, 
    float *delta, float *rho, float *dlam, cblas_int_t *info);

/* Subroutine */ int slaed5_(cblas_int_t *i__, float *d__, float *z__, float *delta, 
    float *rho, float *dlam);

/* Subroutine */ int slaed6_(cblas_int_t *kniter, logical *orgati, float *rho, 
    float *d__, float *z__, float *finit, float *tau, cblas_int_t *info);

/* Subroutine */ int slaed7_(cblas_int_t *icompq, cblas_int_t *n, cblas_int_t *qsiz, 
    cblas_int_t *tlvls, cblas_int_t *curlvl, cblas_int_t *curpbm, float *d__, float *q, 
    cblas_int_t *ldq, cblas_int_t *indxq, float *rho, cblas_int_t *cutpnt, float *
    qstore, cblas_int_t *qptr, cblas_int_t *prmptr, cblas_int_t *perm, cblas_int_t *
    givptr, cblas_int_t *givcol, float *givnum, float *work, cblas_int_t *iwork, 
    cblas_int_t *info);

/* Subroutine */ int slaed8_(cblas_int_t *icompq, cblas_int_t *k, cblas_int_t *n, cblas_int_t 
    *qsiz, float *d__, float *q, cblas_int_t *ldq, cblas_int_t *indxq, float *rho, 
    cblas_int_t *cutpnt, float *z__, float *dlamda, float *q2, cblas_int_t *ldq2, 
    float *w, cblas_int_t *perm, cblas_int_t *givptr, cblas_int_t *givcol, float *
    givnum, cblas_int_t *indxp, cblas_int_t *indx, cblas_int_t *info);

/* Subroutine */ int slaed9_(cblas_int_t *k, cblas_int_t *kstart, cblas_int_t *kstop, 
    cblas_int_t *n, float *d__, float *q, cblas_int_t *ldq, float *rho, float *dlamda,
     float *w, float *s, cblas_int_t *lds, cblas_int_t *info);

/* Subroutine */ int slaeda_(cblas_int_t *n, cblas_int_t *tlvls, cblas_int_t *curlvl, 
    cblas_int_t *curpbm, cblas_int_t *prmptr, cblas_int_t *perm, cblas_int_t *givptr, 
    cblas_int_t *givcol, float *givnum, float *q, cblas_int_t *qptr, float *z__, 
    float *ztemp, cblas_int_t *info);

/* Subroutine */ int slaein_(logical *rightv, logical *noinit, cblas_int_t *n, 
    float *h__, cblas_int_t *ldh, float *wr, float *wi, float *vr, float *vi, float 
    *b, cblas_int_t *ldb, float *work, float *eps3, float *smlnum, float *bignum, 
    cblas_int_t *info);

/* Subroutine */ int slaev2_(float *a, float *b, float *c__, float *rt1, float *
    rt2, float *cs1, float *sn1);

/* Subroutine */ int slaexc_(logical *wantq, cblas_int_t *n, float *t, cblas_int_t *
    ldt, float *q, cblas_int_t *ldq, cblas_int_t *j1, cblas_int_t *n1, cblas_int_t *n2, 
    float *work, cblas_int_t *info);

/* Subroutine */ int slag2_(float *a, cblas_int_t *lda, float *b, cblas_int_t *ldb, 
    float *safmin, float *scale1, float *scale2, float *wr1, float *wr2, float *
    wi);

/* Subroutine */ int slag2d_(cblas_int_t *m, cblas_int_t *n, float *sa, cblas_int_t *ldsa, 
    double *a, cblas_int_t *lda, cblas_int_t *info);

/* Subroutine */ int slags2_(logical *upper, float *a1, float *a2, float *a3, 
    float *b1, float *b2, float *b3, float *csu, float *snu, float *csv, float *
    snv, float *csq, float *snq);

/* Subroutine */ int slagtf_(cblas_int_t *n, float *a, float *lambda, float *b, float 
    *c__, float *tol, float *d__, cblas_int_t *in, cblas_int_t *info);

/* Subroutine */ int slagtm_(const char *trans, cblas_int_t *n, cblas_int_t *nrhs, float *
    alpha, float *dl, float *d__, float *du, float *x, cblas_int_t *ldx, float *
    beta, float *b, cblas_int_t *ldb);

/* Subroutine */ int slagts_(cblas_int_t *job, cblas_int_t *n, float *a, float *b, float 
    *c__, float *d__, cblas_int_t *in, float *y, float *tol, cblas_int_t *info);

/* Subroutine */ int slagv2_(float *a, cblas_int_t *lda, float *b, cblas_int_t *ldb, 
    float *alphar, float *alphai, float *beta, float *csl, float *snl, float *
    csr, float *snr);

/* Subroutine */ int slahqr_(logical *wantt, logical *wantz, cblas_int_t *n, 
    cblas_int_t *ilo, cblas_int_t *ihi, float *h__, cblas_int_t *ldh, float *wr, float *
    wi, cblas_int_t *iloz, cblas_int_t *ihiz, float *z__, cblas_int_t *ldz, cblas_int_t *
    info);

/* Subroutine */ int slahr2_(cblas_int_t *n, cblas_int_t *k, cblas_int_t *nb, float *a, 
    cblas_int_t *lda, float *tau, float *t, cblas_int_t *ldt, float *y, cblas_int_t *ldy);

/* Subroutine */ int slahrd_(cblas_int_t *n, cblas_int_t *k, cblas_int_t *nb, float *a, 
    cblas_int_t *lda, float *tau, float *t, cblas_int_t *ldt, float *y, cblas_int_t *ldy);

/* Subroutine */ int slaic1_(cblas_int_t *job, cblas_int_t *j, float *x, float *sest, 
    float *w, float *gamma, float *sestpr, float *s, float *c__);

/* Subroutine */ int slaln2_(logical *ltrans, cblas_int_t *na, cblas_int_t *nw, float *
    smin, float *ca, float *a, cblas_int_t *lda, float *d1, float *d2, float *b, 
    cblas_int_t *ldb, float *wr, float *wi, float *x, cblas_int_t *ldx, float *scale, 
    float *xnorm, cblas_int_t *info);

/* Subroutine */ int slals0_(cblas_int_t *icompq, cblas_int_t *nl, cblas_int_t *nr, 
    cblas_int_t *sqre, cblas_int_t *nrhs, float *b, cblas_int_t *ldb, float *bx, 
    cblas_int_t *ldbx, cblas_int_t *perm, cblas_int_t *givptr, cblas_int_t *givcol, 
    cblas_int_t *ldgcol, float *givnum, cblas_int_t *ldgnum, float *poles, float *
    difl, float *difr, float *z__, cblas_int_t *k, float *c__, float *s, float *
    work, cblas_int_t *info);

/* Subroutine */ int slalsa_(cblas_int_t *icompq, cblas_int_t *smlsiz, cblas_int_t *n, 
    cblas_int_t *nrhs, float *b, cblas_int_t *ldb, float *bx, cblas_int_t *ldbx, float *
    u, cblas_int_t *ldu, float *vt, cblas_int_t *k, float *difl, float *difr, float *
    z__, float *poles, cblas_int_t *givptr, cblas_int_t *givcol, cblas_int_t *ldgcol, 
    cblas_int_t *perm, float *givnum, float *c__, float *s, float *work, cblas_int_t *
    iwork, cblas_int_t *info);

/* Subroutine */ int slalsd_(const char *uplo, cblas_int_t *smlsiz, cblas_int_t *n, cblas_int_t 
    *nrhs, float *d__, float *e, float *b, cblas_int_t *ldb, float *rcond, 
    cblas_int_t *rank, float *work, cblas_int_t *iwork, cblas_int_t *info);

/* Subroutine */ int slamrg_(cblas_int_t *n1, cblas_int_t *n2, float *a, cblas_int_t *
    strd1, cblas_int_t *strd2, cblas_int_t *index);

/* Subroutine */ int slanv2_(float *a, float *b, float *c__, float *d__, float *
    rt1r, float *rt1i, float *rt2r, float *rt2i, float *cs, float *sn);

/* Subroutine */ int slapll_(cblas_int_t *n, float *x, cblas_int_t *incx, float *y, 
    cblas_int_t *incy, float *ssmin);

/* Subroutine */ int slapmt_(logical *forwrd, cblas_int_t *m, cblas_int_t *n, float *x,
     cblas_int_t *ldx, cblas_int_t *k);

/* Subroutine */ int slaqgb_(cblas_int_t *m, cblas_int_t *n, cblas_int_t *kl, cblas_int_t *ku,
     float *ab, cblas_int_t *ldab, float *r__, float *c__, float *rowcnd, float *
    colcnd, float *amax, const char *equed);

/* Subroutine */ int slaqge_(cblas_int_t *m, cblas_int_t *n, float *a, cblas_int_t *lda, 
    float *r__, float *c__, float *rowcnd, float *colcnd, float *amax, const char *
    equed);

/* Subroutine */ int slaqp2_(cblas_int_t *m, cblas_int_t *n, cblas_int_t *offset, float *a,
     cblas_int_t *lda, cblas_int_t *jpvt, float *tau, float *vn1, float *vn2, float *
    work);

/* Subroutine */ int slaqps_(cblas_int_t *m, cblas_int_t *n, cblas_int_t *offset, cblas_int_t 
    *nb, cblas_int_t *kb, float *a, cblas_int_t *lda, cblas_int_t *jpvt, float *tau, 
    float *vn1, float *vn2, float *auxv, float *f, cblas_int_t *ldf);

/* Subroutine */ int slaqr0_(logical *wantt, logical *wantz, cblas_int_t *n, 
    cblas_int_t *ilo, cblas_int_t *ihi, float *h__, cblas_int_t *ldh, float *wr, float *
    wi, cblas_int_t *iloz, cblas_int_t *ihiz, float *z__, cblas_int_t *ldz, float *work,
     cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int slaqr1_(cblas_int_t *n, float *h__, cblas_int_t *ldh, float *sr1, 
    float *si1, float *sr2, float *si2, float *v);

/* Subroutine */ int slaqr2_(logical *wantt, logical *wantz, cblas_int_t *n, 
    cblas_int_t *ktop, cblas_int_t *kbot, cblas_int_t *nw, float *h__, cblas_int_t *ldh, 
    cblas_int_t *iloz, cblas_int_t *ihiz, float *z__, cblas_int_t *ldz, cblas_int_t *ns, 
    cblas_int_t *nd, float *sr, float *si, float *v, cblas_int_t *ldv, cblas_int_t *nh, 
    float *t, cblas_int_t *ldt, cblas_int_t *nv, float *wv, cblas_int_t *ldwv, float *
    work, cblas_int_t *lwork);

/* Subroutine */ int slaqr3_(logical *wantt, logical *wantz, cblas_int_t *n, 
    cblas_int_t *ktop, cblas_int_t *kbot, cblas_int_t *nw, float *h__, cblas_int_t *ldh, 
    cblas_int_t *iloz, cblas_int_t *ihiz, float *z__, cblas_int_t *ldz, cblas_int_t *ns, 
    cblas_int_t *nd, float *sr, float *si, float *v, cblas_int_t *ldv, cblas_int_t *nh, 
    float *t, cblas_int_t *ldt, cblas_int_t *nv, float *wv, cblas_int_t *ldwv, float *
    work, cblas_int_t *lwork);

/* Subroutine */ int slaqr4_(logical *wantt, logical *wantz, cblas_int_t *n, 
    cblas_int_t *ilo, cblas_int_t *ihi, float *h__, cblas_int_t *ldh, float *wr, float *
    wi, cblas_int_t *iloz, cblas_int_t *ihiz, float *z__, cblas_int_t *ldz, float *work,
     cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int slaqr5_(logical *wantt, logical *wantz, cblas_int_t *kacc22, 
    cblas_int_t *n, cblas_int_t *ktop, cblas_int_t *kbot, cblas_int_t *nshfts, float *sr, 
    float *si, float *h__, cblas_int_t *ldh, cblas_int_t *iloz, cblas_int_t *ihiz, float 
    *z__, cblas_int_t *ldz, float *v, cblas_int_t *ldv, float *u, cblas_int_t *ldu, 
    cblas_int_t *nv, float *wv, cblas_int_t *ldwv, cblas_int_t *nh, float *wh, cblas_int_t *
    ldwh);

/* Subroutine */ int slaqsb_(const char *uplo, cblas_int_t *n, cblas_int_t *kd, float *ab, 
    cblas_int_t *ldab, float *s, float *scond, float *amax, const char *equed   );

/* Subroutine */ int slaqsp_(const char *uplo, cblas_int_t *n, float *ap, float *s, float *
    scond, float *amax, const char *equed);

/* Subroutine */ int slaqsy_(const char *uplo, cblas_int_t *n, float *a, cblas_int_t *lda, 
    float *s, float *scond, float *amax, const char *equed);

/* Subroutine */ int slaqtr_(logical *ltran, logical *lfloat, cblas_int_t *n, float 
    *t, cblas_int_t *ldt, float *b, float *w, float *scale, float *x, float *work, 
    cblas_int_t *info);

/* Subroutine */ int slar1v_(cblas_int_t *n, cblas_int_t *b1, cblas_int_t *bn, float *
    lambda, float *d__, float *l, float *ld, float *lld, float *pivmin, float *
    gaptol, float *z__, logical *wantnc, cblas_int_t *negcnt, float *ztz, float *
    mingma, cblas_int_t *r__, cblas_int_t *isuppz, float *nrminv, float *resid, 
    float *rqcorr, float *work);

/* Subroutine */ int slar2v_(cblas_int_t *n, float *x, float *y, float *z__, cblas_int_t 
    *incx, float *c__, float *s, cblas_int_t *incc);

/* Subroutine */ int slarf_(const char *side, cblas_int_t *m, cblas_int_t *n, float *v, 
    cblas_int_t *incv, float *tau, float *c__, cblas_int_t *ldc, float *work    );

/* Subroutine */ int slarfb_(const char *side, const char *trans, const char *direct, const char *
    storev, cblas_int_t *m, cblas_int_t *n, cblas_int_t *k, float *v, cblas_int_t *ldv, 
    float *t, cblas_int_t *ldt, float *c__, cblas_int_t *ldc, float *work, cblas_int_t *
    ldwork      );

/* Subroutine */ int slarfg_(cblas_int_t *n, float *alpha, float *x, cblas_int_t *incx, 
    float *tau);

/* Subroutine */ int slarft_(const char *direct, const char *storev, cblas_int_t *n, cblas_int_t *
    k, float *v, cblas_int_t *ldv, float *tau, float *t, cblas_int_t *ldt   );

/* Subroutine */ int slarfx_(const char *side, cblas_int_t *m, cblas_int_t *n, float *v, 
    float *tau, float *c__, cblas_int_t *ldc, float *work);

/* Subroutine */ int slargv_(cblas_int_t *n, float *x, cblas_int_t *incx, float *y, 
    cblas_int_t *incy, float *c__, cblas_int_t *incc);

/* Subroutine */ int slarnv_(cblas_int_t *idist, cblas_int_t *iseed, cblas_int_t *n, float 
    *x);

/* Subroutine */ int slarra_(cblas_int_t *n, float *d__, float *e, float *e2, float *
    spltol, float *tnrm, cblas_int_t *nsplit, cblas_int_t *isplit, cblas_int_t *info);

/* Subroutine */ int slarrb_(cblas_int_t *n, float *d__, float *lld, cblas_int_t *
    ifirst, cblas_int_t *ilast, float *rtol1, float *rtol2, cblas_int_t *offset, 
    float *w, float *wgap, float *werr, float *work, cblas_int_t *iwork, float *
    pivmin, float *spdiam, cblas_int_t *twist, cblas_int_t *info);

/* Subroutine */ int slarrc_(const char *jobt, cblas_int_t *n, float *vl, float *vu, float 
    *d__, float *e, float *pivmin, cblas_int_t *eigcnt, cblas_int_t *lcnt, cblas_int_t *
    rcnt, cblas_int_t *info);

/* Subroutine */ int slarrd_(const char *range, const char *order, cblas_int_t *n, float *vl, 
    float *vu, cblas_int_t *il, cblas_int_t *iu, float *gers, float *reltol, float *
    d__, float *e, float *e2, float *pivmin, cblas_int_t *nsplit, cblas_int_t *
    isplit, cblas_int_t *m, float *w, float *werr, float *wl, float *wu, cblas_int_t *
    iblock, cblas_int_t *indexw, float *work, cblas_int_t *iwork, cblas_int_t *info);

/* Subroutine */ int slarre_(const char *range, cblas_int_t *n, float *vl, float *vu, 
    cblas_int_t *il, cblas_int_t *iu, float *d__, float *e, float *e2, float *rtol1, 
    float *rtol2, float *spltol, cblas_int_t *nsplit, cblas_int_t *isplit, cblas_int_t *
    m, float *w, float *werr, float *wgap, cblas_int_t *iblock, cblas_int_t *indexw, 
    float *gers, float *pivmin, float *work, cblas_int_t *iwork, cblas_int_t *info);

/* Subroutine */ int slarrf_(cblas_int_t *n, float *d__, float *l, float *ld, 
    cblas_int_t *clstrt, cblas_int_t *clend, float *w, float *wgap, float *werr, 
    float *spdiam, float *clgapl, float *clgapr, float *pivmin, float *sigma, 
    float *dplus, float *lplus, float *work, cblas_int_t *info);

/* Subroutine */ int slarrj_(cblas_int_t *n, float *d__, float *e2, cblas_int_t *ifirst,
     cblas_int_t *ilast, float *rtol, cblas_int_t *offset, float *w, float *werr, 
    float *work, cblas_int_t *iwork, float *pivmin, float *spdiam, cblas_int_t *info);

/* Subroutine */ int slarrk_(cblas_int_t *n, cblas_int_t *iw, float *gl, float *gu, 
    float *d__, float *e2, float *pivmin, float *reltol, float *w, float *werr, 
    cblas_int_t *info);

/* Subroutine */ int slarrr_(cblas_int_t *n, float *d__, float *e, cblas_int_t *info);

/* Subroutine */ int slarrv_(cblas_int_t *n, float *vl, float *vu, float *d__, float *
    l, float *pivmin, cblas_int_t *isplit, cblas_int_t *m, cblas_int_t *dol, cblas_int_t *
    dou, float *minrgp, float *rtol1, float *rtol2, float *w, float *werr, 
    float *wgap, cblas_int_t *iblock, cblas_int_t *indexw, float *gers, float *z__, 
    cblas_int_t *ldz, cblas_int_t *isuppz, float *work, cblas_int_t *iwork, cblas_int_t *
    info);

/* Subroutine */ int slartg_(float *f, float *g, float *cs, float *sn, float *r__);

/* Subroutine */ int slartv_(cblas_int_t *n, float *x, cblas_int_t *incx, float *y, 
    cblas_int_t *incy, float *c__, float *s, cblas_int_t *incc);

/* Subroutine */ int slaruv_(cblas_int_t *iseed, cblas_int_t *n, float *x);

/* Subroutine */ int slarz_(const char *side, cblas_int_t *m, cblas_int_t *n, cblas_int_t *l, 
    float *v, cblas_int_t *incv, float *tau, float *c__, cblas_int_t *ldc, float *
    work);

/* Subroutine */ int slarzb_(const char *side, const char *trans, const char *direct, const char *
    storev, cblas_int_t *m, cblas_int_t *n, cblas_int_t *k, cblas_int_t *l, float *v, 
    cblas_int_t *ldv, float *t, cblas_int_t *ldt, float *c__, cblas_int_t *ldc, float *
    work, cblas_int_t *ldwork   );

/* Subroutine */ int slarzt_(const char *direct, const char *storev, cblas_int_t *n, cblas_int_t *
    k, float *v, cblas_int_t *ldv, float *tau, float *t, cblas_int_t *ldt   );

/* Subroutine */ int slas2_(float *f, float *g, float *h__, float *ssmin, float *
    ssmax);

/* Subroutine */ int slascl_(const char *type__, cblas_int_t *kl, cblas_int_t *ku, float *
    cfrom, float *cto, cblas_int_t *m, cblas_int_t *n, float *a, cblas_int_t *lda, 
    cblas_int_t *info);

/* Subroutine */ int slasd0_(cblas_int_t *n, cblas_int_t *sqre, float *d__, float *e, 
    float *u, cblas_int_t *ldu, float *vt, cblas_int_t *ldvt, cblas_int_t *smlsiz, 
    cblas_int_t *iwork, float *work, cblas_int_t *info);

/* Subroutine */ int slasd1_(cblas_int_t *nl, cblas_int_t *nr, cblas_int_t *sqre, float *
    d__, float *alpha, float *beta, float *u, cblas_int_t *ldu, float *vt, 
    cblas_int_t *ldvt, cblas_int_t *idxq, cblas_int_t *iwork, float *work, cblas_int_t *
    info);

/* Subroutine */ int slasd2_(cblas_int_t *nl, cblas_int_t *nr, cblas_int_t *sqre, cblas_int_t 
    *k, float *d__, float *z__, float *alpha, float *beta, float *u, cblas_int_t *
    ldu, float *vt, cblas_int_t *ldvt, float *dsigma, float *u2, cblas_int_t *ldu2, 
    float *vt2, cblas_int_t *ldvt2, cblas_int_t *idxp, cblas_int_t *idx, cblas_int_t *idxc,
     cblas_int_t *idxq, cblas_int_t *coltyp, cblas_int_t *info);

/* Subroutine */ int slasd3_(cblas_int_t *nl, cblas_int_t *nr, cblas_int_t *sqre, cblas_int_t 
    *k, float *d__, float *q, cblas_int_t *ldq, float *dsigma, float *u, cblas_int_t *
    ldu, float *u2, cblas_int_t *ldu2, float *vt, cblas_int_t *ldvt, float *vt2, 
    cblas_int_t *ldvt2, cblas_int_t *idxc, cblas_int_t *ctot, float *z__, cblas_int_t *
    info);

/* Subroutine */ int slasd4_(cblas_int_t *n, cblas_int_t *i__, float *d__, float *z__, 
    float *delta, float *rho, float *sigma, float *work, cblas_int_t *info);

/* Subroutine */ int slasd5_(cblas_int_t *i__, float *d__, float *z__, float *delta, 
    float *rho, float *dsigma, float *work);

/* Subroutine */ int slasd6_(cblas_int_t *icompq, cblas_int_t *nl, cblas_int_t *nr, 
    cblas_int_t *sqre, float *d__, float *vf, float *vl, float *alpha, float *beta,
     cblas_int_t *idxq, cblas_int_t *perm, cblas_int_t *givptr, cblas_int_t *givcol, 
    cblas_int_t *ldgcol, float *givnum, cblas_int_t *ldgnum, float *poles, float *
    difl, float *difr, float *z__, cblas_int_t *k, float *c__, float *s, float *
    work, cblas_int_t *iwork, cblas_int_t *info);

/* Subroutine */ int slasd7_(cblas_int_t *icompq, cblas_int_t *nl, cblas_int_t *nr, 
    cblas_int_t *sqre, cblas_int_t *k, float *d__, float *z__, float *zw, float *vf, 
    float *vfw, float *vl, float *vlw, float *alpha, float *beta, float *dsigma,
     cblas_int_t *idx, cblas_int_t *idxp, cblas_int_t *idxq, cblas_int_t *perm, cblas_int_t *
    givptr, cblas_int_t *givcol, cblas_int_t *ldgcol, float *givnum, cblas_int_t *
    ldgnum, float *c__, float *s, cblas_int_t *info);

/* Subroutine */ int slasd8_(cblas_int_t *icompq, cblas_int_t *k, float *d__, float *
    z__, float *vf, float *vl, float *difl, float *difr, cblas_int_t *lddifr, 
    float *dsigma, float *work, cblas_int_t *info);

/* Subroutine */ int slasd9_(cblas_int_t *icompq, cblas_int_t *ldu, cblas_int_t *k, float *
    d__, float *z__, float *vf, float *vl, float *difl, float *difr, float *
    dsigma, float *work, cblas_int_t *info);

/* Subroutine */ int slasda_(cblas_int_t *icompq, cblas_int_t *smlsiz, cblas_int_t *n, 
    cblas_int_t *sqre, float *d__, float *e, float *u, cblas_int_t *ldu, float *vt, 
    cblas_int_t *k, float *difl, float *difr, float *z__, float *poles, cblas_int_t *
    givptr, cblas_int_t *givcol, cblas_int_t *ldgcol, cblas_int_t *perm, float *givnum,
     float *c__, float *s, float *work, cblas_int_t *iwork, cblas_int_t *info);

/* Subroutine */ int slasdq_(const char *uplo, cblas_int_t *sqre, cblas_int_t *n, cblas_int_t *
    ncvt, cblas_int_t *nru, cblas_int_t *ncc, float *d__, float *e, float *vt, 
    cblas_int_t *ldvt, float *u, cblas_int_t *ldu, float *c__, cblas_int_t *ldc, float *
    work, cblas_int_t *info);

/* Subroutine */ int slasdt_(cblas_int_t *n, cblas_int_t *lvl, cblas_int_t *nd, cblas_int_t *
    inode, cblas_int_t *ndiml, cblas_int_t *ndimr, cblas_int_t *msub);

/* Subroutine */ int slaset_(const char *uplo, cblas_int_t *m, cblas_int_t *n, float *alpha, 
    float *beta, float *a, cblas_int_t *lda);

/* Subroutine */ int slasq1_(cblas_int_t *n, float *d__, float *e, float *work, 
    cblas_int_t *info);

/* Subroutine */ int slasq2_(cblas_int_t *n, float *z__, cblas_int_t *info);

/* Subroutine */ int slasq3_(cblas_int_t *i0, cblas_int_t *n0, float *z__, cblas_int_t *pp,
     float *dmin__, float *sigma, float *desig, float *qmax, cblas_int_t *nfail, 
    cblas_int_t *iter, cblas_int_t *ndiv, logical *ieee);

/* Subroutine */ int slasq4_(cblas_int_t *i0, cblas_int_t *n0, float *z__, cblas_int_t *pp,
     cblas_int_t *n0in, float *dmin__, float *dmin1, float *dmin2, float *dn, 
    float *dn1, float *dn2, float *tau, cblas_int_t *ttype);

/* Subroutine */ int slasq5_(cblas_int_t *i0, cblas_int_t *n0, float *z__, cblas_int_t *pp,
     float *tau, float *dmin__, float *dmin1, float *dmin2, float *dn, float *
    dnm1, float *dnm2, logical *ieee);

/* Subroutine */ int slasq6_(cblas_int_t *i0, cblas_int_t *n0, float *z__, cblas_int_t *pp,
     float *dmin__, float *dmin1, float *dmin2, float *dn, float *dnm1, float *
    dnm2);

/* Subroutine */ int slasr_(const char *side, const char *pivot, const char *direct, cblas_int_t *m,
     cblas_int_t *n, float *c__, float *s, float *a, cblas_int_t *lda   );

/* Subroutine */ int slasrt_(const char *id, cblas_int_t *n, float *d__, cblas_int_t *info);

/* Subroutine */ int slassq_(cblas_int_t *n, float *x, cblas_int_t *incx, float *scale, 
    float *sumsq);

/* Subroutine */ int slasv2_(float *f, float *g, float *h__, float *ssmin, float *
    ssmax, float *snr, float *csr, float *snl, float *csl);

/* Subroutine */ int slaswp_(cblas_int_t *n, float *a, cblas_int_t *lda, cblas_int_t *k1, 
    cblas_int_t *k2, cblas_int_t *ipiv, cblas_int_t *incx);

/* Subroutine */ int slasy2_(logical *ltranl, logical *ltranr, cblas_int_t *isgn, 
    cblas_int_t *n1, cblas_int_t *n2, float *tl, cblas_int_t *ldtl, float *tr, cblas_int_t *
    ldtr, float *b, cblas_int_t *ldb, float *scale, float *x, cblas_int_t *ldx, float 
    *xnorm, cblas_int_t *info);

/* Subroutine */ int slasyf_(const char *uplo, cblas_int_t *n, cblas_int_t *nb, cblas_int_t *kb,
     float *a, cblas_int_t *lda, cblas_int_t *ipiv, float *w, cblas_int_t *ldw, cblas_int_t 
    *info);

/* Subroutine */ int slatbs_(const char *uplo, const char *trans, const char *diag, const char *
    normin, cblas_int_t *n, cblas_int_t *kd, float *ab, cblas_int_t *ldab, float *x, 
    float *scale, float *cnorm, cblas_int_t *info);

/* Subroutine */ int slatdf_(cblas_int_t *ijob, cblas_int_t *n, float *z__, cblas_int_t *
    ldz, float *rhs, float *rdsum, float *rdscal, cblas_int_t *ipiv, cblas_int_t *
    jpiv);

/* Subroutine */ int slatps_(const char *uplo, const char *trans, const char *diag, const char *
    normin, cblas_int_t *n, float *ap, float *x, float *scale, float *cnorm, 
    cblas_int_t *info);

/* Subroutine */ int slatrd_(const char *uplo, cblas_int_t *n, cblas_int_t *nb, float *a, 
    cblas_int_t *lda, float *e, float *tau, float *w, cblas_int_t *ldw      );

/* Subroutine */ int slatrs_(const char *uplo, const char *trans, const char *diag, const char *
    normin, cblas_int_t *n, float *a, cblas_int_t *lda, float *x, float *scale, float 
    *cnorm, cblas_int_t *info);

/* Subroutine */ int slatrz_(cblas_int_t *m, cblas_int_t *n, cblas_int_t *l, float *a, 
    cblas_int_t *lda, float *tau, float *work);

/* Subroutine */ int slatzm_(const char *side, cblas_int_t *m, cblas_int_t *n, float *v, 
    cblas_int_t *incv, float *tau, float *c1, float *c2, cblas_int_t *ldc, float *
    work);

/* Subroutine */ int slauu2_(const char *uplo, cblas_int_t *n, float *a, cblas_int_t *lda, 
    cblas_int_t *info);

/* Subroutine */ int slauum_(const char *uplo, cblas_int_t *n, float *a, cblas_int_t *lda, 
    cblas_int_t *info);

/* Subroutine */ int slazq3_(cblas_int_t *i0, cblas_int_t *n0, float *z__, cblas_int_t *pp,
     float *dmin__, float *sigma, float *desig, float *qmax, cblas_int_t *nfail, 
    cblas_int_t *iter, cblas_int_t *ndiv, logical *ieee, cblas_int_t *ttype, float *
    dmin1, float *dmin2, float *dn, float *dn1, float *dn2, float *tau);

/* Subroutine */ int slazq4_(cblas_int_t *i0, cblas_int_t *n0, float *z__, cblas_int_t *pp,
     cblas_int_t *n0in, float *dmin__, float *dmin1, float *dmin2, float *dn, 
    float *dn1, float *dn2, float *tau, cblas_int_t *ttype, float *g);

/* Subroutine */ int sopgtr_(const char *uplo, cblas_int_t *n, float *ap, float *tau, 
    float *q, cblas_int_t *ldq, float *work, cblas_int_t *info);

/* Subroutine */ int sopmtr_(const char *side, const char *uplo, const char *trans, cblas_int_t *m, 
    cblas_int_t *n, float *ap, float *tau, float *c__, cblas_int_t *ldc, float *work, 
    cblas_int_t *info);

/* Subroutine */ int sorg2l_(cblas_int_t *m, cblas_int_t *n, cblas_int_t *k, float *a, 
    cblas_int_t *lda, float *tau, float *work, cblas_int_t *info);

/* Subroutine */ int sorg2r_(cblas_int_t *m, cblas_int_t *n, cblas_int_t *k, float *a, 
    cblas_int_t *lda, float *tau, float *work, cblas_int_t *info);

/* Subroutine */ int sorgbr_(const char *vect, cblas_int_t *m, cblas_int_t *n, cblas_int_t *k, 
    float *a, cblas_int_t *lda, float *tau, float *work, cblas_int_t *lwork, cblas_int_t 
    *info);

/* Subroutine */ int sorghr_(cblas_int_t *n, cblas_int_t *ilo, cblas_int_t *ihi, float *a, 
    cblas_int_t *lda, float *tau, float *work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int sorgl2_(cblas_int_t *m, cblas_int_t *n, cblas_int_t *k, float *a, 
    cblas_int_t *lda, float *tau, float *work, cblas_int_t *info);

/* Subroutine */ int sorglq_(cblas_int_t *m, cblas_int_t *n, cblas_int_t *k, float *a, 
    cblas_int_t *lda, float *tau, float *work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int sorgql_(cblas_int_t *m, cblas_int_t *n, cblas_int_t *k, float *a, 
    cblas_int_t *lda, float *tau, float *work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int sorgqr_(cblas_int_t *m, cblas_int_t *n, cblas_int_t *k, float *a, 
    cblas_int_t *lda, float *tau, float *work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int sorgr2_(cblas_int_t *m, cblas_int_t *n, cblas_int_t *k, float *a, 
    cblas_int_t *lda, float *tau, float *work, cblas_int_t *info);

/* Subroutine */ int sorgrq_(cblas_int_t *m, cblas_int_t *n, cblas_int_t *k, float *a, 
    cblas_int_t *lda, float *tau, float *work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int sorgtr_(const char *uplo, cblas_int_t *n, float *a, cblas_int_t *lda, 
    float *tau, float *work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int sorm2l_(const char *side, const char *trans, cblas_int_t *m, cblas_int_t *n, 
    cblas_int_t *k, float *a, cblas_int_t *lda, float *tau, float *c__, cblas_int_t *ldc,
     float *work, cblas_int_t *info);

/* Subroutine */ int sorm2r_(const char *side, const char *trans, cblas_int_t *m, cblas_int_t *n, 
    cblas_int_t *k, float *a, cblas_int_t *lda, float *tau, float *c__, cblas_int_t *ldc,
     float *work, cblas_int_t *info);

/* Subroutine */ int sormbr_(const char *vect, const char *side, const char *trans, cblas_int_t *m, 
    cblas_int_t *n, cblas_int_t *k, float *a, cblas_int_t *lda, float *tau, float *c__,     
    cblas_int_t *ldc, float *work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int sormhr_(const char *side, const char *trans, cblas_int_t *m, cblas_int_t *n, 
    cblas_int_t *ilo, cblas_int_t *ihi, float *a, cblas_int_t *lda, float *tau, float *     
    c__, cblas_int_t *ldc, float *work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int sorml2_(const char *side, const char *trans, cblas_int_t *m, cblas_int_t *n, 
    cblas_int_t *k, float *a, cblas_int_t *lda, float *tau, float *c__, cblas_int_t *ldc,
     float *work, cblas_int_t *info);

/* Subroutine */ int sormlq_(const char *side, const char *trans, cblas_int_t *m, cblas_int_t *n, 
    cblas_int_t *k, float *a, cblas_int_t *lda, float *tau, float *c__, cblas_int_t *ldc,    
    float *work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int sormql_(const char *side, const char *trans, cblas_int_t *m, cblas_int_t *n, 
    cblas_int_t *k, float *a, cblas_int_t *lda, float *tau, float *c__, cblas_int_t *ldc,    
    float *work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int sormqr_(const char *side, const char *trans, cblas_int_t *m, cblas_int_t *n, 
    cblas_int_t *k, float *a, cblas_int_t *lda, float *tau, float *c__, cblas_int_t *ldc,    
    float *work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int sormr2_(const char *side, const char *trans, cblas_int_t *m, cblas_int_t *n, 
    cblas_int_t *k, float *a, cblas_int_t *lda, float *tau, float *c__, cblas_int_t *ldc,
     float *work, cblas_int_t *info);

/* Subroutine */ int sormr3_(const char *side, const char *trans, cblas_int_t *m, cblas_int_t *n, 
    cblas_int_t *k, cblas_int_t *l, float *a, cblas_int_t *lda, float *tau, float *c__, 
    cblas_int_t *ldc, float *work, cblas_int_t *info);

/* Subroutine */ int sormrq_(const char *side, const char *trans, cblas_int_t *m, cblas_int_t *n, 
    cblas_int_t *k, float *a, cblas_int_t *lda, float *tau, float *c__, cblas_int_t *ldc,    
    float *work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int sormrz_(const char *side, const char *trans, cblas_int_t *m, cblas_int_t *n, 
    cblas_int_t *k, cblas_int_t *l, float *a, cblas_int_t *lda, float *tau, float *c__,     
    cblas_int_t *ldc, float *work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int sormtr_(const char *side, const char *uplo, const char *trans, cblas_int_t *m, 
    cblas_int_t *n, float *a, cblas_int_t *lda, float *tau, float *c__, cblas_int_t *ldc,    
    float *work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int spbcon_(const char *uplo, cblas_int_t *n, cblas_int_t *kd, float *ab, 
    cblas_int_t *ldab, float *anorm, float *rcond, float *work, cblas_int_t *iwork, 
    cblas_int_t *info);

/* Subroutine */ int spbequ_(const char *uplo, cblas_int_t *n, cblas_int_t *kd, float *ab, 
    cblas_int_t *ldab, float *s, float *scond, float *amax, cblas_int_t *info);

/* Subroutine */ int spbrfs_(const char *uplo, cblas_int_t *n, cblas_int_t *kd, cblas_int_t *
    nrhs, float *ab, cblas_int_t *ldab, float *afb, cblas_int_t *ldafb, float *b, 
    cblas_int_t *ldb, float *x, cblas_int_t *ldx, float *ferr, float *berr, float *
    work, cblas_int_t *iwork, cblas_int_t *info);

/* Subroutine */ int spbstf_(const char *uplo, cblas_int_t *n, cblas_int_t *kd, float *ab, 
    cblas_int_t *ldab, cblas_int_t *info);

/* Subroutine */ int spbsv_(const char *uplo, cblas_int_t *n, cblas_int_t *kd, cblas_int_t *
    nrhs, float *ab, cblas_int_t *ldab, float *b, cblas_int_t *ldb, cblas_int_t *info);

/* Subroutine */ int spbsvx_(const char *fact, const char *uplo, cblas_int_t *n, cblas_int_t *kd, 
    cblas_int_t *nrhs, float *ab, cblas_int_t *ldab, float *afb, cblas_int_t *ldafb, 
    const char *equed, float *s, float *b, cblas_int_t *ldb, float *x, cblas_int_t *ldx, 
    float *rcond, float *ferr, float *berr, float *work, cblas_int_t *iwork, 
    cblas_int_t *info);

/* Subroutine */ int spbtf2_(const char *uplo, cblas_int_t *n, cblas_int_t *kd, float *ab, 
    cblas_int_t *ldab, cblas_int_t *info);

/* Subroutine */ int spbtrf_(const char *uplo, cblas_int_t *n, cblas_int_t *kd, float *ab, 
    cblas_int_t *ldab, cblas_int_t *info);

/* Subroutine */ int spbtrs_(const char *uplo, cblas_int_t *n, cblas_int_t *kd, cblas_int_t *
    nrhs, float *ab, cblas_int_t *ldab, float *b, cblas_int_t *ldb, cblas_int_t *info);

/* Subroutine */ int spocon_(const char *uplo, cblas_int_t *n, float *a, cblas_int_t *lda, 
    float *anorm, float *rcond, float *work, cblas_int_t *iwork, cblas_int_t *info);

/* Subroutine */ int spoequ_(cblas_int_t *n, float *a, cblas_int_t *lda, float *s, float 
    *scond, float *amax, cblas_int_t *info);

/* Subroutine */ int sporfs_(const char *uplo, cblas_int_t *n, cblas_int_t *nrhs, float *a, 
    cblas_int_t *lda, float *af, cblas_int_t *ldaf, float *b, cblas_int_t *ldb, float *x,
     cblas_int_t *ldx, float *ferr, float *berr, float *work, cblas_int_t *iwork, 
    cblas_int_t *info);

/* Subroutine */ int sposv_(const char *uplo, cblas_int_t *n, cblas_int_t *nrhs, float *a, 
    cblas_int_t *lda, float *b, cblas_int_t *ldb, cblas_int_t *info);

/* Subroutine */ int sposvx_(const char *fact, const char *uplo, cblas_int_t *n, cblas_int_t *
    nrhs, float *a, cblas_int_t *lda, float *af, cblas_int_t *ldaf, const char *equed, 
    float *s, float *b, cblas_int_t *ldb, float *x, cblas_int_t *ldx, float *rcond, 
    float *ferr, float *berr, float *work, cblas_int_t *iwork, cblas_int_t *info);

/* Subroutine */ int spotf2_(const char *uplo, cblas_int_t *n, float *a, cblas_int_t *lda, 
    cblas_int_t *info);

/* Subroutine */ int spotrf_(const char *uplo, cblas_int_t *n, float *a, cblas_int_t *lda, 
    cblas_int_t *info);

/* Subroutine */ int spotri_(const char *uplo, cblas_int_t *n, float *a, cblas_int_t *lda, 
    cblas_int_t *info);

/* Subroutine */ int spotrs_(const char *uplo, cblas_int_t *n, cblas_int_t *nrhs, float *a, 
    cblas_int_t *lda, float *b, cblas_int_t *ldb, cblas_int_t *info);

/* Subroutine */ int sppcon_(const char *uplo, cblas_int_t *n, float *ap, float *anorm, 
    float *rcond, float *work, cblas_int_t *iwork, cblas_int_t *info);

/* Subroutine */ int sppequ_(const char *uplo, cblas_int_t *n, float *ap, float *s, float *
    scond, float *amax, cblas_int_t *info);

/* Subroutine */ int spprfs_(const char *uplo, cblas_int_t *n, cblas_int_t *nrhs, float *ap, 
    float *afp, float *b, cblas_int_t *ldb, float *x, cblas_int_t *ldx, float *ferr, 
    float *berr, float *work, cblas_int_t *iwork, cblas_int_t *info);

/* Subroutine */ int sppsv_(const char *uplo, cblas_int_t *n, cblas_int_t *nrhs, float *ap, 
    float *b, cblas_int_t *ldb, cblas_int_t *info);

/* Subroutine */ int sppsvx_(const char *fact, const char *uplo, cblas_int_t *n, cblas_int_t *
    nrhs, float *ap, float *afp, const char *equed, float *s, float *b, cblas_int_t *
    ldb, float *x, cblas_int_t *ldx, float *rcond, float *ferr, float *berr, float 
    *work, cblas_int_t *iwork, cblas_int_t *info);

/* Subroutine */ int spptrf_(const char *uplo, cblas_int_t *n, float *ap, cblas_int_t *info);

/* Subroutine */ int spptri_(const char *uplo, cblas_int_t *n, float *ap, cblas_int_t *info);

/* Subroutine */ int spptrs_(const char *uplo, cblas_int_t *n, cblas_int_t *nrhs, float *ap, 
    float *b, cblas_int_t *ldb, cblas_int_t *info);

/* Subroutine */ int sptcon_(cblas_int_t *n, float *d__, float *e, float *anorm, 
    float *rcond, float *work, cblas_int_t *info);

/* Subroutine */ int spteqr_(const char *compz, cblas_int_t *n, float *d__, float *e, 
    float *z__, cblas_int_t *ldz, float *work, cblas_int_t *info);

/* Subroutine */ int sptrfs_(cblas_int_t *n, cblas_int_t *nrhs, float *d__, float *e, 
    float *df, float *ef, float *b, cblas_int_t *ldb, float *x, cblas_int_t *ldx, 
    float *ferr, float *berr, float *work, cblas_int_t *info);

/* Subroutine */ int sptsv_(cblas_int_t *n, cblas_int_t *nrhs, float *d__, float *e, 
    float *b, cblas_int_t *ldb, cblas_int_t *info);

/* Subroutine */ int sptsvx_(const char *fact, cblas_int_t *n, cblas_int_t *nrhs, float *d__,
     float *e, float *df, float *ef, float *b, cblas_int_t *ldb, float *x, cblas_int_t 
    *ldx, float *rcond, float *ferr, float *berr, float *work, cblas_int_t *info);

/* Subroutine */ int spttrf_(cblas_int_t *n, float *d__, float *e, cblas_int_t *info);

/* Subroutine */ int spttrs_(cblas_int_t *n, cblas_int_t *nrhs, float *d__, float *e, 
    float *b, cblas_int_t *ldb, cblas_int_t *info);

/* Subroutine */ int sptts2_(cblas_int_t *n, cblas_int_t *nrhs, float *d__, float *e, 
    float *b, cblas_int_t *ldb);

/* Subroutine */ int srscl_(cblas_int_t *n, float *sa, float *sx, cblas_int_t *incx);

/* Subroutine */ int ssbev_(const char *jobz, const char *uplo, cblas_int_t *n, cblas_int_t *kd, 
    float *ab, cblas_int_t *ldab, float *w, float *z__, cblas_int_t *ldz, float *work,
     cblas_int_t *info);

/* Subroutine */ int ssbevd_(const char *jobz, const char *uplo, cblas_int_t *n, cblas_int_t *kd, 
    float *ab, cblas_int_t *ldab, float *w, float *z__, cblas_int_t *ldz, float *work,
     cblas_int_t *lwork, cblas_int_t *iwork, cblas_int_t *liwork, cblas_int_t *info);

/* Subroutine */ int ssbevx_(const char *jobz, const char *range, const char *uplo, cblas_int_t *n, 
    cblas_int_t *kd, float *ab, cblas_int_t *ldab, float *q, cblas_int_t *ldq, float *vl,
     float *vu, cblas_int_t *il, cblas_int_t *iu, float *abstol, cblas_int_t *m, float *
    w, float *z__, cblas_int_t *ldz, float *work, cblas_int_t *iwork, cblas_int_t *
    ifail, cblas_int_t *info);

/* Subroutine */ int ssbgst_(const char *vect, const char *uplo, cblas_int_t *n, cblas_int_t *ka, 
    cblas_int_t *kb, float *ab, cblas_int_t *ldab, float *bb, cblas_int_t *ldbb, float *
    x, cblas_int_t *ldx, float *work, cblas_int_t *info);

/* Subroutine */ int ssbgv_(const char *jobz, const char *uplo, cblas_int_t *n, cblas_int_t *ka, 
    cblas_int_t *kb, float *ab, cblas_int_t *ldab, float *bb, cblas_int_t *ldbb, float *
    w, float *z__, cblas_int_t *ldz, float *work, cblas_int_t *info);

/* Subroutine */ int ssbgvd_(const char *jobz, const char *uplo, cblas_int_t *n, cblas_int_t *ka, 
    cblas_int_t *kb, float *ab, cblas_int_t *ldab, float *bb, cblas_int_t *ldbb, float *
    w, float *z__, cblas_int_t *ldz, float *work, cblas_int_t *lwork, cblas_int_t *
    iwork, cblas_int_t *liwork, cblas_int_t *info);

/* Subroutine */ int ssbgvx_(const char *jobz, const char *range, const char *uplo, cblas_int_t *n, 
    cblas_int_t *ka, cblas_int_t *kb, float *ab, cblas_int_t *ldab, float *bb, cblas_int_t *
    ldbb, float *q, cblas_int_t *ldq, float *vl, float *vu, cblas_int_t *il, cblas_int_t 
    *iu, float *abstol, cblas_int_t *m, float *w, float *z__, cblas_int_t *ldz, float 
    *work, cblas_int_t *iwork, cblas_int_t *ifail, cblas_int_t *info);

/* Subroutine */ int ssbtrd_(const char *vect, const char *uplo, cblas_int_t *n, cblas_int_t *kd, 
    float *ab, cblas_int_t *ldab, float *d__, float *e, float *q, cblas_int_t *ldq, 
    float *work, cblas_int_t *info);

/* Subroutine */ int sspcon_(const char *uplo, cblas_int_t *n, float *ap, cblas_int_t *ipiv, 
    float *anorm, float *rcond, float *work, cblas_int_t *iwork, cblas_int_t *info);

/* Subroutine */ int sspev_(const char *jobz, const char *uplo, cblas_int_t *n, float *ap, 
    float *w, float *z__, cblas_int_t *ldz, float *work, cblas_int_t *info);

/* Subroutine */ int sspevd_(const char *jobz, const char *uplo, cblas_int_t *n, float *ap, 
    float *w, float *z__, cblas_int_t *ldz, float *work, cblas_int_t *lwork, cblas_int_t 
    *iwork, cblas_int_t *liwork, cblas_int_t *info);

/* Subroutine */ int sspevx_(const char *jobz, const char *range, const char *uplo, cblas_int_t *n, 
    float *ap, float *vl, float *vu, cblas_int_t *il, cblas_int_t *iu, float *abstol, 
    cblas_int_t *m, float *w, float *z__, cblas_int_t *ldz, float *work, cblas_int_t *
    iwork, cblas_int_t *ifail, cblas_int_t *info);

/* Subroutine */ int sspgst_(cblas_int_t *itype, const char *uplo, cblas_int_t *n, float *ap,
     float *bp, cblas_int_t *info);

/* Subroutine */ int sspgv_(cblas_int_t *itype, const char *jobz, const char *uplo, cblas_int_t *
    n, float *ap, float *bp, float *w, float *z__, cblas_int_t *ldz, float *work, 
    cblas_int_t *info);

/* Subroutine */ int sspgvd_(cblas_int_t *itype, const char *jobz, const char *uplo, cblas_int_t *
    n, float *ap, float *bp, float *w, float *z__, cblas_int_t *ldz, float *work, 
    cblas_int_t *lwork, cblas_int_t *iwork, cblas_int_t *liwork, cblas_int_t *info);

/* Subroutine */ int sspgvx_(cblas_int_t *itype, const char *jobz, const char *range, const char *
    uplo, cblas_int_t *n, float *ap, float *bp, float *vl, float *vu, cblas_int_t *il,
     cblas_int_t *iu, float *abstol, cblas_int_t *m, float *w, float *z__, cblas_int_t *
    ldz, float *work, cblas_int_t *iwork, cblas_int_t *ifail, cblas_int_t *info);

/* Subroutine */ int ssprfs_(const char *uplo, cblas_int_t *n, cblas_int_t *nrhs, float *ap, 
    float *afp, cblas_int_t *ipiv, float *b, cblas_int_t *ldb, float *x, cblas_int_t *
    ldx, float *ferr, float *berr, float *work, cblas_int_t *iwork, cblas_int_t *
    info);

/* Subroutine */ int sspsv_(const char *uplo, cblas_int_t *n, cblas_int_t *nrhs, float *ap, 
    cblas_int_t *ipiv, float *b, cblas_int_t *ldb, cblas_int_t *info);

/* Subroutine */ int sspsvx_(const char *fact, const char *uplo, cblas_int_t *n, cblas_int_t *
    nrhs, float *ap, float *afp, cblas_int_t *ipiv, float *b, cblas_int_t *ldb, float 
    *x, cblas_int_t *ldx, float *rcond, float *ferr, float *berr, float *work, 
    cblas_int_t *iwork, cblas_int_t *info);

/* Subroutine */ int ssptrd_(const char *uplo, cblas_int_t *n, float *ap, float *d__, 
    float *e, float *tau, cblas_int_t *info);

/* Subroutine */ int ssptrf_(const char *uplo, cblas_int_t *n, float *ap, cblas_int_t *ipiv, 
    cblas_int_t *info);

/* Subroutine */ int ssptri_(const char *uplo, cblas_int_t *n, float *ap, cblas_int_t *ipiv, 
    float *work, cblas_int_t *info);

/* Subroutine */ int ssptrs_(const char *uplo, cblas_int_t *n, cblas_int_t *nrhs, float *ap, 
    cblas_int_t *ipiv, float *b, cblas_int_t *ldb, cblas_int_t *info);

/* Subroutine */ int sstebz_(const char *range, const char *order, cblas_int_t *n, float *vl, 
    float *vu, cblas_int_t *il, cblas_int_t *iu, float *abstol, float *d__, float *e, 
    cblas_int_t *m, cblas_int_t *nsplit, float *w, cblas_int_t *iblock, cblas_int_t *
    isplit, float *work, cblas_int_t *iwork, cblas_int_t *info);

/* Subroutine */ int sstedc_(const char *compz, cblas_int_t *n, float *d__, float *e, 
    float *z__, cblas_int_t *ldz, float *work, cblas_int_t *lwork, cblas_int_t *iwork, 
    cblas_int_t *liwork, cblas_int_t *info);

/* Subroutine */ int sstegr_(const char *jobz, const char *range, cblas_int_t *n, float *d__, 
    float *e, float *vl, float *vu, cblas_int_t *il, cblas_int_t *iu, float *abstol, 
    cblas_int_t *m, float *w, float *z__, cblas_int_t *ldz, cblas_int_t *isuppz, float *
    work, cblas_int_t *lwork, cblas_int_t *iwork, cblas_int_t *liwork, cblas_int_t *info);

/* Subroutine */ int sstein_(cblas_int_t *n, float *d__, float *e, cblas_int_t *m, float 
    *w, cblas_int_t *iblock, cblas_int_t *isplit, float *z__, cblas_int_t *ldz, float *
    work, cblas_int_t *iwork, cblas_int_t *ifail, cblas_int_t *info);

/* Subroutine */ int sstemr_(const char *jobz, const char *range, cblas_int_t *n, float *d__, 
    float *e, float *vl, float *vu, cblas_int_t *il, cblas_int_t *iu, cblas_int_t *m, 
    float *w, float *z__, cblas_int_t *ldz, cblas_int_t *nzc, cblas_int_t *isuppz, 
    logical *tryrac, float *work, cblas_int_t *lwork, cblas_int_t *iwork, cblas_int_t *
    liwork, cblas_int_t *info);

/* Subroutine */ int ssteqr_(const char *compz, cblas_int_t *n, float *d__, float *e, 
    float *z__, cblas_int_t *ldz, float *work, cblas_int_t *info);

/* Subroutine */ int ssterf_(cblas_int_t *n, float *d__, float *e, cblas_int_t *info);

/* Subroutine */ int sstev_(const char *jobz, cblas_int_t *n, float *d__, float *e, float *
    z__, cblas_int_t *ldz, float *work, cblas_int_t *info);

/* Subroutine */ int sstevd_(const char *jobz, cblas_int_t *n, float *d__, float *e, float 
    *z__, cblas_int_t *ldz, float *work, cblas_int_t *lwork, cblas_int_t *iwork, 
    cblas_int_t *liwork, cblas_int_t *info);

/* Subroutine */ int sstevr_(const char *jobz, const char *range, cblas_int_t *n, float *d__, 
    float *e, float *vl, float *vu, cblas_int_t *il, cblas_int_t *iu, float *abstol, 
    cblas_int_t *m, float *w, float *z__, cblas_int_t *ldz, cblas_int_t *isuppz, float *
    work, cblas_int_t *lwork, cblas_int_t *iwork, cblas_int_t *liwork, cblas_int_t *info);

/* Subroutine */ int sstevx_(const char *jobz, const char *range, cblas_int_t *n, float *d__, 
    float *e, float *vl, float *vu, cblas_int_t *il, cblas_int_t *iu, float *abstol, 
    cblas_int_t *m, float *w, float *z__, cblas_int_t *ldz, float *work, cblas_int_t *
    iwork, cblas_int_t *ifail, cblas_int_t *info);

/* Subroutine */ int ssycon_(const char *uplo, cblas_int_t *n, float *a, cblas_int_t *lda, 
    cblas_int_t *ipiv, float *anorm, float *rcond, float *work, cblas_int_t *iwork, 
    cblas_int_t *info);

/* Subroutine */ int ssyev_(const char *jobz, const char *uplo, cblas_int_t *n, float *a, 
    cblas_int_t *lda, float *w, float *work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int ssyevd_(const char *jobz, const char *uplo, cblas_int_t *n, float *a, 
    cblas_int_t *lda, float *w, float *work, cblas_int_t *lwork, cblas_int_t *iwork, 
    cblas_int_t *liwork, cblas_int_t *info);

/* Subroutine */ int ssyevr_(const char *jobz, const char *range, const char *uplo, cblas_int_t *n, 
    float *a, cblas_int_t *lda, float *vl, float *vu, cblas_int_t *il, cblas_int_t *iu, 
    float *abstol, cblas_int_t *m, float *w, float *z__, cblas_int_t *ldz, cblas_int_t *
    isuppz, float *work, cblas_int_t *lwork, cblas_int_t *iwork, cblas_int_t *liwork, 
    cblas_int_t *info);

/* Subroutine */ int ssyevx_(const char *jobz, const char *range, const char *uplo, cblas_int_t *n, 
    float *a, cblas_int_t *lda, float *vl, float *vu, cblas_int_t *il, cblas_int_t *iu, 
    float *abstol, cblas_int_t *m, float *w, float *z__, cblas_int_t *ldz, float *
    work, cblas_int_t *lwork, cblas_int_t *iwork, cblas_int_t *ifail, cblas_int_t *info);

/* Subroutine */ int ssygs2_(cblas_int_t *itype, const char *uplo, cblas_int_t *n, float *a, 
    cblas_int_t *lda, float *b, cblas_int_t *ldb, cblas_int_t *info);

/* Subroutine */ int ssygst_(cblas_int_t *itype, const char *uplo, cblas_int_t *n, float *a, 
    cblas_int_t *lda, float *b, cblas_int_t *ldb, cblas_int_t *info);

/* Subroutine */ int ssygv_(cblas_int_t *itype, const char *jobz, const char *uplo, cblas_int_t *
    n, float *a, cblas_int_t *lda, float *b, cblas_int_t *ldb, float *w, float *work, 
    cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int ssygvd_(cblas_int_t *itype, const char *jobz, const char *uplo, cblas_int_t *
    n, float *a, cblas_int_t *lda, float *b, cblas_int_t *ldb, float *w, float *work, 
    cblas_int_t *lwork, cblas_int_t *iwork, cblas_int_t *liwork, cblas_int_t *info);

/* Subroutine */ int ssygvx_(cblas_int_t *itype, const char *jobz, const char *range, const char *
    uplo, cblas_int_t *n, float *a, cblas_int_t *lda, float *b, cblas_int_t *ldb, float *
    vl, float *vu, cblas_int_t *il, cblas_int_t *iu, float *abstol, cblas_int_t *m, 
    float *w, float *z__, cblas_int_t *ldz, float *work, cblas_int_t *lwork, cblas_int_t    
    *iwork, cblas_int_t *ifail, cblas_int_t *info);

/* Subroutine */ int ssyrfs_(const char *uplo, cblas_int_t *n, cblas_int_t *nrhs, float *a, 
    cblas_int_t *lda, float *af, cblas_int_t *ldaf, cblas_int_t *ipiv, float *b, 
    cblas_int_t *ldb, float *x, cblas_int_t *ldx, float *ferr, float *berr, float *
    work, cblas_int_t *iwork, cblas_int_t *info);

/* Subroutine */ int ssysv_(const char *uplo, cblas_int_t *n, cblas_int_t *nrhs, float *a, 
    cblas_int_t *lda, cblas_int_t *ipiv, float *b, cblas_int_t *ldb, float *work, 
    cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int ssysvx_(const char *fact, const char *uplo, cblas_int_t *n, cblas_int_t *
    nrhs, float *a, cblas_int_t *lda, float *af, cblas_int_t *ldaf, cblas_int_t *ipiv, 
    float *b, cblas_int_t *ldb, float *x, cblas_int_t *ldx, float *rcond, float *ferr,
     float *berr, float *work, cblas_int_t *lwork, cblas_int_t *iwork, cblas_int_t *
    info);

/* Subroutine */ int ssytd2_(const char *uplo, cblas_int_t *n, float *a, cblas_int_t *lda, 
    float *d__, float *e, float *tau, cblas_int_t *info);

/* Subroutine */ int ssytf2_(const char *uplo, cblas_int_t *n, float *a, cblas_int_t *lda, 
    cblas_int_t *ipiv, cblas_int_t *info);

/* Subroutine */ int ssytrd_(const char *uplo, cblas_int_t *n, float *a, cblas_int_t *lda, 
    float *d__, float *e, float *tau, float *work, cblas_int_t *lwork, cblas_int_t *
    info);

/* Subroutine */ int ssytrf_(const char *uplo, cblas_int_t *n, float *a, cblas_int_t *lda,      
    cblas_int_t *ipiv, float *work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int ssytri_(const char *uplo, cblas_int_t *n, float *a, cblas_int_t *lda, 
    cblas_int_t *ipiv, float *work, cblas_int_t *info);

/* Subroutine */ int ssytrs_(const char *uplo, cblas_int_t *n, cblas_int_t *nrhs, float *a, 
    cblas_int_t *lda, cblas_int_t *ipiv, float *b, cblas_int_t *ldb, cblas_int_t *info);

/* Subroutine */ int stbcon_(const char *norm, const char *uplo, const char *diag, cblas_int_t *n, 
    cblas_int_t *kd, float *ab, cblas_int_t *ldab, float *rcond, float *work, 
    cblas_int_t *iwork, cblas_int_t *info);

/* Subroutine */ int stbrfs_(const char *uplo, const char *trans, const char *diag, cblas_int_t *n, 
    cblas_int_t *kd, cblas_int_t *nrhs, float *ab, cblas_int_t *ldab, float *b, cblas_int_t 
    *ldb, float *x, cblas_int_t *ldx, float *ferr, float *berr, float *work, 
    cblas_int_t *iwork, cblas_int_t *info);

/* Subroutine */ int stbtrs_(const char *uplo, const char *trans, const char *diag, cblas_int_t *n, 
    cblas_int_t *kd, cblas_int_t *nrhs, float *ab, cblas_int_t *ldab, float *b, cblas_int_t 
    *ldb, cblas_int_t *info);

/* Subroutine */ int stgevc_(const char *side, const char *howmny, logical *select, 
    cblas_int_t *n, float *s, cblas_int_t *lds, float *p, cblas_int_t *ldp, float *vl, 
    cblas_int_t *ldvl, float *vr, cblas_int_t *ldvr, cblas_int_t *mm, cblas_int_t *m, float 
    *work, cblas_int_t *info);

/* Subroutine */ int stgex2_(logical *wantq, logical *wantz, cblas_int_t *n, float 
    *a, cblas_int_t *lda, float *b, cblas_int_t *ldb, float *q, cblas_int_t *ldq, float *
    z__, cblas_int_t *ldz, cblas_int_t *j1, cblas_int_t *n1, cblas_int_t *n2, float *work, 
    cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int stgexc_(logical *wantq, logical *wantz, cblas_int_t *n, float 
    *a, cblas_int_t *lda, float *b, cblas_int_t *ldb, float *q, cblas_int_t *ldq, float *
    z__, cblas_int_t *ldz, cblas_int_t *ifst, cblas_int_t *ilst, float *work, cblas_int_t *
    lwork, cblas_int_t *info);

/* Subroutine */ int stgsen_(cblas_int_t *ijob, logical *wantq, logical *wantz, 
    logical *select, cblas_int_t *n, float *a, cblas_int_t *lda, float *b, cblas_int_t *
    ldb, float *alphar, float *alphai, float *beta, float *q, cblas_int_t *ldq, 
    float *z__, cblas_int_t *ldz, cblas_int_t *m, float *pl, float *pr, float *dif, 
    float *work, cblas_int_t *lwork, cblas_int_t *iwork, cblas_int_t *liwork, cblas_int_t *
    info);

/* Subroutine */ int stgsja_(const char *jobu, const char *jobv, const char *jobq, cblas_int_t *m, 
    cblas_int_t *p, cblas_int_t *n, cblas_int_t *k, cblas_int_t *l, float *a, cblas_int_t *lda,
     float *b, cblas_int_t *ldb, float *tola, float *tolb, float *alpha, float *
    beta, float *u, cblas_int_t *ldu, float *v, cblas_int_t *ldv, float *q, cblas_int_t *
    ldq, float *work, cblas_int_t *ncycle, cblas_int_t *info);

/* Subroutine */ int stgsna_(const char *job, const char *howmny, logical *select, 
    cblas_int_t *n, float *a, cblas_int_t *lda, float *b, cblas_int_t *ldb, float *vl, 
    cblas_int_t *ldvl, float *vr, cblas_int_t *ldvr, float *s, float *dif, cblas_int_t *
    mm, cblas_int_t *m, float *work, cblas_int_t *lwork, cblas_int_t *iwork, cblas_int_t *
    info);

/* Subroutine */ int stgsy2_(const char *trans, cblas_int_t *ijob, cblas_int_t *m, cblas_int_t *
    n, float *a, cblas_int_t *lda, float *b, cblas_int_t *ldb, float *c__, cblas_int_t *
    ldc, float *d__, cblas_int_t *ldd, float *e, cblas_int_t *lde, float *f, cblas_int_t 
    *ldf, float *scale, float *rdsum, float *rdscal, cblas_int_t *iwork, cblas_int_t 
    *pq, cblas_int_t *info);

/* Subroutine */ int stgsyl_(const char *trans, cblas_int_t *ijob, cblas_int_t *m, cblas_int_t *
    n, float *a, cblas_int_t *lda, float *b, cblas_int_t *ldb, float *c__, cblas_int_t *
    ldc, float *d__, cblas_int_t *ldd, float *e, cblas_int_t *lde, float *f, cblas_int_t 
    *ldf, float *scale, float *dif, float *work, cblas_int_t *lwork, cblas_int_t *
    iwork, cblas_int_t *info);

/* Subroutine */ int stpcon_(const char *norm, const char *uplo, const char *diag, cblas_int_t *n, 
    float *ap, float *rcond, float *work, cblas_int_t *iwork, cblas_int_t *info);

/* Subroutine */ int stprfs_(const char *uplo, const char *trans, const char *diag, cblas_int_t *n, 
    cblas_int_t *nrhs, float *ap, float *b, cblas_int_t *ldb, float *x, cblas_int_t *ldx,
     float *ferr, float *berr, float *work, cblas_int_t *iwork, cblas_int_t *info);

/* Subroutine */ int stptri_(const char *uplo, const char *diag, cblas_int_t *n, float *ap, 
    cblas_int_t *info);

/* Subroutine */ int stptrs_(const char *uplo, const char *trans, const char *diag, cblas_int_t *n, 
    cblas_int_t *nrhs, float *ap, float *b, cblas_int_t *ldb, cblas_int_t *info);

/* Subroutine */ int strcon_(const char *norm, const char *uplo, const char *diag, cblas_int_t *n, 
    float *a, cblas_int_t *lda, float *rcond, float *work, cblas_int_t *iwork, 
    cblas_int_t *info);

/* Subroutine */ int strevc_(const char *side, const char *howmny, logical *select, 
    cblas_int_t *n, float *t, cblas_int_t *ldt, float *vl, cblas_int_t *ldvl, float *vr, 
    cblas_int_t *ldvr, cblas_int_t *mm, cblas_int_t *m, float *work, cblas_int_t *info);

/* Subroutine */ int strexc_(const char *compq, cblas_int_t *n, float *t, cblas_int_t *ldt, 
    float *q, cblas_int_t *ldq, cblas_int_t *ifst, cblas_int_t *ilst, float *work, 
    cblas_int_t *info);

/* Subroutine */ int strrfs_(const char *uplo, const char *trans, const char *diag, cblas_int_t *n, 
    cblas_int_t *nrhs, float *a, cblas_int_t *lda, float *b, cblas_int_t *ldb, float *x, 
    cblas_int_t *ldx, float *ferr, float *berr, float *work, cblas_int_t *iwork, 
    cblas_int_t *info);

/* Subroutine */ int strsen_(const char *job, const char *compq, logical *select, cblas_int_t 
    *n, float *t, cblas_int_t *ldt, float *q, cblas_int_t *ldq, float *wr, float *wi, 
    cblas_int_t *m, float *s, float *sep, float *work, cblas_int_t *lwork, cblas_int_t *
    iwork, cblas_int_t *liwork, cblas_int_t *info);

/* Subroutine */ int strsna_(const char *job, const char *howmny, logical *select, 
    cblas_int_t *n, float *t, cblas_int_t *ldt, float *vl, cblas_int_t *ldvl, float *vr, 
    cblas_int_t *ldvr, float *s, float *sep, cblas_int_t *mm, cblas_int_t *m, float *
    work, cblas_int_t *ldwork, cblas_int_t *iwork, cblas_int_t *info);

/* Subroutine */ int strsyl_(const char *trana, const char *tranb, cblas_int_t *isgn, cblas_int_t 
    *m, cblas_int_t *n, float *a, cblas_int_t *lda, float *b, cblas_int_t *ldb, float *
    c__, cblas_int_t *ldc, float *scale, cblas_int_t *info);

/* Subroutine */ int strti2_(const char *uplo, const char *diag, cblas_int_t *n, float *a, 
    cblas_int_t *lda, cblas_int_t *info);

/* Subroutine */ int strtri_(const char *uplo, const char *diag, cblas_int_t *n, float *a, 
    cblas_int_t *lda, cblas_int_t *info);

/* Subroutine */ int strtrs_(const char *uplo, const char *trans, const char *diag, cblas_int_t *n, 
    cblas_int_t *nrhs, float *a, cblas_int_t *lda, float *b, cblas_int_t *ldb, cblas_int_t *
    info);

/* Subroutine */ int stzrqf_(cblas_int_t *m, cblas_int_t *n, float *a, cblas_int_t *lda, 
    float *tau, cblas_int_t *info);

/* Subroutine */ int stzrzf_(cblas_int_t *m, cblas_int_t *n, float *a, cblas_int_t *lda, 
    float *tau, float *work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int xerbla_(const char *srname, cblas_int_t *info);

/* Subroutine */ int zbdsqr_(const char *uplo, cblas_int_t *n, cblas_int_t *ncvt, cblas_int_t *
    nru, cblas_int_t *ncc, double *d__, double *e, doublecomplex *vt, 
    cblas_int_t *ldvt, doublecomplex *u, cblas_int_t *ldu, doublecomplex *c__, 
    cblas_int_t *ldc, double *rwork, cblas_int_t *info);

/* Subroutine */ int zcgesv_(cblas_int_t *n, cblas_int_t *nrhs, doublecomplex *a, 
    cblas_int_t *lda, cblas_int_t *ipiv, doublecomplex *b, cblas_int_t *ldb, 
    doublecomplex *x, cblas_int_t *ldx, doublecomplex *work, complex *swork, 
    cblas_int_t *iter, cblas_int_t *info);

/* Subroutine */ int zdrscl_(cblas_int_t *n, double *sa, doublecomplex *sx, 
    cblas_int_t *incx);

/* Subroutine */ int zgbbrd_(const char *vect, cblas_int_t *m, cblas_int_t *n, cblas_int_t *ncc,
     cblas_int_t *kl, cblas_int_t *ku, doublecomplex *ab, cblas_int_t *ldab, 
    double *d__, double *e, doublecomplex *q, cblas_int_t *ldq, 
    doublecomplex *pt, cblas_int_t *ldpt, doublecomplex *c__, cblas_int_t *ldc, 
    doublecomplex *work, double *rwork, cblas_int_t *info);

/* Subroutine */ int zgbcon_(const char *norm, cblas_int_t *n, cblas_int_t *kl, cblas_int_t *ku,
     doublecomplex *ab, cblas_int_t *ldab, cblas_int_t *ipiv, double *anorm, 
    double *rcond, doublecomplex *work, double *rwork, cblas_int_t *
    info);

/* Subroutine */ int zgbequ_(cblas_int_t *m, cblas_int_t *n, cblas_int_t *kl, cblas_int_t *ku,
     doublecomplex *ab, cblas_int_t *ldab, double *r__, double *c__, 
    double *rowcnd, double *colcnd, double *amax, cblas_int_t *
    info);

/* Subroutine */ int zgbrfs_(const char *trans, cblas_int_t *n, cblas_int_t *kl, cblas_int_t *
    ku, cblas_int_t *nrhs, doublecomplex *ab, cblas_int_t *ldab, doublecomplex *
    afb, cblas_int_t *ldafb, cblas_int_t *ipiv, doublecomplex *b, cblas_int_t *ldb, 
    doublecomplex *x, cblas_int_t *ldx, double *ferr, double *berr, 
    doublecomplex *work, double *rwork, cblas_int_t *info);

/* Subroutine */ int zgbsv_(cblas_int_t *n, cblas_int_t *kl, cblas_int_t *ku, cblas_int_t *
    nrhs, doublecomplex *ab, cblas_int_t *ldab, cblas_int_t *ipiv, doublecomplex *
    b, cblas_int_t *ldb, cblas_int_t *info);

/* Subroutine */ int zgbsvx_(const char *fact, const char *trans, cblas_int_t *n, cblas_int_t *kl,
     cblas_int_t *ku, cblas_int_t *nrhs, doublecomplex *ab, cblas_int_t *ldab, 
    doublecomplex *afb, cblas_int_t *ldafb, cblas_int_t *ipiv, const char *equed, 
    double *r__, double *c__, doublecomplex *b, cblas_int_t *ldb, 
    doublecomplex *x, cblas_int_t *ldx, double *rcond, double *ferr, 
    double *berr, doublecomplex *work, double *rwork, cblas_int_t *
    info);

/* Subroutine */ int zgbtf2_(cblas_int_t *m, cblas_int_t *n, cblas_int_t *kl, cblas_int_t *ku,
     doublecomplex *ab, cblas_int_t *ldab, cblas_int_t *ipiv, cblas_int_t *info);

/* Subroutine */ int zgbtrf_(cblas_int_t *m, cblas_int_t *n, cblas_int_t *kl, cblas_int_t *ku,
     doublecomplex *ab, cblas_int_t *ldab, cblas_int_t *ipiv, cblas_int_t *info);

/* Subroutine */ int zgbtrs_(const char *trans, cblas_int_t *n, cblas_int_t *kl, cblas_int_t *
    ku, cblas_int_t *nrhs, doublecomplex *ab, cblas_int_t *ldab, cblas_int_t *ipiv, 
    doublecomplex *b, cblas_int_t *ldb, cblas_int_t *info);

/* Subroutine */ int zgebak_(const char *job, const char *side, cblas_int_t *n, cblas_int_t *ilo, 
    cblas_int_t *ihi, double *scale, cblas_int_t *m, doublecomplex *v, 
    cblas_int_t *ldv, cblas_int_t *info);

/* Subroutine */ int zgebal_(const char *job, cblas_int_t *n, doublecomplex *a, cblas_int_t 
    *lda, cblas_int_t *ilo, cblas_int_t *ihi, double *scale, cblas_int_t *info);

/* Subroutine */ int zgebd2_(cblas_int_t *m, cblas_int_t *n, doublecomplex *a, 
    cblas_int_t *lda, double *d__, double *e, doublecomplex *tauq, 
    doublecomplex *taup, doublecomplex *work, cblas_int_t *info);

/* Subroutine */ int zgebrd_(cblas_int_t *m, cblas_int_t *n, doublecomplex *a, 
    cblas_int_t *lda, double *d__, double *e, doublecomplex *tauq, 
    doublecomplex *taup, doublecomplex *work, cblas_int_t *lwork, cblas_int_t *
    info);

/* Subroutine */ int zgecon_(const char *norm, cblas_int_t *n, doublecomplex *a, 
    cblas_int_t *lda, double *anorm, double *rcond, doublecomplex *
    work, double *rwork, cblas_int_t *info);

/* Subroutine */ int zgeequ_(cblas_int_t *m, cblas_int_t *n, doublecomplex *a, 
    cblas_int_t *lda, double *r__, double *c__, double *rowcnd, 
    double *colcnd, double *amax, cblas_int_t *info);

/* Subroutine */ int zgees_(const char *jobvs, const char *sort, L_fp select, cblas_int_t *n, 
    doublecomplex *a, cblas_int_t *lda, cblas_int_t *sdim, doublecomplex *w, 
    doublecomplex *vs, cblas_int_t *ldvs, doublecomplex *work, cblas_int_t *lwork,
     double *rwork, logical *bwork, cblas_int_t *info);

/* Subroutine */ int zgeesx_(const char *jobvs, const char *sort, L_fp select, const char *
    sense, cblas_int_t *n, doublecomplex *a, cblas_int_t *lda, cblas_int_t *sdim, 
    doublecomplex *w, doublecomplex *vs, cblas_int_t *ldvs, double *
    rconde, double *rcondv, doublecomplex *work, cblas_int_t *lwork, 
    double *rwork, logical *bwork, cblas_int_t *info);

/* Subroutine */ int zgeev_(const char *jobvl, const char *jobvr, cblas_int_t *n, 
    doublecomplex *a, cblas_int_t *lda, doublecomplex *w, doublecomplex *vl, 
    cblas_int_t *ldvl, doublecomplex *vr, cblas_int_t *ldvr, doublecomplex *work, 
    cblas_int_t *lwork, double *rwork, cblas_int_t *info);

/* Subroutine */ int zgeevx_(const char *balanc, const char *jobvl, const char *jobvr, const char *
    sense, cblas_int_t *n, doublecomplex *a, cblas_int_t *lda, doublecomplex *w, 
    doublecomplex *vl, cblas_int_t *ldvl, doublecomplex *vr, cblas_int_t *ldvr, 
    cblas_int_t *ilo, cblas_int_t *ihi, double *scale, double *abnrm, 
    double *rconde, double *rcondv, doublecomplex *work, cblas_int_t *  
    lwork, double *rwork, cblas_int_t *info);

/* Subroutine */ int zgegs_(const char *jobvsl, const char *jobvsr, cblas_int_t *n, 
    doublecomplex *a, cblas_int_t *lda, doublecomplex *b, cblas_int_t *ldb, 
    doublecomplex *alpha, doublecomplex *beta, doublecomplex *vsl, 
    cblas_int_t *ldvsl, doublecomplex *vsr, cblas_int_t *ldvsr, doublecomplex *     
    work, cblas_int_t *lwork, double *rwork, cblas_int_t *info);

/* Subroutine */ int zgegv_(const char *jobvl, const char *jobvr, cblas_int_t *n, 
    doublecomplex *a, cblas_int_t *lda, doublecomplex *b, cblas_int_t *ldb, 
    doublecomplex *alpha, doublecomplex *beta, doublecomplex *vl, cblas_int_t 
    *ldvl, doublecomplex *vr, cblas_int_t *ldvr, doublecomplex *work, cblas_int_t 
    *lwork, double *rwork, cblas_int_t *info);

/* Subroutine */ int zgehd2_(cblas_int_t *n, cblas_int_t *ilo, cblas_int_t *ihi, 
    doublecomplex *a, cblas_int_t *lda, doublecomplex *tau, doublecomplex *
    work, cblas_int_t *info);

/* Subroutine */ int zgehrd_(cblas_int_t *n, cblas_int_t *ilo, cblas_int_t *ihi, 
    doublecomplex *a, cblas_int_t *lda, doublecomplex *tau, doublecomplex *
    work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int zgelq2_(cblas_int_t *m, cblas_int_t *n, doublecomplex *a, 
    cblas_int_t *lda, doublecomplex *tau, doublecomplex *work, cblas_int_t *info);

/* Subroutine */ int zgelqf_(cblas_int_t *m, cblas_int_t *n, doublecomplex *a, 
    cblas_int_t *lda, doublecomplex *tau, doublecomplex *work, cblas_int_t *lwork,
     cblas_int_t *info);

/* Subroutine */ int zgels_(const char *trans, cblas_int_t *m, cblas_int_t *n, cblas_int_t *
    nrhs, doublecomplex *a, cblas_int_t *lda, doublecomplex *b, cblas_int_t *ldb, 
    doublecomplex *work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int zgelsd_(cblas_int_t *m, cblas_int_t *n, cblas_int_t *nrhs, 
    doublecomplex *a, cblas_int_t *lda, doublecomplex *b, cblas_int_t *ldb, 
    double *s, double *rcond, cblas_int_t *rank, doublecomplex *work, 
    cblas_int_t *lwork, double *rwork, cblas_int_t *iwork, cblas_int_t *info);

/* Subroutine */ int zgelss_(cblas_int_t *m, cblas_int_t *n, cblas_int_t *nrhs, 
    doublecomplex *a, cblas_int_t *lda, doublecomplex *b, cblas_int_t *ldb, 
    double *s, double *rcond, cblas_int_t *rank, doublecomplex *work, 
    cblas_int_t *lwork, double *rwork, cblas_int_t *info);

/* Subroutine */ int zgelsx_(cblas_int_t *m, cblas_int_t *n, cblas_int_t *nrhs, 
    doublecomplex *a, cblas_int_t *lda, doublecomplex *b, cblas_int_t *ldb, 
    cblas_int_t *jpvt, double *rcond, cblas_int_t *rank, doublecomplex *work, 
    double *rwork, cblas_int_t *info);

/* Subroutine */ int zgelsy_(cblas_int_t *m, cblas_int_t *n, cblas_int_t *nrhs, 
    doublecomplex *a, cblas_int_t *lda, doublecomplex *b, cblas_int_t *ldb, 
    cblas_int_t *jpvt, double *rcond, cblas_int_t *rank, doublecomplex *work, 
    cblas_int_t *lwork, double *rwork, cblas_int_t *info);

/* Subroutine */ int zgeql2_(cblas_int_t *m, cblas_int_t *n, doublecomplex *a, 
    cblas_int_t *lda, doublecomplex *tau, doublecomplex *work, cblas_int_t *info);

/* Subroutine */ int zgeqlf_(cblas_int_t *m, cblas_int_t *n, doublecomplex *a, 
    cblas_int_t *lda, doublecomplex *tau, doublecomplex *work, cblas_int_t *lwork,
     cblas_int_t *info);

/* Subroutine */ int zgeqp3_(cblas_int_t *m, cblas_int_t *n, doublecomplex *a, 
    cblas_int_t *lda, cblas_int_t *jpvt, doublecomplex *tau, doublecomplex *work, 
    cblas_int_t *lwork, double *rwork, cblas_int_t *info);

/* Subroutine */ int zgeqpf_(cblas_int_t *m, cblas_int_t *n, doublecomplex *a, 
    cblas_int_t *lda, cblas_int_t *jpvt, doublecomplex *tau, doublecomplex *work, 
    double *rwork, cblas_int_t *info);

/* Subroutine */ int zgeqr2_(cblas_int_t *m, cblas_int_t *n, doublecomplex *a, 
    cblas_int_t *lda, doublecomplex *tau, doublecomplex *work, cblas_int_t *info);

/* Subroutine */ int zgeqrf_(cblas_int_t *m, cblas_int_t *n, doublecomplex *a, 
    cblas_int_t *lda, doublecomplex *tau, doublecomplex *work, cblas_int_t *lwork,
     cblas_int_t *info);

/* Subroutine */ int zgerfs_(const char *trans, cblas_int_t *n, cblas_int_t *nrhs, 
    doublecomplex *a, cblas_int_t *lda, doublecomplex *af, cblas_int_t *ldaf, 
    cblas_int_t *ipiv, doublecomplex *b, cblas_int_t *ldb, doublecomplex *x, 
    cblas_int_t *ldx, double *ferr, double *berr, doublecomplex *work,
     double *rwork, cblas_int_t *info);

/* Subroutine */ int zgerq2_(cblas_int_t *m, cblas_int_t *n, doublecomplex *a, 
    cblas_int_t *lda, doublecomplex *tau, doublecomplex *work, cblas_int_t *info);

/* Subroutine */ int zgerqf_(cblas_int_t *m, cblas_int_t *n, doublecomplex *a, 
    cblas_int_t *lda, doublecomplex *tau, doublecomplex *work, cblas_int_t *lwork,
     cblas_int_t *info);

/* Subroutine */ int zgesc2_(cblas_int_t *n, doublecomplex *a, cblas_int_t *lda, 
    doublecomplex *rhs, cblas_int_t *ipiv, cblas_int_t *jpiv, double *scale);

/* Subroutine */ int zgesdd_(const char *jobz, cblas_int_t *m, cblas_int_t *n, 
    doublecomplex *a, cblas_int_t *lda, double *s, doublecomplex *u, 
    cblas_int_t *ldu, doublecomplex *vt, cblas_int_t *ldvt, doublecomplex *work, 
    cblas_int_t *lwork, double *rwork, cblas_int_t *iwork, cblas_int_t *info);

/* Subroutine */ int zgesv_(cblas_int_t *n, cblas_int_t *nrhs, doublecomplex *a, 
    cblas_int_t *lda, cblas_int_t *ipiv, doublecomplex *b, cblas_int_t *ldb, cblas_int_t *
    info);

/* Subroutine */ int zgesvd_(const char *jobu, const char *jobvt, cblas_int_t *m, cblas_int_t *n, 
    doublecomplex *a, cblas_int_t *lda, double *s, doublecomplex *u, 
    cblas_int_t *ldu, doublecomplex *vt, cblas_int_t *ldvt, doublecomplex *work, 
    cblas_int_t *lwork, double *rwork, cblas_int_t *info);

/* Subroutine */ int zgesvx_(const char *fact, const char *trans, cblas_int_t *n, cblas_int_t *
    nrhs, doublecomplex *a, cblas_int_t *lda, doublecomplex *af, cblas_int_t *
    ldaf, cblas_int_t *ipiv, const char *equed, double *r__, double *c__, 
    doublecomplex *b, cblas_int_t *ldb, doublecomplex *x, cblas_int_t *ldx, 
    double *rcond, double *ferr, double *berr, doublecomplex *
    work, double *rwork, cblas_int_t *info);

/* Subroutine */ int zgetc2_(cblas_int_t *n, doublecomplex *a, cblas_int_t *lda, 
    cblas_int_t *ipiv, cblas_int_t *jpiv, cblas_int_t *info);

/* Subroutine */ int zgetf2_(cblas_int_t *m, cblas_int_t *n, doublecomplex *a, 
    cblas_int_t *lda, cblas_int_t *ipiv, cblas_int_t *info);

/* Subroutine */ int zgetrf_(cblas_int_t *m, cblas_int_t *n, doublecomplex *a, 
    cblas_int_t *lda, cblas_int_t *ipiv, cblas_int_t *info);

/* Subroutine */ int zgetri_(cblas_int_t *n, doublecomplex *a, cblas_int_t *lda, 
    cblas_int_t *ipiv, doublecomplex *work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int zgetrs_(const char *trans, cblas_int_t *n, cblas_int_t *nrhs, 
    doublecomplex *a, cblas_int_t *lda, cblas_int_t *ipiv, doublecomplex *b, 
    cblas_int_t *ldb, cblas_int_t *info);

/* Subroutine */ int zggbak_(const char *job, const char *side, cblas_int_t *n, cblas_int_t *ilo, 
    cblas_int_t *ihi, double *lscale, double *rscale, cblas_int_t *m, 
    doublecomplex *v, cblas_int_t *ldv, cblas_int_t *info);

/* Subroutine */ int zggbal_(const char *job, cblas_int_t *n, doublecomplex *a, cblas_int_t 
    *lda, doublecomplex *b, cblas_int_t *ldb, cblas_int_t *ilo, cblas_int_t *ihi, 
    double *lscale, double *rscale, double *work, cblas_int_t *
    info);

/* Subroutine */ int zgges_(const char *jobvsl, const char *jobvsr, const char *sort, L_fp 
    selctg, cblas_int_t *n, doublecomplex *a, cblas_int_t *lda, doublecomplex *b, 
    cblas_int_t *ldb, cblas_int_t *sdim, doublecomplex *alpha, doublecomplex *
    beta, doublecomplex *vsl, cblas_int_t *ldvsl, doublecomplex *vsr, cblas_int_t 
    *ldvsr, doublecomplex *work, cblas_int_t *lwork, double *rwork, 
    logical *bwork, cblas_int_t *info);

/* Subroutine */ int zggesx_(const char *jobvsl, const char *jobvsr, const char *sort, L_fp 
    selctg, const char *sense, cblas_int_t *n, doublecomplex *a, cblas_int_t *lda, 
    doublecomplex *b, cblas_int_t *ldb, cblas_int_t *sdim, doublecomplex *alpha, 
    doublecomplex *beta, doublecomplex *vsl, cblas_int_t *ldvsl, 
    doublecomplex *vsr, cblas_int_t *ldvsr, double *rconde, double *
    rcondv, doublecomplex *work, cblas_int_t *lwork, double *rwork, 
    cblas_int_t *iwork, cblas_int_t *liwork, logical *bwork, cblas_int_t *info);

/* Subroutine */ int zggev_(const char *jobvl, const char *jobvr, cblas_int_t *n, 
    doublecomplex *a, cblas_int_t *lda, doublecomplex *b, cblas_int_t *ldb, 
    doublecomplex *alpha, doublecomplex *beta, doublecomplex *vl, cblas_int_t 
    *ldvl, doublecomplex *vr, cblas_int_t *ldvr, doublecomplex *work, cblas_int_t
    *lwork, double *rwork, cblas_int_t *info);

/* Subroutine */ int zggevx_(const char *balanc, const char *jobvl, const char *jobvr, const char *
    sense, cblas_int_t *n, doublecomplex *a, cblas_int_t *lda, doublecomplex *b, 
    cblas_int_t *ldb, doublecomplex *alpha, doublecomplex *beta, 
    doublecomplex *vl, cblas_int_t *ldvl, doublecomplex *vr, cblas_int_t *ldvr, 
    cblas_int_t *ilo, cblas_int_t *ihi, double *lscale, double *rscale, 
    double *abnrm, double *bbnrm, double *rconde, double *
    rcondv, doublecomplex *work, cblas_int_t *lwork, double *rwork, 
    cblas_int_t *iwork, logical *bwork, cblas_int_t *info);

/* Subroutine */ int zggglm_(cblas_int_t *n, cblas_int_t *m, cblas_int_t *p, 
    doublecomplex *a, cblas_int_t *lda, doublecomplex *b, cblas_int_t *ldb, 
    doublecomplex *d__, doublecomplex *x, doublecomplex *y, doublecomplex 
    *work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int zgghrd_(const char *compq, const char *compz, cblas_int_t *n, cblas_int_t *
    ilo, cblas_int_t *ihi, doublecomplex *a, cblas_int_t *lda, doublecomplex *b, 
    cblas_int_t *ldb, doublecomplex *q, cblas_int_t *ldq, doublecomplex *z__, 
    cblas_int_t *ldz, cblas_int_t *info);

/* Subroutine */ int zgglse_(cblas_int_t *m, cblas_int_t *n, cblas_int_t *p, 
    doublecomplex *a, cblas_int_t *lda, doublecomplex *b, cblas_int_t *ldb, 
    doublecomplex *c__, doublecomplex *d__, doublecomplex *x, 
    doublecomplex *work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int zggqrf_(cblas_int_t *n, cblas_int_t *m, cblas_int_t *p, 
    doublecomplex *a, cblas_int_t *lda, doublecomplex *taua, doublecomplex *b,
     cblas_int_t *ldb, doublecomplex *taub, doublecomplex *work, cblas_int_t *
    lwork, cblas_int_t *info);

/* Subroutine */ int zggrqf_(cblas_int_t *m, cblas_int_t *p, cblas_int_t *n, 
    doublecomplex *a, cblas_int_t *lda, doublecomplex *taua, doublecomplex *b,
     cblas_int_t *ldb, doublecomplex *taub, doublecomplex *work, cblas_int_t *
    lwork, cblas_int_t *info);

/* Subroutine */ int zggsvd_(const char *jobu, const char *jobv, const char *jobq, cblas_int_t *m, 
    cblas_int_t *n, cblas_int_t *p, cblas_int_t *k, cblas_int_t *l, doublecomplex *a, 
    cblas_int_t *lda, doublecomplex *b, cblas_int_t *ldb, double *alpha, 
    double *beta, doublecomplex *u, cblas_int_t *ldu, doublecomplex *v, 
    cblas_int_t *ldv, doublecomplex *q, cblas_int_t *ldq, doublecomplex *work, 
    double *rwork, cblas_int_t *iwork, cblas_int_t *info);

/* Subroutine */ int zggsvp_(const char *jobu, const char *jobv, const char *jobq, cblas_int_t *m, 
    cblas_int_t *p, cblas_int_t *n, doublecomplex *a, cblas_int_t *lda, doublecomplex 
    *b, cblas_int_t *ldb, double *tola, double *tolb, cblas_int_t *k, 
    cblas_int_t *l, doublecomplex *u, cblas_int_t *ldu, doublecomplex *v, cblas_int_t 
    *ldv, doublecomplex *q, cblas_int_t *ldq, cblas_int_t *iwork, double *
    rwork, doublecomplex *tau, doublecomplex *work, cblas_int_t *info);

/* Subroutine */ int zgtcon_(const char *norm, cblas_int_t *n, doublecomplex *dl, 
    doublecomplex *d__, doublecomplex *du, doublecomplex *du2, cblas_int_t *
    ipiv, double *anorm, double *rcond, doublecomplex *work, 
    cblas_int_t *info);

/* Subroutine */ int zgtrfs_(const char *trans, cblas_int_t *n, cblas_int_t *nrhs, 
    doublecomplex *dl, doublecomplex *d__, doublecomplex *du, 
    doublecomplex *dlf, doublecomplex *df, doublecomplex *duf, 
    doublecomplex *du2, cblas_int_t *ipiv, doublecomplex *b, cblas_int_t *ldb, 
    doublecomplex *x, cblas_int_t *ldx, double *ferr, double *berr, 
    doublecomplex *work, double *rwork, cblas_int_t *info);

/* Subroutine */ int zgtsv_(cblas_int_t *n, cblas_int_t *nrhs, doublecomplex *dl, 
    doublecomplex *d__, doublecomplex *du, doublecomplex *b, cblas_int_t *ldb,
     cblas_int_t *info);

/* Subroutine */ int zgtsvx_(const char *fact, const char *trans, cblas_int_t *n, cblas_int_t *
    nrhs, doublecomplex *dl, doublecomplex *d__, doublecomplex *du, 
    doublecomplex *dlf, doublecomplex *df, doublecomplex *duf, 
    doublecomplex *du2, cblas_int_t *ipiv, doublecomplex *b, cblas_int_t *ldb, 
    doublecomplex *x, cblas_int_t *ldx, double *rcond, double *ferr, 
    double *berr, doublecomplex *work, double *rwork, cblas_int_t *
    info);

/* Subroutine */ int zgttrf_(cblas_int_t *n, doublecomplex *dl, doublecomplex *
    d__, doublecomplex *du, doublecomplex *du2, cblas_int_t *ipiv, cblas_int_t *
    info);

/* Subroutine */ int zgttrs_(const char *trans, cblas_int_t *n, cblas_int_t *nrhs, 
    doublecomplex *dl, doublecomplex *d__, doublecomplex *du, 
    doublecomplex *du2, cblas_int_t *ipiv, doublecomplex *b, cblas_int_t *ldb, 
    cblas_int_t *info);

/* Subroutine */ int zgtts2_(cblas_int_t *itrans, cblas_int_t *n, cblas_int_t *nrhs, 
    doublecomplex *dl, doublecomplex *d__, doublecomplex *du, 
    doublecomplex *du2, cblas_int_t *ipiv, doublecomplex *b, cblas_int_t *ldb);

/* Subroutine */ int zhbev_(const char *jobz, const char *uplo, cblas_int_t *n, cblas_int_t *kd, 
    doublecomplex *ab, cblas_int_t *ldab, double *w, doublecomplex *z__, 
    cblas_int_t *ldz, doublecomplex *work, double *rwork, cblas_int_t *info);

/* Subroutine */ int zhbevd_(const char *jobz, const char *uplo, cblas_int_t *n, cblas_int_t *kd, 
    doublecomplex *ab, cblas_int_t *ldab, double *w, doublecomplex *z__, 
    cblas_int_t *ldz, doublecomplex *work, cblas_int_t *lwork, double *rwork, 
    cblas_int_t *lrwork, cblas_int_t *iwork, cblas_int_t *liwork, cblas_int_t *info);

/* Subroutine */ int zhbevx_(const char *jobz, const char *range, const char *uplo, cblas_int_t *n, 
    cblas_int_t *kd, doublecomplex *ab, cblas_int_t *ldab, doublecomplex *q, 
    cblas_int_t *ldq, double *vl, double *vu, cblas_int_t *il, cblas_int_t *
    iu, double *abstol, cblas_int_t *m, double *w, doublecomplex *z__,
     cblas_int_t *ldz, doublecomplex *work, double *rwork, cblas_int_t *iwork,
     cblas_int_t *ifail, cblas_int_t *info);

/* Subroutine */ int zhbgst_(const char *vect, const char *uplo, cblas_int_t *n, cblas_int_t *ka, 
    cblas_int_t *kb, doublecomplex *ab, cblas_int_t *ldab, doublecomplex *bb, 
    cblas_int_t *ldbb, doublecomplex *x, cblas_int_t *ldx, doublecomplex *work, 
    double *rwork, cblas_int_t *info);

/* Subroutine */ int zhbgv_(const char *jobz, const char *uplo, cblas_int_t *n, cblas_int_t *ka, 
    cblas_int_t *kb, doublecomplex *ab, cblas_int_t *ldab, doublecomplex *bb, 
    cblas_int_t *ldbb, double *w, doublecomplex *z__, cblas_int_t *ldz, 
    doublecomplex *work, double *rwork, cblas_int_t *info);

/* Subroutine */ int zhbgvd_(const char *jobz, const char *uplo, cblas_int_t *n, cblas_int_t *ka, 
    cblas_int_t *kb, doublecomplex *ab, cblas_int_t *ldab, doublecomplex *bb, 
    cblas_int_t *ldbb, double *w, doublecomplex *z__, cblas_int_t *ldz, 
    doublecomplex *work, cblas_int_t *lwork, double *rwork, cblas_int_t *
    lrwork, cblas_int_t *iwork, cblas_int_t *liwork, cblas_int_t *info);

/* Subroutine */ int zhbgvx_(const char *jobz, const char *range, const char *uplo, cblas_int_t *n, 
    cblas_int_t *ka, cblas_int_t *kb, doublecomplex *ab, cblas_int_t *ldab, 
    doublecomplex *bb, cblas_int_t *ldbb, doublecomplex *q, cblas_int_t *ldq, 
    double *vl, double *vu, cblas_int_t *il, cblas_int_t *iu, double *
    abstol, cblas_int_t *m, double *w, doublecomplex *z__, cblas_int_t *ldz, 
    doublecomplex *work, double *rwork, cblas_int_t *iwork, cblas_int_t *
    ifail, cblas_int_t *info);

/* Subroutine */ int zhbtrd_(const char *vect, const char *uplo, cblas_int_t *n, cblas_int_t *kd, 
    doublecomplex *ab, cblas_int_t *ldab, double *d__, double *e, 
    doublecomplex *q, cblas_int_t *ldq, doublecomplex *work, cblas_int_t *info);

/* Subroutine */ int zhecon_(const char *uplo, cblas_int_t *n, doublecomplex *a, 
    cblas_int_t *lda, cblas_int_t *ipiv, double *anorm, double *rcond, 
    doublecomplex *work, cblas_int_t *info);

/* Subroutine */ int zheev_(const char *jobz, const char *uplo, cblas_int_t *n, doublecomplex 
    *a, cblas_int_t *lda, double *w, doublecomplex *work, cblas_int_t *lwork, 
    double *rwork, cblas_int_t *info);

/* Subroutine */ int zheevd_(const char *jobz, const char *uplo, cblas_int_t *n, 
    doublecomplex *a, cblas_int_t *lda, double *w, doublecomplex *work, 
    cblas_int_t *lwork, double *rwork, cblas_int_t *lrwork, cblas_int_t *iwork, 
    cblas_int_t *liwork, cblas_int_t *info);

/* Subroutine */ int zheevr_(const char *jobz, const char *range, const char *uplo, cblas_int_t *n, 
    doublecomplex *a, cblas_int_t *lda, double *vl, double *vu, 
    cblas_int_t *il, cblas_int_t *iu, double *abstol, cblas_int_t *m, double *
    w, doublecomplex *z__, cblas_int_t *ldz, cblas_int_t *isuppz, doublecomplex *
    work, cblas_int_t *lwork, double *rwork, cblas_int_t *lrwork, cblas_int_t *
    iwork, cblas_int_t *liwork, cblas_int_t *info);

/* Subroutine */ int zheevx_(const char *jobz, const char *range, const char *uplo, cblas_int_t *n, 
    doublecomplex *a, cblas_int_t *lda, double *vl, double *vu, 
    cblas_int_t *il, cblas_int_t *iu, double *abstol, cblas_int_t *m, double *
    w, doublecomplex *z__, cblas_int_t *ldz, doublecomplex *work, cblas_int_t *
    lwork, double *rwork, cblas_int_t *iwork, cblas_int_t *ifail, cblas_int_t *
    info);

/* Subroutine */ int zhegs2_(cblas_int_t *itype, const char *uplo, cblas_int_t *n, 
    doublecomplex *a, cblas_int_t *lda, doublecomplex *b, cblas_int_t *ldb, 
    cblas_int_t *info);

/* Subroutine */ int zhegst_(cblas_int_t *itype, const char *uplo, cblas_int_t *n, 
    doublecomplex *a, cblas_int_t *lda, doublecomplex *b, cblas_int_t *ldb, 
    cblas_int_t *info);

/* Subroutine */ int zhegv_(cblas_int_t *itype, const char *jobz, const char *uplo, cblas_int_t *
    n, doublecomplex *a, cblas_int_t *lda, doublecomplex *b, cblas_int_t *ldb, 
    double *w, doublecomplex *work, cblas_int_t *lwork, double *rwork,
     cblas_int_t *info);

/* Subroutine */ int zhegvd_(cblas_int_t *itype, const char *jobz, const char *uplo, cblas_int_t *
    n, doublecomplex *a, cblas_int_t *lda, doublecomplex *b, cblas_int_t *ldb, 
    double *w, doublecomplex *work, cblas_int_t *lwork, double *rwork,
     cblas_int_t *lrwork, cblas_int_t *iwork, cblas_int_t *liwork, cblas_int_t *info);

/* Subroutine */ int zhegvx_(cblas_int_t *itype, const char *jobz, const char *range, const char *
    uplo, cblas_int_t *n, doublecomplex *a, cblas_int_t *lda, doublecomplex *b, 
    cblas_int_t *ldb, double *vl, double *vu, cblas_int_t *il, cblas_int_t *
    iu, double *abstol, cblas_int_t *m, double *w, doublecomplex *z__,
     cblas_int_t *ldz, doublecomplex *work, cblas_int_t *lwork, double *rwork,
     cblas_int_t *iwork, cblas_int_t *ifail, cblas_int_t *info);

/* Subroutine */ int zherfs_(const char *uplo, cblas_int_t *n, cblas_int_t *nrhs, 
    doublecomplex *a, cblas_int_t *lda, doublecomplex *af, cblas_int_t *ldaf, 
    cblas_int_t *ipiv, doublecomplex *b, cblas_int_t *ldb, doublecomplex *x, 
    cblas_int_t *ldx, double *ferr, double *berr, doublecomplex *work,
     double *rwork, cblas_int_t *info);

/* Subroutine */ int zhesv_(const char *uplo, cblas_int_t *n, cblas_int_t *nrhs, 
    doublecomplex *a, cblas_int_t *lda, cblas_int_t *ipiv, doublecomplex *b, 
    cblas_int_t *ldb, doublecomplex *work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int zhesvx_(const char *fact, const char *uplo, cblas_int_t *n, cblas_int_t *
    nrhs, doublecomplex *a, cblas_int_t *lda, doublecomplex *af, cblas_int_t *
    ldaf, cblas_int_t *ipiv, doublecomplex *b, cblas_int_t *ldb, doublecomplex *x,
     cblas_int_t *ldx, double *rcond, double *ferr, double *berr, 
    doublecomplex *work, cblas_int_t *lwork, double *rwork, cblas_int_t *info);

/* Subroutine */ int zhetd2_(const char *uplo, cblas_int_t *n, doublecomplex *a, 
    cblas_int_t *lda, double *d__, double *e, doublecomplex *tau, 
    cblas_int_t *info);

/* Subroutine */ int zhetf2_(const char *uplo, cblas_int_t *n, doublecomplex *a, 
    cblas_int_t *lda, cblas_int_t *ipiv, cblas_int_t *info);

/* Subroutine */ int zhetrd_(const char *uplo, cblas_int_t *n, doublecomplex *a, 
    cblas_int_t *lda, double *d__, double *e, doublecomplex *tau, 
    doublecomplex *work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int zhetrf_(const char *uplo, cblas_int_t *n, doublecomplex *a, 
    cblas_int_t *lda, cblas_int_t *ipiv, doublecomplex *work, cblas_int_t *lwork, 
    cblas_int_t *info);

/* Subroutine */ int zhetri_(const char *uplo, cblas_int_t *n, doublecomplex *a, 
    cblas_int_t *lda, cblas_int_t *ipiv, doublecomplex *work, cblas_int_t *info);

/* Subroutine */ int zhetrs_(const char *uplo, cblas_int_t *n, cblas_int_t *nrhs, 
    doublecomplex *a, cblas_int_t *lda, cblas_int_t *ipiv, doublecomplex *b, 
    cblas_int_t *ldb, cblas_int_t *info);

/* Subroutine */ int zhgeqz_(const char *job, const char *compq, const char *compz, cblas_int_t *n, 
    cblas_int_t *ilo, cblas_int_t *ihi, doublecomplex *h__, cblas_int_t *ldh, 
    doublecomplex *t, cblas_int_t *ldt, doublecomplex *alpha, doublecomplex *
    beta, doublecomplex *q, cblas_int_t *ldq, doublecomplex *z__, cblas_int_t *
    ldz, doublecomplex *work, cblas_int_t *lwork, double *rwork, cblas_int_t *
    info);

/* Subroutine */ int zhpcon_(const char *uplo, cblas_int_t *n, doublecomplex *ap, 
    cblas_int_t *ipiv, double *anorm, double *rcond, doublecomplex *
    work, cblas_int_t *info);

/* Subroutine */ int zhpev_(const char *jobz, const char *uplo, cblas_int_t *n, doublecomplex 
    *ap, double *w, doublecomplex *z__, cblas_int_t *ldz, doublecomplex *
    work, double *rwork, cblas_int_t *info);

/* Subroutine */ int zhpevd_(const char *jobz, const char *uplo, cblas_int_t *n, 
    doublecomplex *ap, double *w, doublecomplex *z__, cblas_int_t *ldz, 
    doublecomplex *work, cblas_int_t *lwork, double *rwork, cblas_int_t *
    lrwork, cblas_int_t *iwork, cblas_int_t *liwork, cblas_int_t *info);

/* Subroutine */ int zhpevx_(const char *jobz, const char *range, const char *uplo, cblas_int_t *n, 
    doublecomplex *ap, double *vl, double *vu, cblas_int_t *il, 
    cblas_int_t *iu, double *abstol, cblas_int_t *m, double *w, 
    doublecomplex *z__, cblas_int_t *ldz, doublecomplex *work, double *
    rwork, cblas_int_t *iwork, cblas_int_t *ifail, cblas_int_t *info);

/* Subroutine */ int zhpgst_(cblas_int_t *itype, const char *uplo, cblas_int_t *n, 
    doublecomplex *ap, doublecomplex *bp, cblas_int_t *info);

/* Subroutine */ int zhpgv_(cblas_int_t *itype, const char *jobz, const char *uplo, cblas_int_t *
    n, doublecomplex *ap, doublecomplex *bp, double *w, doublecomplex 
    *z__, cblas_int_t *ldz, doublecomplex *work, double *rwork, cblas_int_t *
    info);

/* Subroutine */ int zhpgvd_(cblas_int_t *itype, const char *jobz, const char *uplo, cblas_int_t *
    n, doublecomplex *ap, doublecomplex *bp, double *w, doublecomplex 
    *z__, cblas_int_t *ldz, doublecomplex *work, cblas_int_t *lwork, double *
    rwork, cblas_int_t *lrwork, cblas_int_t *iwork, cblas_int_t *liwork, cblas_int_t *
    info);

/* Subroutine */ int zhpgvx_(cblas_int_t *itype, const char *jobz, const char *range, const char *
    uplo, cblas_int_t *n, doublecomplex *ap, doublecomplex *bp, double *
    vl, double *vu, cblas_int_t *il, cblas_int_t *iu, double *abstol, 
    cblas_int_t *m, double *w, doublecomplex *z__, cblas_int_t *ldz, 
    doublecomplex *work, double *rwork, cblas_int_t *iwork, cblas_int_t *
    ifail, cblas_int_t *info);

/* Subroutine */ int zhprfs_(const char *uplo, cblas_int_t *n, cblas_int_t *nrhs, 
    doublecomplex *ap, doublecomplex *afp, cblas_int_t *ipiv, doublecomplex *
    b, cblas_int_t *ldb, doublecomplex *x, cblas_int_t *ldx, double *ferr, 
    double *berr, doublecomplex *work, double *rwork, cblas_int_t *
    info);

/* Subroutine */ int zhpsv_(const char *uplo, cblas_int_t *n, cblas_int_t *nrhs, 
    doublecomplex *ap, cblas_int_t *ipiv, doublecomplex *b, cblas_int_t *ldb, 
    cblas_int_t *info);

/* Subroutine */ int zhpsvx_(const char *fact, const char *uplo, cblas_int_t *n, cblas_int_t *
    nrhs, doublecomplex *ap, doublecomplex *afp, cblas_int_t *ipiv, 
    doublecomplex *b, cblas_int_t *ldb, doublecomplex *x, cblas_int_t *ldx, 
    double *rcond, double *ferr, double *berr, doublecomplex *
    work, double *rwork, cblas_int_t *info);

/* Subroutine */ int zhptrd_(const char *uplo, cblas_int_t *n, doublecomplex *ap, 
    double *d__, double *e, doublecomplex *tau, cblas_int_t *info);

/* Subroutine */ int zhptrf_(const char *uplo, cblas_int_t *n, doublecomplex *ap, 
    cblas_int_t *ipiv, cblas_int_t *info);

/* Subroutine */ int zhptri_(const char *uplo, cblas_int_t *n, doublecomplex *ap, 
    cblas_int_t *ipiv, doublecomplex *work, cblas_int_t *info);

/* Subroutine */ int zhptrs_(const char *uplo, cblas_int_t *n, cblas_int_t *nrhs, 
    doublecomplex *ap, cblas_int_t *ipiv, doublecomplex *b, cblas_int_t *ldb, 
    cblas_int_t *info);

/* Subroutine */ int zhsein_(const char *side, const char *eigsrc, const char *initv, logical *
    select, cblas_int_t *n, doublecomplex *h__, cblas_int_t *ldh, doublecomplex *
    w, doublecomplex *vl, cblas_int_t *ldvl, doublecomplex *vr, cblas_int_t *ldvr,
     cblas_int_t *mm, cblas_int_t *m, doublecomplex *work, double *rwork, 
    cblas_int_t *ifaill, cblas_int_t *ifailr, cblas_int_t *info);

/* Subroutine */ int zhseqr_(const char *job, const char *compz, cblas_int_t *n, cblas_int_t *ilo,
     cblas_int_t *ihi, doublecomplex *h__, cblas_int_t *ldh, doublecomplex *w, 
    doublecomplex *z__, cblas_int_t *ldz, doublecomplex *work, cblas_int_t *lwork,
     cblas_int_t *info);

/* Subroutine */ int zlabrd_(cblas_int_t *m, cblas_int_t *n, cblas_int_t *nb, 
    doublecomplex *a, cblas_int_t *lda, double *d__, double *e, 
    doublecomplex *tauq, doublecomplex *taup, doublecomplex *x, cblas_int_t *
    ldx, doublecomplex *y, cblas_int_t *ldy);

/* Subroutine */ int zlacgv_(cblas_int_t *n, doublecomplex *x, cblas_int_t *incx);

/* Subroutine */ int zlacn2_(cblas_int_t *n, doublecomplex *v, doublecomplex *x, 
    double *est, cblas_int_t *kase, cblas_int_t *isave);

/* Subroutine */ int zlacon_(cblas_int_t *n, doublecomplex *v, doublecomplex *x, 
    double *est, cblas_int_t *kase);

/* Subroutine */ int zlacp2_(const char *uplo, cblas_int_t *m, cblas_int_t *n, double *
    a, cblas_int_t *lda, doublecomplex *b, cblas_int_t *ldb);

/* Subroutine */ int zlacpy_(const char *uplo, cblas_int_t *m, cblas_int_t *n, 
    doublecomplex *a, cblas_int_t *lda, doublecomplex *b, cblas_int_t *ldb);

/* Subroutine */ int zlacrm_(cblas_int_t *m, cblas_int_t *n, doublecomplex *a, 
    cblas_int_t *lda, double *b, cblas_int_t *ldb, doublecomplex *c__, 
    cblas_int_t *ldc, double *rwork);

/* Subroutine */ int zlacrt_(cblas_int_t *n, doublecomplex *cx, cblas_int_t *incx, 
    doublecomplex *cy, cblas_int_t *incy, doublecomplex *c__, doublecomplex *
    s);

/* Double Complex */ VOID zladiv_(doublecomplex * ret_val, doublecomplex *x, 
    doublecomplex *y);

/* Subroutine */ int zlaed0_(cblas_int_t *qsiz, cblas_int_t *n, double *d__, 
    double *e, doublecomplex *q, cblas_int_t *ldq, doublecomplex *qstore, 
    cblas_int_t *ldqs, double *rwork, cblas_int_t *iwork, cblas_int_t *info);

/* Subroutine */ int zlaed7_(cblas_int_t *n, cblas_int_t *cutpnt, cblas_int_t *qsiz, 
    cblas_int_t *tlvls, cblas_int_t *curlvl, cblas_int_t *curpbm, double *d__, 
    doublecomplex *q, cblas_int_t *ldq, double *rho, cblas_int_t *indxq, 
    double *qstore, cblas_int_t *qptr, cblas_int_t *prmptr, cblas_int_t *perm, 
    cblas_int_t *givptr, cblas_int_t *givcol, double *givnum, doublecomplex *
    work, double *rwork, cblas_int_t *iwork, cblas_int_t *info);

/* Subroutine */ int zlaed8_(cblas_int_t *k, cblas_int_t *n, cblas_int_t *qsiz, 
    doublecomplex *q, cblas_int_t *ldq, double *d__, double *rho, 
    cblas_int_t *cutpnt, double *z__, double *dlamda, doublecomplex *
    q2, cblas_int_t *ldq2, double *w, cblas_int_t *indxp, cblas_int_t *indx, 
    cblas_int_t *indxq, cblas_int_t *perm, cblas_int_t *givptr, cblas_int_t *givcol, 
    double *givnum, cblas_int_t *info);

/* Subroutine */ int zlaein_(logical *rightv, logical *noinit, cblas_int_t *n, 
    doublecomplex *h__, cblas_int_t *ldh, doublecomplex *w, doublecomplex *v, 
    doublecomplex *b, cblas_int_t *ldb, double *rwork, double *eps3, 
    double *smlnum, cblas_int_t *info);

/* Subroutine */ int zlaesy_(doublecomplex *a, doublecomplex *b, 
    doublecomplex *c__, doublecomplex *rt1, doublecomplex *rt2, 
    doublecomplex *evscal, doublecomplex *cs1, doublecomplex *sn1);

/* Subroutine */ int zlaev2_(doublecomplex *a, doublecomplex *b, 
    doublecomplex *c__, double *rt1, double *rt2, double *cs1,
     doublecomplex *sn1);

/* Subroutine */ int zlag2c_(cblas_int_t *m, cblas_int_t *n, doublecomplex *a, 
    cblas_int_t *lda, complex *sa, cblas_int_t *ldsa, cblas_int_t *info);

/* Subroutine */ int zlags2_(logical *upper, double *a1, doublecomplex *
    a2, double *a3, double *b1, doublecomplex *b2, double *b3,
     double *csu, doublecomplex *snu, double *csv, doublecomplex *
    snv, double *csq, doublecomplex *snq);

/* Subroutine */ int zlagtm_(const char *trans, cblas_int_t *n, cblas_int_t *nrhs, 
    double *alpha, doublecomplex *dl, doublecomplex *d__, 
    doublecomplex *du, doublecomplex *x, cblas_int_t *ldx, double *beta, 
    doublecomplex *b, cblas_int_t *ldb);

/* Subroutine */ int zlahef_(const char *uplo, cblas_int_t *n, cblas_int_t *nb, cblas_int_t *kb,
     doublecomplex *a, cblas_int_t *lda, cblas_int_t *ipiv, doublecomplex *w, 
    cblas_int_t *ldw, cblas_int_t *info);

/* Subroutine */ int zlahqr_(logical *wantt, logical *wantz, cblas_int_t *n, 
    cblas_int_t *ilo, cblas_int_t *ihi, doublecomplex *h__, cblas_int_t *ldh, 
    doublecomplex *w, cblas_int_t *iloz, cblas_int_t *ihiz, doublecomplex *z__, 
    cblas_int_t *ldz, cblas_int_t *info);

/* Subroutine */ int zlahr2_(cblas_int_t *n, cblas_int_t *k, cblas_int_t *nb, 
    doublecomplex *a, cblas_int_t *lda, doublecomplex *tau, doublecomplex *t, 
    cblas_int_t *ldt, doublecomplex *y, cblas_int_t *ldy);

/* Subroutine */ int zlahrd_(cblas_int_t *n, cblas_int_t *k, cblas_int_t *nb, 
    doublecomplex *a, cblas_int_t *lda, doublecomplex *tau, doublecomplex *t, 
    cblas_int_t *ldt, doublecomplex *y, cblas_int_t *ldy);

/* Subroutine */ int zlaic1_(cblas_int_t *job, cblas_int_t *j, doublecomplex *x, 
    double *sest, doublecomplex *w, doublecomplex *gamma, double *
    sestpr, doublecomplex *s, doublecomplex *c__);

/* Subroutine */ int zlals0_(cblas_int_t *icompq, cblas_int_t *nl, cblas_int_t *nr, 
    cblas_int_t *sqre, cblas_int_t *nrhs, doublecomplex *b, cblas_int_t *ldb, 
    doublecomplex *bx, cblas_int_t *ldbx, cblas_int_t *perm, cblas_int_t *givptr, 
    cblas_int_t *givcol, cblas_int_t *ldgcol, double *givnum, cblas_int_t *ldgnum,
     double *poles, double *difl, double *difr, double *
    z__, cblas_int_t *k, double *c__, double *s, double *rwork, 
    cblas_int_t *info);

/* Subroutine */ int zlalsa_(cblas_int_t *icompq, cblas_int_t *smlsiz, cblas_int_t *n, 
    cblas_int_t *nrhs, doublecomplex *b, cblas_int_t *ldb, doublecomplex *bx, 
    cblas_int_t *ldbx, double *u, cblas_int_t *ldu, double *vt, cblas_int_t *
    k, double *difl, double *difr, double *z__, double *
    poles, cblas_int_t *givptr, cblas_int_t *givcol, cblas_int_t *ldgcol, cblas_int_t *
    perm, double *givnum, double *c__, double *s, double *
    rwork, cblas_int_t *iwork, cblas_int_t *info);

/* Subroutine */ int zlalsd_(const char *uplo, cblas_int_t *smlsiz, cblas_int_t *n, cblas_int_t 
    *nrhs, double *d__, double *e, doublecomplex *b, cblas_int_t *ldb,
     double *rcond, cblas_int_t *rank, doublecomplex *work, double *
    rwork, cblas_int_t *iwork, cblas_int_t *info);

/* Subroutine */ int zlapll_(cblas_int_t *n, doublecomplex *x, cblas_int_t *incx, 
    doublecomplex *y, cblas_int_t *incy, double *ssmin);

/* Subroutine */ int zlapmt_(logical *forwrd, cblas_int_t *m, cblas_int_t *n, 
    doublecomplex *x, cblas_int_t *ldx, cblas_int_t *k);

/* Subroutine */ int zlaqgb_(cblas_int_t *m, cblas_int_t *n, cblas_int_t *kl, cblas_int_t *ku,
     doublecomplex *ab, cblas_int_t *ldab, double *r__, double *c__, 
    double *rowcnd, double *colcnd, double *amax, const char *equed);

/* Subroutine */ int zlaqge_(cblas_int_t *m, cblas_int_t *n, doublecomplex *a, 
    cblas_int_t *lda, double *r__, double *c__, double *rowcnd, 
    double *colcnd, double *amax, const char *equed);

/* Subroutine */ int zlaqhb_(const char *uplo, cblas_int_t *n, cblas_int_t *kd, 
    doublecomplex *ab, cblas_int_t *ldab, double *s, double *scond, 
    double *amax, const char *equed);

/* Subroutine */ int zlaqhe_(const char *uplo, cblas_int_t *n, doublecomplex *a, 
    cblas_int_t *lda, double *s, double *scond, double *amax, 
    const char *equed);

/* Subroutine */ int zlaqhp_(const char *uplo, cblas_int_t *n, doublecomplex *ap, 
    double *s, double *scond, double *amax, const char *equed);

/* Subroutine */ int zlaqp2_(cblas_int_t *m, cblas_int_t *n, cblas_int_t *offset, 
    doublecomplex *a, cblas_int_t *lda, cblas_int_t *jpvt, doublecomplex *tau, 
    double *vn1, double *vn2, doublecomplex *work);

/* Subroutine */ int zlaqps_(cblas_int_t *m, cblas_int_t *n, cblas_int_t *offset, cblas_int_t 
    *nb, cblas_int_t *kb, doublecomplex *a, cblas_int_t *lda, cblas_int_t *jpvt, 
    doublecomplex *tau, double *vn1, double *vn2, doublecomplex *
    auxv, doublecomplex *f, cblas_int_t *ldf);

/* Subroutine */ int zlaqr0_(logical *wantt, logical *wantz, cblas_int_t *n, 
    cblas_int_t *ilo, cblas_int_t *ihi, doublecomplex *h__, cblas_int_t *ldh, 
    doublecomplex *w, cblas_int_t *iloz, cblas_int_t *ihiz, doublecomplex *z__, 
    cblas_int_t *ldz, doublecomplex *work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int zlaqr1_(cblas_int_t *n, doublecomplex *h__, cblas_int_t *ldh, 
    doublecomplex *s1, doublecomplex *s2, doublecomplex *v);

/* Subroutine */ int zlaqr2_(logical *wantt, logical *wantz, cblas_int_t *n, 
    cblas_int_t *ktop, cblas_int_t *kbot, cblas_int_t *nw, doublecomplex *h__, 
    cblas_int_t *ldh, cblas_int_t *iloz, cblas_int_t *ihiz, doublecomplex *z__, 
    cblas_int_t *ldz, cblas_int_t *ns, cblas_int_t *nd, doublecomplex *sh, 
    doublecomplex *v, cblas_int_t *ldv, cblas_int_t *nh, doublecomplex *t, 
    cblas_int_t *ldt, cblas_int_t *nv, doublecomplex *wv, cblas_int_t *ldwv, 
    doublecomplex *work, cblas_int_t *lwork);

/* Subroutine */ int zlaqr3_(logical *wantt, logical *wantz, cblas_int_t *n, 
    cblas_int_t *ktop, cblas_int_t *kbot, cblas_int_t *nw, doublecomplex *h__, 
    cblas_int_t *ldh, cblas_int_t *iloz, cblas_int_t *ihiz, doublecomplex *z__, 
    cblas_int_t *ldz, cblas_int_t *ns, cblas_int_t *nd, doublecomplex *sh, 
    doublecomplex *v, cblas_int_t *ldv, cblas_int_t *nh, doublecomplex *t, 
    cblas_int_t *ldt, cblas_int_t *nv, doublecomplex *wv, cblas_int_t *ldwv, 
    doublecomplex *work, cblas_int_t *lwork);

/* Subroutine */ int zlaqr4_(logical *wantt, logical *wantz, cblas_int_t *n, 
    cblas_int_t *ilo, cblas_int_t *ihi, doublecomplex *h__, cblas_int_t *ldh, 
    doublecomplex *w, cblas_int_t *iloz, cblas_int_t *ihiz, doublecomplex *z__, 
    cblas_int_t *ldz, doublecomplex *work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int zlaqr5_(logical *wantt, logical *wantz, cblas_int_t *kacc22, 
    cblas_int_t *n, cblas_int_t *ktop, cblas_int_t *kbot, cblas_int_t *nshfts, 
    doublecomplex *s, doublecomplex *h__, cblas_int_t *ldh, cblas_int_t *iloz, 
    cblas_int_t *ihiz, doublecomplex *z__, cblas_int_t *ldz, doublecomplex *v, 
    cblas_int_t *ldv, doublecomplex *u, cblas_int_t *ldu, cblas_int_t *nv, 
    doublecomplex *wv, cblas_int_t *ldwv, cblas_int_t *nh, doublecomplex *wh, 
    cblas_int_t *ldwh);

/* Subroutine */ int zlaqsb_(const char *uplo, cblas_int_t *n, cblas_int_t *kd, 
    doublecomplex *ab, cblas_int_t *ldab, double *s, double *scond, 
    double *amax, const char *equed);

/* Subroutine */ int zlaqsp_(const char *uplo, cblas_int_t *n, doublecomplex *ap, 
    double *s, double *scond, double *amax, const char *equed);

/* Subroutine */ int zlaqsy_(const char *uplo, cblas_int_t *n, doublecomplex *a, 
    cblas_int_t *lda, double *s, double *scond, double *amax, 
    const char *equed);

/* Subroutine */ int zlar1v_(cblas_int_t *n, cblas_int_t *b1, cblas_int_t *bn, double 
    *lambda, double *d__, double *l, double *ld, double *
    lld, double *pivmin, double *gaptol, doublecomplex *z__, 
    logical *wantnc, cblas_int_t *negcnt, double *ztz, double *mingma,
     cblas_int_t *r__, cblas_int_t *isuppz, double *nrminv, double *resid,
     double *rqcorr, double *work);

/* Subroutine */ int zlar2v_(cblas_int_t *n, doublecomplex *x, doublecomplex *y, 
    doublecomplex *z__, cblas_int_t *incx, double *c__, doublecomplex *s, 
    cblas_int_t *incc);

/* Subroutine */ int zlarcm_(cblas_int_t *m, cblas_int_t *n, double *a, cblas_int_t *
    lda, doublecomplex *b, cblas_int_t *ldb, doublecomplex *c__, cblas_int_t *ldc,
     double *rwork);

/* Subroutine */ int zlarf_(const char *side, cblas_int_t *m, cblas_int_t *n, doublecomplex 
    *v, cblas_int_t *incv, doublecomplex *tau, doublecomplex *c__, cblas_int_t *
    ldc, doublecomplex *work);

/* Subroutine */ int zlarfb_(const char *side, const char *trans, const char *direct, const char *
    storev, cblas_int_t *m, cblas_int_t *n, cblas_int_t *k, doublecomplex *v, cblas_int_t 
    *ldv, doublecomplex *t, cblas_int_t *ldt, doublecomplex *c__, cblas_int_t *
    ldc, doublecomplex *work, cblas_int_t *ldwork   );

/* Subroutine */ int zlarfg_(cblas_int_t *n, doublecomplex *alpha, doublecomplex *
    x, cblas_int_t *incx, doublecomplex *tau);

/* Subroutine */ int zlarft_(const char *direct, const char *storev, cblas_int_t *n, cblas_int_t *
    k, doublecomplex *v, cblas_int_t *ldv, doublecomplex *tau, doublecomplex *
    t, cblas_int_t *ldt);

/* Subroutine */ int zlarfx_(const char *side, cblas_int_t *m, cblas_int_t *n, 
    doublecomplex *v, doublecomplex *tau, doublecomplex *c__, cblas_int_t *
    ldc, doublecomplex *work);

/* Subroutine */ int zlargv_(cblas_int_t *n, doublecomplex *x, cblas_int_t *incx, 
    doublecomplex *y, cblas_int_t *incy, double *c__, cblas_int_t *incc);

/* Subroutine */ int zlarnv_(cblas_int_t *idist, cblas_int_t *iseed, cblas_int_t *n, 
    doublecomplex *x);

/* Subroutine */ int zlarrv_(cblas_int_t *n, double *vl, double *vu, 
    double *d__, double *l, double *pivmin, cblas_int_t *isplit, 
    cblas_int_t *m, cblas_int_t *dol, cblas_int_t *dou, double *minrgp, 
    double *rtol1, double *rtol2, double *w, double *werr,
     double *wgap, cblas_int_t *iblock, cblas_int_t *indexw, double *gers,
     doublecomplex *z__, cblas_int_t *ldz, cblas_int_t *isuppz, double *work, 
    cblas_int_t *iwork, cblas_int_t *info);

/* Subroutine */ int zlartg_(doublecomplex *f, doublecomplex *g, double *
    cs, doublecomplex *sn, doublecomplex *r__);

/* Subroutine */ int zlartv_(cblas_int_t *n, doublecomplex *x, cblas_int_t *incx, 
    doublecomplex *y, cblas_int_t *incy, double *c__, doublecomplex *s, 
    cblas_int_t *incc);

/* Subroutine */ int zlarz_(const char *side, cblas_int_t *m, cblas_int_t *n, cblas_int_t *l, 
    doublecomplex *v, cblas_int_t *incv, doublecomplex *tau, doublecomplex *
    c__, cblas_int_t *ldc, doublecomplex *work);

/* Subroutine */ int zlarzb_(const char *side, const char *trans, const char *direct, const char *
    storev, cblas_int_t *m, cblas_int_t *n, cblas_int_t *k, cblas_int_t *l, doublecomplex 
    *v, cblas_int_t *ldv, doublecomplex *t, cblas_int_t *ldt, doublecomplex *c__, 
    cblas_int_t *ldc, doublecomplex *work, cblas_int_t *ldwork);

/* Subroutine */ int zlarzt_(const char *direct, const char *storev, cblas_int_t *n, cblas_int_t *
    k, doublecomplex *v, cblas_int_t *ldv, doublecomplex *tau, doublecomplex *
    t, cblas_int_t *ldt);

/* Subroutine */ int zlascl_(const char *type__, cblas_int_t *kl, cblas_int_t *ku, 
    double *cfrom, double *cto, cblas_int_t *m, cblas_int_t *n, 
    doublecomplex *a, cblas_int_t *lda, cblas_int_t *info);

/* Subroutine */ int zlaset_(const char *uplo, cblas_int_t *m, cblas_int_t *n, 
    doublecomplex *alpha, doublecomplex *beta, doublecomplex *a, cblas_int_t *
    lda);

/* Subroutine */ int zlasr_(const char *side, const char *pivot, const char *direct, cblas_int_t *m,
     cblas_int_t *n, double *c__, double *s, doublecomplex *a, 
    cblas_int_t *lda);

/* Subroutine */ int zlassq_(cblas_int_t *n, doublecomplex *x, cblas_int_t *incx, 
    double *scale, double *sumsq);

/* Subroutine */ int zlaswp_(cblas_int_t *n, doublecomplex *a, cblas_int_t *lda, 
    cblas_int_t *k1, cblas_int_t *k2, cblas_int_t *ipiv, cblas_int_t *incx);

/* Subroutine */ int zlasyf_(const char *uplo, cblas_int_t *n, cblas_int_t *nb, cblas_int_t *kb,
     doublecomplex *a, cblas_int_t *lda, cblas_int_t *ipiv, doublecomplex *w, 
    cblas_int_t *ldw, cblas_int_t *info);

/* Subroutine */ int zlatbs_(const char *uplo, const char *trans, const char *diag, const char *
    normin, cblas_int_t *n, cblas_int_t *kd, doublecomplex *ab, cblas_int_t *ldab, 
    doublecomplex *x, double *scale, double *cnorm, cblas_int_t *info);

/* Subroutine */ int zlatdf_(cblas_int_t *ijob, cblas_int_t *n, doublecomplex *z__, 
    cblas_int_t *ldz, doublecomplex *rhs, double *rdsum, double *
    rdscal, cblas_int_t *ipiv, cblas_int_t *jpiv);

/* Subroutine */ int zlatps_(const char *uplo, const char *trans, const char *diag, const char *
    normin, cblas_int_t *n, doublecomplex *ap, doublecomplex *x, double *
    scale, double *cnorm, cblas_int_t *info);

/* Subroutine */ int zlatrd_(const char *uplo, cblas_int_t *n, cblas_int_t *nb, 
    doublecomplex *a, cblas_int_t *lda, double *e, doublecomplex *tau, 
    doublecomplex *w, cblas_int_t *ldw);

/* Subroutine */ int zlatrs_(const char *uplo, const char *trans, const char *diag, const char *
    normin, cblas_int_t *n, doublecomplex *a, cblas_int_t *lda, doublecomplex *x, 
    double *scale, double *cnorm, cblas_int_t *info);

/* Subroutine */ int zlatrz_(cblas_int_t *m, cblas_int_t *n, cblas_int_t *l, 
    doublecomplex *a, cblas_int_t *lda, doublecomplex *tau, doublecomplex *
    work);

/* Subroutine */ int zlatzm_(const char *side, cblas_int_t *m, cblas_int_t *n, 
    doublecomplex *v, cblas_int_t *incv, doublecomplex *tau, doublecomplex *
    c1, doublecomplex *c2, cblas_int_t *ldc, doublecomplex *work    );

/* Subroutine */ int zlauu2_(const char *uplo, cblas_int_t *n, doublecomplex *a, 
    cblas_int_t *lda, cblas_int_t *info);

/* Subroutine */ int zlauum_(const char *uplo, cblas_int_t *n, doublecomplex *a, 
    cblas_int_t *lda, cblas_int_t *info);

/* Subroutine */ int zpbcon_(const char *uplo, cblas_int_t *n, cblas_int_t *kd, 
    doublecomplex *ab, cblas_int_t *ldab, double *anorm, double *
    rcond, doublecomplex *work, double *rwork, cblas_int_t *info);

/* Subroutine */ int zpbequ_(const char *uplo, cblas_int_t *n, cblas_int_t *kd, 
    doublecomplex *ab, cblas_int_t *ldab, double *s, double *scond, 
    double *amax, cblas_int_t *info);

/* Subroutine */ int zpbrfs_(const char *uplo, cblas_int_t *n, cblas_int_t *kd, cblas_int_t *
    nrhs, doublecomplex *ab, cblas_int_t *ldab, doublecomplex *afb, cblas_int_t *
    ldafb, doublecomplex *b, cblas_int_t *ldb, doublecomplex *x, cblas_int_t *ldx,
     double *ferr, double *berr, doublecomplex *work, double *
    rwork, cblas_int_t *info);

/* Subroutine */ int zpbstf_(const char *uplo, cblas_int_t *n, cblas_int_t *kd, 
    doublecomplex *ab, cblas_int_t *ldab, cblas_int_t *info);

/* Subroutine */ int zpbsv_(const char *uplo, cblas_int_t *n, cblas_int_t *kd, cblas_int_t *
    nrhs, doublecomplex *ab, cblas_int_t *ldab, doublecomplex *b, cblas_int_t *
    ldb, cblas_int_t *info);

/* Subroutine */ int zpbsvx_(const char *fact, const char *uplo, cblas_int_t *n, cblas_int_t *kd, 
    cblas_int_t *nrhs, doublecomplex *ab, cblas_int_t *ldab, doublecomplex *afb, 
    cblas_int_t *ldafb, const char *equed, double *s, doublecomplex *b, cblas_int_t 
    *ldb, doublecomplex *x, cblas_int_t *ldx, double *rcond, double *
    ferr, double *berr, doublecomplex *work, double *rwork, 
    cblas_int_t *info);

/* Subroutine */ int zpbtf2_(const char *uplo, cblas_int_t *n, cblas_int_t *kd, 
    doublecomplex *ab, cblas_int_t *ldab, cblas_int_t *info);

/* Subroutine */ int zpbtrf_(const char *uplo, cblas_int_t *n, cblas_int_t *kd, 
    doublecomplex *ab, cblas_int_t *ldab, cblas_int_t *info);

/* Subroutine */ int zpbtrs_(const char *uplo, cblas_int_t *n, cblas_int_t *kd, cblas_int_t *
    nrhs, doublecomplex *ab, cblas_int_t *ldab, doublecomplex *b, cblas_int_t *
    ldb, cblas_int_t *info);

/* Subroutine */ int zpocon_(const char *uplo, cblas_int_t *n, doublecomplex *a, 
    cblas_int_t *lda, double *anorm, double *rcond, doublecomplex *
    work, double *rwork, cblas_int_t *info);

/* Subroutine */ int zpoequ_(cblas_int_t *n, doublecomplex *a, cblas_int_t *lda, 
    double *s, double *scond, double *amax, cblas_int_t *info);

/* Subroutine */ int zporfs_(const char *uplo, cblas_int_t *n, cblas_int_t *nrhs, 
    doublecomplex *a, cblas_int_t *lda, doublecomplex *af, cblas_int_t *ldaf, 
    doublecomplex *b, cblas_int_t *ldb, doublecomplex *x, cblas_int_t *ldx, 
    double *ferr, double *berr, doublecomplex *work, double *
    rwork, cblas_int_t *info);

/* Subroutine */ int zposv_(const char *uplo, cblas_int_t *n, cblas_int_t *nrhs, 
    doublecomplex *a, cblas_int_t *lda, doublecomplex *b, cblas_int_t *ldb, 
    cblas_int_t *info);

/* Subroutine */ int zposvx_(const char *fact, const char *uplo, cblas_int_t *n, cblas_int_t *
    nrhs, doublecomplex *a, cblas_int_t *lda, doublecomplex *af, cblas_int_t *
    ldaf, const char *equed, double *s, doublecomplex *b, cblas_int_t *ldb, 
    doublecomplex *x, cblas_int_t *ldx, double *rcond, double *ferr, 
    double *berr, doublecomplex *work, double *rwork, cblas_int_t *
    info);

/* Subroutine */ int zpotf2_(const char *uplo, cblas_int_t *n, doublecomplex *a, 
    cblas_int_t *lda, cblas_int_t *info);

/* Subroutine */ int zpotrf_(const char *uplo, cblas_int_t *n, doublecomplex *a, 
    cblas_int_t *lda, cblas_int_t *info);

/* Subroutine */ int zpotri_(const char *uplo, cblas_int_t *n, doublecomplex *a, 
    cblas_int_t *lda, cblas_int_t *info);

/* Subroutine */ int zpotrs_(const char *uplo, cblas_int_t *n, cblas_int_t *nrhs, 
    doublecomplex *a, cblas_int_t *lda, doublecomplex *b, cblas_int_t *ldb, 
    cblas_int_t *info);

/* Subroutine */ int zppcon_(const char *uplo, cblas_int_t *n, doublecomplex *ap, 
    double *anorm, double *rcond, doublecomplex *work, double 
    *rwork, cblas_int_t *info);

/* Subroutine */ int zppequ_(const char *uplo, cblas_int_t *n, doublecomplex *ap, 
    double *s, double *scond, double *amax, cblas_int_t *info);

/* Subroutine */ int zpprfs_(const char *uplo, cblas_int_t *n, cblas_int_t *nrhs, 
    doublecomplex *ap, doublecomplex *afp, doublecomplex *b, cblas_int_t *ldb,
     doublecomplex *x, cblas_int_t *ldx, double *ferr, double *berr, 
    doublecomplex *work, double *rwork, cblas_int_t *info);

/* Subroutine */ int zppsv_(const char *uplo, cblas_int_t *n, cblas_int_t *nrhs, 
    doublecomplex *ap, doublecomplex *b, cblas_int_t *ldb, cblas_int_t *info);

/* Subroutine */ int zppsvx_(const char *fact, const char *uplo, cblas_int_t *n, cblas_int_t *
    nrhs, doublecomplex *ap, doublecomplex *afp, const char *equed, double *
    s, doublecomplex *b, cblas_int_t *ldb, doublecomplex *x, cblas_int_t *ldx, 
    double *rcond, double *ferr, double *berr, doublecomplex *
    work, double *rwork, cblas_int_t *info);

/* Subroutine */ int zpptrf_(const char *uplo, cblas_int_t *n, doublecomplex *ap, 
    cblas_int_t *info);

/* Subroutine */ int zpptri_(const char *uplo, cblas_int_t *n, doublecomplex *ap, 
    cblas_int_t *info);

/* Subroutine */ int zpptrs_(const char *uplo, cblas_int_t *n, cblas_int_t *nrhs, 
    doublecomplex *ap, doublecomplex *b, cblas_int_t *ldb, cblas_int_t *info);

/* Subroutine */ int zptcon_(cblas_int_t *n, double *d__, doublecomplex *e, 
    double *anorm, double *rcond, double *rwork, cblas_int_t *
    info);

/* Subroutine */ int zpteqr_(const char *compz, cblas_int_t *n, double *d__, 
    double *e, doublecomplex *z__, cblas_int_t *ldz, double *work, 
    cblas_int_t *info);

/* Subroutine */ int zptrfs_(const char *uplo, cblas_int_t *n, cblas_int_t *nrhs, 
    double *d__, doublecomplex *e, double *df, doublecomplex *ef, 
    doublecomplex *b, cblas_int_t *ldb, doublecomplex *x, cblas_int_t *ldx, 
    double *ferr, double *berr, doublecomplex *work, double *
    rwork, cblas_int_t *info);

/* Subroutine */ int zptsv_(cblas_int_t *n, cblas_int_t *nrhs, double *d__, 
    doublecomplex *e, doublecomplex *b, cblas_int_t *ldb, cblas_int_t *info);

/* Subroutine */ int zptsvx_(const char *fact, cblas_int_t *n, cblas_int_t *nrhs, 
    double *d__, doublecomplex *e, double *df, doublecomplex *ef, 
    doublecomplex *b, cblas_int_t *ldb, doublecomplex *x, cblas_int_t *ldx, 
    double *rcond, double *ferr, double *berr, doublecomplex *
    work, double *rwork, cblas_int_t *info);

/* Subroutine */ int zpttrf_(cblas_int_t *n, double *d__, doublecomplex *e, 
    cblas_int_t *info);

/* Subroutine */ int zpttrs_(const char *uplo, cblas_int_t *n, cblas_int_t *nrhs, 
    double *d__, doublecomplex *e, doublecomplex *b, cblas_int_t *ldb, 
    cblas_int_t *info);

/* Subroutine */ int zptts2_(cblas_int_t *iuplo, cblas_int_t *n, cblas_int_t *nrhs, 
    double *d__, doublecomplex *e, doublecomplex *b, cblas_int_t *ldb);

/* Subroutine */ int zrot_(cblas_int_t *n, doublecomplex *cx, cblas_int_t *incx, 
    doublecomplex *cy, cblas_int_t *incy, double *c__, doublecomplex *s);

/* Subroutine */ int zspcon_(const char *uplo, cblas_int_t *n, doublecomplex *ap, 
    cblas_int_t *ipiv, double *anorm, double *rcond, doublecomplex *
    work, cblas_int_t *info);

/* Subroutine */ int zspmv_(const char *uplo, cblas_int_t *n, doublecomplex *alpha, 
    doublecomplex *ap, doublecomplex *x, cblas_int_t *incx, doublecomplex *
    beta, doublecomplex *y, cblas_int_t *incy);

/* Subroutine */ int zspr_(const char *uplo, cblas_int_t *n, doublecomplex *alpha, 
    doublecomplex *x, cblas_int_t *incx, doublecomplex *ap);

/* Subroutine */ int zsprfs_(const char *uplo, cblas_int_t *n, cblas_int_t *nrhs, 
    doublecomplex *ap, doublecomplex *afp, cblas_int_t *ipiv, doublecomplex *
    b, cblas_int_t *ldb, doublecomplex *x, cblas_int_t *ldx, double *ferr, 
    double *berr, doublecomplex *work, double *rwork, cblas_int_t *
    info);

/* Subroutine */ int zspsv_(const char *uplo, cblas_int_t *n, cblas_int_t *nrhs, 
    doublecomplex *ap, cblas_int_t *ipiv, doublecomplex *b, cblas_int_t *ldb, 
    cblas_int_t *info);

/* Subroutine */ int zspsvx_(const char *fact, const char *uplo, cblas_int_t *n, cblas_int_t *
    nrhs, doublecomplex *ap, doublecomplex *afp, cblas_int_t *ipiv, 
    doublecomplex *b, cblas_int_t *ldb, doublecomplex *x, cblas_int_t *ldx, 
    double *rcond, double *ferr, double *berr, doublecomplex *
    work, double *rwork, cblas_int_t *info);

/* Subroutine */ int zsptrf_(const char *uplo, cblas_int_t *n, doublecomplex *ap, 
    cblas_int_t *ipiv, cblas_int_t *info);

/* Subroutine */ int zsptri_(const char *uplo, cblas_int_t *n, doublecomplex *ap, 
    cblas_int_t *ipiv, doublecomplex *work, cblas_int_t *info);

/* Subroutine */ int zsptrs_(const char *uplo, cblas_int_t *n, cblas_int_t *nrhs, 
    doublecomplex *ap, cblas_int_t *ipiv, doublecomplex *b, cblas_int_t *ldb, 
    cblas_int_t *info);

/* Subroutine */ int zstedc_(const char *compz, cblas_int_t *n, double *d__, 
    double *e, doublecomplex *z__, cblas_int_t *ldz, doublecomplex *work, 
    cblas_int_t *lwork, double *rwork, cblas_int_t *lrwork, cblas_int_t *iwork, 
    cblas_int_t *liwork, cblas_int_t *info);

/* Subroutine */ int zstegr_(const char *jobz, const char *range, cblas_int_t *n, double *
    d__, double *e, double *vl, double *vu, cblas_int_t *il, 
    cblas_int_t *iu, double *abstol, cblas_int_t *m, double *w, 
    doublecomplex *z__, cblas_int_t *ldz, cblas_int_t *isuppz, double *work, 
    cblas_int_t *lwork, cblas_int_t *iwork, cblas_int_t *liwork, cblas_int_t *info);

/* Subroutine */ int zstein_(cblas_int_t *n, double *d__, double *e, 
    cblas_int_t *m, double *w, cblas_int_t *iblock, cblas_int_t *isplit, 
    doublecomplex *z__, cblas_int_t *ldz, double *work, cblas_int_t *iwork, 
    cblas_int_t *ifail, cblas_int_t *info);

/* Subroutine */ int zstemr_(const char *jobz, const char *range, cblas_int_t *n, double *
    d__, double *e, double *vl, double *vu, cblas_int_t *il, 
    cblas_int_t *iu, cblas_int_t *m, double *w, doublecomplex *z__, cblas_int_t *
    ldz, cblas_int_t *nzc, cblas_int_t *isuppz, logical *tryrac, double *work,
     cblas_int_t *lwork, cblas_int_t *iwork, cblas_int_t *liwork, cblas_int_t *info);

/* Subroutine */ int zsteqr_(const char *compz, cblas_int_t *n, double *d__, 
    double *e, doublecomplex *z__, cblas_int_t *ldz, double *work, 
    cblas_int_t *info);

/* Subroutine */ int zsycon_(const char *uplo, cblas_int_t *n, doublecomplex *a, 
    cblas_int_t *lda, cblas_int_t *ipiv, double *anorm, double *rcond, 
    doublecomplex *work, cblas_int_t *info);

/* Subroutine */ int zsymv_(const char *uplo, cblas_int_t *n, doublecomplex *alpha, 
    doublecomplex *a, cblas_int_t *lda, doublecomplex *x, cblas_int_t *incx, 
    doublecomplex *beta, doublecomplex *y, cblas_int_t *incy);

/* Subroutine */ int zsyr_(const char *uplo, cblas_int_t *n, doublecomplex *alpha, 
    doublecomplex *x, cblas_int_t *incx, doublecomplex *a, cblas_int_t *lda);

/* Subroutine */ int zsyrfs_(const char *uplo, cblas_int_t *n, cblas_int_t *nrhs, 
    doublecomplex *a, cblas_int_t *lda, doublecomplex *af, cblas_int_t *ldaf, 
    cblas_int_t *ipiv, doublecomplex *b, cblas_int_t *ldb, doublecomplex *x, 
    cblas_int_t *ldx, double *ferr, double *berr, doublecomplex *work,
     double *rwork, cblas_int_t *info);

/* Subroutine */ int zsysv_(const char *uplo, cblas_int_t *n, cblas_int_t *nrhs, 
    doublecomplex *a, cblas_int_t *lda, cblas_int_t *ipiv, doublecomplex *b, 
    cblas_int_t *ldb, doublecomplex *work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int zsysvx_(const char *fact, const char *uplo, cblas_int_t *n, cblas_int_t *
    nrhs, doublecomplex *a, cblas_int_t *lda, doublecomplex *af, cblas_int_t *
    ldaf, cblas_int_t *ipiv, doublecomplex *b, cblas_int_t *ldb, doublecomplex *x,
     cblas_int_t *ldx, double *rcond, double *ferr, double *berr, 
    doublecomplex *work, cblas_int_t *lwork, double *rwork, cblas_int_t *info);

/* Subroutine */ int zsytf2_(const char *uplo, cblas_int_t *n, doublecomplex *a, 
    cblas_int_t *lda, cblas_int_t *ipiv, cblas_int_t *info);

/* Subroutine */ int zsytrf_(const char *uplo, cblas_int_t *n, doublecomplex *a, 
    cblas_int_t *lda, cblas_int_t *ipiv, doublecomplex *work, cblas_int_t *lwork, 
    cblas_int_t *info);

/* Subroutine */ int zsytri_(const char *uplo, cblas_int_t *n, doublecomplex *a, 
    cblas_int_t *lda, cblas_int_t *ipiv, doublecomplex *work, cblas_int_t *info);

/* Subroutine */ int zsytrs_(const char *uplo, cblas_int_t *n, cblas_int_t *nrhs, 
    doublecomplex *a, cblas_int_t *lda, cblas_int_t *ipiv, doublecomplex *b, 
    cblas_int_t *ldb, cblas_int_t *info);

/* Subroutine */ int ztbcon_(const char *norm, const char *uplo, const char *diag, cblas_int_t *n, 
    cblas_int_t *kd, doublecomplex *ab, cblas_int_t *ldab, double *rcond, 
    doublecomplex *work, double *rwork, cblas_int_t *info);

/* Subroutine */ int ztbrfs_(const char *uplo, const char *trans, const char *diag, cblas_int_t *n, 
    cblas_int_t *kd, cblas_int_t *nrhs, doublecomplex *ab, cblas_int_t *ldab, 
    doublecomplex *b, cblas_int_t *ldb, doublecomplex *x, cblas_int_t *ldx, 
    double *ferr, double *berr, doublecomplex *work, double *
    rwork, cblas_int_t *info);

/* Subroutine */ int ztbtrs_(const char *uplo, const char *trans, const char *diag, cblas_int_t *n, 
    cblas_int_t *kd, cblas_int_t *nrhs, doublecomplex *ab, cblas_int_t *ldab, 
    doublecomplex *b, cblas_int_t *ldb, cblas_int_t *info);

/* Subroutine */ int ztgevc_(const char *side, const char *howmny, logical *select, 
    cblas_int_t *n, doublecomplex *s, cblas_int_t *lds, doublecomplex *p, cblas_int_t 
    *ldp, doublecomplex *vl, cblas_int_t *ldvl, doublecomplex *vr, cblas_int_t *
    ldvr, cblas_int_t *mm, cblas_int_t *m, doublecomplex *work, double *rwork,
     cblas_int_t *info);

/* Subroutine */ int ztgex2_(logical *wantq, logical *wantz, cblas_int_t *n, 
    doublecomplex *a, cblas_int_t *lda, doublecomplex *b, cblas_int_t *ldb, 
    doublecomplex *q, cblas_int_t *ldq, doublecomplex *z__, cblas_int_t *ldz, 
    cblas_int_t *j1, cblas_int_t *info);

/* Subroutine */ int ztgexc_(logical *wantq, logical *wantz, cblas_int_t *n, 
    doublecomplex *a, cblas_int_t *lda, doublecomplex *b, cblas_int_t *ldb, 
    doublecomplex *q, cblas_int_t *ldq, doublecomplex *z__, cblas_int_t *ldz, 
    cblas_int_t *ifst, cblas_int_t *ilst, cblas_int_t *info);

/* Subroutine */ int ztgsen_(cblas_int_t *ijob, logical *wantq, logical *wantz, 
    logical *select, cblas_int_t *n, doublecomplex *a, cblas_int_t *lda, 
    doublecomplex *b, cblas_int_t *ldb, doublecomplex *alpha, doublecomplex *
    beta, doublecomplex *q, cblas_int_t *ldq, doublecomplex *z__, cblas_int_t *
    ldz, cblas_int_t *m, double *pl, double *pr, double *dif, 
    doublecomplex *work, cblas_int_t *lwork, cblas_int_t *iwork, cblas_int_t *liwork, 
    cblas_int_t *info);

/* Subroutine */ int ztgsja_(const char *jobu, const char *jobv, const char *jobq, cblas_int_t *m, 
    cblas_int_t *p, cblas_int_t *n, cblas_int_t *k, cblas_int_t *l, doublecomplex *a, 
    cblas_int_t *lda, doublecomplex *b, cblas_int_t *ldb, double *tola, 
    double *tolb, double *alpha, double *beta, doublecomplex *
    u, cblas_int_t *ldu, doublecomplex *v, cblas_int_t *ldv, doublecomplex *q, 
    cblas_int_t *ldq, doublecomplex *work, cblas_int_t *ncycle, cblas_int_t *info);

/* Subroutine */ int ztgsna_(const char *job, const char *howmny, logical *select, 
    cblas_int_t *n, doublecomplex *a, cblas_int_t *lda, doublecomplex *b, cblas_int_t 
    *ldb, doublecomplex *vl, cblas_int_t *ldvl, doublecomplex *vr, cblas_int_t *
    ldvr, double *s, double *dif, cblas_int_t *mm, cblas_int_t *m, 
    doublecomplex *work, cblas_int_t *lwork, cblas_int_t *iwork, cblas_int_t *info);

/* Subroutine */ int ztgsy2_(const char *trans, cblas_int_t *ijob, cblas_int_t *m, cblas_int_t *
    n, doublecomplex *a, cblas_int_t *lda, doublecomplex *b, cblas_int_t *ldb, 
    doublecomplex *c__, cblas_int_t *ldc, doublecomplex *d__, cblas_int_t *ldd, 
    doublecomplex *e, cblas_int_t *lde, doublecomplex *f, cblas_int_t *ldf, 
    double *scale, double *rdsum, double *rdscal, cblas_int_t *
    info);

/* Subroutine */ int ztgsyl_(const char *trans, cblas_int_t *ijob, cblas_int_t *m, cblas_int_t *
    n, doublecomplex *a, cblas_int_t *lda, doublecomplex *b, cblas_int_t *ldb, 
    doublecomplex *c__, cblas_int_t *ldc, doublecomplex *d__, cblas_int_t *ldd, 
    doublecomplex *e, cblas_int_t *lde, doublecomplex *f, cblas_int_t *ldf, 
    double *scale, double *dif, doublecomplex *work, cblas_int_t *
    lwork, cblas_int_t *iwork, cblas_int_t *info);

/* Subroutine */ int ztpcon_(const char *norm, const char *uplo, const char *diag, cblas_int_t *n, 
    doublecomplex *ap, double *rcond, doublecomplex *work, double 
    *rwork, cblas_int_t *info);

/* Subroutine */ int ztprfs_(const char *uplo, const char *trans, const char *diag, cblas_int_t *n, 
    cblas_int_t *nrhs, doublecomplex *ap, doublecomplex *b, cblas_int_t *ldb, 
    doublecomplex *x, cblas_int_t *ldx, double *ferr, double *berr, 
    doublecomplex *work, double *rwork, cblas_int_t *info);

/* Subroutine */ int ztptri_(const char *uplo, const char *diag, cblas_int_t *n, 
    doublecomplex *ap, cblas_int_t *info);

/* Subroutine */ int ztptrs_(const char *uplo, const char *trans, const char *diag, cblas_int_t *n, 
    cblas_int_t *nrhs, doublecomplex *ap, doublecomplex *b, cblas_int_t *ldb, 
    cblas_int_t *info);

/* Subroutine */ int ztrcon_(const char *norm, const char *uplo, const char *diag, cblas_int_t *n, 
    doublecomplex *a, cblas_int_t *lda, double *rcond, doublecomplex *
    work, double *rwork, cblas_int_t *info);

/* Subroutine */ int ztrevc_(const char *side, const char *howmny, logical *select, 
    cblas_int_t *n, doublecomplex *t, cblas_int_t *ldt, doublecomplex *vl, 
    cblas_int_t *ldvl, doublecomplex *vr, cblas_int_t *ldvr, cblas_int_t *mm, cblas_int_t 
    *m, doublecomplex *work, double *rwork, cblas_int_t *info);

/* Subroutine */ int ztrexc_(const char *compq, cblas_int_t *n, doublecomplex *t, 
    cblas_int_t *ldt, doublecomplex *q, cblas_int_t *ldq, cblas_int_t *ifst, cblas_int_t *
    ilst, cblas_int_t *info);

/* Subroutine */ int ztrrfs_(const char *uplo, const char *trans, const char *diag, cblas_int_t *n, 
    cblas_int_t *nrhs, doublecomplex *a, cblas_int_t *lda, doublecomplex *b, 
    cblas_int_t *ldb, doublecomplex *x, cblas_int_t *ldx, double *ferr, 
    double *berr, doublecomplex *work, double *rwork, cblas_int_t *
    info);

/* Subroutine */ int ztrsen_(const char *job, const char *compq, logical *select, cblas_int_t 
    *n, doublecomplex *t, cblas_int_t *ldt, doublecomplex *q, cblas_int_t *ldq, 
    doublecomplex *w, cblas_int_t *m, double *s, double *sep, 
    doublecomplex *work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int ztrsna_(const char *job, const char *howmny, logical *select, 
    cblas_int_t *n, doublecomplex *t, cblas_int_t *ldt, doublecomplex *vl, 
    cblas_int_t *ldvl, doublecomplex *vr, cblas_int_t *ldvr, double *s, 
    double *sep, cblas_int_t *mm, cblas_int_t *m, doublecomplex *work, 
    cblas_int_t *ldwork, double *rwork, cblas_int_t *info);

/* Subroutine */ int ztrsyl_(const char *trana, const char *tranb, cblas_int_t *isgn, cblas_int_t 
    *m, cblas_int_t *n, doublecomplex *a, cblas_int_t *lda, doublecomplex *b, 
    cblas_int_t *ldb, doublecomplex *c__, cblas_int_t *ldc, double *scale, 
    cblas_int_t *info);

/* Subroutine */ int ztrti2_(const char *uplo, const char *diag, cblas_int_t *n, 
    doublecomplex *a, cblas_int_t *lda, cblas_int_t *info);

/* Subroutine */ int ztrtri_(const char *uplo, const char *diag, cblas_int_t *n, 
    doublecomplex *a, cblas_int_t *lda, cblas_int_t *info);

/* Subroutine */ int ztrtrs_(const char *uplo, const char *trans, const char *diag, cblas_int_t *n, 
    cblas_int_t *nrhs, doublecomplex *a, cblas_int_t *lda, doublecomplex *b, 
    cblas_int_t *ldb, cblas_int_t *info);

/* Subroutine */ int ztzrqf_(cblas_int_t *m, cblas_int_t *n, doublecomplex *a, 
    cblas_int_t *lda, doublecomplex *tau, cblas_int_t *info);

/* Subroutine */ int ztzrzf_(cblas_int_t *m, cblas_int_t *n, doublecomplex *a, 
    cblas_int_t *lda, doublecomplex *tau, doublecomplex *work, cblas_int_t *lwork,
     cblas_int_t *info);

/* Subroutine */ int zung2l_(cblas_int_t *m, cblas_int_t *n, cblas_int_t *k, 
    doublecomplex *a, cblas_int_t *lda, doublecomplex *tau, doublecomplex *
    work, cblas_int_t *info);

/* Subroutine */ int zung2r_(cblas_int_t *m, cblas_int_t *n, cblas_int_t *k, 
    doublecomplex *a, cblas_int_t *lda, doublecomplex *tau, doublecomplex *
    work, cblas_int_t *info);

/* Subroutine */ int zungbr_(const char *vect, cblas_int_t *m, cblas_int_t *n, cblas_int_t *k, 
    doublecomplex *a, cblas_int_t *lda, doublecomplex *tau, doublecomplex *
    work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int zunghr_(cblas_int_t *n, cblas_int_t *ilo, cblas_int_t *ihi, 
    doublecomplex *a, cblas_int_t *lda, doublecomplex *tau, doublecomplex *
    work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int zungl2_(cblas_int_t *m, cblas_int_t *n, cblas_int_t *k, 
    doublecomplex *a, cblas_int_t *lda, doublecomplex *tau, doublecomplex *
    work, cblas_int_t *info);

/* Subroutine */ int zunglq_(cblas_int_t *m, cblas_int_t *n, cblas_int_t *k, 
    doublecomplex *a, cblas_int_t *lda, doublecomplex *tau, doublecomplex *
    work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int zungql_(cblas_int_t *m, cblas_int_t *n, cblas_int_t *k, 
    doublecomplex *a, cblas_int_t *lda, doublecomplex *tau, doublecomplex *
    work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int zungqr_(cblas_int_t *m, cblas_int_t *n, cblas_int_t *k, 
    doublecomplex *a, cblas_int_t *lda, doublecomplex *tau, doublecomplex *
    work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int zungr2_(cblas_int_t *m, cblas_int_t *n, cblas_int_t *k, 
    doublecomplex *a, cblas_int_t *lda, doublecomplex *tau, doublecomplex *
    work, cblas_int_t *info);

/* Subroutine */ int zungrq_(cblas_int_t *m, cblas_int_t *n, cblas_int_t *k, 
    doublecomplex *a, cblas_int_t *lda, doublecomplex *tau, doublecomplex *
    work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int zungtr_(const char *uplo, cblas_int_t *n, doublecomplex *a, 
    cblas_int_t *lda, doublecomplex *tau, doublecomplex *work, cblas_int_t *lwork,
     cblas_int_t *info);

/* Subroutine */ int zunm2l_(const char *side, const char *trans, cblas_int_t *m, cblas_int_t *n, 
    cblas_int_t *k, doublecomplex *a, cblas_int_t *lda, doublecomplex *tau, 
    doublecomplex *c__, cblas_int_t *ldc, doublecomplex *work, cblas_int_t *info);

/* Subroutine */ int zunm2r_(const char *side, const char *trans, cblas_int_t *m, cblas_int_t *n, 
    cblas_int_t *k, doublecomplex *a, cblas_int_t *lda, doublecomplex *tau, 
    doublecomplex *c__, cblas_int_t *ldc, doublecomplex *work, cblas_int_t *info);

/* Subroutine */ int zunmbr_(const char *vect, const char *side, const char *trans, cblas_int_t *m, 
    cblas_int_t *n, cblas_int_t *k, doublecomplex *a, cblas_int_t *lda, doublecomplex 
    *tau, doublecomplex *c__, cblas_int_t *ldc, doublecomplex *work, cblas_int_t *  
    lwork, cblas_int_t *info);

/* Subroutine */ int zunmhr_(const char *side, const char *trans, cblas_int_t *m, cblas_int_t *n, 
    cblas_int_t *ilo, cblas_int_t *ihi, doublecomplex *a, cblas_int_t *lda, 
    doublecomplex *tau, doublecomplex *c__, cblas_int_t *ldc, doublecomplex *   
    work, cblas_int_t *lwork, cblas_int_t *info);

/* Subroutine */ int zunml2_(const char *side, const char *trans, cblas_int_t *m, cblas_int_t *n, 
    cblas_int_t *k, doublecomplex *a, cblas_int_t *lda, doublecomplex *tau, 
    doublecomplex *c__, cblas_int_t *ldc, doublecomplex *work, cblas_int_t *info);

/* Subroutine */ int zunmlq_(const char *side, const char *trans, cblas_int_t *m, cblas_int_t *n, 
    cblas_int_t *k, doublecomplex *a, cblas_int_t *lda, doublecomplex *tau, 
    doublecomplex *c__, cblas_int_t *ldc, doublecomplex *work, cblas_int_t *lwork,
     cblas_int_t *info);

/* Subroutine */ int zunmql_(const char *side, const char *trans, cblas_int_t *m, cblas_int_t *n, 
    cblas_int_t *k, doublecomplex *a, cblas_int_t *lda, doublecomplex *tau, 
    doublecomplex *c__, cblas_int_t *ldc, doublecomplex *work, cblas_int_t *lwork,
     cblas_int_t *info);

/* Subroutine */ int zunmqr_(const char *side, const char *trans, cblas_int_t *m, cblas_int_t *n, 
    cblas_int_t *k, doublecomplex *a, cblas_int_t *lda, doublecomplex *tau, 
    doublecomplex *c__, cblas_int_t *ldc, doublecomplex *work, cblas_int_t *lwork,
     cblas_int_t *info);

/* Subroutine */ int zunmr2_(const char *side, const char *trans, cblas_int_t *m, cblas_int_t *n, 
    cblas_int_t *k, doublecomplex *a, cblas_int_t *lda, doublecomplex *tau, 
    doublecomplex *c__, cblas_int_t *ldc, doublecomplex *work, cblas_int_t *info);

/* Subroutine */ int zunmr3_(const char *side, const char *trans, cblas_int_t *m, cblas_int_t *n, 
    cblas_int_t *k, cblas_int_t *l, doublecomplex *a, cblas_int_t *lda, doublecomplex 
    *tau, doublecomplex *c__, cblas_int_t *ldc, doublecomplex *work, cblas_int_t *
    info);

/* Subroutine */ int zunmrq_(const char *side, const char *trans, cblas_int_t *m, cblas_int_t *n, 
    cblas_int_t *k, doublecomplex *a, cblas_int_t *lda, doublecomplex *tau, 
    doublecomplex *c__, cblas_int_t *ldc, doublecomplex *work, cblas_int_t *lwork,
     cblas_int_t *info);

/* Subroutine */ int zunmrz_(const char *side, const char *trans, cblas_int_t *m, cblas_int_t *n, 
    cblas_int_t *k, cblas_int_t *l, doublecomplex *a, cblas_int_t *lda, doublecomplex 
    *tau, doublecomplex *c__, cblas_int_t *ldc, doublecomplex *work, cblas_int_t *
    lwork, cblas_int_t *info);

/* Subroutine */ int zunmtr_(const char *side, const char *uplo, const char *trans, cblas_int_t *m, 
    cblas_int_t *n, doublecomplex *a, cblas_int_t *lda, doublecomplex *tau, 
    doublecomplex *c__, cblas_int_t *ldc, doublecomplex *work, cblas_int_t *lwork,
     cblas_int_t *info);

/* Subroutine */ int zupgtr_(const char *uplo, cblas_int_t *n, doublecomplex *ap, 
    doublecomplex *tau, doublecomplex *q, cblas_int_t *ldq, doublecomplex *
    work, cblas_int_t *info);

/* Subroutine */ int zupmtr_(const char *side, const char *uplo, const char *trans, cblas_int_t *m, 
    cblas_int_t *n, doublecomplex *ap, doublecomplex *tau, doublecomplex *c__,
     cblas_int_t *ldc, doublecomplex *work, cblas_int_t *info);

#ifdef __cplusplus
}
#endif

#endif // _GLUE_CLAPACK_H_INCLUDE_

