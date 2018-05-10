/** 
 * @file lapack.h
 * @brief      F77 LAPACK partial header file
 * @author     Christian Kirches
 * @date       2009 Jun 1
 */

#pragma once
#ifndef GLUE_LAPACK_H_INCLUDED_
#define GLUE_LAPACK_H_INCLUDED_

#ifdef __cplusplus
extern "C" {
#endif

#include "util_types.h"
#include "blas.h"

/** DGEQRF: blocked or unblocked QR decomposition
 */
extern void dgeqrf_ (
    const blas_int_t *m,
    const blas_int_t *n,
    double     *A,
    const blas_int_t *lda,
    double     *Tau,
    double     *work,
    const blas_int_t *lwork,
    blas_int_t       *info
);

/** DGERQF: blocked or unblocked RQ decomposition
 */
extern void dgerqf_ (
    const blas_int_t *m,
    const blas_int_t *n,
    double     *A,
    const blas_int_t *lda,
    double     *Tau,
    double     *work,
    const blas_int_t *lwork,
    blas_int_t       *info
);

/** DGEQRF: blocked or unblocked QL decomposition
 */
extern void dgeqlf_ (
    const blas_int_t *m,
    const blas_int_t *n,
    double     *A,
    const blas_int_t *lda,
    double     *Tau,
    double     *work,
    const blas_int_t *lwork,
    blas_int_t       *info
);

/** DGERQF: blocked or unblocked LQ decomposition
 */
extern void dgelqf_ (
    const blas_int_t *m,
    const blas_int_t *n,
    double     *A,
    const blas_int_t *lda,
    double     *Tau,
    double     *work,
    const blas_int_t *lwork,
    blas_int_t       *info
);

/** DORGQR: compute unitary factor Q from DGEQRF outputs
 */
extern void dorgqr_ (
    const blas_int_t   *m,
    const blas_int_t   *n,
    const blas_int_t   *k,
    double       *A,
    const blas_int_t   *lda,
    const double *Tau,
    double       *work,
    const blas_int_t   *lwork,
    blas_int_t         *info
);

/** DORGRQ: compute unitary factor Q from DGERQF outputs
 */
extern void dorgrq_ (
    const blas_int_t   *m,
    const blas_int_t   *n,
    const blas_int_t   *k,
    double       *A,
    const blas_int_t   *lda,
    const double *Tau,
    double       *work,
    const blas_int_t   *lwork,
    blas_int_t         *info
);

/** DORGQL: compute unitary factor Q from DGEQLF outputs
 */
extern void dorgql_ (
    const blas_int_t   *m,
    const blas_int_t   *n,
    const blas_int_t   *k,
    double       *A,
    const blas_int_t   *lda,
    const double *Tau,
    double       *work,
    const blas_int_t   *lwork,
    blas_int_t         *info
);

/** DORGLQ: compute unitary factor Q from DGELQF outputs
 */
extern void dorglq_ (
    const blas_int_t   *m,
    const blas_int_t   *n,
    const blas_int_t   *k,
    double       *A,
    const blas_int_t   *lda,
    const double *Tau,
    double       *work,
    const blas_int_t   *lwork,
    blas_int_t         *info
);

/** DORMQR: multiply with unitary factor Q from DGEQRF outputs
 */
extern void dormqr_ (
    const char   *side,
    const char   *trans,
    const blas_int_t   *m,
    const blas_int_t   *n,
    const blas_int_t   *k,
    const double *A,
    const blas_int_t   *lda,
    const double *Tau,
    double       *C,
    const blas_int_t   *ldc,
    double       *work,
    const blas_int_t   *lwork,
    blas_int_t         *info,
    blas_int_t          side_len,
    blas_int_t          trans_len
);

/** DORMRQ: multiply with unitary factor Q from DGERQF outputs
 */
extern void dormrq_ (
    const char   *side,
    const char   *trans,
    const blas_int_t   *m,
    const blas_int_t   *n,
    const blas_int_t   *k,
    const double *A,
    const blas_int_t   *lda,
    const double *Tau,
    double       *C,
    const blas_int_t   *ldc,
    double       *work,
    const blas_int_t   *lwork,
    blas_int_t         *info,
    blas_int_t          side_len,
    blas_int_t          trans_len
);

/** DPOTRF: Cholesky decomposition
 */
extern void dpotrf_ (
    const char *uplo,
    const blas_int_t *n,
    double     *A,
    const blas_int_t *lda,
    blas_int_t       *info,
    blas_int_t        uplo_len
);

/** DTRTRS: Backsolve with a triangular matrix
 */
extern void dtrtrs_ (
    const char   *uplo,
    const char   *trans,
    const char   *diag,
    const blas_int_t   *n,
    const blas_int_t   *nrhs,
    const double *A,
    const blas_int_t   *lda,
    double       *B,
    const blas_int_t   *ldb,
    blas_int_t         *info,
    blas_int_t          uplo_len,
    blas_int_t          trans_len,
    blas_int_t          diag_len
);


/** DGESV: Solve a linear equation system using a pivoting LU decomposition
 */
extern void dgesv_ (
    const blas_int_t *n,
    const blas_int_t *nrhs,
    double     *A,
    const blas_int_t *lda,
    blas_int_t       *ipiv,
    double     *B,
    const blas_int_t *ldb,
    blas_int_t       *info
);

/** DSYSV: Solve a linear symmetric indefinite equation system using a
 *        pivoting LBL^T decomposition
 */
extern void dsysv_ (
    const char *uplo,
    const blas_int_t *n,
    const blas_int_t *nrhs,
    double     *A,
    const blas_int_t *lda,
    blas_int_t       *ipiv,
    double     *B,
    const blas_int_t *ldb,
    double     *work,
    blas_int_t       *lwork,
    blas_int_t       *info,
    blas_int_t        uplo_len
);

/** DSYTRF: Symmetric indefinite decomposition using Bunch-Kaufman diagonal
 *         pivoting
 */
extern void dsytrf_ (
    const char *uplo,
    const blas_int_t *n,
    double     *A,
    const blas_int_t *lda,
    blas_int_t       *ipiv,
    double     *work,
    blas_int_t       *lwork,
    blas_int_t       *info,
    blas_int_t        uplo_len
);

/** DSYTRS: Solve system A*X=B using decomposition from DSYTRF
 */
extern void dsytrs_ (
    const char *uplo,
    const blas_int_t *n,
    const blas_int_t *nrhs,
    double     *A,
    const blas_int_t *lda,
    blas_int_t       *ipiv,
    double     *b,
    const blas_int_t *ldb,
    blas_int_t       *info,
    blas_int_t        uplo_len
);

/** DPOCON: Condition number estimate from Cholesky factor DPOTRF
 */
extern void dpocon_ (
    const char   *uplo,
    const blas_int_t   *n,
    const double *A,
    const blas_int_t   *lda,
    const double *Anorm,
    double       *rcond,
    double       *work,
    blas_int_t         *iwork,
    blas_int_t         *info,
    blas_int_t          uplo_len
);

/** DGBTRF: Band LU decomposition
 */
extern void dgbtrf_ (
    const blas_int_t *m,
    const blas_int_t *n,
    const blas_int_t *kl,
    const blas_int_t *ku,
    double     *ab,
    const blas_int_t *ldab,
    blas_int_t       *ipiv,
    blas_int_t       *info    
);
 
/** DGBTRS: Band LU backsolve
 */
extern void dgbtrs_ (
    const char   *trans,
    const blas_int_t   *n,
    const blas_int_t   *kl,
    const blas_int_t   *ku, 
    const blas_int_t   *nrhs,
    const double *ab,
    const blas_int_t   *ldab,
    const blas_int_t   *ipiv,
    double       *b,
    const blas_int_t   *ldb,
    blas_int_t         *info,
    blas_int_t          trans_len
);

#ifdef __cplusplus
}
#endif

#endif // GLUE_LAPACK_H_INCLUDED_
