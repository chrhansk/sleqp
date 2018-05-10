/** @file anylapack.h
 * @author     Christian Kirches
 * @author     $LastChangedBy: ckirches $
 * @date       2009 Feb 6
 * @date       $LastChangedDate: 2010-02-08 15:51:05 +0100 (Mon, 08 Feb 2010) $
 *
 * $Id: anylapack.h 867 2010-02-08 14:51:05Z ckirches $
 *
 * @brief      LAPACK wrapper transparent to either C LAPACK or F77 LAPACK
 *
 *             Provides access to a subset of LAPACK only: Matrix decompositions
 *             (QR and Cholesky) and related utility functions
 *
 *             This file is part of proprietary software owned by the Simulation
 *             & Optimization Group, Interdisciplinary Center for Scientific
 *             Computing, Ruprecht-Karls-University of Heidelberg, Germany.
 *             Copyright (C) 2009. All Rights Reserved.
 */

#ifndef ANY_LAPACK_H_
#define ANY_LAPACK_H_


////////////////////////////////////////////////////////////////////////////////
// includes //
////////////////////////////////////////////////////////////////////////////////

#include "anyblas.h"


////////////////////////////////////////////////////////////////////////////////
// defines //
////////////////////////////////////////////////////////////////////////////////

#ifndef LAPACK_BLOCKSIZE

/** Maximum size for unblocked decompositions. Larger matrices use blocked
 *  decompositions, which are faster if a tuned DGEMM routine is available.
 *  The value 64 is recommended. */
#define LAPACK_BLOCKSIZE 64

#endif


////////////////////////////////////////////////////////////////////////////////
// LAPACK selection //
////////////////////////////////////////////////////////////////////////////////

#if !defined(HAVE_CLAPACK) && !defined(HAVE_F77LAPACK)
    #error "Fatal error: Neither HAVE_CLAPACK nor HAVE_F77LAPACK set !"
#endif

#ifdef HAVE_CLAPACK

    #include "clapack.h"

    #define CLAPACKNAME(name) clapack_##name

    /** Selector for lower or upper triangular matrices
     */
    enum lapack_UpLo {
        lapackLower   = CblasLower,
        lapackUpper   = CblasUpper
    };

    /** Selector for transposed or non-transposed matrices
     */
    enum lapack_Trans {
        lapackNoTrans = CblasNoTrans,
        lapackTrans   = CblasTrans
    };

    /** Selector for left-hand or right-hand side multiplication
     */
    enum lapack_Side {
        lapackLeft    = CblasLeft,
        lapackRight   = CblasRight
    };

    /** Selector for unit or non-unit diagonals of matrices
     */
    enum lapack_Diag {
        lapackNonUnit = CblasNonUnit,
        lapackUnit    = CblasUnit
    };

#endif // HAVE_CLAPACK

#ifdef HAVE_F77LAPACK

    #include "f77lapack.h"

    #ifndef F77NAME
        #define F77NAME(name)   name##_
    #endif // F77NAME

    /** Selector for lower or upper triangular matrices
     */
    enum lapack_UpLo {
        lapackLower   = 'L',
        lapackUpper   = 'U'
    };

    /** Selector for transposed or non-transposed matrices
     */
    enum lapack_Trans {
        lapackNoTrans = 'N',
        lapackTrans   = 'T'
    };

    /** Selector for left-hand or right-hand side multiplication
     */
    enum lapack_Side {
        lapackLeft    = 'L',
        lapackRight   = 'R'
    };

    /** Selector for unit or non-unit diagonals of matrices
     */
    enum lapack_Diag {
        lapackNonUnit = 'N',
        lapackUnit    = 'U'
    };

#endif // HAVE_F77LAPACK



////////////////////////////////////////////////////////////////////////////////
// functions //
////////////////////////////////////////////////////////////////////////////////

/** \brief Compute a QR decomposition of A
 *
 *  On call, A holds the matrix to be decomposed. On return, the upper triangle
 *  of A including the diagonal holds the factor R. The first k columns of the
 *  lower triangle of A hold the elementary reflectors \f$v_i\f$, where
 *  \f$k=\min(m,n)\f$. The vector tau holds the scalar factors \f$\tau_i\f$ of
 *  the elementary reflectors.
 *
 *  The factor Q is then represented as
 *  \f$Q:=\prod_{i=1}^k\left(I-\tau_i v_i v_i^T\right)\f$
 *  where \f$v_i:=\left(0,\ldots,0, A_{i,i+1},\ldots, A_{i,m}\right)\f$
 *
 *  \param m   Number of rows of the matrix A
 *  \param n   Number of columns of the matrix A
 *  \param A   Pointer to the upper left element of the matrix A
 *  \param lda Leading dimension of the matrix A
 *  \param tau Pointer to the first element of the scalar factors vector tau
 *  \return    0 upon success, < 0 indicates the number of the invalid argument
 */
INLINE long
lapack_dgeqrf (
    const long  m,
    const long  n,
    double     *A,
    const long  lda,
    double     *tau
)
{
    // local variables
    long    info  = 0;
    long    lWork = (m > 0) ? (m * LAPACK_BLOCKSIZE) : 1l;
    double *pWork = (double *) alloca (lWork * sizeof (double));

    F77NAME(dgeqrf)(&m, &n, A, &lda, tau, pWork, &lWork, &info);
    return info;
}


////////////////////////////////////////////////////////////////////////////////

/** \brief Compute an RQ decomposition of A.
 *  \see   ::lapack_dgeqrf
 */
INLINE long
lapack_dgerqf (
    const long  m,
    const long  n,
    double     *A,
    const long  lda,
    double     *tau
)
{
    // local variables
    long    info  = 0;
    long    lWork = (m > 0) ? (m * LAPACK_BLOCKSIZE) : 1l;
    double *pWork = (double *) alloca (lWork * sizeof (double));

    F77NAME(dgerqf)(&m, &n, A, &lda, tau, pWork, &lWork, &info);
    return info;
}


////////////////////////////////////////////////////////////////////////////////

/** \brief Compute an QL decomposition of A.
 *  \see   ::lapack_dgeqrf
 */
INLINE long
lapack_dgeqlf (
    const long  m,
    const long  n,
    double     *A,
    const long  lda,
    double     *tau
)
{
    // local variables
    long    info  = 0;
    long    lWork = (m > 0) ? (m * LAPACK_BLOCKSIZE) : 1l;
    double *pWork = (double *) alloca (lWork * sizeof (double));

    F77NAME(dgeqlf)(&m, &n, A, &lda, tau, pWork, &lWork, &info);
    return info;
}


////////////////////////////////////////////////////////////////////////////////

/** \brief Compute an LQ decomposition of A.
 *  \see   ::lapack_dgeqrf
 */
INLINE long
lapack_dgelqf (
    const long  m,
    const long  n,
    double     *A,
    const long  lda,
    double     *tau
)
{
    // local variables
    long    info  = 0;
    long    lWork = (m > 0) ? (m * LAPACK_BLOCKSIZE) : 1l;
    double *pWork = (double *) alloca (lWork * sizeof (double));

    F77NAME(dgelqf)(&m, &n, A, &lda, tau, pWork, &lWork, &info);
    return info;
}


////////////////////////////////////////////////////////////////////////////////

/** \brief Explicitly compute the unitary matrix Q after a QR decomposition.
 *
 *  On call A, holds the upper triangle R and the elementary reflectors v
 *  as obtained from a call to ::lapack_dgeqrf. On return, A holds the QR
 *  decomposition factor Q which has been computed from
 *  \f$Q:=\prod_{i=1}^k\left(I-\tau_i v_i v_i^T\right)\f$
 *
 *  \param m   Number of rows of the matrix A
 *  \param n   Number of columns of the matrix A
 *  \param A   Pointer to the upper left element of the matrix A
 *             as obtained from a call to ::lapack_dgeqrf
 *  \param lda Leading dimension of the matrix A
 *  \param tau Pointer to the first element of the scalar factors vector tau
 *             as obtained from a call to ::lapack_dgeqrf
 *  \return    0 upon success, < 0 indicates the number of the invalid argument
 */
INLINE long
lapack_dorgqr (
    const long    m,
    const long    n,
    const long    k,
    double       *A,
    const long    lda,
    const double *tau
)
{
    // local variables
    long    info  = 0;
    long    lWork = (n > 0) ? (n * LAPACK_BLOCKSIZE) : 1l;
    double *pWork = (double *) alloca (lWork * sizeof (double));

    F77NAME(dorgqr)(&m, &n, &k, A, &lda, tau, pWork, &lWork, &info);
    return info;
}


////////////////////////////////////////////////////////////////////////////////

/** \brief Explicitly compute the unitary matrix Q after an RQ decomposition.
 *  \see   ::lapack_dorgqr
 */
INLINE long
lapack_dorgrq (
    const long    m,
    const long    n,
    const long    k,
    double       *A,
    const long    lda,
    const double *tau
)
{
    // local variables
    long    info  = 0;
    long    lWork = (n > 0) ? (n * LAPACK_BLOCKSIZE) : 1l;
    double *pWork = (double *) alloca (lWork * sizeof (double));

    F77NAME(dorgrq)(&m, &n, &k, A, &lda, tau, pWork, &lWork, &info);
    return info;
}


////////////////////////////////////////////////////////////////////////////////

/** \brief Explicitly compute the unitary matrix Q after an QL decomposition.
 *  \see   ::lapack_dorgqr
 */
INLINE long
lapack_dorgql (
    const long    m,
    const long    n,
    const long    k,
    double       *A,
    const long    lda,
    const double *tau
)
{
    // local variables
    long    info  = 0;
    long    lWork = (n > 0) ? (n * LAPACK_BLOCKSIZE) : 1l;
    double *pWork = (double *) alloca (lWork * sizeof (double));

    F77NAME(dorgql)(&m, &n, &k, A, &lda, tau, pWork, &lWork, &info);
    return info;
}


////////////////////////////////////////////////////////////////////////////////

/** \brief Explicitly compute the unitary matrix Q after an LQ decomposition.
 *  \see   ::lapack_dorgqr
 */
INLINE long
lapack_dorglq (
    const long    m,
    const long    n,
    const long    k,
    double       *A,
    const long    lda,
    const double *tau
)
{
    // local variables
    long    info  = 0;
    long    lWork = (n > 0) ? (n * LAPACK_BLOCKSIZE) : 1l;
    double *pWork = (double *) alloca (lWork * sizeof (double));

    F77NAME(dorglq)(&m, &n, &k, A, &lda, tau, pWork, &lWork, &info);
    return info;
}


////////////////////////////////////////////////////////////////////////////////

/** \brief Multiply a matrix with the unitary matrix Q obtained from a QR
 *         decomposition
 *
 *  On call, C holds the matrix onto which Q is to be multiplied. On return,
 *  C holds the product \f$CQ\f$, \f$CQ^T\f$, \f$QC\f$, or \f$Q^TC\f$.
 *
 *  \param side  Whether to apply the matrix from the left or from the right
 *  \param trans Whether to apply \f$Q\f$ or \f$Q^T\f$
 *  \param m     Number of rows of the matrix A
 *  \param n     Number of columns of the matrix A
 *  \param A     Pointer to the upper left element of the matrix A
 *               as obtained from a call to ::lapack_dgeqrf
 *  \param lda   Leading dimension of the matrix A
 *  \param tau   Pointer to the first element of the scalar factors vector tau
 *               as obtained from a call to ::lapack_dgeqrf
 *  \param C     Pointer to the upper left element of the matrix C onto which Q
 *               is to be multiplied.
 *  \param ldc   Leading dimension of the matrix C.
 *  \return      0 upon success, < 0 indicates the number of the invalid argument.
 */
INLINE long
lapack_dormqr (
    const enum lapack_Side   side,
    const enum lapack_Trans  trans,
    const long               m,
    const long               n,
    const long               k,
    const double            *A,
    const long               lda,
    const double            *tau,
    double                  *C,
    const long               ldc
)
{
    // local variables
    long    info  = 0;
    long    l     = (side == lapackLeft) ? n : m;
    long    lWork = (l > 0) ? (l * LAPACK_BLOCKSIZE) : 1l;
    double *pWork = (double *) alloca (lWork * sizeof (double));

    F77NAME(dormqr)((const char *)&side, (const char *)&trans, &m, &n, &k,
                A, &lda, tau, C, &ldc, pWork, &lWork, &info, 1, 1);
    return info;
}


////////////////////////////////////////////////////////////////////////////////

INLINE long
lapack_dormrq (
    const enum lapack_Side   side,
    const enum lapack_Trans  trans,
    const long               m,
    const long               n,
    const long               k,
    const double            *A,
    const long               lda,
    const double            *tau,
    double                  *C,
    const long               ldc
)
{
    // local variables
    long    info  = 0;
    long    l     = (side == lapackLeft) ? n : m;
    long    lWork = (l > 0) ? (l * LAPACK_BLOCKSIZE) : 1l;
    double *pWork = (double *) alloca (lWork * sizeof (double));

    F77NAME(dormrq)((const char *)&side, (const char *)&trans, &m, &n, &k,
                A, &lda, tau, C, &ldc, pWork, &lWork, &info, 1, 1);
    return info;
}


////////////////////////////////////////////////////////////////////////////////

/** \brief Computes the Cholesky (\f$LL^T\f$ or \f$R^TR\f$) decomposition of a
 *         symmetric positive definite matrix
 *
 *  Upon call, A holds the upper or lower triangle of the symmetric matrix to
 *  be decomposed. Upon return, A holds the triangular factor R if the upper
 *  triangle of A was given, and the triangular factor L if the lower triangle
 *  of A was given.
 *  \param uplo Whether to compute the lower triangular factor L from the lower
 *              triangle of A, or the upper triangular factor R from the upper
 *              triangle of A.
 *  \param n    Number of rows and columns of the square matrix A
 *  \param A    Pointer to the upper left element of the matrix A
 *  \param lda  Leading dimension of the matrix A
 *  \return     0 upon success,
 *              < 0 indicates the number of the invalid argument
 *              > 0 indicates the diagonal position of the minor that is not
 *                  positive definite
 */
INLINE long
lapack_dpotrf (
    const enum lapack_UpLo  uplo,
    const long              n,
    double                 *A,
    const long              lda
)
{
    // local variables
    long info = 0;

    F77NAME(dpotrf)((const char *) &uplo, &n, A, &lda, &info, 1);
    return info;
}


////////////////////////////////////////////////////////////////////////////////

/** \brief Solves a general triangular system of linear equations
 *
 *  Upon call, B holds the linear equation system's right hand side. On return,
 *  B holds the solution \f$X=A^ {-1}B\f$.
 *
 *  \param uplo   Whether A is an upper or lower triangular matrix
 *  \param trans  Whether to solve \f$AX=B\f$ or \f$A^TX=B\f$
 *  \param diag   Whether or not A has a unit diagonal
 *  \param n      Number of rows and columns of A, number of rows of X and B
 *  \param nrhs   Number of columns of the matrices X and B
 *  \param A      Pointer to the upper left element of the matrix A
 *  \param lda    Leading dimension of the matrix A
 *  \param B      Pointer to the upper left element of the matrix B
 *  \param ldb    leading dimension of the matrix B
 *  \return       0 upon success
 *                <0 indicates number of invalid argument
 *                >0 indicates position of the first zero diagonal element of A
 */
INLINE long
lapack_dtrtrs (
    const enum lapack_UpLo   uplo,
    const enum lapack_Trans  trans,
    const enum lapack_Diag   diag,
    const long               n,
    const long               nrhs,
    const double            *A,
    const long               lda,
    double                  *B,
    const long               ldb
)
{
    // local variables
    long info = 0;

    F77NAME(dtrtrs)((const char *)&uplo, (const char *)&trans,
                (const char *)&diag, &n, &nrhs, A, &lda, B, &ldb, &info,
                1, 1, 1);
    return info;
}



////////////////////////////////////////////////////////////////////////////////

/** \brief Compute a condition number estimate from a Cholesky factor of a
 *         symmetric positive definite matrix.
 *  \param uplo   Whether A is an upper or lower triangular matrix
 *  \param n      Number of rows and columns of A, number of rows of X and B
 *  \param A      Pointer to the upper left element of the matrix A
 *  \param lda    Leading dimension of the matrix A
 *  \param rcond  Pointer to the output condition number estimate.
 *  \return       0 upon success
 *                <0 indicates number of invalid argument
 */
INLINE long
lapack_dpocon (
    const enum lapack_UpLo  uplo,
    const long              n,
    const double           *A,
    const long              lda,
    double                 *rcond
)
{
    // local variables
    double *work  = (double *) alloca (3 * n * sizeof (double));
    long   *iwork = (long *) alloca (n * sizeof (long));
    double anorm  = 0.0;
    long   info   = 0;
    long   ii, jj;

    // compute 1-norm of symmetric matrix A
    if (uplo == lapackUpper) {
        for (jj = 0; jj < n; ++jj) {    // column j
            double colnorm = 0.0;
            for (ii = 0; ii <= jj; ++ii)    // row i
                colnorm += fabs (A[ii+lda*jj]);
            for (ii = jj+1; ii < n; ++ii)
                colnorm += fabs (A[jj*lda+ii]);
            if (anorm < colnorm)
                anorm = colnorm;
        }
    } else {    // uplo == lapackLower
        for (jj = 0; jj < n; ++jj) {    // column j
            double colnorm = 0.0;
            for (ii = 0; ii <= jj; ++ii)    // row i
                colnorm += fabs (A[jj+lda*ii]);
            for (ii = jj+1; ii < n; ++ii)
                colnorm += fabs (A[ii*lda+jj]);
            if (anorm < colnorm)
                anorm = colnorm;
        }
    }

    // compute condition number estimate
    F77NAME(dpocon)((const char *)&uplo, &n, A, &lda, &anorm, rcond, work,
                    iwork, &info, 1);
    return info;
}


////////////////////////////////////////////////////////////////////////////////

/** \brief Compute the LU decomposition of a banded matrix */
INLINE long
lapack_dgbtrf (
    const long  m,
    const long  n,
    const long  kl,
    const long  ku,
    double     *ab,
    const long  ldab,
    long       *ipiv
) {
    // local variables
    long info = 0;
        
    F77NAME (dgbtrf)(&m, &n, &kl, &ku, ab, &ldab, ipiv, &info);
    return info;
}


////////////////////////////////////////////////////////////////////////////////

/** \brief Backsolve with a LU decomposition of a banded matrix */
INLINE long
lapack_dgbtrs (
    const enum lapack_Trans  trans,
    const long               n,
    const long               kl,
    const long               ku,    
    const long               nrhs,
    const double            *ab,
    const long               ldab,
    const long              *ipiv,
    double                  *b,
    const long               ldb      
) {
    // local variables
    long info = 0;    
    
    F77NAME (dgbtrs)((const char *) &trans, &n, &kl, &ku, &nrhs, ab, &ldab,
                     ipiv, b, &ldb, &info, 1);
    return info;
} 


////////////////////////////////////////////////////////////////////////////////

#endif /* ANY_LAPACK_H_ */

/* end of file */
