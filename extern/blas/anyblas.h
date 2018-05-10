/** 
 * @file anyblas.h
 * @author     Christian Kirches
 * @author     $LastChangedBy: ckirches $
 * @date       2008 Oct 29
 * @date       $LastChangedDate: 2010-02-11 19:21:25 +0100 (Thu, 11 Feb 2010) $
 *
 * $Id: anyblas.h 871 2010-02-11 18:21:25Z ckirches $
 *
 * @brief      BLAS wrapper transparent to either C BLAS or F77 BLAS
 *
 *             This file is part of proprietary software owned by the Simulation
 *             & Optimization Group, Interdisciplinary Center for Scientific
 *             Computing, Ruprecht-Karls-University of Heidelberg, Germany.
 *             Copyright (C) 2008--2009. All Rights Reserved.
 */

#ifndef ANY_BLAS_H_
#define ANY_BLAS_H_


////////////////////////////////////////////////////////////////////////////////
// BLAS selection //
////////////////////////////////////////////////////////////////////////////////

#if !defined(HAVE_CBLAS) && !defined(HAVE_F77BLAS)
    #error "Fatal error: Neither HAVE_CBLAS nor HAVE_F77BLAS is set."
#endif

#ifdef HAVE_CBLAS

    #include "cblas.h"

    #define CBLASNAME(name) cblas_##name

    /** Selector for lower or upper triangular matrices
     */
    enum blas_UpLo {
        blasLower   = CblasLower,  /**< Matrix is lower triangular */
        blasUpper   = CblasUpper   /**< Matrix is upper triangular */
    };

    /** Selector for transposed or non-transposed matrices
     */
    enum blas_Trans {
        blasNoTrans = CblasNoTrans,  /**< Use untransposed matrix */
        blasTrans   = CblasTrans     /**< Use transposed matrix */
    };

    /** Selector for left-hand or right-hand side multiplication
     */
    enum blas_Side {
        blasLeft    = CblasLeft,  /** Multiply from the left */
        blasRight   = CblasRight  /** Multiply from the right */
    };

    /** Selector for unit or non-unit diagonals of matrices
     */
    enum blas_Diag {
        blasNonUnit = CblasNonUnit, /** Matrix has arbitrary diagonal */
        blasUnit    = CblasUnit     /** Matrix has a unit diagonal */
    };

#endif // HAVE_CBLAS

#if defined(HAVE_F77BLAS) && !defined(HAVE_CBLAS)

    #include "f77blas.h"

    /** Selector for lower or upper triangular matrices
     */
    enum blas_UpLo {
        blasLower   = 'L',  /**< Matrix is lower triangular */
        blasUpper   = 'U'   /**< Matrix is upper triangular */
    };

    /** Selector for transposed or non-transposed matrices
     */
    enum blas_Trans {
        blasNoTrans = 'N',  /**< Use untransposed matrix */
        blasTrans   = 'T'   /**< Use transposed matrix */
    };

    /** Selector for left-hand or right-hand side multiplication
     */
    enum blas_Side {
        blasLeft    = 'L', /** Multiply from the left */
        blasRight   = 'R'  /** Multiply from the right */
    };

    /** Selector for unit or non-unit diagonals of matrices
     */
    enum blas_Diag {
        blasNonUnit = 'N', /** Matrix has arbitrary diagonal */
        blasUnit    = 'U'  /** Matrix has a unit diagonal */
    };

#endif // HAVE_F77BLAS && !HAVE_CBLAS


////////////////////////////////////////////////////////////////////////////////
// CBLAS style inline functions //
////////////////////////////////////////////////////////////////////////////////

/** \brief Find the index of the vector element with largest magnitude, i.e,
 *         compute \f$i:=\text{argmax}\{x_i~|~i=1,\ldots,n\}\f$
 *  \param n  Length of the vector x
 *  \param x  Pointer to the first element of the vector x
 *  \param ix Stride of the pointer \a x
 *  \return   Index of the element with largest magnitude. This is a base-1
 *            index, i.e., the first element has index 1, the last element has
 *            index n. A result of 0 indicates an error.
 *  \warning  This call always uses the Fortran 77 base-1 indices even if the C
 *            BLAS interface is used. While it would be more convenient to make
 *            blas_idamax always use base-0 indices, it would be a major hassle
 *            to adapt QPOPT to that.
 */
INLINE long
blas_idamax (
    const long    n,
    const double *x,
    const long    ix
)
{
#ifdef HAVE_CBLAS
    return CBLASNAME(idamax) (n, x, ix);
#else
    return F77NAME(idamax) (&n, x, &ix);
#endif // HAVE_CBLAS
}


////////////////////////////////////////////////////////////////////////////////

/** \brief Copy one vector onto another, i.e.,
 *         compute \f$y:=x\f$
 *  \param n  Length of the vectors x and y
 *  \param x  Pointer to the first element of the vector x
 *  \param ix Stride of the pointer \a x
 *  \param y  Pointer to the first element of the vector y
 *  \param iy Stride of the pointer \a y
 *  \return   Always zero.
 */
INLINE long
blas_dcopy (
    const long    n,
    const double *x,
    const long    ix,
    double       *y,
    const long    iy
)
{
#ifdef HAVE_CBLAS
    CBLASNAME(dcopy) (n, x, ix, y, iy);
#else
    F77NAME(dcopy) (&n, x, &ix, y, &iy);
#endif // HAVE_CBLAS
    return 0;
}


////////////////////////////////////////////////////////////////////////////////

/** \brief Swap two vectors x and y, i.e., compute \f$x\longleftrightarrow y\f$
 *  \param n  Length of the vectors x and y
 *  \param x  Pointer to the first element of the vector x
 *  \param ix Stride of the pointer \a x
 *  \param y  Pointer to the first element of the vector y
 *  \param iy Stride of the pointer \a y
 *  \return   Always zero.
 */
INLINE long
blas_dswap (
    const long    n,
    double       *x,
    const long    ix,
    double       *y,
    const long    iy
)
{
#ifdef HAVE_CBLAS
    CBLASNAME(dswap) (n, x, ix, y, iy);
#else
    F77NAME(dswap) (&n, x, &ix, y, &iy);
#endif // HAVE_CBLAS
    return 0;
}


////////////////////////////////////////////////////////////////////////////////

/** \brief Compute a vector multiplication-addition \f$y := \alpha x + y \f$
 *  \param n     Length of the vectors x and y
 *  \param alpha Scalar multiplier for the  vector x
 *  \param x     Pointer to the first element of the vector x
 *  \param ix    Stride of the pointer \a x
 *  \param y     Pointer to the first element of the vector y
 *  \param iy    Stride of the pointer \a y
 *  \return      Always zero.
 */
INLINE long
blas_daxpy (
    const long    n,
    const double  alpha,
    const double *x,
    const long    ix,
    double       *y,
    const long    iy
)
{
#ifdef HAVE_CBLAS
    CBLASNAME(daxpy) (n, alpha, x, ix, y, iy);
#else
    F77NAME(daxpy) (&n, &alpha, x, &ix, y, &iy);
#endif // HAVE_CBLAS
    return 0;
}


////////////////////////////////////////////////////////////////////////////////

/** \brief Compute a vector multiplication \f$x := \alpha x \f$
 *  \param n     Length of the vector x
 *  \param alpha Scalar multiplier for the vector x
 *  \param x     Pointer to the first element of the vector x
 *  \param ix    Stride of the pointer \a x
 *  \return      Always zero.
 */
INLINE long
blas_dscal (
    const long    n,
    const double  alpha,
    double       *x,
    const long    ix
)
{
#ifdef HAVE_CBLAS
    CBLASNAME(dscal) (n, alpha, x, ix);
#else
    F77NAME(dscal) (&n, &alpha, x, &ix);
#endif // HAVE_CBLAS
    return 0;
}


////////////////////////////////////////////////////////////////////////////////

/** \brief Compute a dot product \f$r:=x^Ty\f$
 *  \param n     Length of the vectors x and y
 *  \param x     Pointer to the first element of the vector x
 *  \param ix    Stride of the pointer \a x
 *  \param y     Pointer to the first element of the vector y
 *  \param iy    Stride of the pointer \a y
 *  \return      Dot product r of the two vectors.
 */
INLINE double
blas_ddot (
    const long    n,
    const double *x,
    const long    ix,
    const double *y,
    const long    iy
)
{
#ifdef HAVE_CBLAS
    return CBLASNAME(ddot) (n, x, ix, y, iy);
#else
    return F77NAME(ddot) (&n, x, &ix, y, &iy);
#endif // HAVE_CBLAS
}


////////////////////////////////////////////////////////////////////////////////

/** \brief Compute \f$||x||_2=(x^Tx)^{\frac{1}{2}}\f$
 *  \param n     Length of the vector x
 *  \param x     Pointer to the first element of the vector x
 *  \param ix    Stride of the pointer \a x
 *  \return      Euclidean length \f$||x||_2\f$ of the vector x
 */
INLINE double
blas_dnrm2 (
    const long    n,
    const double *x,
    const long    ix
)
{
#ifdef HAVE_CBLAS
    return CBLASNAME(dnrm2) (n, x, ix);
#else
    return F77NAME(dnrm2) (&n, x, &ix);
#endif // HAVE_CBLAS
}


////////////////////////////////////////////////////////////////////////////////

/** \brief Compute \f$||x||_1=\sum_{i}|x_i|\f$
 *  \param n     Length of the vector x
 *  \param x     Pointer to the first element of the vector x
 *  \param ix    Stride of the pointer \a x
 *  \return      1-norm \f$||x||_1\f$ of the vector x
 */
INLINE double
blas_dasum (
    const long    n,
    const double *x,
    const long    ix
)
{
#ifdef HAVE_CBLAS
    return CBLASNAME(dasum) (n, x, ix);
#else
    return F77NAME(dasum) (&n, x, &ix);
#endif // HAVE_CBLAS
}


////////////////////////////////////////////////////////////////////////////////

INLINE long
blas_drotg (
    double *a,
    double *b,
    double *c,
    double *s
) {
#ifdef HAVE_CBLAS
    CBLASNAME(drotg) (a, b, c, s);
#else
    F77NAME(drotg) (a, b, c, s);
#endif // HAVE_CBLAS
    return 0;
}


////////////////////////////////////////////////////////////////////////////////

INLINE long
blas_drot (
    const long    n,
    double       *x,
    const long    ix,
    double       *y,
    const long    iy,
    const double *c,
    const double *s
) {
#ifdef HAVE_CBLAS
    CBLASNAME(drot) (n, x, ix, y, iy, c, s);
#else
    F77NAME(drot) (&n, x, &ix, y, &iy, c, s);
#endif // HAVE_CBLAS
    return 0;
}


////////////////////////////////////////////////////////////////////////////////

/** \brief Compute \f$y := \alpha \text{op}(A) x + \beta y \f$
 *         where \f$\text{op}(A)=A\f$ or \f$\text{op}(A)=A^T\f$
 *  \param trans  'N' to use \f$\text{op}(A)=A\f$,
 *                'T' to use \f$\text{op}(A)=A^T\f$
 *  \param m      Number of rows in A
 *  \param n      Number of columns in A
 *  \param alpha  Scalar multiplier for A
 *  \param a      Pointer to the upper left element of A
 *  \param lda    Leading dimension of A
 *  \param x      Pointer to the first element of x
 *  \param ix     Stride of the pointer \a x
 *  \param beta   Scalar multiplier for \a y
 *  \param y      Pointer to the first element of y
 *  \param iy     Stride of the pointer \a y
 *  \return       Always zero.
 */
INLINE long
blas_dgemv (
    const enum blas_Trans  trans,
    const long             m,
    const long             n,
    const double           alpha,
    const double          *a,
    const long             lda,
    const double          *x,
    const long             ix,
    const double           beta,
    double                *y,
    const long             iy
)
{
#ifdef HAVE_CBLAS
    CBLASNAME (dgemv) (CblasColMajor, trans, m, n, alpha, a, lda, x, ix,
               beta, y, iy);
#else
    F77NAME(dgemv) ((const char *) &trans, &m, &n, &alpha, a, &lda, x, &ix,
            &beta, y, &iy, 1);
#endif // HAVE_CBLAS
    return 0;
}



////////////////////////////////////////////////////////////////////////////////

/** \brief Compute \f$y := \text{op}(A) x\f$ where A is a triangular matrix
 *         and \f$\text{op}(A)=A\f$ or \f$\text{op}(A)=A^T\f$
 *  \param uplo  Whether A is an upper or lower triangular matrix
 *  \param trans Whether to use \f$\text{op}(A)=A\f$ or \f$\text{op}(A)=A^T\f$
 *  \param diag  Whether A has a unit or non-uint diagonal
 *  \param n     Number of elements in \a x
 *  \param a     Pointer to the upper left element of A
 *  \param lda   Leading dimension of A
 *  \param x     Pointer to the first element of x
 *  \param incx  Stride of the pointer \a x
 *  \return      Always zero.
 */
INLINE long
blas_dtrmv (
    const enum blas_UpLo  uplo,
    const enum blas_Trans trans,
    const enum blas_Diag  diag,
    const long            n,
    const double         *a,
    const long            lda,
    double               *x,
    const long            incx
)
{
#ifdef HAVE_CBLAS
    CBLASNAME(dtrmv) (CblasColMajor, uplo, trans, diag, n, a, lda, x, incx);
#else
    F77NAME(dtrmv) ((const char *) &uplo, (const char * )&trans,
            (const char *) &diag, &n, a, &lda, x, &incx, 1, 1, 1);
#endif // HAVE_CBLAS
    return 0;
}


////////////////////////////////////////////////////////////////////////////////

/** \brief Compute \f$y := \alpha A x + \beta y\f$
 *         where A is a symmetric matrix whose upper or lower triangle is given
 *  \param uplo   Whether the upper or lower triangle of A (symmetric) is given
 *  \param n      Length of the vector y
 *  \param alpha  Scalar factor \f$\alpha\f$ for the product Ax
 *  \param a      Pointer to the first element of the matrix A
 *  \param lda    Leading dimension of A
 *  \param x      Pointer to the first element of the vector x
 *  \param ix     Stride of the pointer \a x
 *  \param beta   Scalar factor \f$\beta\f$ for the vector y
 *  \param y      Pointer to the first element of the vector y
 *  \param iy     Stride of the pointer \a y
 *  \return       Always zero.
 */
INLINE long
blas_dsymv (
    const enum blas_UpLo  uplo,
    const long            n,
    const double          alpha,
    const double         *a,
    const long            lda,
    const double         *x,
    const long            ix,
    const double          beta,
    double               *y,
    const long            iy
)
{
#ifdef HAVE_CBLAS
    CBLASNAME(dsymv) (CblasColMajor, uplo, n, alpha, a, lda, x, ix, beta, y,
              iy);
#else
    F77NAME(dsymv) ((const char *) &uplo, &n, &alpha, a, &lda, x, &ix, &beta,
                y, &iy, 1);
#endif // HAVE_CBLAS
    return 0;
}

////////////////////////////////////////////////////////////////////////////////

/** \brief Solve \f$\text{op}(A)x=b\f$ where A is an upper or lower
 *         triangular matrix with unit or non-unit diagonal, and
 *         \f$\text{op}(A)=A\f$ or \f$\text{op}(A)=A^T\f$
 *  \param uplo   Whether A is upper or lower triangular
 *  \param transa Whether \f$\text{op}(A)=A\f$ or \f$\text{op}(A)=A^T\f$
 *  \param diag   Whether A has a unit or non-unit diagonal
 *  \param n      Length of the vectors x and b
 *  \param a      Pointer to the upper left element of the matrix A
 *  \param lda    Leading dimension of the matrix A
 *  \param x      Pointer to the first element of the vectors x and b
 *  \param ix     Stride of the pointer \a x
 *  \return       Always zero.
 */
INLINE long
blas_dtrsv (
    const enum blas_UpLo   uplo,
    const enum blas_Trans  transa,
    const enum blas_Diag   diag,
    const long             n,
    const double          *a,
    const long             lda,
    double                *x,
    const long             ix
)
{
#ifdef HAVE_CBLAS
    CBLASNAME(dtrsv) (CblasColMajor, uplo, transa, diag, n, a, lda, x, ix);
#else
    F77NAME(dtrsv) ((const char *) &uplo, (const char *) &transa,
            (const char * )&diag, &n, a, &lda, x, &ix, 1, 1, 1);
#endif // HAVE_CBLAS
    return 0;
}


////////////////////////////////////////////////////////////////////////////////

/** \brief Compute \f$A := \alpha x y^T\f$
 *  \param m      Number of rows in A and length of the vector x
 *  \param n      Number of columns in A and length of the vector y
 *  \param alpha  Scalar factor \f$\alpha\f$ for the dyadic product \f$xy^T\f$
 *  \param x      Pointer to the first element of the vector x
 *  \param ix     Stride of the pointer \a x
 *  \param y      Pointer to the first element of the vector y
 *  \param iy     Stride of the pointer \a y
 *  \param a      Pointer to the upper left element of the matrix A
 *  \param lda    Leading dimension of the matrix A
 *  \return       Always zero.
 */
INLINE long
blas_dger (
    const long    m,
    const long    n,
    const double  alpha,
    const double *x,
    const long    ix,
    const double *y,
    const long    iy,
    double       *a,
    const long    lda
)
{
#ifdef HAVE_CBLAS
    CBLASNAME(dger) (CblasColMajor, m, n, alpha, x, ix, y, iy, a, lda);
#else
    F77NAME(dger) (&m, &n, &alpha, x, &ix, y, &iy, a, &lda);
#endif // HAVE_CBLAS
    return 0;
}


////////////////////////////////////////////////////////////////////////////////

/** \brief Compute \f$A := \alpha x x^T\f$ where A is a symmetric matrix whose
 *         upper or lower triangle is given
 *  \param uplo   Whether the upper or lower triangle of the matrix A is given
 *  \param n      Number of rows and columns in A and length of the vector x
 *  \param alpha  Scalar factor \f$\alpha\f$ for the dyadic product \f$xx^T\f$
 *  \param x      Pointer to the first element of the vector x
 *  \param ix     Stride of the pointer \a x
 *  \param a      Pointer to the upper left element of the matrix A
 *  \param lda    Leading dimension of the matrix A
 *  \return       Always zero.
 */
INLINE long
blas_dsyr (
    const enum blas_UpLo  uplo,
    const long            n,
    const double          alpha,
    const double         *x,
    const long            ix,
    double               *a,
    const long            lda
)
{
#ifdef HAVE_CBLAS
    CBLASNAME(dsyr) (CblasColMajor, uplo, n, alpha, x, ix, a, lda);
#else
    F77NAME(dsyr) ((const char *) &uplo, &n, &alpha, x, &ix, a, &lda, 1);
#endif // HAVE_CBLAS
    return 0;
}


////////////////////////////////////////////////////////////////////////////////

/** \brief Compute \f$C := \alpha \text{op}(A)\text{op}(B) + \beta C\f$
 *         where \f$\text{op}(A)=A\f$ or \f$\text{op}(A)=A^T\f$
 *  \param transa Whether \f$\text{op}(A)=A\f$ or \f$\text{op}(A)=A^T\f$
 *  \param transb Whether \f$\text{op}(B)=B\f$ or \f$\text{op}(B)=B^T\f$
 *  \param m      Number of rows in C and in \f$\text{op}(A)\f$
 *  \param n      Number of columns in C and in \f$\text{op}(B)\f$
 *  \param k      Number of columns in \f$\text{op}(A)\f$
 *                and of rows in \f$\text{op}(B)\f$
 *  \param alpha  Scalar factor \f$\alpha\f$ for the product
 *                \f$\text{op}(A)\text{op}(B)\f$
 *  \param a      Pointer to the upper left element of the matrix A
 *  \param lda    Leading dimension of the matrix A
 *  \param b      Pointer to the upper left element of the matrix B
 *  \param ldb    Leading dimension of the matrix B
 *  \param beta   Scalar factor \f$\beta\f$ for the matrix C
 *  \param c      Pointer to the upper left element of the matrix C
 *  \param ldc    Leading dimension of the matrix C
 *  \return       Always zero
 */
INLINE long
blas_dgemm (
    const enum blas_Trans  transa,
    const enum blas_Trans  transb,
    const long             m,
    const long             n,
    const long             k,
    const double           alpha,
    const double          *a,
    const long             lda,
    const double          *b,
    const long             ldb,
    const double           beta,
    double                *c,
    const long             ldc
)
{
#ifdef HAVE_CBLAS
    CBLASNAME(dgemm) (CblasColMajor, transa, transb, m, n, k, alpha, a, lda,
              b, ldb, beta, c, ldc);
#else
    F77NAME(dgemm) ((const char *) &transa, (const char *) &transb, &m, &n,
                &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc, 1, 1);
#endif // HAVE_CBLAS
    return 0;
}


////////////////////////////////////////////////////////////////////////////////

/** \brief Compute \f$C := \alpha\text{op}(A,B) + \beta C\f$
 *  \param side   Whether \f$\text{op}(A,B)=AB\f$ or \f$\text{op}(A,B)=BA\f$
 *  \param uplo   Whether the upper or lower triangle of A is given
 *  \param m      Number of rows in the matrix C
 *  \param n      Number of columns in the matrix C
 *  \param alpha  Scalar factor \f$\alpha\f$ for the product \f$\text{op}(A,B)\f$
 *  \param a      Pointer to the upper left element of the matrix A
 *  \param lda    Leading dimension of the matrix A
 *  \param b      Pointer to the upper left element of the matrix B
 *  \param ldb    Leading dimension of the matrix B
 *  \param beta   Scalar factor \f$\beta\f$ for the matrix C
 *  \param c      Pointer to the upper left element of the matrix C
 *  \param ldc    Leading dimension of the matrix C
 *  \return       Always zero
 */
INLINE long
blas_dsymm (
    const enum blas_Side  side,
    const enum blas_UpLo  uplo,
    const long            m,
    const long            n,
    const double          alpha,
    const double         *a,
    const long            lda,
    const double         *b,
    const long            ldb,
    const double          beta,
    double               *c,
    const long            ldc
)
{
#ifdef HAVE_CBLAS
    CBLASNAME(dsymm) (CblasColMajor, side, uplo, m, n, alpha, a, lda, b,
              ldb, beta, c, ldc);
#else
    F77NAME(dsymm) ((const char *) &side, (const char * )&uplo, &m, &n,
             &alpha, a, &lda, b, &ldb, &beta, c, &ldc, 1, 1);
#endif // HAVE_CBLAS
    return 0;
}


////////////////////////////////////////////////////////////////////////////////

/** \brief Compute \f$C := \alpha\text{op}(A)\text{op(A)}^T + \beta C\f$
 *         where C is a symmetric matrix whose upper or lower triangle is given
 *  \param uplo   Whether the upper or lower triangle of A is given
 *  \param trans  Whether \f$\text{op}(A)=A\f$ or \f$\text{op}(A)=A^T\f$
 *  \param n      Number of rows and columns in the matrix C
 *  \param k      Number of columns in \f$\text{op}(A)\f$
 *  \param alpha  Scalar factor \f$\alpha\f$ for the product
 *                \f$\text{op}(A)\text{op}(A)^T\f$
 *  \param a      Pointer to the upper left element of the matrix A
 *  \param lda    Leading dimension of the matrix A
 *  \param beta   Scalar factor \f$\beta\f$ for the matrix C
 *  \param c      Pointer to the upper left element of the matrix C
 *  \param ldc    Leading dimension of the matrix C
 *  \return       Always zero
 */
INLINE long
blas_dsyrk (
    const enum blas_UpLo  uplo,
    const enum blas_Trans trans,
    const long            n,
    const long            k,
    const double          alpha,
    const double         *a,
    const long            lda,
    const double          beta,
    double               *c,
    const long            ldc
)
{
#ifdef HAVE_CBLAS
    CBLASNAME(dsyrk) (CblasColMajor, uplo, trans, n, k, alpha, a, lda, beta,
              c, ldc);
#else
    F77NAME(dsyrk) ((const char *) &uplo, (const char *) &trans, &n, &k,
             &alpha, a, &lda, &beta, c, &ldc, 1, 1);
#endif // HAVE_CBLAS
    return 0;
}


////////////////////////////////////////////////////////////////////////////////

/** \brief Solve \f$\text{op}(A)X=\alpha B\f$ or \f$X\text{op}A=\alpha B\f$
 *         where A is an upper or lower triangular matrix with unit or non-unit
 *         diagonal, and\f$\text{op}(A)=A\f$ or \f$\text{op}(A)=A^T\f$
 *  \param side   Side of X on which \f$\text{op}(A)\f$ appears
 *  \param uplo   Whether A is upper or lower triangular
 *  \param transa Whether \f$\text{op}(A)=A\f$ or \f$\text{op}(A)=A^T\f$
 *  \param diag   Whether A has a unit or non-unit diagonal
 *  \param m      Number of rows and columns of the matrix A
 *  \param n      Number of columns of the matrices X and B
 *  \param a      Pointer to the upper left element of the matrix A
 *  \param alpha  Scalar multiplier \f$\alpha\f$ for right hand side B
 *  \param lda    Leading dimension of the matrix A
 *  \param b      Pointer to the upper left element of the matrices X and B
 *  \param ldb    Leading dimension of the matrices X and B
 *  \return       Always zero.
 */
INLINE long
blas_dtrsm (
    const enum blas_Side  side,
    const enum blas_UpLo  uplo,
    const enum blas_Trans trans,
    const enum blas_Diag  diag,
    const long            m,
    const long            n,
    const double          alpha,
    const double         *a,
    const long            lda,
    double               *b,
    const long            ldb
)
{
#ifdef HAVE_CBLAS
    CBLASNAME(dtrsm) (CblasColMajor, side, uplo, trans, diag, m, n, alpha,
                  a, lda, b, ldb);
#else
    F77NAME(dtrsm) ((const char *) &side, (const char *) &uplo,
                (const char *) &trans, (const char *) &diag, &m, &n,
                &alpha, a, &lda, b, &ldb, 1, 1, 1, 1);
#endif // HAVE_CBLAS
    return 0;
}


////////////////////////////////////////////////////////////////////////////////

#endif /* ANY_BLAS_H_ */

/* end of file */
