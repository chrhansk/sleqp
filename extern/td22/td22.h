/** \file   hsl_td22.h
 *  \author Christian Kirches
 *  \author $LastChangedBy: ckirches $
 *  \date   2010 Mar 2
 *  \date   $LastChangedDate: 2009-04-29 17:20:12 +0200 (Mi, 29 Apr 2009) $
 *
 *  $Id: anylapack.h 623 2009-04-29 15:20:12Z ckirches $
 *
 *  \brief  Harwell Subroutine Library (HSL) TD22 wrapper for C
 *
 *  All documentation is abridged and taken from the HSL 2007 package
 *  specification documents for TD22.
 *
 *  References:
 *  -# Curtis, Powell, Reid. J. Inst. Math. Applics. 13 (1974), pp. 117--120
 *  -# Curtis, Reid. J. Inst. Math. Applics. 13 (1974), pp. 121--126
 */

#pragma once
#ifndef _GLUE_TD22_H_INCLUDED_
#define _GLUE_TD22_H_INCLUDED_

#ifdef __cplusplus
extern "C" {
#endif

#include "util_types.h"

/**
 * @brief      Sets default values for the components of the arrays that hold
 *             control parameters.
 *
 *             Normally the user will call TD22I/TD22ID prior to any call to
 *             TD22A/AD. If non- default values for any of the control
 *             parameters are required, they should be set immediately after the
 *             call to TD22I/ID.
 *
 * @param      icntl  an INTEGER array of length 5 that need not be set by the
 *                    user. On return it contains default values.
 * @param      cntl   a double PRECISION array of length 7 that need not be set
 *                    by the user. On return it contains default values.
 * @param      keep   an INTEGER array of length 20 that need not be set by the
 *                    user. It is used by TD22 as private workspace and must not
 *                    be altered by the user.
 * @param      rkeep  a double PRECISION array of length 10 that need not be
 *                    set by the user. It is used by TD22 as private workspace
 *                    and must not be altered by the user.
 *
 *             icntl(1) specifies the unit number to be used to output error
 *             messages. A negative number will suppress output. It has a
 *             default value of 6. icntl(2) has default value 2. This may be
 *             reset by the user to 1 to economise in the number of function
 *             calls. In that case, no estimation of errors is made, no
 *             adjustments are made to step sizes, and the less accurate
 *             one-sided difference approximation is used in place of the
 *             two-sided one.
 *
 *             info(1)  is used as an error indicator. It is set to zero when
 *             there are no errors otherwise it is set to a nonzero value. At
 *             present, only one error is possible, i.e. when @a ia is found to
 *             be too small by td22bd or td22cd; then info(1) is returned set to
 *             1, info(2) set to a suggested value for @a ia and an error
 *             message is output on unit icntl(1). For td22bd the calculation
 *             may be continued by copying the data in @a irn to an array of
 *             size at least info(2) and recalling the subroutine with the
 *             larger array, its size in @a ia, and @a iflag set to 2. For
 *             td22cd, the calculation will be successful if repeated with a new
 *             array @a irn of size info(2). info(2)  is used to return a
 *             supplementary value associated with the error indicated in
 *             info(1) info(3)  is only set by td22ad and returns the number of
 *             groups used.
 */
extern void td22id_ (
    fint_t   *icntl,
    double *cntl,
    fint_t   *keep,
    double *rkeep
);

/**
 * @brief      Calculates the approximate jacobian
 *
 * @param      m      an INTEGER variable set by the user to the number of
 *                    functions
 *             @f$f_i(x)\f$, that is the number of rows in the jacobian matrix.
 *             It is not altered by the routine.
 * @param      n      an INTEGER variable set by the user to the number of
 *                    variables
 *             @f$x_i\f$, that is the number of columns in the jacobian matrix.
 *             It is not altered by the routine.
 * @param      irn    an INTEGER array of size of the number of nonzeros in the
 *                    jacobian matrix. It must be set by the user to hold the
 *                    row indices of the nonzeros, stored by columns.
 * @param      ip     an INTEGER array of size @a n+1 that must be set by the
 *                    user so that @a ip(j) is the position in @a irn of the
 *                    start of column j in the jacobian matrix, j=1,...,n, and
 *                    @a ip(n+1) is the position of the first unused location in
 *                    @a irn. It is not altered by the subroutine. Fortran
 *                    base-1 indices are used !
 * @param      h      is a double PRECISION array of size @a n required by
 *                    td22ad to specify the step lengths to be used when
 *                    estimating the derivatives by finite differences. The step
 *                    lengths are adjusted by td22ad to give good accuracy and
 *                    are left set on the assumption that td22ad is likely to be
 *                    called later with similar values for the variables x(j),
 *                    j=1,...,n, in which case no resetting is necessary.
 * @param      x      is a double PRECISION array of size @a n that must be set
 *                    initially by the user to the point
 *             @f$x_j\f$, j=1,...,n, at which the approximate jacobian is
 *             required. It is altered by the subroutine between intermediate
 *             calls and is finally restored to its initial value.
 * @param      f      is a double PRECISION array of size @a m whose elements
 *                    must be set by the user initially and befoe every
 *                    intermediate call to
 *             @f$f_i(x)\f$, i=1,...,m. It is not altered by the subroutine,
 *             except on the final return when it is restored to its initial
 *             value.
 * @param      hmax   is a double PRECISION array of size @a n set by the user
 *                    to upper bounds for
 *             @f$h_j\f$, j=1,...,n, but if hmax(1)<0 the all bounds are taken
 *             as |hmax(1)| and the size may be one. Lower bounds for
 *             @f$h_j\f$ are taken as eps1 times the corresponding upper bound,
 *             where eps1 specifies the relative machine precision, see cntl(5).
 * @param      a      is a double PRECISION array of size the number of nonzero
 *                    derivatives into which the nonzero derivatives are placed.
 * @param      ig     is an INTEGER work array if size at least 2*n+1. ig(1)
 *                    should be set to zero before a first entry to td22ad for a
 *                    particular sparsity structure; during this entry
 *                    ig(1),...,ig(2*n+1) are found and are required for
 *                    subsequent entries to td22ad for this sparsity structure.
 * @param      w      is a double PRECISION array of size m+2*n used as
 *                    workspace.
 * @param      iflag  is an INTEGER variabke which must be set by the user to
 *                    the value 1 on the first entry to the subroutine. On
 *                    return, the value 0 indicates that the calculation is
 *                    complete and positive values indicate that the subroutine
 *                    should be recalled after calculating
 *             @f$f_i(x)\f$, i=1,...,m, and placing them in @a f without
 *             changing other arguments.
 * @param      icntl  is an INTEGER array of length 5 that contains control
 *                    parameters. Default values for the elements are set by
 *                    td22id. This argument is not altered by the subroutine.
 * @param      cntl   is a double PRECISION array of length 7 that contains
 *                    control parameters. Defaukt values for the elements are
 *                    set by td22id. This argument is not altered by the
 *                    subroutine.
 * @param      info   is an INTEGER array of length 5 that need not be set by
 *                    the user. On return, it holds information on the
 *                    calculation.
 * @param      keep   is an INTEGER array of length 20 used by td22 for private
 *                    workspace. It must be initialized by calling td22id and
 *                    must not be altered by the user.
 * @param      rkeep  is a double PRECISION array of length 10 used by td22 for
 *                    private workspace. It must be initialized by calling
 *                    td22id and must not be altered by the user.
 */
extern void td22ad_ (
    fint_t   *m,
    fint_t   *n,
    fint_t   *irn,
    fint_t   *ip,
    double *h,
    double *x,
    double *f,
    double *hmax,
    double *a,
    fint_t   *ig,
    double *w,
    fint_t   *iflag,
    fint_t   *icntl,
    double *cntl,
    fint_t   *info,
    fint_t   *keep,
    double *rkeep
);

/**
 * @brief      Construct the sparsity pattern using function values.
 *
 * @param      m      an INTEGER variable set by the user to the number of
 *                    functions
 *             @f$f_i(x)\f$, that is the number of rows in the jacobian matrix.
 *             It is not altered by the routine.
 * @param      n      an INTEGER variable set by the user to the number of
 *                    variables
 *             @f$x_i\f$, that is the number of columns in the jacobian matrix.
 *             It is not altered by the routine.
 * @param      irn    Need not be set on entry and on final return holds the
 *                    sparsity structure of the jacobian matrix.
 * @param      ip     Need not be set on entry and on final return holds the
 *                    sparsity structure of the jacobian matrix.
 * @param      h      is a double PRECISION array of size @a n required by
 *                    td22ad to specify the step lengths to be used when
 *                    altering the components of @a x to find the dependencies.
 *                    It is not altered by the subroutine.
 * @param      x      is a double PRECISION array of size @a n that must be set
 *                    initially by the user to the point
 *             @f$x_j\f$, j=1,...,n, at which the approximate jacobian is
 *             required. It is altered by the subroutine between intermediate
 *             calls and is finally restored to its initial value.
 * @param      f      is a double PRECISION array of size @a m whose elements
 *                    must be set by the user initially and befoe every
 *                    intermediate call to
 *             @f$f_i(x)\f$, i=1,...,m. It is not altered by the subroutine.
 * @param      w      is a double PRECISION array of size m used as workspace.
 * @param      ia     is an INTEGER variable specifying the size of the arrays
 *                    @a a and @a irn. It is not altered by the routine.
 * @param      iflag  is an INTEGER variabke which must be set by the user to
 *                    the value 1 on the first entry to the subroutine. On
 *                    return, the value 0 indicates that the calculation is
 *                    complete and positive values indicate that the subroutine
 *                    should be recalled after calculating
 *             @f$f_i(x)\f$, i=1,...,m, and placing them in @a f without
 *             changing other arguments. On an error, @a iflag has a negative
 *             value.
 * @param      icntl  is an INTEGER array of length 5 that contains control
 *                    parameters. Default values for the elements are set by
 *                    td22id. This argument is not altered by the subroutine.
 * @param      info   is an INTEGER array of length 5 that need not be set by
 *                    the user. On return, it holds information on the
 *                    calculation.
 * @param      keep   is an INTEGER array of length 20 used by td22 for private
 *                    workspace. It must be initialized by calling td22id and
 *                    must not be altered by the user.
 * @param      rkeep  is a double PRECISION array of length 10 used by td22 for
 *                    private workspace. It must be initialized by calling
 *                    td22id and must not be altered by the user.
 */
extern void td22bd_ (
    fint_t   *m,
    fint_t   *n,
    fint_t   *irn,
    fint_t   *ip,
    double *h,
    double *x,
    double *f,
    double *w,
    fint_t   *ia,
    fint_t   *iflag,
    fint_t   *icntl,
    fint_t   *info,
    fint_t   *keep,
    double *rkeep
);

/**
 * @brief      Construct the sparsity pattern of a band matrix.
 *
 * @param      n      is an INTEGER variable set by the user to the number of
 *                    variables
 *             @f$x_j\f$, that is the number of columns in the jacobian matrix.
 *             It is not altered by the subroutine.
 * @param      irn    an INTEGER array of size of the number of nonzeros in the
 *                    jacobian matrix. It must be set by the user to hold the
 *                    row indices of the nonzeros, stored by columns.
 * @param      ip     an INTEGER array of size @a n+1 that must be set by the
 *                    user so that @a ip(j) is the position in @a irn of the
 *                    start of column j in the jacobian matrix, j=1,...,n, and
 *                    @a ip(n+1) is the position of the first unused location in
 *                    @a irn. It is not altered by the subroutine. Fortran
 *                    base-1 indices are used !
 * @param      ia     is an INTEGER variable specifying the size of the arrays
 *                    @a a and @a irn. It is not altered by the routine.
 * @param      mbd    is an INTEGER variable which must be set by the user to
 *                    the semi-bandwidth of the jacobian, so that
 * @param      icntl  is an INTEGER array of length 5 that contains control
 *                    parameters. Default values for the elements are set by
 *                    td22id. This argument is not altered by the subroutine.
 * @param      info   is an INTEGER array of length 5 that need not be set by
 *                    the user. On return, it holds information on the
 *                    calculation.
 * @f$
 * @par tial f_i/\partial x_j\=0\f$ if
 *             @f$|i-j\geq mbd\f$. It is not altered by the subroutine.
 */
extern void td22cd_ (
    fint_t *n,
    fint_t *irn,
    fint_t *ip,
    fint_t *ia,
    fint_t *mbd,
    fint_t *icntl,
    fint_t *info
);

#ifdef __cplusplus
}
#endif

#endif // _GLUE_TD22_H_INCLUDED_
