/** @file ma57.h
 *  @brief Header file for MA57
 *  @author Christian Kirches
 */

#pragma once
#ifndef _GLUE_MA57_H_INCLUDED_
#define _GLUE_MA57_H_INCLUDED_

#ifdef __cplusplus
extern "C" {
#endif

#include "util_types.h"

/**
 * @brief      Sets default values for the components of the arrays that hold
 *             control parameters.
 *
 *             Normally the user will call MA57I/ID prior to any call to
 *             MA57A/AD. If non- default values for any of the control
 *             parameters are required, they should be set immediately after the
 *             call to MA57I/ID.
 *
 * @param[in]  <unnamed>  { parameter_description }
 * @param      cntl   a double PRECISION array of length 5 that need not be set by
 *                    the user. On return it contains default values.
 * @param      icntl  an INTEGER array of length 20 that need not be set by the
 *                    user. On return it contains default values.
 */
extern void ma57id_ (
    double *cntl,
    fint_t  *icntl
);

/**
 * @brief      Accepts the pattern of A and chooses pivots for Gaussian
 *             elimination.
 *
 *             Uses a selection criterion to preserve sparsity. It subsequently
 *             constructs subsidiary information for the actual factorization by
 *             MA57B/BD. The user may provide the pivot sequence; in which case
 *             only the necessary information for MA57B/BD will be generated.
 *
 * @param[in]  <unnamed>  { parameter_description }
 * @param      n      an INTEGER variable that must be set by the user to the order
 *                    n of the matrix A. It is not altered by the subroutine.
 *                    Restriction: @a n >= 1.
 * @param      ne     an INTEGER variable that must be set by the user to the
 *                    number of entries input in IRN and JCN. It is not altered by
 *                    the subroutine. Restriction: @a ne >= 0.
 * @param      irn    INTEGER array of length @a ne.
 * @param      jcn    INTEGER array of length @a ne. The user must set them so that
 *                    each diagonal entry a(ii) is represented by IRN(k)=i and
 *                    JCN(k)=i and each pair of off-diagonal entries a(ij) and
 *                    a(ji) is represented by IRN(k)=i and JCN(k)=j or by IRN(k)=j
 *                    and JCN(k)=i. Entries (on or off the diagonal) that are known
 *                    to be zero can be excluded. Multiple entries are permitted.
 *                    If IRN(k) or JCN(k) are less than 1 or greater than N the
 *                    entry is ignored. These arrays are not altered by any of the
 *                    calls to the MA57 package. They must be preserved by the user
 *                    between this call and a call to MA57D/DD for the same matrix.
 * @param      lkeep  an INTEGER variable that must be set by the user to the
 *                    length of array @a keep. It might be more efficient to
 *                    allocate more than the minimum required, say about N to 2*N
 *                    more space. Restriction: @a lkeep >= 5*n+ne+max(n,ne)+42.
 * @param      keep   an INTEGER array of length LKEEP. It need not be set by the
 *                    user and must be preserved between a call to MA57A/AD and
 *                    subsequent calls to MA57B/BD. If the user wishes to input the
 *                    pivot sequence, the position of variable i in the pivot order
 *                    should be placed in KEEP(I), I = 1, 2,..., N and ICNTL(6)
 *                    should be set to 1. The subroutine may replace the given
 *                    order by another that gives the same fill-in pattern and
 *                    virtually identical numerical results.
 * @param      iwork  an INTEGER array of length 5*N that is used as workspace.
 * @param      icntl  an INTEGER array of length 20 that contains control
 *                    parameters and must be set by the user. Default values for
 *                    the components may be set by a call to MA57I/ID.
 * @param      info   an INTEGER array of length 40 that need not be set by the
 *                    user. On return from MA57A/AD, a non-negative value for
 *                    INFO(1) indicates that the subroutine has performed
 *                    successfully.
 * @param      rinfo  a double PRECISION array of length 20 that need not be set by
 *                    the user. This array supplies information on the execution of
 *                    the subroutine.
 * @see        @a jcn for further details.
 */
extern void ma57ad_ (
    const fint_t *n,
    const fint_t *ne,
    const fint_t *irn,
    const fint_t *jcn,
    fint_t       *lkeep,
    fint_t       *keep,
    fint_t       *iwork,
    fint_t       *icntl,
    fint_t       *info,
    double      *rinfo
);

/** \brief Factorizes a matrix A
 *
 * Uses the information from a previous call to MA57A/AD. The actual pivot
 * sequence used may differ from that of MA57A/AD if A is not definite.
 *
 * \param n     is an INTEGER variable that must be set by the user to the order
 *              n of the matrix A. It must be unchanged since the call to MA57A/
 *              AD and is not altered by the subroutine. Restriction: \a n >= 1.
 * \param ne    is an INTEGER variable that must be set by the user to the
 *              number of entries in the matrix A. It must be unchanged since
 *              the call to MA57A/AD and is not altered by the subroutine.
 *              Restriction: \a ne >= 0.
 * \param a     is a double PRECISION array of length NE that must be set by the
 *              user so that A(k) holds the value of the diagonal entry or pair
 *              of off-diagonal entries whose indices were held in IRN(k) and
 *              JCN(k) on entry to MA57A/AD, for k = 1, 2,..., NE. Multiple
 *              entries are summed and any that correspond to an IRN(k) or
 *              JCN(k) value that was out of range are ignored. The array A is
 *              not altered by the subroutine and must be preserved by the user
 *              between this call and a call to MA57D/DD for the same matrix.
 * \param fact  is a double PRECISION array of length LFACT. It need not be set
 *              by the user and, on exit, will hold the entries of the factors
 *              of the matrix A. These entries in FACT must be preserved by the
 *              user between calls to this subroutine and subsequent calls to
 *              MA57C/CD or MA57D/DD. If MA57B/BD is called with ICNTL(7) set to
 *              4, the entries LFACT-N+1 to LFACT of FACT will hold the changes
 *              made to the pivot entries by the matrix modification. If MA57B/
 *              BD is being called after a call to MA57E/ED, then FACT must be
 *              the same array as the array NEWFAC as output from MA57E/ED.
 * \param lfact is an INTEGER variable that must be set by the user to the
 *              length of array FACT. It must be at least as great as INFO(9) as
 *              output from MA57A/AD (see Section 2.2). It is advisable to allow
 *              a slightly greater value because numerical pivoting might
 *              increase storage requirements. In extreme cases, for example
 *              when there are many zeros on the diagonal, a value of more than
 *              twice INFO(9) might be needed. It is not altered by the
 *              subroutine.
 * \param ifact is an INTEGER array of length LIFACT. It need not be set by the
 *              user and, on exit, holds integer indexing information on the
 *              matrix factors. It must be preserved by the user between calls
 *              to this subroutine and subsequent calls to MA57C/CD or MA57D/DD.
 *              If MA57B/BD is being called after a call to MA57E/ED, then IFACT
 *              must be the same array as the array NEWIFC as output from MA57E/
 *              ED.
 * \param lifact is an INTEGER variable that must be set by the user to the
 *               length of array IFACT. It must be at least as great as INFO(10)
 *               as output from MA57A/AD (see Section 2.2). A slightly greater
 *               value is recommended because numerical pivoting may increase
 *               storage requirements. It is not altered by the subroutine.
 * \param lkeep is an INTEGER variable that must be set by the user to the
 *              length of array KEEP. Restriction: LKEEP >= 5*n+ne+max(n,ne)+42.
 * \param keep  is an INTEGER array of length LKEEP that must be unchanged since
 *              the call to MA57A/AD. It is not altered by MA57B/BD.
 * \param iwork is an INTEGER array of length N. It is used as workspace by the
 *              subroutine.
 * \param icntl is an INTEGER array of length 20 that contains control
 *              parameters and must be set by the user. Default values for the
 *              components may be set by a call to MA57I/ID.
 * \param cntl  is a double PRECISION array of length 5 that contains control
 *              parameters and must be set by the user. Default values for the
 *              components may be set by a call to MA57I/ID.
 * \param info  is an INTEGER array of length 40 that need not be set by the
 *              user. On return from MA57B/BD, a non-negative value for INFO(1)
 *              indicates that the subroutine has performed successfully.
 * \param rinfo  a double PRECISION array of length 20 that need not be set by
 *               the user. This array supplies information on the execution of
 *               MA57B/BD.
 */
extern void ma57bd_ (
    const fint_t  *n,
    const fint_t  *ne,
    const double *a,
    double       *fact,
    const fint_t  *lfact,
    fint_t        *ifact,
    const fint_t  *lifact,
    fint_t        *lkeep,
    const fint_t  *keep,
    fint_t        *iwork,
    fint_t        *icntl,
    double       *cntl,
    fint_t        *info,
    double       *rinfo
);

/** \brief Solve a system of equations Ax=b (AX=B).
 *
 * Uses the factors generated by MA57B/BD to solve a system of equations Ax=b
 * (AX=B) using iterative refinement and (optionally) returning estimates of the
 * error.
 *
 * \param job    must be set by the user to determine what action is desired by
 *               the user. Values of JOB and their effect are: If ICNTL(9)>1
 *               then JOB=0 if no estimate of solution in X and JOB=2 if
 *               estimate of solution in X. If ICNTL(9)=1, then:
 *                 0: Solve Ax=b, calculate residual r=b-Ax and exit. (Note that
 *                    MA57C/CD should be used if solution without residual is
 *                    required)
 *                 1: Solve Ax=b, calculate residual r=b-Ax, solve A(dx)=r,
 *                    update solution and exit.
 *               If JOB > 1, an estimate of the solution must be input in X.
 *                 2: Calculate residual r=b-Ax, solve A(dx)=r, update solution
 *                    and exit.
 *               If JOB > 2, the residual for this estimate must also be input.
 *                 3: Solve A(dx)=r, update solution and exit.
 *                 4: Solve A(dx)=r, update solution, calculate residual for new
 *                    solution and exit.
 *               Restriction: 0 <= JOB <= 4 and JOB = 0 or 2 if ICNTL(9) > 1.
 * \param n      is an INTEGER variable that must be set by the user to the
 *               order n of the matrix A. It must be unchanged since the call to
 *               MA57B/BD and is not altered by the subroutine.
 *               Restriction: N >= 1.
 * \param ne     is an INTEGER variable that must be set by the user to the
 *               number of entries in the matrix A. It must be unchanged since
 *               the call to MA57A/AD and is not altered by the subroutine.
 *               Restriction: NE >= 0.
 * \param a      is a double PRECISION array of length NE that must be unchanged
 *               since the call to MA57B/BD. It is not altered by the
 *               subroutine.
 * \param irn    INTEGER array of length NE. See \a jcn for further details.
 * \param jcn    INTEGER array of length NE that must be unchanged since the
 *               call to MA57A/AD. They are not altered by the subroutine.
 * \param fact   is a double PRECISION array of length LFACT that must be
 *               unchanged since the call to MA57B/BD. It is not altered by the
 *               subroutine.
 * \param lfact  an INTEGER variable that must be set by the user to the length
 *               of array FACT. It is not altered by the subroutine.
 * \param ifact  an INTEGER array of length LIFACT that must be unchanged since
 *               the call to MA57B/BD. It is not altered by the subroutine.
 * \param lifact an INTEGER variable that must be set by the user to the length
 *               of array IFACT. It is not altered by the subroutine.
 * \param rhs    a double PRECISION array of of length N that must be set by the
 *               user to the right-hand side for the equation being solved. It
 *               is not altered by the subroutine.
 * \param x      a double PRECISION array of of length N. If JOB>=2, it must be
 *               set by the user to an estimated solution. Otherwise, it need
 *               not be set by the user. On exit, it holds the improved solution
 *               vector.
 * \param resid  a double PRECISION array of length N. If JOB>2, it must be set
 *               on entry to the value of the residual for the current solution
 *               estimate held in X. Otherwise, it need not be set by the user
 *               on entry. If JOB=0 or 4 or if ICNTL(9)>1, on exit it will hold
 *               the residual vector for the equations being solved. If 1<=JOB
 *               <=3, then RESID will hold on exit the last correction vector
 *               added to the solution X.
 * \param work   a double PRECISION work array. If ICNTL(9)=1, it must be of
 *               length at least N. If ICNTL(9)>1, it must be of length at least
 *               3*N. If ICNTL(9)>1 and ICNTL(10)>0, then WORK must be of length
 *               at least 4*N.
 * \param iwork  an INTEGER of length N that is used as a work array if ICNTL(1)
 *               >9. It is not accessed if ICNTL(9)=1.
 * \param icntl  an INTEGER array of length 20 that contains control parameters
 *               and must be set by the user. Default values for the components
 *               may be set by a call to MA57I/ID.
 * \param cntl   a double PRECISION array of length 5 that contains control
 *               parameters and must be set by the user. Default values for the
 *               components may be set by a call to MA57I/ID.
 * \param info   an INTEGER array of length 40 that need not be set by the user.
 *               On return from MA57D/DD, a non-negative value for INFO(1)
 *               indicates that the subroutine has performed successfully.
 * \param rinfo  a double PRECISION array of length 20 that need not be set by
 *               the user. This array supplies information on the execution of
 *               MA57D/DD.
 */
extern void ma57dd_ (
    const fint_t  *job,
    const fint_t  *n,
    const fint_   *ne,
    const double *a,
    const fint_t  *irn,
    const fint_t  *jcn,
    const double *fact,
    const fint_t  *lfact,
    const fint_t  *ifact,
    const fint_t  *lifact,
    const double *rhs,
    double       *x,
    double       *resid,
    double       *work,
    const fint_t  *iwork,
    fint_t        *icntl,
    double       *cntl,
    fint_t        *info,
    double       *rinfo
);

#ifdef __cplusplus
}
#endif

#endif // _GLUE_MA57_H_INCLUDED_
