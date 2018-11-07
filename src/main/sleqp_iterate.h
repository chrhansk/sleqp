#ifndef SLEQP_ITERATE_H
#define SLEQP_ITERATE_H

/**
 * @file sleqp_iterate.h
 * @brief Definition of iterate.
 **/


#include "sleqp_active_set.h"
#include "sleqp_problem.h"

#include "sparse/sleqp_sparse_matrix.h"
#include "sparse/sleqp_sparse_vec.h"

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct SleqpIterate
  {
    /**
     * The current point. Has dimension = num_variables.
     **/
    SleqpSparseVec* x;

    /**
     * The current function value
     **/
    double func_val;

    /**
     * The current function gradient. Has dimension = num_variables.
     **/
    SleqpSparseVec* func_grad;

    /**
     * The current constraint values. Has dimension = num_constraints.
     **/
    SleqpSparseVec* cons_val;

    /**
     * The Jacobian of the constraitns at the current iterate.
     * Has num_constraints many rows, num_variables many columns.
     */
    SleqpSparseMatrix* cons_jac;

    /**
     * The current active set.
     **/
    SleqpActiveSet* active_set;

    /**
     * The dual values of the constraints. Has dimension = num_constraints.
     */
    SleqpSparseVec* cons_dual;

    /**
     * The dual values of the variable bounds. Has dimension = num_variables.
     */
    SleqpSparseVec* vars_dual;

  } SleqpIterate;

  SLEQP_RETCODE sleqp_iterate_create(SleqpIterate** star,
                                     SleqpProblem* problem,
                                     SleqpSparseVec* x);

  SLEQP_RETCODE sleqp_iterate_active_set_contains(SleqpIterate* iterate,
                                                  SleqpProblem* problem,
                                                  SleqpSparseVec* direction,
                                                  double eps,
                                                  double* cache,
                                                  bool* contained);

  double sleqp_iterate_slackness_residuum(SleqpIterate* iterate,
                                          SleqpProblem* problem);

  SLEQP_RETCODE sleqp_iterate_free(SleqpIterate** star);


#ifdef __cplusplus
}
#endif

#endif /* SLEQP_ITERATE_H */
