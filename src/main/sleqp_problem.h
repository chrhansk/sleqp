#ifndef SLEQP_PROBLEM_H
#define SLEQP_PROBLEM_H

#include "sleqp_func.h"
#include "sleqp_types.h"
#include "sleqp_params.h"

#include "sparse/sleqp_sparse_matrix.h"
#include "sparse/sleqp_sparse_vec.h"

#ifdef __cplusplus
extern "C" {
#endif


  /**
   * Main problem definion.
   **/
  typedef struct SleqpProblem
  {
    /**
     * The function is vector-valued
     * containing both the objective
     * and the constraints.
     **/
    SleqpFunc* func;

    /**
     * lower bounds on variables.
     **/
    SleqpSparseVec* var_lb;

    /**
     * upper bounds on variables.
     **/
    SleqpSparseVec* var_ub;

    /**
     * lower bounds on constraints.
     **/
    SleqpSparseVec* cons_lb;

    /**
     * upper bounds on constraints.
     **/
    SleqpSparseVec* cons_ub;

    /**
     * number of variables.
     **/
    int num_variables;

    /**
     * number of constraints.
     **/
    int num_constraints;

  } SleqpProblem;


  /**
   * Creates a new problem, which will be saved
   * in the newly allocated pointer.
   *
   * The function pointer is kept as a reference,
   * The upper / lower bound vectors are copied.
   *
   **/
  SLEQP_RETCODE sleqp_problem_create(SleqpProblem** star,
                                     SleqpFunc* func,
                                     SleqpParams* params,
                                     SleqpSparseVec* var_lb,
                                     SleqpSparseVec* var_ub,
                                     SleqpSparseVec* cons_lb,
                                     SleqpSparseVec* cons_ub);

  /**
   * Frees a previously created problem.
   **/
  SLEQP_RETCODE sleqp_problem_free(SleqpProblem** star);


#ifdef __cplusplus
}
#endif

#endif /* SLEQP_PROBLEM_H */
