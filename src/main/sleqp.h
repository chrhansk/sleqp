#ifndef SLEQP_H
#define SLEQP_H

#include "sleqp_func.h"
#include "sleqp_types.h"

#include "sparse/sleqp_sparse.h"

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
    size_t num_variables;

    /**
     * number of constraints.
     **/
    size_t num_constraints;

  } SleqpProblem;

  typedef struct SleqpSolver SleqpSolver;


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
                                     SleqpSparseVec* var_lb,
                                     SleqpSparseVec* var_ub,
                                     SleqpSparseVec* cons_lb,
                                     SleqpSparseVec* cons_ub);


  typedef struct SleqpIterate SleqpIterate;

  SLEQP_RETCODE sleqp_set_and_evaluate(SleqpProblem* problem,
                                       SleqpIterate* iterate);

  /**
   * Frees a problem created previously.
   **/
  SLEQP_RETCODE sleqp_problem_free(SleqpProblem** star);


  SLEQP_RETCODE sleqp_solver_create(SleqpSolver** star,
                                    SleqpProblem* problem,
                                    SleqpSparseVec* x);

  SLEQP_RETCODE sleqp_solve(SleqpSolver* solver);

  SLEQP_RETCODE sleqp_solver_free(SleqpSolver** star);


#ifdef __cplusplus
}
#endif

#endif /* SLEQP_H */
