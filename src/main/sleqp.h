#ifndef SLEQP_H
#define SLEQP_H

#include "sleqp_func.h"
#include "sleqp_types.h"

#include "sparse/sleqp_sparse.h"

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct SleqpProblem
  {
    SleqpFunc* func;

    SleqpSparseVec* var_lb;
    SleqpSparseVec* var_ub;
    SleqpSparseVec* cons_lb;
    SleqpSparseVec* cons_ub;

    size_t num_variables, num_constraints;

  } SleqpProblem;

  typedef struct SleqpSolver SleqpSolver;


  SLEQP_RETCODE sleqp_problem_create(SleqpProblem** star,
                                     SleqpFunc* func,
                                     SleqpSparseVec* var_lb,
                                     SleqpSparseVec* var_ub,
                                     SleqpSparseVec* cons_lb,
                                     SleqpSparseVec* cons_ub);

  /*
  SLEQP_RETCODE sleqp_problem_solve(Sleqp* sleqp,
                                    SleqpSparseVec* x);
  */

  SLEQP_RETCODE sleqp_problem_free(SleqpProblem** star);


  SLEQP_RETCODE sleqp_solver_create(SleqpSolver** star,
                                    SleqpProblem* problem);

  SLEQP_RETCODE sleqp_solve(SleqpSolver* solver,
                            SleqpSparseVec* x);

  SLEQP_RETCODE sleqp_solver_free(SleqpSolver** star);


#ifdef __cplusplus
}
#endif

#endif /* SLEQP_H */
