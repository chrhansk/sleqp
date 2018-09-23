#ifndef SLEQP_H
#define SLEQP_H

#include "sleqp_func.h"
#include "sleqp_types.h"
#include "sleqp_problem.h"

#include "sparse/sleqp_sparse.h"

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct SleqpSolver SleqpSolver;


  typedef struct SleqpIterate SleqpIterate;

  SLEQP_RETCODE sleqp_set_and_evaluate(SleqpProblem* problem,
                                       SleqpIterate* iterate);


  SLEQP_RETCODE sleqp_solver_create(SleqpSolver** star,
                                    SleqpProblem* problem,
                                    SleqpSparseVec* x);

  SLEQP_RETCODE sleqp_solve(SleqpSolver* solver);

  SLEQP_RETCODE sleqp_solver_free(SleqpSolver** star);


#ifdef __cplusplus
}
#endif

#endif /* SLEQP_H */
