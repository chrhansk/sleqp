#ifndef SLEQP_H
#define SLEQP_H

#include "sleqp_func.h"
#include "sleqp_iterate.h"
#include "sleqp_problem.h"
#include "sleqp_types.h"

#include "sparse/sleqp_sparse_matrix.h"
#include "sparse/sleqp_sparse_vec.h"

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct SleqpSolver SleqpSolver;

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
