#ifndef SLEQP_ITERATE_H
#define SLEQP_ITERATE_H

#include "sleqp_types.h"

#include "sleqp.h"
#include "sparse/sleqp_sparse.h"

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct SleqpIterate
  {
    SleqpSparseVec* x;

    double func_val;
    SleqpSparseVec* func_grad;

    SleqpSparseVec* cons_val;
    SleqpSparseMatrix* cons_jac;

    SleqpSparseVec* cons_dual;
    SleqpSparseVec* vars_dual;

  } SleqpIterate;

  SLEQP_RETCODE sleqp_iterate_create(SleqpIterate** star,
                                     SleqpProblem* problem,
                                     SleqpSparseVec* x);

  SLEQP_RETCODE sleqp_iterate_free(SleqpIterate** star);


#ifdef __cplusplus
}
#endif

#endif /* SLEQP_ITERATE_H */
