#ifndef SLEQP_STANDARD_AUG_JAC_H
#define SLEQP_STANDARD_AUG_JAC_H

#include "aug_jac.h"
#include "problem.h"
#include "iterate.h"
#include "params.h"
#include "sparse/sparse_factorization.h"

#ifdef __cplusplus
extern "C" {
#endif

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_standard_aug_jac_create(SleqpAugJac** star,
                                              SleqpProblem* problem,
                                              SleqpParams* params,
                                              SleqpSparseFactorization* sparse_factorization);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_STANDARD_AUG_JAC_H */
