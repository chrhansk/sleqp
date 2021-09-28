#ifndef SLEQP_UNCONSTRAINED_AUG_JAC_H
#define SLEQP_UNCONSTRAINED_AUG_JAC_H

#include "aug_jac.h"

#ifdef __cplusplus
extern "C" {
#endif

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_unconstrained_aug_jac_create(SleqpAugJac** star,
                                                   SleqpProblem* problem);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_UNCONSTRAINED_AUG_JAC_H */
