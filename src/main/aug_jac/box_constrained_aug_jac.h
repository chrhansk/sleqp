#ifndef SLEQP_BOX_CONSTRAINED_AUG_JAC_H
#define SLEQP_BOX_CONSTRAINED_AUG_JAC_H

#include "aug_jac.h"

#ifdef __cplusplus
extern "C" {
#endif

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_box_constrained_aug_jac_create(SleqpAugJac** star,
                                                     SleqpProblem* problem);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_BOX_CONSTRAINED_AUG_JAC_H */
