#ifndef SLEQP_BOX_CONSTRAINED_AUG_JAC_H
#define SLEQP_BOX_CONSTRAINED_AUG_JAC_H

#include "aug_jac.h"

SLEQP_NODISCARD
SLEQP_RETCODE sleqp_box_constrained_aug_jac_create(SleqpAugJac** star,
                                                   SleqpProblem* problem);

#endif /* SLEQP_BOX_CONSTRAINED_AUG_JAC_H */
