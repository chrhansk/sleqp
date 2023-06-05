#ifndef SLEQP_REDUCED_AUG_JAC_H
#define SLEQP_REDUCED_AUG_JAC_H

#include "aug_jac.h"
#include "fact/fact.h"
#include "iterate.h"
#include "problem.h"

/**
 * Create a reduced augmented Jacobian operating
 * on the system \f$ A_W A_W^{T} \f$
 *
 **/
SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_reduced_aug_jac_create(SleqpAugJac** star,
                             SleqpProblem* problem,
                             SleqpSettings* settings,
                             SleqpFact* factorization);

#endif /* SLEQP_REDUCED_AUG_JAC_H */
