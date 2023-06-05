#ifndef SLEQP_STANDARD_AUG_JAC_H
#define SLEQP_STANDARD_AUG_JAC_H

#include "aug_jac.h"
#include "fact/fact.h"
#include "iterate.h"
#include "problem.h"

/**
 * Create a standard augmented Jacobian operating
 * on the system
 *
 * \f[
 * \pmatrix{
 * I & A_W^T \\
 * A_W & 0
 * }
 * \f]
 *
 **/
SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_standard_aug_jac_create(SleqpAugJac** star,
                              SleqpProblem* problem,
                              SleqpSettings* settings,
                              SleqpFact* factorization);

#endif /* SLEQP_STANDARD_AUG_JAC_H */
