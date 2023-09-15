#ifndef SLEQP_DIRECT_AUG_JAC_H
#define SLEQP_DIRECT_AUG_JAC_H

#include "aug_jac.h"
#include "fact/fact_qr.h"

/**
 * Works on the augmented Jacobian
 * by QR factorizing \f$ A_W \f$
 *
 **/
SLEQP_WARNUNUSED
SLEQP_RETCODE
sleqp_direct_aug_jac_create(SleqpAugJac** star,
                            SleqpProblem* problem,
                            SleqpSettings* settings,
                            SleqpFactQR* fact);

#endif /* SLEQP_DIRECT_AUG_JAC_H */
