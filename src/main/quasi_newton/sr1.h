#ifndef SLEQP_SR1_H
#define SLEQP_SR1_H

/**
 * @file sr1.h
 * @brief Defintion of SR1 method.
 **/

#include "func.h"
#include "iterate.h"
#include "options.h"
#include "params.h"
#include "timer.h"

#include "quasi_newton_types.h"

#ifdef __cplusplus
extern "C" {
#endif

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_sr1_create(SleqpQuasiNewton** star,
                                 SleqpFunc* func,
                                 SleqpParams* params,
                                 SleqpOptions* options);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_SR1_H */
