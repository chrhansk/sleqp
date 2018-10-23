#ifndef SLEQP_ACTIVE_SET_H
#define SLEQP_ACTIVE_SET_H

/**
 * @file sleqp_active_set.h
 * @brief Definition of active sets.
 **/

#include "sleqp_func.h"
#include "sleqp_problem.h"
#include "sleqp_types.h"

#ifdef __cplusplus
extern "C" {
#endif

  enum SLEQP_Active_State {SLEQP_INACTIVE = 0,
                           SLEQP_ACTIVE_LOWER = 1,
                           SLEQP_ACTIVE_UPPER = 2};

  typedef enum SLEQP_Active_State SLEQP_ACTIVE_STATE;

  typedef struct SleqpActiveSet SleqpActiveSet;

  SLEQP_RETCODE sleqp_active_set_create(SleqpActiveSet** star,
                                        SleqpProblem* problem);

  SLEQP_RETCODE sleqp_active_set_reset(SleqpActiveSet* active_set);

  SLEQP_ACTIVE_STATE* sleqp_active_set_var_states(SleqpActiveSet* active_set);
  SLEQP_ACTIVE_STATE* sleqp_active_set_cons_states(SleqpActiveSet* active_set);

  SLEQP_RETCODE sleqp_active_set_free(SleqpActiveSet** star);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_ACTIVE_SET_H */
