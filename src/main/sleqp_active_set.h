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

  typedef struct SleqpActiveSet SleqpActiveSet;

  SLEQP_RETCODE sleqp_active_set_create(SleqpActiveSet** star,
                                        SleqpProblem* problem);

  SLEQP_RETCODE sleqp_active_set_reset(SleqpActiveSet* active_set);

  SLEQP_RETCODE sleqp_active_set_add_variable(SleqpActiveSet* active_set,
                                              int index,
                                              SLEQP_ACTIVE_STATE state);

  SLEQP_RETCODE sleqp_active_set_add_constraint(SleqpActiveSet* active_set,
                                                int index,
                                                SLEQP_ACTIVE_STATE state);

  int sleqp_active_set_get_constraint_index(SleqpActiveSet* active_set,
                                            int index);

  int sleqp_active_set_get_variable_index(SleqpActiveSet* active_set,
                                          int index);

  int sleqp_active_set_get_content(SleqpActiveSet* active_set,
                                   int index);

  SLEQP_ACTIVE_STATE sleqp_active_set_get_variable_state(SleqpActiveSet* active_set,
                                                         int index);

  SLEQP_ACTIVE_STATE sleqp_active_set_get_constraint_state(SleqpActiveSet* active_set,
                                                           int index);

  int sleqp_active_set_num_active_vars(SleqpActiveSet* active_set);

  int sleqp_active_set_num_active_conss(SleqpActiveSet* active_set);

  int sleqp_active_set_size(SleqpActiveSet* active_set);

  SLEQP_RETCODE sleqp_active_set_free(SleqpActiveSet** star);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_ACTIVE_SET_H */
