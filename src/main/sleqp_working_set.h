#ifndef SLEQP_WORKING_SET_H
#define SLEQP_WORKING_SET_H

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

  typedef struct SleqpWorkingSet SleqpWorkingSet;

  SLEQP_RETCODE sleqp_working_set_create(SleqpWorkingSet** star,
                                         SleqpProblem* problem);

  SLEQP_RETCODE sleqp_working_set_reset(SleqpWorkingSet* working_set);

  SLEQP_RETCODE sleqp_working_set_add_variable(SleqpWorkingSet* working_set,
                                               int index,
                                               SLEQP_ACTIVE_STATE state);

  SLEQP_RETCODE sleqp_working_set_add_constraint(SleqpWorkingSet* working_set,
                                                 int index,
                                                 SLEQP_ACTIVE_STATE state);

  int sleqp_working_set_get_constraint_index(SleqpWorkingSet* working_set,
                                             int index);

  int sleqp_working_set_get_variable_index(SleqpWorkingSet* working_set,
                                           int index);

  int sleqp_working_set_get_content(SleqpWorkingSet* working_set,
                                    int index);

  SLEQP_ACTIVE_STATE sleqp_working_set_get_variable_state(SleqpWorkingSet* working_set,
                                                          int index);

  SLEQP_ACTIVE_STATE sleqp_working_set_get_constraint_state(SleqpWorkingSet* working_set,
                                                            int index);

  int sleqp_working_set_num_active_vars(SleqpWorkingSet* working_set);

  int sleqp_working_set_num_active_cons(SleqpWorkingSet* working_set);

  int sleqp_working_set_size(SleqpWorkingSet* working_set);

  SLEQP_RETCODE sleqp_working_set_fprintf(SleqpWorkingSet* working_set,
                                          FILE* output);

  SLEQP_RETCODE sleqp_working_set_copy(SleqpWorkingSet* source,
                                       SleqpWorkingSet* target);

  SLEQP_RETCODE sleqp_working_set_free(SleqpWorkingSet** star);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_WORKING_SET_H */
