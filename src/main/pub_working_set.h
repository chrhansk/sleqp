#ifndef SLEQP_PUB_WORKING_SET_H
#define SLEQP_PUB_WORKING_SET_H

/**
 * @file pub_working_set.h
 * @brief Definition of working sets.
 **/

#include "sleqp/export.h"

#include "sleqp/pub_problem.h"

typedef struct SleqpWorkingSet SleqpWorkingSet;

/**
 * Creates a new working set
 *
 * @param[out] star            A pointer to the working set to be created
 * @param[int] problem         The underlying problem
 **/
SLEQP_EXPORT SLEQP_NODISCARD SLEQP_RETCODE
sleqp_working_set_create(SleqpWorkingSet** star, SleqpProblem* problem);

/**
 * Returns the state of the given variable with respect to the given working set
 *
 * @param[in]  working_set           The working set
 * @param[in]  index                 The variable index
 **/
SLEQP_EXPORT SLEQP_ACTIVE_STATE
sleqp_working_set_var_state(const SleqpWorkingSet* working_set, int index);

/**
 * Returns the state of the given constraint with respect to the given working
 * set
 *
 * @param[in]  working_set           The working set
 * @param[in]  index                 The constraint index
 **/
SLEQP_EXPORT SLEQP_ACTIVE_STATE
sleqp_working_set_cons_state(const SleqpWorkingSet* working_set, int index);

/**
 * Returns an array containing the states of all variables with respect to the
 * given working set
 **/
SLEQP_EXPORT const SLEQP_ACTIVE_STATE*
sleqp_working_set_var_states(const SleqpWorkingSet* working_set);

/**
 * Returns an array containing the states of all constraints with respect to the
 * given working set
 **/
SLEQP_EXPORT const SLEQP_ACTIVE_STATE*
sleqp_working_set_cons_states(const SleqpWorkingSet* working_set);

SLEQP_EXPORT SLEQP_ACTIVE_STATE
sleqp_working_set_state(const SleqpWorkingSet* working_set,
                        bool constraint,
                        int index);

/**
 * Returns the problem underling the given working set
 *
 * @param[in]  working_set           The working set
 **/
SLEQP_EXPORT SleqpProblem*
sleqp_working_set_problem(const SleqpWorkingSet* working_set);

/**
 * Returns the number of variables contained in the given working set
 *
 * @param[in]  working_set           The working set
 **/
SLEQP_EXPORT int
sleqp_working_set_num_active_vars(const SleqpWorkingSet* working_set);

/**
 * Returns the number of constraints contained in the given working set
 *
 * @param[in]  working_set           The working set
 **/
SLEQP_EXPORT int
sleqp_working_set_num_active_cons(const SleqpWorkingSet* working_set);

/**
 * Returns the size of the given set, i.e. the number of contained variables
 *plus the number of contained constraints
 *
 * @param[in]  working_set           The working set
 **/
SLEQP_EXPORT int
sleqp_working_set_size(const SleqpWorkingSet* working_set);

SLEQP_EXPORT SLEQP_NODISCARD SLEQP_RETCODE
sleqp_working_set_capture(SleqpWorkingSet* working_set);

SLEQP_EXPORT SLEQP_NODISCARD SLEQP_RETCODE
sleqp_working_set_release(SleqpWorkingSet** star);

#endif /* SLEQP_PUB_WORKING_SET_H */
