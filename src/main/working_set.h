#ifndef SLEQP_WORKING_SET_H
#define SLEQP_WORKING_SET_H

#include "pub_working_set.h"

/**
 * Returns whether the given working set is *valid*, i.e., whether
 * - all interal indices are consistent
 * - the internally stored sizes are consistent
 * - the intenal indices are consistent with variables states
 *
 * @param[in]  working_set           The working set
 **/
bool
sleqp_working_set_valid(const SleqpWorkingSet* working_set);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_working_set_supports_cons_dual(const SleqpWorkingSet* working_set,
                                     SleqpVec* cons_dual,
                                     bool* supports);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_working_set_supports_vars_dual(const SleqpWorkingSet* working_set,
                                     SleqpVec* vars_dual,
                                     bool* supports);

/**
 * Prints the given working set to the given file
 *
 * @param[in]  vec     A pointer to the vector
 * @param[in]  output  A pointer to an output `FILE*`
 **/
SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_working_set_fprintf(const SleqpWorkingSet* working_set, FILE* output);

/**
 * Copies one working set to another
 *
 * @param[in]  source  A pointer to the copy source
 * @param[in]  target  A pointer to the copy target
 **/
SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_working_set_copy(const SleqpWorkingSet* source, SleqpWorkingSet* target);

/**
 * Resets this working set by removing all variables and constraints from it
 *
 * @param[in]  working_set           The working set
 **/
SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_working_set_reset(SleqpWorkingSet* working_set);

SLEQP_NODISCARD
bool
sleqp_working_set_eq(SleqpWorkingSet* first, SleqpWorkingSet* second);

/**
 * Adds a variable to the given working set
 *
 * @param[in]  working_set           The working set
 * @param[in]  index                 The variable index
 * @param[in]  state                 The variable state
 *
 * @note The `state` must not be @ref SLEQP_INACTIVE
 * @note Variables must be added before constraints
 **/
SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_working_set_add_var(SleqpWorkingSet* working_set,
                          int index,
                          SLEQP_ACTIVE_STATE state);

/**
 * Adds a constraint to the given working set
 *
 * @param[in]  working_set           The working set
 * @param[in]  index                 The constraint index
 * @param[in]  state                 The constraint state
 *
 * @note The `state` must not be @ref SLEQP_INACTIVE
 * @note Variables must be added before constraints
 **/
SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_working_set_add_cons(SleqpWorkingSet* working_set,
                           int index,
                           SLEQP_ACTIVE_STATE state);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_working_set_add(SleqpWorkingSet* working_set,
                      int index,
                      bool constraint,
                      SLEQP_ACTIVE_STATE state);

/**
 * Returns the index of the given constraint with respect to the given
 * working set, or @ref SLEQP_NONE if the constraint is not contained
 * in the working set
 *
 * @param[in]  working_set           The working set
 * @param[in]  index                 The constraint index
 **/
int
sleqp_working_set_cons_index(const SleqpWorkingSet* working_set, int index);

/**
 * Returns the index of the given variable with respect to the given
 * working set, or @ref SLEQP_NONE if the variable is not contained
 * in the working set
 *
 * @param[in]  working_set           The working set
 * @param[in]  index                 The variable index
 **/
int
sleqp_working_set_var_index(const SleqpWorkingSet* working_set, int index);

/**
 * Returns the content of the working set at the given working set index,
 * which must be at least zero and less than @ref sleqp_working_set_size.
 * If the returned index is less than the number of variables, it corresponds
 * to a variable. Otherwise, index - num_variables corresponds to a constraint.
 *
 * @param[in]  working_set           The working set
 * @param[in]  index                 The working set index
 **/
int
sleqp_working_set_content(const SleqpWorkingSet* working_set, int index);

#endif /* SLEQP_WORKING_SET_H */
