#ifndef SLEQP_WORKING_SET_H
#define SLEQP_WORKING_SET_H

/**
 * @file sleqp_working_set.h
 * @brief Definition of working sets.
 **/

#include "sleqp_export.h"
#include "sleqp_func.h"
#include "sleqp_problem.h"
#include "sleqp_types.h"

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct SleqpWorkingSet SleqpWorkingSet;

  /**
   * Creates a new working set
   *
   * @param[out] star            A pointer to the working set to be created
   * @param[int] problem         The underlying problem
   **/
  SLEQP_EXPORT SLEQP_RETCODE sleqp_working_set_create(SleqpWorkingSet** star,
                                                      SleqpProblem* problem);

  /**
   * Resets this working set by removing all variables and constraints from it
   *
   * @param[in]  working_set           The working set
   **/
  SLEQP_EXPORT SLEQP_RETCODE sleqp_working_set_reset(SleqpWorkingSet* working_set);

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
  SLEQP_EXPORT SLEQP_RETCODE sleqp_working_set_add_variable(SleqpWorkingSet* working_set,
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
  SLEQP_EXPORT SLEQP_RETCODE sleqp_working_set_add_constraint(SleqpWorkingSet* working_set,
                                                              int index,
                                                              SLEQP_ACTIVE_STATE state);

  /**
   * Returns the index of the given constraint with respect to the given
   * working set, or @ref SLEQP_NONE if the constraint is not contained
   * in the working set
   *
   * @param[in]  working_set           The working set
   * @param[in]  index                 The constraint index
   **/
  SLEQP_EXPORT int sleqp_working_set_get_constraint_index(const SleqpWorkingSet* working_set,
                                                          int index);

  /**
   * Returns the index of the given variable with respect to the given
   * working set, or @ref SLEQP_NONE if the variable is not contained
   * in the working set
   *
   * @param[in]  working_set           The working set
   * @param[in]  index                 The variable index
   **/
  SLEQP_EXPORT int sleqp_working_set_get_variable_index(const SleqpWorkingSet* working_set,
                                                        int index);

  /**
   * Returns the content of the working set at the given working set index,
   SLEQP_EXPORT * which must be at least zero and less than @ref sleqp_working_set_size
   *
   * @param[in]  working_set           The working set
   * @param[in]  index                 The working set index
   **/
  SLEQP_EXPORT int sleqp_working_set_get_content(const SleqpWorkingSet* working_set,
                                                 int index);

  /**
   * Returns the state of the given variable with respect to the given working set
   *
   * @param[in]  working_set           The working set
   * @param[in]  index                 The variable index
   **/
  SLEQP_EXPORT SLEQP_ACTIVE_STATE sleqp_working_set_get_variable_state(const SleqpWorkingSet* working_set,
                                                                       int index);

  /**
   * Returns the state of the given constraint with respect to the given working set
   *
   * @param[in]  working_set           The working set
   * @param[in]  index                 The constraint index
   **/
  SLEQP_EXPORT SLEQP_ACTIVE_STATE sleqp_working_set_get_constraint_state(const SleqpWorkingSet* working_set,
                                                                         int index);

  /**
   * Returns the problem underling the given working set
   *
   * @param[in]  working_set           The working set
   **/
  SLEQP_EXPORT SleqpProblem* sleqp_working_set_get_problem(const SleqpWorkingSet* working_set);

  /**
   * Returns the number of variables contained in the given working set
   *
   * @param[in]  working_set           The working set
   **/
  SLEQP_EXPORT int sleqp_working_set_num_active_vars(const SleqpWorkingSet* working_set);

  /**
   * Returns the number of constraints contained in the given working set
   *
   * @param[in]  working_set           The working set
   **/
  SLEQP_EXPORT int sleqp_working_set_num_active_cons(const SleqpWorkingSet* working_set);

  /**
   * Returns the size of the given set, i.e. the number of contained variables plus
   * the number of contained constraints
   *
   * @param[in]  working_set           The working set
   **/
  SLEQP_EXPORT int sleqp_working_set_size(const SleqpWorkingSet* working_set);

  /**
   * Returns whether the given working set is *valid*, i.e., whether
   * - all interal indices are consistent
   * - the internally stored sizes are consistent
   * - the intenal indices are consistent with variables states
   *
   * @param[in]  working_set           The working set
   **/
  SLEQP_EXPORT bool sleqp_working_set_valid(const SleqpWorkingSet* working_set);

  /**
   * Prints the given working set to the given file
   *
   * @param[in]  vec     A pointer to the vector
   * @param[in]  output  A pointer to an output `FILE*`
   **/
  SLEQP_EXPORT SLEQP_RETCODE sleqp_working_set_fprintf(const SleqpWorkingSet* working_set,
                                                       FILE* output);

  /**
   * Copies one working set to another
   *
   * @param[in]  source  A pointer to the copy source
   * @param[in]  target  A pointer to the copy target
   **/
  SLEQP_EXPORT SLEQP_RETCODE sleqp_working_set_copy(const SleqpWorkingSet* source,
                                                    SleqpWorkingSet* target);

  SLEQP_EXPORT SLEQP_RETCODE sleqp_working_set_capture(SleqpWorkingSet* working_set);

  SLEQP_EXPORT SLEQP_RETCODE sleqp_working_set_release(SleqpWorkingSet** star);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_WORKING_SET_H */
