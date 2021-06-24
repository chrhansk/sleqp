#ifndef SLEQP_PREPROCESSING_STATE_H
#define SLEQP_PREPROCESSING_STATE_H

#include "sleqp_problem.h"

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct
  {
    int constraint;
    int variable;
    double factor;
    double var_lb;
    double var_ub;

    enum
    {
      SLEQP_BOUND_STATE_LOWER = (1 << 1),
      SLEQP_BOUND_STATE_UPPER = (1 << 2),
      SLEQP_BOUND_STATE_BOTH  = (SLEQP_BOUND_STATE_LOWER | SLEQP_BOUND_STATE_UPPER)
    } state;

  } SleqpConvertedBound;

  typedef struct
  {
    enum
    {
      SLEQP_VAR_UNCHANGED,
      SLEQP_VAR_BOUNDFIXED,
    } state;

    double value;

  } SleqpVariableState;

  typedef struct
  {
    enum
    {
      SLEQP_CONS_UNCHANGED,
      SLEQP_CONS_REDUNDANT,
      SLEQP_CONS_BOUNDCONVERTED
    } state;

    int bound;

  } SleqpConstraintState;

  typedef enum
  {
    SleqpLowerBound = (1 << 1),
    SleqpUpperBound = (1 << 2),
    SleqpBothBounds = (SleqpLowerBound | SleqpUpperBound)

  } SleqpBoundState;


  typedef struct SleqpPreprocessingState SleqpPreprocessingState;

  SLEQP_RETCODE sleqp_preprocessing_state_create(SleqpPreprocessingState** star,
                                                 SleqpProblem* problem);

  SLEQP_RETCODE sleqp_preprocessing_state_reset(SleqpPreprocessingState* state);

  SLEQP_RETCODE sleqp_preprocessing_state_convert_linear_constraint_to_bound(SleqpPreprocessingState* state,
                                                                             int constraint,
                                                                             int variable,
                                                                             double factor,
                                                                             double var_lb,
                                                                             double var_ub,
                                                                             SleqpBoundState bound_state);

  SLEQP_RETCODE sleqp_preprocessing_state_remove_linear_constraint(SleqpPreprocessingState* state,
                                                                   int constraint);

  SLEQP_RETCODE sleqp_preprocessing_state_fix_variable(SleqpPreprocessingState* state,
                                                       int variable,
                                                       double value);

  SLEQP_RETCODE sleqp_preprocessing_state_converted_bounds(SleqpPreprocessingState* state,
                                                           SleqpConvertedBound** star,
                                                           int* num_converted_bounds);

  SleqpVariableState* sleqp_preprocessing_state_variable_states(const SleqpPreprocessingState* state);

  SleqpConstraintState* sleqp_preprocessing_state_linear_constraint_states(const SleqpPreprocessingState* state);

  int sleqp_preprocessing_state_num_fixed_variables(const SleqpPreprocessingState* state);

  int sleqp_preprocessing_state_num_removed_linear_constraints(const SleqpPreprocessingState* state);

  SleqpProblem* sleqp_preprocessing_state_get_problem(const SleqpPreprocessingState* state);

  SLEQP_RETCODE sleqp_preprocessing_state_capture(SleqpPreprocessingState* state);

  SLEQP_RETCODE sleqp_preprocessing_state_release(SleqpPreprocessingState** star);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_PREPROCESSING_STATE_H */
