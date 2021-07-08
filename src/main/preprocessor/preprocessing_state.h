#ifndef SLEQP_PREPROCESSING_STATE_H
#define SLEQP_PREPROCESSING_STATE_H

#include "problem.h"

#ifdef __cplusplus
extern "C" {
#endif

  typedef enum
  {
    SLEQP_LOWER_BOUND = (1 << 1),
    SLEQP_UPPER_BOUND = (1 << 2),
    SLEQP_BOTH_BOUNDS = (SLEQP_LOWER_BOUND | SLEQP_UPPER_BOUND)
  } SleqpBoundState;

  typedef struct
  {
    int constraint;
    int variable;
    double factor;
    double var_lb;
    double var_ub;

    SleqpBoundState state;

  } SleqpConvertedBound;

  typedef struct
  {
    int constraint;

    int* variables;
    double* factors;
    int num_variables;

    SleqpBoundState state;

  } SleqpForcingConstraint;

  typedef struct
  {
    enum
    {
      SLEQP_VAR_UNCHANGED,
      SLEQP_VAR_BOUND_FIXED,
      SLEQP_VAR_FORCING_FIXED,
      SLEQP_VAR_FIXED = (SLEQP_VAR_BOUND_FIXED | SLEQP_VAR_FORCING_FIXED)
    } state;

    double value;

  } SleqpVariableState;

  typedef struct
  {
    enum
    {
      SLEQP_CONS_UNCHANGED,
      SLEQP_CONS_REDUNDANT,
      SLEQP_CONS_BOUNDCONVERTED,
      SLEQP_CONS_FORCING,
    } state;

    int bound;

  } SleqpConstraintState;


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

  SLEQP_RETCODE sleqp_preprocessing_state_add_forcing_constraint(SleqpPreprocessingState* state,
                                                                 int constraint,
                                                                 SleqpBoundState bound_state,
                                                                 double* var_lb,
                                                                 double* var_ub);

  SLEQP_RETCODE sleqp_preprocessing_state_remove_linear_constraint(SleqpPreprocessingState* state,
                                                                   int constraint);

  SLEQP_RETCODE sleqp_preprocessing_state_fix_variable_to_bounds(SleqpPreprocessingState* state,
                                                       int variable,
                                                       double value);

  SLEQP_RETCODE sleqp_preprocessing_state_converted_bounds(SleqpPreprocessingState* state,
                                                           SleqpConvertedBound** star,
                                                           int* num_converted_bounds);

  SLEQP_RETCODE sleqp_preprocessing_state_forcing_constraints(SleqpPreprocessingState* state,
                                                              SleqpForcingConstraint** star,
                                                              int* num_forcing_constraints);

  SleqpVariableState* sleqp_preprocessing_state_variable_states(const SleqpPreprocessingState* state);

  SleqpConstraintState* sleqp_preprocessing_state_linear_constraint_states(const SleqpPreprocessingState* state);

  int sleqp_preprocessing_state_num_fixed_variables(const SleqpPreprocessingState* state);

  int sleqp_preprocessing_state_num_removed_linear_constraints(const SleqpPreprocessingState* state);

  SLEQP_RETCODE sleqp_preprocessing_state_flush(SleqpPreprocessingState* state);

  SLEQP_RETCODE sleqp_preprocessing_state_fixed_variables(SleqpPreprocessingState* state,
                                                          int* num_fixed_vars,
                                                          int** fixed_var_indices,
                                                          double** fixed_var_values);

  SLEQP_RETCODE sleqp_preprocessing_state_removed_linear_constraints(SleqpPreprocessingState* state,
                                                                     int* num_removed_cons,
                                                                     int** removed_cons_indices);

  SleqpProblem* sleqp_preprocessing_state_get_problem(const SleqpPreprocessingState* state);

  SLEQP_RETCODE sleqp_preprocessing_state_capture(SleqpPreprocessingState* state);

  SLEQP_RETCODE sleqp_preprocessing_state_release(SleqpPreprocessingState** star);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_PREPROCESSING_STATE_H */
