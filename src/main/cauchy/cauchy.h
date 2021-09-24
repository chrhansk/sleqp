#ifndef SLEQP_CAUCHY_H
#define SLEQP_CAUCHY_H

#include "cauchy_types.h"

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct SleqpCauchy SleqpCauchy;

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_cauchy_create(SleqpCauchy** cauchy,
                                    SleqpCauchyCallbacks* callbacks,
                                    void* cauchy_data);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_cauchy_set_iterate(SleqpCauchy* cauchy,
                                         SleqpIterate* iterate,
                                         double trust_radius);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_cauchy_set_trust_radius(SleqpCauchy* cauchy,
                                              double trust_radius);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_cauchy_solve(SleqpCauchy* cauchy,
                                   SleqpSparseVec* gradient,
                                   double penalty,
                                   SLEQP_CAUCHY_OBJECTIVE_TYPE objective_type);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_cauchy_get_objective_value(SleqpCauchy* cauchy,
                                                 double* objective_value);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_cauchy_get_working_set(SleqpCauchy* cauchy,
                                             SleqpIterate* iterate);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_cauchy_get_direction(SleqpCauchy* cauchy,
                                           SleqpSparseVec* direction);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_cauchy_locally_infeasible(SleqpCauchy* cauchy,
                                                bool* locally_infeasible);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_cauchy_get_dual_estimation(SleqpCauchy* cauchy,
                                                 SleqpIterate* iterate);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_cauchy_get_violation(SleqpCauchy* cauchy,
                                           double* violation);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_cauchy_get_basis_condition(SleqpCauchy* cauchy,
                                                 bool* exact,
                                                 double* condition);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_cauchy_capture(SleqpCauchy* cauchy);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_cauchy_release(SleqpCauchy** star);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_CAUCHY_H */
