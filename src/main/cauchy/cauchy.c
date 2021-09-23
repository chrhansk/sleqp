#include "cauchy.h"

#include "mem.h"

struct SleqpCauchy
{
  int refcount;
  SleqpCauchyCallbacks callbacks;
  void* cauchy_data;
};

SLEQP_RETCODE sleqp_cauchy_create(SleqpCauchy** star,
                                  SleqpCauchyCallbacks* callbacks,
                                  void* cauchy_data)
{
  SLEQP_CALL(sleqp_malloc(star));

  SleqpCauchy* cauchy = *star;

  *cauchy = (SleqpCauchy){0};

  cauchy->refcount = 1;

  cauchy->callbacks = *callbacks;

  cauchy->cauchy_data = cauchy_data;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_cauchy_set_iterate(SleqpCauchy* cauchy,
                                       SleqpIterate* iterate,
                                       double trust_radius)
{
  return cauchy->callbacks.set_iterate(iterate,
                                       trust_radius,
                                       cauchy->cauchy_data);
}

SLEQP_RETCODE sleqp_cauchy_set_trust_radius(SleqpCauchy* cauchy,
                                            double trust_radius)
{
  return cauchy->callbacks.set_trust_radius(trust_radius,
                                            cauchy->cauchy_data);
}

SLEQP_RETCODE sleqp_cauchy_solve(SleqpCauchy* cauchy,
                                 SleqpSparseVec* gradient,
                                 double penalty,
                                 SLEQP_CAUCHY_OBJECTIVE_TYPE objective_type)
{
  return cauchy->callbacks.solve(gradient,
                                 penalty,
                                 objective_type,
                                 cauchy->cauchy_data);
}

SLEQP_RETCODE sleqp_cauchy_get_objective_value(SleqpCauchy* cauchy,
                                               double* objective_value)
{
  return cauchy->callbacks.get_objective_value(objective_value,
                                               cauchy->cauchy_data);
}

SLEQP_RETCODE sleqp_cauchy_get_working_set(SleqpCauchy* cauchy,
                                           SleqpIterate* iterate)
{
  return cauchy->callbacks.get_working_set(iterate,
                                           cauchy->cauchy_data);
}

SLEQP_RETCODE sleqp_cauchy_get_direction(SleqpCauchy* cauchy,
                                         SleqpSparseVec* direction)
{
  return cauchy->callbacks.get_direction(direction,
                                         cauchy->cauchy_data);
}

SLEQP_RETCODE sleqp_cauchy_locally_infeasible(SleqpCauchy* cauchy,
                                              bool* locally_infeasible)
{
  return cauchy->callbacks.locally_infeasible(locally_infeasible,
                                              cauchy->cauchy_data);
}

SLEQP_RETCODE sleqp_cauchy_get_dual_estimation(SleqpCauchy* cauchy,
                                               SleqpIterate* iterate)
{
  return cauchy->callbacks.get_dual_estimation(iterate,
                                               cauchy->cauchy_data);
}

SLEQP_RETCODE sleqp_cauchy_get_violation(SleqpCauchy* cauchy,
                                         double* violation)
{
  return cauchy->callbacks.get_violation(violation,
                                         cauchy->cauchy_data);
}

SLEQP_RETCODE sleqp_cauchy_get_basis_condition(SleqpCauchy* cauchy,
                                               bool* exact,
                                               double* condition)
{
  return cauchy->callbacks.get_basis_condition(exact,
                                               condition,
                                               cauchy->cauchy_data);
}

SLEQP_RETCODE sleqp_cauchy_capture(SleqpCauchy* cauchy)
{
  ++cauchy->refcount;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE cauchy_free(SleqpCauchy** star)
{
  SleqpCauchy* cauchy = *star;

  SLEQP_CALL(cauchy->callbacks.free(cauchy->cauchy_data));

  sleqp_free(star);

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_cauchy_release(SleqpCauchy** star)
{
  SleqpCauchy* cauchy = *star;

  if(!cauchy)
  {
    return SLEQP_OKAY;
  }

  if(--cauchy->refcount == 0)
  {
    SLEQP_CALL(cauchy_free(star));
  }

  *star = NULL;

  return SLEQP_OKAY;
}
