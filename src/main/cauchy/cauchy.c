#include "cauchy.h"

#include "cmp.h"
#include "feas.h"
#include "mem.h"
#include "pub_iterate.h"

struct SleqpCauchy
{
  int refcount;
  double trust_radius;
  SleqpCauchyCallbacks callbacks;
  void* cauchy_data;
};

SLEQP_RETCODE
sleqp_cauchy_create(SleqpCauchy** star,
                    SleqpCauchyCallbacks* callbacks,
                    void* cauchy_data)
{
  SLEQP_CALL(sleqp_malloc(star));

  SleqpCauchy* cauchy = *star;

  *cauchy = (SleqpCauchy){0};

  cauchy->refcount     = 1;
  cauchy->trust_radius = SLEQP_NONE;

  cauchy->callbacks = *callbacks;

  cauchy->cauchy_data = cauchy_data;

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_cauchy_set_iterate(SleqpCauchy* cauchy,
                         SleqpIterate* iterate,
                         double trust_radius)
{
  cauchy->trust_radius = trust_radius;

  return cauchy->callbacks.set_iterate(iterate,
                                       trust_radius,
                                       cauchy->cauchy_data);
}

SLEQP_RETCODE
sleqp_cauchy_set_trust_radius(SleqpCauchy* cauchy, double trust_radius)
{
  cauchy->trust_radius = trust_radius;

  return cauchy->callbacks.set_trust_radius(trust_radius, cauchy->cauchy_data);
}

SLEQP_RETCODE
sleqp_cauchy_solve(SleqpCauchy* cauchy,
                   SleqpVec* gradient,
                   double penalty,
                   SLEQP_CAUCHY_OBJTYPE objective_type)
{
  assert(cauchy->trust_radius != SLEQP_NONE);

  return cauchy->callbacks.solve(gradient,
                                 penalty,
                                 objective_type,
                                 cauchy->cauchy_data);
}

SLEQP_RETCODE
sleqp_cauchy_obj_val(SleqpCauchy* cauchy, double* objective_value)
{
  return cauchy->callbacks.obj_val(objective_value, cauchy->cauchy_data);
}

SLEQP_RETCODE
sleqp_cauchy_working_set(SleqpCauchy* cauchy, SleqpIterate* iterate)
{
  return cauchy->callbacks.working_set(iterate, cauchy->cauchy_data);
}

SLEQP_RETCODE
sleqp_cauchy_lp_step(SleqpCauchy* cauchy, SleqpVec* direction)
{
  return cauchy->callbacks.lp_step(direction, cauchy->cauchy_data);
}

SLEQP_RETCODE
sleqp_cauchy_locally_infeasible(SleqpCauchy* cauchy, bool* locally_infeasible)
{
  return cauchy->callbacks.locally_infeasible(locally_infeasible,
                                              cauchy->cauchy_data);
}

SLEQP_RETCODE
sleqp_cauchy_estimate_duals(SleqpCauchy* cauchy,
                            const SleqpWorkingSet* working_set,
                            SleqpVec* cons_dual,
                            SleqpVec* vars_dual)
{
  return cauchy->callbacks.estimate_duals(working_set,
                                          cons_dual,
                                          vars_dual,
                                          cauchy->cauchy_data);
}

SLEQP_RETCODE
sleqp_cauchy_violation(SleqpCauchy* cauchy, double* violation)
{
  return cauchy->callbacks.violation(violation, cauchy->cauchy_data);
}

SLEQP_RETCODE
sleqp_cauchy_set_time_limit(SleqpCauchy* cauchy, double time_limit)
{
  return cauchy->callbacks.set_time_limit(time_limit, cauchy->cauchy_data);
}

SLEQP_RETCODE
sleqp_cauchy_basis_condition(SleqpCauchy* cauchy,
                             bool* exact,
                             double* condition)
{
  return cauchy->callbacks.basis_condition(exact,
                                           condition,
                                           cauchy->cauchy_data);
}

SLEQP_RETCODE
sleqp_cauchy_print_stats(SleqpCauchy* cauchy, double total_elapsed)
{
  return cauchy->callbacks.print_stats(total_elapsed, cauchy->cauchy_data);
}

SLEQP_RETCODE
sleqp_cauchy_compute_criticality_bound(SleqpCauchy* cauchy,
                                       double merit_value,
                                       double* criticality_bound)
{
  double objective_value;

  SLEQP_CALL(sleqp_cauchy_obj_val(cauchy, &objective_value));

  const double reduction = merit_value - objective_value;

  *criticality_bound = reduction / SLEQP_MIN(cauchy->trust_radius, 1.);

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_cauchy_capture(SleqpCauchy* cauchy)
{
  ++cauchy->refcount;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
cauchy_free(SleqpCauchy** star)
{
  SleqpCauchy* cauchy = *star;

  SLEQP_CALL(cauchy->callbacks.free(cauchy->cauchy_data));

  sleqp_free(star);

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_cauchy_release(SleqpCauchy** star)
{
  SleqpCauchy* cauchy = *star;

  if (!cauchy)
  {
    return SLEQP_OKAY;
  }

  if (--cauchy->refcount == 0)
  {
    SLEQP_CALL(cauchy_free(star));
  }

  *star = NULL;

  return SLEQP_OKAY;
}
