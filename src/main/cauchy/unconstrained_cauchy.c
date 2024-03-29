#include "unconstrained_cauchy.h"

#include "mem.h"
#include "working_set.h"

typedef struct
{
  SleqpProblem* problem;
  SleqpSettings* settings;

  SleqpIterate* iterate;
  double trust_radius;

  SleqpVec* direction;
  double objective;

} CauchyData;

static SLEQP_RETCODE
unconstrained_cauchy_set_iterate(SleqpIterate* iterate,
                                 double trust_radius,
                                 void* data)
{
  CauchyData* cauchy_data = (CauchyData*)data;

  SLEQP_CALL(sleqp_iterate_release(&cauchy_data->iterate));

  SLEQP_CALL(sleqp_iterate_capture(iterate));

  cauchy_data->iterate = iterate;

  cauchy_data->trust_radius = trust_radius;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
unconstrained_cauchy_set_trust_radius(double trust_radius, void* data)
{
  CauchyData* cauchy_data = (CauchyData*)data;

  cauchy_data->trust_radius = trust_radius;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
unconstrained_cauchy_solve(SleqpVec* gradient,
                           double penalty,
                           SLEQP_CAUCHY_OBJTYPE objective_type,
                           void* data)
{
  CauchyData* cauchy_data = (CauchyData*)data;

  SleqpProblem* problem = cauchy_data->problem;
  SleqpIterate* iterate = cauchy_data->iterate;

  SleqpVec* direction = cauchy_data->direction;

  double* objective = &(cauchy_data->objective);

  (*objective) = sleqp_iterate_obj_val(iterate);

  const SleqpVec* grad = sleqp_iterate_obj_grad(iterate);

  const int num_variables   = sleqp_problem_num_vars(problem);
  const double trust_radius = cauchy_data->trust_radius;

  SLEQP_CALL(sleqp_vec_clear(direction));

  int k_g = 0;

  for (int j = 0; j < num_variables; ++j)
  {
    while (k_g < grad->nnz && grad->indices[k_g] < j)
    {
      ++k_g;
    }

    const double g_val
      = (k_g < grad->nnz && grad->indices[k_g] == j) ? grad->data[k_g] : 0.;
    const double ng_val = -g_val;

    if (ng_val >= 0.)
    {
      SLEQP_CALL(sleqp_vec_push(direction, j, trust_radius));

      (*objective) += trust_radius * g_val;
    }
    else
    {
      SLEQP_CALL(sleqp_vec_push(direction, j, -trust_radius));

      (*objective) += (-trust_radius) * g_val;
    }
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
unconstrained_cauchy_obj_val(double* objective_value, void* data)
{
  CauchyData* cauchy_data = (CauchyData*)data;

  *objective_value = cauchy_data->objective;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
unconstrained_cauchy_working_set(SleqpIterate* iterate, void* data)
{
  SleqpWorkingSet* working_set = sleqp_iterate_working_set(iterate);

  SLEQP_CALL(sleqp_working_set_reset(working_set));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
unconstrained_cauchy_lp_step(SleqpVec* direction, void* data)
{
  CauchyData* cauchy_data = (CauchyData*)data;

  SLEQP_CALL(sleqp_vec_copy(cauchy_data->direction, direction));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
unconstrained_cauchy_locally_infeasible(bool* locally_infeasible, void* data)
{
  *locally_infeasible = false;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
unconstrained_cauchy_estimate_duals(const SleqpWorkingSet* working_set,
                                    SleqpVec* cons_dual,
                                    SleqpVec* vars_dual,
                                    void* data)
{
  SLEQP_CALL(sleqp_vec_clear(vars_dual));
  SLEQP_CALL(sleqp_vec_clear(cons_dual));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
unconstrained_cauchy_violation(double* violation, void* data)
{
  *violation = 0.;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
unconstrained_cauchy_set_time_limit(double time_limit, void* data)
{
  return SLEQP_OKAY;
}

static SLEQP_RETCODE
unconstrained_cauchy_basis_condition(bool* exact, double* condition, void* data)
{
  *condition = 1.;
  *exact     = true;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
unconstrained_cauchy_print_stats(double total_elapsed, void* data)
{
  return SLEQP_OKAY;
}

static SLEQP_RETCODE
unconstrained_cauchy_free(void* data)
{
  CauchyData* cauchy_data = (CauchyData*)data;

  SLEQP_CALL(sleqp_vec_free(&cauchy_data->direction));

  SLEQP_CALL(sleqp_iterate_release(&cauchy_data->iterate));

  SLEQP_CALL(sleqp_settings_release(&cauchy_data->settings));

  SLEQP_CALL(sleqp_problem_release(&cauchy_data->problem));

  sleqp_free(&cauchy_data);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
cauchy_data_create(CauchyData** star,
                   SleqpProblem* problem,
                   SleqpSettings* settings)
{
  SLEQP_CALL(sleqp_malloc(star));

  CauchyData* cauchy_data = *star;

  *cauchy_data = (CauchyData){0};

  SLEQP_CALL(sleqp_problem_capture(problem));

  cauchy_data->problem = problem;

  SLEQP_CALL(sleqp_settings_capture(settings));

  cauchy_data->settings = settings;

  const int num_variables = sleqp_problem_num_vars(problem);

  SLEQP_CALL(sleqp_vec_create_full(&cauchy_data->direction, num_variables));

  cauchy_data->trust_radius = SLEQP_NONE;

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_unconstrained_cauchy_create(SleqpCauchy** star,
                                  SleqpProblem* problem,
                                  SleqpSettings* settings)
{
  CauchyData* cauchy_data;

  SLEQP_CALL(cauchy_data_create(&cauchy_data, problem, settings));

  SleqpCauchyCallbacks callbacks
    = {.set_iterate        = unconstrained_cauchy_set_iterate,
       .set_trust_radius   = unconstrained_cauchy_set_trust_radius,
       .solve              = unconstrained_cauchy_solve,
       .obj_val            = unconstrained_cauchy_obj_val,
       .working_set        = unconstrained_cauchy_working_set,
       .lp_step            = unconstrained_cauchy_lp_step,
       .locally_infeasible = unconstrained_cauchy_locally_infeasible,
       .estimate_duals     = unconstrained_cauchy_estimate_duals,
       .violation          = unconstrained_cauchy_violation,
       .set_time_limit     = unconstrained_cauchy_set_time_limit,
       .basis_condition    = unconstrained_cauchy_basis_condition,
       .print_stats        = unconstrained_cauchy_print_stats,
       .free               = unconstrained_cauchy_free};

  SLEQP_CALL(sleqp_cauchy_create(star, &callbacks, (void*)cauchy_data));

  return SLEQP_OKAY;
}
