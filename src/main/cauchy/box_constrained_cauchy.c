#include "box_constrained_cauchy.h"

#include "cmp.h"
#include "mem.h"
#include "working_set.h"

typedef struct
{
  SleqpProblem* problem;
  SleqpParams* params;

  SleqpIterate* iterate;
  double trust_radius;

  SleqpVec* lower_diff;
  SleqpVec* upper_diff;

  SleqpVec* direction;
  SleqpVec* duals;

  SLEQP_ACTIVE_STATE* var_states;

  bool* fixed_vars;

  double objective;

} CauchyData;

static SLEQP_RETCODE
compute_fixed_vars(CauchyData* cauchy_data)
{
  SleqpProblem* problem = cauchy_data->problem;

  const int num_variables = sleqp_problem_num_vars(problem);

  const SleqpVec* var_lb = sleqp_problem_vars_lb(problem);
  const SleqpVec* var_ub = sleqp_problem_vars_ub(problem);

  int k_l = 0, k_u = 0;

  for (int j = 0; j < num_variables; ++j)
  {
    while (k_l < var_lb->nnz && var_lb->indices[k_l] < j)
    {
      ++k_l;
    }

    while (k_u < var_ub->nnz && var_ub->indices[k_u] < j)
    {
      ++k_u;
    }

    const double l_val = (k_l < var_lb->nnz && var_lb->indices[k_l] == j)
                           ? var_lb->data[k_l]
                           : 0.;
    const double u_val = (k_u < var_ub->nnz && var_ub->indices[k_u] == j)
                           ? var_ub->data[k_u]
                           : 0.;

    if (sleqp_is_finite(l_val) && sleqp_is_finite(u_val) && l_val == u_val)
    {
      cauchy_data->fixed_vars[j] = true;
    }
    else
    {
      cauchy_data->fixed_vars[j] = false;
    }
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
compute_diffs(CauchyData* cauchy_data)
{
  SleqpProblem* problem = cauchy_data->problem;
  SleqpIterate* iterate = cauchy_data->iterate;

  const SleqpVec* var_lb = sleqp_problem_vars_lb(problem);
  const SleqpVec* var_ub = sleqp_problem_vars_ub(problem);
  SleqpVec* primal       = sleqp_iterate_primal(iterate);

  assert(sleqp_vec_is_boxed(primal, var_lb, var_ub));

  const double zero_eps
    = sleqp_params_value(cauchy_data->params, SLEQP_PARAM_ZERO_EPS);

  SLEQP_CALL(sleqp_vec_add_scaled(var_lb,
                                  primal,
                                  1.,
                                  -1.,
                                  zero_eps,
                                  cauchy_data->lower_diff));

  SLEQP_CALL(sleqp_vec_add_scaled(var_ub,
                                  primal,
                                  1.,
                                  -1.,
                                  zero_eps,
                                  cauchy_data->upper_diff));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
box_constrained_cauchy_set_iterate(SleqpIterate* iterate,
                                   double trust_radius,
                                   void* data)
{
  CauchyData* cauchy_data = (CauchyData*)data;

  SLEQP_CALL(sleqp_iterate_release(&cauchy_data->iterate));

  SLEQP_CALL(sleqp_iterate_capture(iterate));

  cauchy_data->iterate = iterate;

  cauchy_data->trust_radius = trust_radius;

  SLEQP_CALL(compute_diffs(cauchy_data));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
box_constrained_cauchy_set_trust_radius(double trust_radius, void* data)
{
  CauchyData* cauchy_data = (CauchyData*)data;

  cauchy_data->trust_radius = trust_radius;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
box_constrained_cauchy_solve(SleqpVec* gradient,
                             double penalty,
                             SLEQP_CAUCHY_OBJTYPE objective_type,
                             void* data)
{
  CauchyData* cauchy_data = (CauchyData*)data;
  SleqpProblem* problem   = cauchy_data->problem;
  SleqpIterate* iterate   = cauchy_data->iterate;

  const int num_variables   = sleqp_problem_num_vars(problem);
  const double trust_radius = cauchy_data->trust_radius;

  SLEQP_CALL(sleqp_vec_clear(cauchy_data->direction));
  SLEQP_CALL(sleqp_vec_clear(cauchy_data->duals));

  const SleqpVec* ld   = cauchy_data->lower_diff;
  const SleqpVec* ud   = cauchy_data->upper_diff;
  const SleqpVec* grad = sleqp_iterate_obj_grad(iterate);

  SleqpVec* direction = cauchy_data->direction;
  SleqpVec* duals     = cauchy_data->duals;

  double* objective = &(cauchy_data->objective);
  bool* fixed_vars  = cauchy_data->fixed_vars;

  (*objective) = sleqp_iterate_obj_val(iterate);

  int k_l = 0, k_u = 0, k_g = 0;

  for (int j = 0; j < num_variables; ++j)
  {
    while (k_l < ld->nnz && ld->indices[k_l] < j)
    {
      ++k_l;
    }

    while (k_u < ud->nnz && ud->indices[k_u] < j)
    {
      ++k_u;
    }

    while (k_g < grad->nnz && grad->indices[k_g] < j)
    {
      ++k_g;
    }

    const double l_val
      = (k_l < ld->nnz && ld->indices[k_l] == j) ? ld->data[k_l] : 0.;
    const double u_val
      = (k_u < ud->nnz && ud->indices[k_u] == j) ? ud->data[k_u] : 0.;
    const double g_val
      = (k_g < grad->nnz && grad->indices[k_g] == j) ? grad->data[k_g] : 0.;
    const double ng_val = -g_val;

    if (ng_val >= 0.)
    {
      if (trust_radius < u_val)
      {
        cauchy_data->var_states[j] = SLEQP_INACTIVE;

        SLEQP_CALL(sleqp_vec_push(direction, j, trust_radius));

        (*objective) += trust_radius * g_val;
      }
      else
      {
        cauchy_data->var_states[j]
          = (fixed_vars[j]) ? SLEQP_ACTIVE_BOTH : SLEQP_ACTIVE_UPPER;

        SLEQP_CALL(sleqp_vec_push(direction, j, u_val));

        SLEQP_CALL(sleqp_vec_push(duals, j, ng_val));

        (*objective) += u_val * g_val;
      }
    }
    else
    {
      if (l_val < -trust_radius)
      {
        cauchy_data->var_states[j] = SLEQP_INACTIVE;

        SLEQP_CALL(sleqp_vec_push(direction, j, -trust_radius));

        (*objective) += (-trust_radius) * g_val;
      }
      else
      {
        cauchy_data->var_states[j]
          = (fixed_vars[j]) ? SLEQP_ACTIVE_BOTH : SLEQP_ACTIVE_LOWER;

        SLEQP_CALL(sleqp_vec_push(direction, j, l_val));

        SLEQP_CALL(sleqp_vec_push(duals, j, ng_val));

        (*objective) += l_val * g_val;
      }
    }
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
box_constrained_cauchy_obj_val(double* objective_value, void* data)
{
  CauchyData* cauchy_data = (CauchyData*)data;

  *objective_value = cauchy_data->objective;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
box_constrained_cauchy_working_set(SleqpIterate* iterate, void* data)
{
  CauchyData* cauchy_data      = (CauchyData*)data;
  SleqpProblem* problem        = cauchy_data->problem;
  SleqpWorkingSet* working_set = sleqp_iterate_working_set(iterate);

  SLEQP_CALL(sleqp_working_set_reset(working_set));

  const int num_variables = sleqp_problem_num_vars(problem);

  for (int j = 0; j < num_variables; ++j)
  {
    SLEQP_ACTIVE_STATE state = cauchy_data->var_states[j];

    if (state != SLEQP_INACTIVE)
    {
      SLEQP_CALL(sleqp_working_set_add_var(working_set, j, state));
    }
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
box_constrained_cauchy_lp_step(SleqpVec* direction, void* data)
{
  CauchyData* cauchy_data = (CauchyData*)data;

  SLEQP_CALL(sleqp_vec_copy(cauchy_data->direction, direction));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
box_constrained_cauchy_locally_infeasible(bool* locally_infeasible, void* data)
{
  *locally_infeasible = false;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
box_constrained_cauchy_estimate_duals(const SleqpWorkingSet* working_set,
                                      SleqpVec* cons_dual,
                                      SleqpVec* vars_dual,
                                      void* data)
{
  CauchyData* cauchy_data = (CauchyData*)data;

  SLEQP_CALL(sleqp_vec_copy(cauchy_data->duals, vars_dual));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
box_constrained_cauchy_violation(double* violation, void* data)
{
  *violation = 0.;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
box_constrained_cauchy_set_time_limit(double time_limit, void* data)
{
  return SLEQP_OKAY;
}

static SLEQP_RETCODE
box_constrained_cauchy_basis_condition(bool* exact,
                                       double* condition,
                                       void* data)
{
  *condition = 1.;
  *exact     = true;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
box_constrained_cauchy_print_stats(double total_elapsed, void* data)
{
  return SLEQP_OKAY;
}

static SLEQP_RETCODE
box_constrained_cauchy_free(void* data)
{
  CauchyData* cauchy_data = (CauchyData*)data;

  sleqp_free(&cauchy_data->fixed_vars);

  sleqp_free(&cauchy_data->var_states);

  SLEQP_CALL(sleqp_vec_free(&cauchy_data->duals));

  SLEQP_CALL(sleqp_vec_free(&cauchy_data->direction));

  SLEQP_CALL(sleqp_vec_free(&cauchy_data->upper_diff));

  SLEQP_CALL(sleqp_vec_free(&cauchy_data->lower_diff));

  SLEQP_CALL(sleqp_iterate_release(&cauchy_data->iterate));

  SLEQP_CALL(sleqp_params_release(&cauchy_data->params));

  SLEQP_CALL(sleqp_problem_release(&cauchy_data->problem));

  sleqp_free(&cauchy_data);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
cauchy_data_create(CauchyData** star,
                   SleqpProblem* problem,
                   SleqpParams* params)
{
  SLEQP_CALL(sleqp_malloc(star));

  CauchyData* cauchy_data = *star;

  assert(sleqp_problem_num_cons(problem) == 0);

  const int num_variables = sleqp_problem_num_vars(problem);

  *cauchy_data = (CauchyData){0};

  SLEQP_CALL(sleqp_problem_capture(problem));

  cauchy_data->problem = problem;

  SLEQP_CALL(sleqp_params_capture(params));

  cauchy_data->params = params;

  cauchy_data->trust_radius = SLEQP_NONE;

  SLEQP_CALL(sleqp_vec_create_empty(&cauchy_data->lower_diff, num_variables));

  SLEQP_CALL(sleqp_vec_create_empty(&cauchy_data->upper_diff, num_variables));

  SLEQP_CALL(sleqp_vec_create_full(&cauchy_data->direction, num_variables));

  SLEQP_CALL(sleqp_vec_create_full(&cauchy_data->duals, num_variables));

  SLEQP_CALL(sleqp_alloc_array(&cauchy_data->var_states, num_variables));

  SLEQP_CALL(sleqp_alloc_array(&cauchy_data->fixed_vars, num_variables));

  SLEQP_CALL(compute_fixed_vars(cauchy_data));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_box_constrained_cauchy_create(SleqpCauchy** star,
                                    SleqpProblem* problem,
                                    SleqpParams* params)
{
  CauchyData* cauchy_data;

  SLEQP_CALL(cauchy_data_create(&cauchy_data, problem, params));

  SleqpCauchyCallbacks callbacks
    = {.set_iterate        = box_constrained_cauchy_set_iterate,
       .set_trust_radius   = box_constrained_cauchy_set_trust_radius,
       .solve              = box_constrained_cauchy_solve,
       .obj_val            = box_constrained_cauchy_obj_val,
       .working_set        = box_constrained_cauchy_working_set,
       .lp_step            = box_constrained_cauchy_lp_step,
       .locally_infeasible = box_constrained_cauchy_locally_infeasible,
       .estimate_duals     = box_constrained_cauchy_estimate_duals,
       .violation          = box_constrained_cauchy_violation,
       .set_time_limit     = box_constrained_cauchy_set_time_limit,
       .basis_condition    = box_constrained_cauchy_basis_condition,
       .print_stats        = box_constrained_cauchy_print_stats,
       .free               = box_constrained_cauchy_free};

  SLEQP_CALL(sleqp_cauchy_create(star, &callbacks, (void*)cauchy_data));

  return SLEQP_OKAY;
}
