#include "box_constrained_cauchy.h"

#include "cmp.h"
#include "mem.h"
#include "working_set.h"

typedef struct {
  SleqpProblem* problem;
  SleqpParams* params;

  SleqpIterate* iterate;
  double trust_radius;

  SleqpSparseVec* lower_diff;
  SleqpSparseVec* upper_diff;

  SleqpSparseVec* direction;
  SleqpSparseVec* duals;

  SLEQP_ACTIVE_STATE* var_states;

  bool* fixed_vars;

  double objective;

} CauchyData;

static SLEQP_RETCODE
compute_fixed_vars(CauchyData* cauchy_data)
{
  SleqpProblem* problem = cauchy_data->problem;

  const int num_variables = sleqp_problem_num_variables(problem);

  const SleqpSparseVec* var_lb = sleqp_problem_var_lb(problem);
  const SleqpSparseVec* var_ub = sleqp_problem_var_ub(problem);

  int k_l = 0, k_u = 0;

  for(int j = 0; j < num_variables; ++j)
  {
    while(k_l < var_lb->nnz && var_lb->indices[k_l] < j)
    {
      ++k_l;
    }

    while(k_u < var_ub->nnz && var_ub->indices[k_u] < j)
    {
      ++k_u;
    }

    const double l_val = (k_l < var_lb->nnz && var_lb->indices[k_l] == j) ? var_lb->data[k_l] : 0.;
    const double u_val = (k_u < var_ub->nnz && var_ub->indices[k_u] == j) ? var_ub->data[k_u] : 0.;

    if(sleqp_is_finite(l_val) &&
       sleqp_is_finite(u_val) &&
       l_val == u_val)
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

  const SleqpSparseVec* var_lb = sleqp_problem_var_lb(problem);
  const SleqpSparseVec* var_ub = sleqp_problem_var_ub(problem);
  SleqpSparseVec* primal = sleqp_iterate_get_primal(iterate);

  assert(sleqp_sparse_vector_is_boxed(primal, var_lb, var_ub));

  const double zero_eps = sleqp_params_get(cauchy_data->params,
                                           SLEQP_PARAM_ZERO_EPS);

  SLEQP_CALL(sleqp_sparse_vector_add_scaled(var_lb,
                                            primal,
                                            1.,
                                            -1.,
                                            zero_eps,
                                            cauchy_data->lower_diff));

  SLEQP_CALL(sleqp_sparse_vector_add_scaled(var_ub,
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
  CauchyData* cauchy_data = (CauchyData* ) data;

  SLEQP_CALL(sleqp_iterate_release(&cauchy_data->iterate));

  SLEQP_CALL(sleqp_iterate_capture(iterate));

  cauchy_data->iterate = iterate;

  cauchy_data->trust_radius = trust_radius;

  SLEQP_CALL(compute_diffs(cauchy_data));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
box_constrained_cauchy_set_trust_radius(double trust_radius,
                                        void* data)
{
  CauchyData* cauchy_data = (CauchyData* ) data;

  cauchy_data->trust_radius = trust_radius;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
box_constrained_cauchy_solve(SleqpSparseVec* gradient,
                             double penalty,
                             SLEQP_CAUCHY_OBJECTIVE_TYPE objective_type,
                             void* data)
{
  CauchyData* cauchy_data = (CauchyData* ) data;
  SleqpProblem* problem = cauchy_data->problem;
  SleqpIterate* iterate = cauchy_data->iterate;

  const int num_variables = sleqp_problem_num_variables(problem);
  const double trust_radius = cauchy_data->trust_radius;

  SLEQP_CALL(sleqp_sparse_vector_clear(cauchy_data->direction));
  SLEQP_CALL(sleqp_sparse_vector_clear(cauchy_data->duals));

  const SleqpSparseVec* ld = cauchy_data->lower_diff;
  const SleqpSparseVec* ud = cauchy_data->upper_diff;
  const SleqpSparseVec* grad = sleqp_iterate_get_func_grad(iterate);

  SleqpSparseVec* direction = cauchy_data->direction;
  SleqpSparseVec* duals = cauchy_data->duals;

  double* objective = &(cauchy_data->objective);
  bool* fixed_vars = cauchy_data->fixed_vars;

  (*objective) = sleqp_iterate_get_func_val(iterate);

  int k_l = 0, k_u = 0, k_g = 0;

  for(int j = 0; j < num_variables; ++j)
  {
    while(k_l < ld->nnz && ld->indices[k_l] < j)
    {
      ++k_l;
    }

    while(k_u < ud->nnz && ud->indices[k_u] < j)
    {
      ++k_u;
    }

    while(k_g < grad->nnz && grad->indices[k_g] < j)
    {
      ++k_g;
    }

    const double l_val = (k_l < ld->nnz && ld->indices[k_l] == j) ? ld->data[k_l] : 0.;
    const double u_val = (k_u < ud->nnz && ud->indices[k_u] == j) ? ud->data[k_u] : 0.;
    const double g_val = (k_g < grad->nnz && grad->indices[k_g] == j) ? grad->data[k_g] : 0.;
    const double ng_val = -g_val;

    if(ng_val >= 0.)
    {
      if(trust_radius < u_val)
      {
        cauchy_data->var_states[j] = SLEQP_INACTIVE;

        SLEQP_CALL(sleqp_sparse_vector_push(direction,
                                            j,
                                            trust_radius));

        (*objective) += trust_radius * g_val;
      }
      else
      {
        cauchy_data->var_states[j] = (fixed_vars[j]) ? SLEQP_ACTIVE_BOTH : SLEQP_ACTIVE_UPPER;

        SLEQP_CALL(sleqp_sparse_vector_push(direction,
                                            j,
                                            u_val));

        SLEQP_CALL(sleqp_sparse_vector_push(duals,
                                            j,
                                            ng_val));

        (*objective) += u_val * g_val;
      }

    }
    else
    {
      if(l_val < -trust_radius)
      {
        cauchy_data->var_states[j] = SLEQP_INACTIVE;

        SLEQP_CALL(sleqp_sparse_vector_push(direction,
                                            j,
                                            -trust_radius));

        (*objective) += (-trust_radius) * g_val;
      }
      else
      {
        cauchy_data->var_states[j] = (fixed_vars[j]) ? SLEQP_ACTIVE_BOTH : SLEQP_ACTIVE_LOWER;

        SLEQP_CALL(sleqp_sparse_vector_push(direction,
                                            j,
                                            l_val));

        SLEQP_CALL(sleqp_sparse_vector_push(duals,
                                            j,
                                            ng_val));

        (*objective) += l_val * g_val;
      }

    }
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
box_constrained_cauchy_get_objective_value(double* objective_value,
                                           void* data)
{
  CauchyData* cauchy_data = (CauchyData* ) data;

  *objective_value = cauchy_data->objective;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
box_constrained_cauchy_get_working_set(SleqpIterate* iterate,
                                       void* data)
{
  CauchyData* cauchy_data = (CauchyData* ) data;
  SleqpProblem* problem = cauchy_data->problem;
  SleqpWorkingSet* working_set = sleqp_iterate_get_working_set(iterate);

  SLEQP_CALL(sleqp_working_set_reset(working_set));

  const int num_variables = sleqp_problem_num_variables(problem);

  for(int j = 0; j < num_variables; ++j)
  {
    SLEQP_ACTIVE_STATE state = cauchy_data->var_states[j];

    if(state != SLEQP_INACTIVE)
    {
      SLEQP_CALL(sleqp_working_set_add_variable(working_set,
                                                j,
                                                state));
    }

  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
box_constrained_cauchy_get_direction(SleqpSparseVec* direction,
                                     void* data)
{
  CauchyData* cauchy_data = (CauchyData* ) data;

  SLEQP_CALL(sleqp_sparse_vector_copy(cauchy_data->direction,
                                      direction));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
box_constrained_cauchy_locally_infeasible(bool* locally_infeasible,
                                          void* data)
{
  *locally_infeasible = false;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
box_constrained_cauchy_get_dual_estimation(SleqpIterate* iterate,
                                           void* data)
{
  CauchyData* cauchy_data = (CauchyData*) data;

  SLEQP_CALL(sleqp_sparse_vector_copy(cauchy_data->duals,
                                      sleqp_iterate_get_vars_dual(iterate)));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
box_constrained_cauchy_get_violation(double* violation,
                                     void* data)
{
  *violation = 0.;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
box_constrained_cauchy_get_basis_condition(bool* exact,
                                           double* condition,
                                           void* data)
{
  *condition = 1.;
  *exact = true;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
box_constrained_cauchy_free(void* data)
{
  CauchyData* cauchy_data = (CauchyData* ) data;

  sleqp_free(&cauchy_data->fixed_vars);

  sleqp_free(&cauchy_data->var_states);

  SLEQP_CALL(sleqp_sparse_vector_free(&cauchy_data->duals));

  SLEQP_CALL(sleqp_sparse_vector_free(&cauchy_data->direction));

  SLEQP_CALL(sleqp_sparse_vector_free(&cauchy_data->upper_diff));

  SLEQP_CALL(sleqp_sparse_vector_free(&cauchy_data->lower_diff));

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

  assert(sleqp_problem_num_constraints(problem) == 0);

  const int num_variables = sleqp_problem_num_variables(problem);

  *cauchy_data = (CauchyData) {0};

  SLEQP_CALL(sleqp_problem_capture(problem));

  cauchy_data->problem = problem;

  SLEQP_CALL(sleqp_params_capture(params));

  cauchy_data->params = params;

  cauchy_data->trust_radius = SLEQP_NONE;

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&cauchy_data->lower_diff,
                                              num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&cauchy_data->upper_diff,
                                              num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_full(&cauchy_data->direction,
                                             num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_full(&cauchy_data->duals,
                                             num_variables));

  SLEQP_CALL(sleqp_alloc_array(&cauchy_data->var_states,
                               num_variables));

  SLEQP_CALL(sleqp_alloc_array(&cauchy_data->fixed_vars,
                               num_variables));

  SLEQP_CALL(compute_fixed_vars(cauchy_data));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_box_constrained_cauchy_create(SleqpCauchy** star,
                                                  SleqpProblem* problem,
                                                  SleqpParams* params)
{
  CauchyData* cauchy_data;

  SLEQP_CALL(cauchy_data_create(&cauchy_data, problem, params));

  SleqpCauchyCallbacks callbacks = {
    .set_iterate         = box_constrained_cauchy_set_iterate,
    .set_trust_radius    = box_constrained_cauchy_set_trust_radius,
    .solve               = box_constrained_cauchy_solve,
    .get_objective_value = box_constrained_cauchy_get_objective_value,
    .get_working_set     = box_constrained_cauchy_get_working_set,
    .get_direction       = box_constrained_cauchy_get_direction,
    .locally_infeasible  = box_constrained_cauchy_locally_infeasible,
    .get_dual_estimation = box_constrained_cauchy_get_dual_estimation,
    .get_violation       = box_constrained_cauchy_get_violation,
    .get_basis_condition = box_constrained_cauchy_get_basis_condition,
    .free                = box_constrained_cauchy_free
  };

  SLEQP_CALL(sleqp_cauchy_create(star,
                                 &callbacks,
                                 (void*) cauchy_data));

  return SLEQP_OKAY;
}
