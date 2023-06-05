#include "standard_cauchy.h"

#include <math.h>

#include "cmp.h"
#include "error.h"
#include "fail.h"
#include "mem.h"
#include "merit.h"
#include "problem.h"
#include "pub_log.h"
#include "timer.h"

typedef enum
{
  NONE        = 0,
  BASE_STATS  = (1 << 0),
  PRIMAL_VALS = (1 << 1),
  DUAL_VALS   = (1 << 2),
  ALL         = (BASE_STATS | PRIMAL_VALS | DUAL_VALS),
} Components;

typedef enum
{
  SLACK_BASIS = SLEQP_CAUCHY_NUM_OBJTYPES
} _;

typedef struct
{
  SleqpProblem* problem;
  SleqpSettings* settings;

  int num_lp_variables;
  int num_lp_constraints;

  double trust_radius;
  SleqpIterate* iterate;

  double time_limit;

  SLEQP_BASESTAT* var_stats;
  SLEQP_BASESTAT* cons_stats;

  bool has_basis[SLEQP_CAUCHY_NUM_OBJTYPES];
  SLEQP_CAUCHY_OBJTYPE current_objective;

  bool first_solve;
  bool has_coefficients;
  bool use_reduced_interface;

  Components dirty_components;

  SleqpLPi* default_interface;
  SleqpLPi* reduced_interface;

  SLEQP_BASESTAT* reduced_cons_stats;

  double* objective;
  double* cons_lb;
  double* cons_ub;
  double* vars_lb;
  double* vars_ub;

  double* primal_values;

  double* vars_dual;
  double* cons_dual;

} CauchyData;

static SLEQP_RETCODE
create_and_set_slack_basis(CauchyData* cauchy_data)
{
  SleqpProblem* problem = cauchy_data->problem;

  SLEQP_BASESTAT* var_stats  = cauchy_data->var_stats;
  SLEQP_BASESTAT* cons_stats = cauchy_data->cons_stats;

  const int num_vars = sleqp_problem_num_vars(problem);
  const int num_cons = sleqp_problem_num_cons(problem);

  SLEQP_BASESTAT* lower_slack_stats = var_stats + num_vars;
  SLEQP_BASESTAT* upper_slack_stats = lower_slack_stats + num_cons;

  double* cons_lb = cauchy_data->cons_lb;
  double* cons_ub = cauchy_data->cons_ub;

  // Set direction variables to lower, i.e., add
  // d_j = 0 to row basis
  for (int j = 0; j < num_vars; ++j)
  {
    var_stats[j] = SLEQP_BASESTAT_LOWER;
  }

  for (int i = 0; i < num_cons; ++i)
  {
    lower_slack_stats[i] = SLEQP_BASESTAT_LOWER;
    upper_slack_stats[i] = SLEQP_BASESTAT_LOWER;

    if (cons_lb[i] > 0.)
    {
      // lower slack (with entry +1) raised from
      // lower (=0) to value > 0
      // row value will be lb(i), achieved by
      // lower slack variable
      lower_slack_stats[i] = SLEQP_BASESTAT_BASIC;
      cons_stats[i]        = SLEQP_BASESTAT_LOWER;
    }
    else if (cons_ub[i] < 0.)
    {
      // same as above, just with opposite sign
      upper_slack_stats[i] = SLEQP_BASESTAT_BASIC;
      cons_stats[i]        = SLEQP_BASESTAT_UPPER;
    }
    else
    {
      // Both slacks set to zero, row
      // set to basic (row value will be zero within [lb(i), ub(i)])
      cons_stats[i] = SLEQP_BASESTAT_BASIC;
    }
  }

  SLEQP_CALL(sleqp_lpi_set_basis(cauchy_data->default_interface,
                                 SLACK_BASIS,
                                 var_stats,
                                 cons_stats));

  SLEQP_CALL(
    sleqp_lpi_restore_basis(cauchy_data->default_interface, SLACK_BASIS));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
cauchy_data_create(CauchyData** star,
                   SleqpProblem* problem,
                   SleqpSettings* settings)
{
  SLEQP_CALL(sleqp_malloc(star));

  const int num_variables = sleqp_problem_num_vars(problem);

  const int num_constraints = sleqp_problem_num_cons(problem);

  CauchyData* data = *star;

  *data = (CauchyData){0};

  data->problem = problem;
  SLEQP_CALL(sleqp_problem_capture(data->problem));

  SLEQP_CALL(sleqp_settings_capture(settings));
  data->settings = settings;

  data->num_lp_variables   = num_variables + 2 * num_constraints;
  data->num_lp_constraints = num_constraints;

  data->trust_radius = SLEQP_NONE;
  data->time_limit   = SLEQP_NONE;

  SLEQP_CALL(sleqp_alloc_array(&data->var_stats, data->num_lp_variables));
  SLEQP_CALL(sleqp_alloc_array(&data->cons_stats, data->num_lp_constraints));

  data->has_basis[SLEQP_CAUCHY_OBJTYPE_DEFAULT] = false;
  data->has_basis[SLEQP_CAUCHY_OBJTYPE_FEAS]    = false;
  data->has_basis[SLEQP_CAUCHY_OBJTYPE_MIXED]   = false;
  data->current_objective                       = SLEQP_NONE;

  data->first_solve      = true;
  data->has_coefficients = false;

  SLEQP_CALL(sleqp_lpi_create_default(&data->default_interface,
                                      data->num_lp_variables,
                                      data->num_lp_constraints,
                                      settings));

  SLEQP_CALL(sleqp_alloc_array(&data->objective, data->num_lp_variables));

  SLEQP_CALL(sleqp_alloc_array(&data->cons_lb, data->num_lp_constraints));
  SLEQP_CALL(sleqp_alloc_array(&data->cons_ub, data->num_lp_constraints));

  SLEQP_CALL(sleqp_alloc_array(&data->vars_lb, data->num_lp_variables));
  SLEQP_CALL(sleqp_alloc_array(&data->vars_ub, data->num_lp_variables));

  SLEQP_CALL(sleqp_alloc_array(&data->primal_values, data->num_lp_variables));

  const int dual_size = data->num_lp_constraints + data->num_lp_variables;

  SLEQP_CALL(sleqp_alloc_array(&data->vars_dual, dual_size));
  data->cons_dual = data->vars_dual + data->num_lp_variables;

  const double inf = sleqp_infinity();

  for (int j = num_variables; j < data->num_lp_variables; ++j)
  {
    data->vars_lb[j] = 0;
    data->vars_ub[j] = inf;
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
append_identities(SleqpMat* cons_jac, int num_variables, int num_constraints)
{
  assert(num_constraints == sleqp_mat_num_rows(cons_jac));
  assert(num_variables == sleqp_mat_num_cols(cons_jac));

  const int nnz      = sleqp_mat_nnz(cons_jac);
  const int num_rows = sleqp_mat_num_rows(cons_jac);
  const int num_cols = sleqp_mat_num_cols(cons_jac);

  /*
   * Reserve a litte more so we can add the two
   * identity matrices afterwards
   */
  SLEQP_CALL(sleqp_mat_reserve(cons_jac, nnz + 2 * num_constraints));

  SLEQP_CALL(
    sleqp_mat_resize(cons_jac, num_rows, num_cols + 2 * num_constraints));

  // append the +I
  for (int i = 0; i < num_constraints; ++i)
  {
    const int con = num_variables + i;

    SLEQP_CALL(sleqp_mat_push_col(cons_jac, con));

    SLEQP_CALL(sleqp_mat_push(cons_jac, i, con, 1.));
  }

  // append the -I
  for (int i = 0; i < num_constraints; ++i)
  {
    const int con = num_variables + num_constraints + i;

    SLEQP_CALL(sleqp_mat_push_col(cons_jac, con));

    SLEQP_CALL(sleqp_mat_push(cons_jac, i, con, -1.));
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
remove_identities(SleqpMat* cons_jac, int num_variables, int num_constraints)
{
  const int num_rows = sleqp_mat_num_rows(cons_jac);
  const int num_cols = sleqp_mat_num_cols(cons_jac);

  assert(num_constraints == num_rows);
  assert(num_variables + 2 * num_constraints == num_cols);

  SLEQP_CALL(sleqp_mat_resize(cons_jac, num_constraints, num_variables));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
create_cons_bounds(CauchyData* cauchy_data,
                   SleqpIterate* iterate,
                   int num_variables,
                   int num_constraints)
{
  int k_c = 0, k_lb = 0, k_ub = 0;

  const SleqpVec* lb = sleqp_problem_cons_lb(cauchy_data->problem);
  const SleqpVec* ub = sleqp_problem_cons_ub(cauchy_data->problem);

  SleqpVec* val = sleqp_iterate_cons_val(iterate);

  const double inf = sleqp_infinity();

  for (int i = 0; i < num_constraints; ++i)
  {
    while (k_c < val->nnz && val->indices[k_c] < i)
    {
      ++k_c;
    }

    while (k_lb < lb->nnz && lb->indices[k_lb] < i)
    {
      ++k_lb;
    }

    while (k_ub < ub->nnz && ub->indices[k_ub] < i)
    {
      ++k_ub;
    }

    const bool valid_ub = (k_ub < ub->nnz && ub->indices[k_ub] == i);
    const bool valid_lb = (k_lb < lb->nnz && lb->indices[k_lb] == i);
    const bool valid_c  = (k_c < val->nnz && val->indices[k_c] == i);

    const double ubval = valid_ub ? ub->data[k_ub] : 0;
    const double lbval = valid_lb ? lb->data[k_lb] : 0;
    const double cval  = valid_c ? val->data[k_c] : 0;

    assert(!sleqp_is_infinite(lbval));
    assert(!sleqp_is_infinite(-ubval));
    assert(sleqp_is_finite(cval));

    if (sleqp_is_infinite(ubval))
    {
      cauchy_data->cons_ub[i] = inf;
    }
    else
    {
      cauchy_data->cons_ub[i] = ubval - cval;
    }

    if (sleqp_is_infinite(-lbval))
    {
      cauchy_data->cons_lb[i] = -inf;
    }
    else
    {
      cauchy_data->cons_lb[i] = lbval - cval;
    }

    assert(cauchy_data->cons_lb[i] <= cauchy_data->cons_ub[i]);
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
create_var_bounds(CauchyData* cauchy_data,
                  SleqpIterate* iterate,
                  int num_variables,
                  int num_constraints)
{
  SleqpVec* x           = sleqp_iterate_primal(iterate);
  SleqpProblem* problem = cauchy_data->problem;

  const SleqpVec* lb = sleqp_problem_vars_lb(problem);
  const SleqpVec* ub = sleqp_problem_vars_ub(problem);

  const double trust_radius = cauchy_data->trust_radius;

  assert(trust_radius > 0.);

  int k_x = 0, k_lb = 0, k_ub = 0;

  for (int j = 0; j < num_variables; ++j)
  {
    while (k_x < x->nnz && x->indices[k_x] < j)
    {
      ++k_x;
    }

    while (k_lb < lb->nnz && lb->indices[k_lb] < j)
    {
      ++k_lb;
    }

    while (k_ub < ub->nnz && ub->indices[k_ub] < j)
    {
      ++k_ub;
    }

    const bool valid_x  = (k_x < x->nnz) && (j == x->indices[k_x]);
    const bool valid_ub = (k_ub < ub->nnz) && (j == ub->indices[k_ub]);
    const bool valid_lb = (k_lb < lb->nnz) && (j == lb->indices[k_lb]);

    const double ubval = valid_ub ? ub->data[k_ub] : 0.;
    const double lbval = valid_lb ? lb->data[k_lb] : 0.;
    const double xval  = valid_x ? x->data[k_x] : 0.;

    assert(!sleqp_is_infinite(lbval));
    assert(!sleqp_is_infinite(-ubval));
    assert(sleqp_is_finite(xval));

    if (sleqp_is_infinite(ubval))
    {
      cauchy_data->vars_ub[j] = trust_radius;
    }
    else
    {
      cauchy_data->vars_ub[j] = SLEQP_MIN(ubval - xval, trust_radius);
    }

    if (sleqp_is_infinite(-lbval))
    {
      cauchy_data->vars_lb[j] = -trust_radius;
    }
    else
    {
      cauchy_data->vars_lb[j] = SLEQP_MAX(lbval - xval, -trust_radius);
    }

    assert(cauchy_data->vars_lb[j] <= cauchy_data->vars_ub[j]);
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
create_default_objective(CauchyData* cauchy_data,
                         SleqpVec* gradient,
                         double penalty)
{
  SleqpProblem* problem = cauchy_data->problem;

  const int num_variables = sleqp_problem_num_vars(problem);

  const int num_lp_variables = cauchy_data->num_lp_variables;

  for (int j = 0; j < num_variables; ++j)
  {
    cauchy_data->objective[j] = 0.;
  }

  for (int j = num_variables; j < num_lp_variables; ++j)
  {
    cauchy_data->objective[j] = penalty;
  }

  if (gradient)
  {
    assert(gradient->dim == num_variables);

    for (int k = 0; k < gradient->nnz; ++k)
    {
      cauchy_data->objective[gradient->indices[k]] = gradient->data[k];
    }
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
set_default_coefficients(CauchyData* cauchy_data, SleqpMat* cons_jac)
{
  SleqpProblem* problem = cauchy_data->problem;

  const int num_variables   = sleqp_mat_num_cols(cons_jac);
  const int num_constraints = sleqp_mat_num_rows(cons_jac);

  bool fixed_jacobian = !(sleqp_problem_has_nonlinear_cons(problem));

  if (fixed_jacobian && cauchy_data->has_coefficients)
  {
    return SLEQP_OKAY;
  }

  SLEQP_CALL(append_identities(cons_jac, num_variables, num_constraints));

  assert(sleqp_mat_is_valid(cons_jac));

  SLEQP_CALL(sleqp_lpi_set_coeffs(cauchy_data->default_interface, cons_jac));

  SLEQP_CALL(remove_identities(cons_jac, num_variables, num_constraints));

  assert(sleqp_mat_is_valid(cons_jac));

  cauchy_data->has_coefficients = true;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
standard_cauchy_set_iterate(SleqpIterate* iterate,
                            double trust_radius,
                            void* data)
{
  CauchyData* cauchy_data = (CauchyData*)data;
  SleqpMat* cons_jac      = sleqp_iterate_cons_jac(iterate);

  const int num_variables   = sleqp_mat_num_cols(cons_jac);
  const int num_constraints = sleqp_mat_num_rows(cons_jac);

  assert(trust_radius > 0.);

  cauchy_data->trust_radius = trust_radius;

  SLEQP_CALL(sleqp_iterate_release(&cauchy_data->iterate));

  cauchy_data->iterate = iterate;

  SLEQP_CALL(sleqp_iterate_capture(cauchy_data->iterate));

  SLEQP_CALL(
    create_var_bounds(cauchy_data, iterate, num_variables, num_constraints));

  SLEQP_CALL(
    create_cons_bounds(cauchy_data, iterate, num_variables, num_constraints));

  SLEQP_CALL(sleqp_lpi_set_bounds(cauchy_data->default_interface,
                                  cauchy_data->cons_lb,
                                  cauchy_data->cons_ub,
                                  cauchy_data->vars_lb,
                                  cauchy_data->vars_ub));

  SLEQP_CALL(set_default_coefficients(cauchy_data, cons_jac));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
standard_cauchy_set_trust_radius(double trust_radius, void* data)
{
  CauchyData* cauchy_data = (CauchyData*)data;

  SleqpMat* cons_jac = sleqp_iterate_cons_jac(cauchy_data->iterate);

  const int num_variables   = sleqp_mat_num_cols(cons_jac);
  const int num_constraints = sleqp_mat_num_rows(cons_jac);

  cauchy_data->trust_radius = trust_radius;

  SLEQP_CALL(create_var_bounds(cauchy_data,
                               cauchy_data->iterate,
                               num_variables,
                               num_constraints));

  SLEQP_CALL(sleqp_lpi_set_bounds(cauchy_data->default_interface,
                                  cauchy_data->cons_lb,
                                  cauchy_data->cons_ub,
                                  cauchy_data->vars_lb,
                                  cauchy_data->vars_ub));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
cauchy_fetch_components(CauchyData* cauchy_data, Components components)
{
  SleqpLPi* solver = cauchy_data->default_interface;

  if (cauchy_data->use_reduced_interface)
  {
    solver = cauchy_data->reduced_interface;
  }

  if ((cauchy_data->dirty_components & BASE_STATS) && (components & BASE_STATS))
  {
    SLEQP_CALL(sleqp_lpi_vars_stats(solver, cauchy_data->var_stats));

    SLEQP_CALL(sleqp_lpi_cons_stats(solver, cauchy_data->cons_stats));

    cauchy_data->dirty_components &= ~(BASE_STATS);
  }

  if ((cauchy_data->dirty_components & PRIMAL_VALS)
      && (components & PRIMAL_VALS))
  {
    SLEQP_CALL(sleqp_lpi_primal_sol(solver, NULL, cauchy_data->primal_values));

    cauchy_data->dirty_components &= ~(PRIMAL_VALS);
  }

  if ((cauchy_data->dirty_components & DUAL_VALS) && (components & DUAL_VALS))
  {

    SLEQP_CALL(sleqp_lpi_dual_sol(solver,
                                  cauchy_data->vars_dual,
                                  cauchy_data->cons_dual));

    cauchy_data->dirty_components &= ~(DUAL_VALS);
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
switch_to_reduced_problem(CauchyData* cauchy_data)
{
  SleqpProblem* problem = cauchy_data->problem;
  SleqpIterate* iterate = cauchy_data->iterate;

  const int num_variables   = sleqp_problem_num_vars(problem);
  const int num_constraints = sleqp_problem_num_cons(problem);

  SleqpMat* cons_jac = sleqp_iterate_cons_jac(iterate);

  if (!cauchy_data->reduced_interface)
  {
    SLEQP_CALL(sleqp_lpi_create_default(&cauchy_data->reduced_interface,
                                        num_variables,
                                        num_constraints,
                                        cauchy_data->settings));
  }

  if (!cauchy_data->reduced_cons_stats)
  {
    SLEQP_CALL(
      sleqp_alloc_array(&cauchy_data->reduced_cons_stats, num_constraints));
  }

  SLEQP_BASESTAT* lower_slack_stats = cauchy_data->var_stats + num_variables;
  SLEQP_BASESTAT* upper_slack_stats = lower_slack_stats + num_constraints;

  for (int i = 0; i < num_constraints; ++i)
  {
    const SLEQP_BASESTAT cons_stat = cauchy_data->cons_stats[i];

    cauchy_data->reduced_cons_stats[i] = SLEQP_BASESTAT_BASIC;

    if (cons_stat == SLEQP_BASESTAT_BASIC)
    {
      continue;
    }

    bool zero_slack = lower_slack_stats[i] == SLEQP_BASESTAT_LOWER
                      && upper_slack_stats[i] == SLEQP_BASESTAT_LOWER;

    if (zero_slack)
    {
      cauchy_data->reduced_cons_stats[i] = cons_stat;
    }
  }

  SLEQP_CALL(cauchy_fetch_components(cauchy_data, BASE_STATS | PRIMAL_VALS));

  SLEQP_CALL(sleqp_lpi_set_objective(cauchy_data->reduced_interface,
                                     cauchy_data->objective));

  const double* lower_slack_values = cauchy_data->primal_values + num_variables;
  const double* upper_slack_values = lower_slack_values + num_constraints;

  for (int i = 0; i < num_constraints; ++i)
  {
    const double slack_diff = lower_slack_values[i] - upper_slack_values[i];

    cauchy_data->cons_lb[i] += slack_diff;
    cauchy_data->cons_ub[i] += slack_diff;
  }

  // adapt cons bounds
  SLEQP_CALL(sleqp_lpi_set_bounds(cauchy_data->reduced_interface,
                                  cauchy_data->cons_lb,
                                  cauchy_data->cons_ub,
                                  cauchy_data->vars_lb,
                                  cauchy_data->vars_ub));

  SLEQP_CALL(sleqp_lpi_set_coeffs(cauchy_data->reduced_interface, cons_jac));

  SLEQP_CALL(sleqp_lpi_set_time_limit(cauchy_data->reduced_interface,
                                      cauchy_data->time_limit));

  SLEQP_CALL(sleqp_lpi_set_basis(cauchy_data->reduced_interface,
                                 0,
                                 cauchy_data->var_stats,
                                 cauchy_data->reduced_cons_stats));

  SLEQP_CALL(sleqp_lpi_restore_basis(cauchy_data->reduced_interface, 0));

  SLEQP_CALL(sleqp_lpi_solve(cauchy_data->reduced_interface));

  SLEQP_LP_STATUS status = sleqp_lpi_status(cauchy_data->reduced_interface);

  if (status != SLEQP_LP_STATUS_OPTIMAL)
  {
    sleqp_raise(SLEQP_INTERNAL_ERROR, "Invalid LP status: %d", status);
  }

  cauchy_data->use_reduced_interface = true;
  cauchy_data->dirty_components      = ALL;

  return SLEQP_OKAY;
}

static SLEQP_ACTIVE_STATE
cauchy_cons_state(SLEQP_BASESTAT cons_stat,
                  SLEQP_BASESTAT lower_slack_stat,
                  SLEQP_BASESTAT upper_slack_stat,
                  bool equal_bounds)
{
  assert((lower_slack_stat == SLEQP_BASESTAT_LOWER)
         || (upper_slack_stat == SLEQP_BASESTAT_LOWER));

  if ((cons_stat == SLEQP_BASESTAT_BASIC) || (cons_stat == SLEQP_BASESTAT_ZERO))
  {
    return SLEQP_INACTIVE;
  }

  const bool zero_slacks = (lower_slack_stat == SLEQP_BASESTAT_LOWER)
                           && (upper_slack_stat == SLEQP_BASESTAT_LOWER);

  if (zero_slacks)
  {
    if (equal_bounds)
    {
      return SLEQP_ACTIVE_BOTH;
    }
    if (cons_stat == SLEQP_BASESTAT_LOWER)
    {
      return SLEQP_ACTIVE_LOWER;
    }
    if (cons_stat == SLEQP_BASESTAT_UPPER)
    {
      return SLEQP_ACTIVE_UPPER;
    }
  }

  return SLEQP_INACTIVE;
}

static SLEQP_RETCODE
needs_reduced_resolve(CauchyData* cauchy_data, bool* resolve)
{
  assert(!(cauchy_data->use_reduced_interface));

  SleqpProblem* problem = cauchy_data->problem;

  const int num_variables   = sleqp_problem_num_vars(problem);
  const int num_constraints = sleqp_problem_num_cons(problem);

  const double eps = sleqp_settings_real_value(cauchy_data->settings, SLEQP_SETTINGS_REAL_EPS);

  *resolve = false;

  SLEQP_CALL(
    cauchy_fetch_components(cauchy_data, BASE_STATS | PRIMAL_VALS | DUAL_VALS));

  SLEQP_BASESTAT* lower_slack_stats = cauchy_data->var_stats + num_variables;
  SLEQP_BASESTAT* upper_slack_stats = lower_slack_stats + num_constraints;

  const double* lower_slack_values = cauchy_data->primal_values + num_variables;
  const double* upper_slack_values = lower_slack_values + num_constraints;

  const double* cons_dual = cauchy_data->cons_dual;

  const SleqpVec* lb = sleqp_problem_cons_lb(problem);
  const SleqpVec* ub = sleqp_problem_cons_ub(problem);

  int k_lb = 0, k_ub = 0;

  bool needs_resolve = false;
  bool feasible      = true;

  for (int i = 0; i < num_constraints; ++i)
  {
    while (k_lb < lb->nnz && lb->indices[k_lb] < i)
    {
      ++k_lb;
    }

    while (k_ub < ub->nnz && ub->indices[k_ub] < i)
    {
      ++k_ub;
    }

    const bool valid_ub = (k_ub < ub->nnz) && (i == ub->indices[k_ub]);
    const bool valid_lb = (k_lb < lb->nnz) && (i == lb->indices[k_lb]);

    const double lb_val = valid_lb ? lb->data[k_lb] : 0.;
    const double ub_val = valid_ub ? ub->data[k_ub] : 0.;

    const bool equal_bounds = (sleqp_is_eq(lb_val, ub_val, eps));

    SLEQP_ACTIVE_STATE cons_stat = cauchy_cons_state(cauchy_data->cons_stats[i],
                                                     lower_slack_stats[i],
                                                     upper_slack_stats[i],
                                                     equal_bounds);

    if (cons_stat != SLEQP_INACTIVE)
    {
      continue;
    }

    const double lower_slack = lower_slack_values[i];
    const double upper_slack = upper_slack_values[i];

    const bool zero_slacks  = (lower_slack == 0) && (upper_slack == 0.);
    const bool nonzero_dual = (cons_dual[i] != 0.);

    if (!zero_slacks)
    {
      feasible = false;
      break;
    }

    // A constraint which is tight in the linearization with
    // a nonzero dual variable in the LP is *not* included
    // in the working set.
    if (zero_slacks && nonzero_dual)
    {
      needs_resolve = true;
    }
  }

  *resolve = (feasible && needs_resolve);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
check_basis(CauchyData* cauchy_data, bool* valid_basis)
{
  (*valid_basis) = true;

  int basis_size = 0;

  SLEQP_CALL(cauchy_fetch_components(cauchy_data, BASE_STATS));

  for (int j = 0; j < cauchy_data->num_lp_variables; ++j)
  {
    SLEQP_BASESTAT var_stat = cauchy_data->var_stats[j];

    if (var_stat != SLEQP_BASESTAT_BASIC)
    {
      ++basis_size;
    }
  }

  for (int i = 0; i < cauchy_data->num_lp_constraints; ++i)
  {
    SLEQP_BASESTAT cons_stat = cauchy_data->cons_stats[i];

    if (cons_stat != SLEQP_BASESTAT_BASIC)
    {
      ++basis_size;
    }
  }

  (*valid_basis) = (basis_size == cauchy_data->num_lp_variables);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
restore_basis(CauchyData* cauchy_data, SLEQP_CAUCHY_OBJTYPE objective_type)
{
  if (cauchy_data->current_objective != objective_type
      && cauchy_data->has_basis[objective_type])
  {
    SLEQP_CALL(
      sleqp_lpi_restore_basis(cauchy_data->default_interface, objective_type));
  }
  else if (cauchy_data->first_solve)
  {
    sleqp_log_debug("Using slack basis for initial solve");
    SLEQP_CALL(create_and_set_slack_basis(cauchy_data));
    cauchy_data->first_solve = false;
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
standard_cauchy_solve(SleqpVec* gradient,
                      double penalty,
                      SLEQP_CAUCHY_OBJTYPE objective_type,
                      void* data)
{
  CauchyData* cauchy_data = (CauchyData*)data;

  SLEQP_CALL(create_default_objective(cauchy_data, gradient, penalty));

  SLEQP_CALL(sleqp_lpi_set_objective(cauchy_data->default_interface,
                                     cauchy_data->objective));

  bool warm_start
    = sleqp_settings_bool_value(cauchy_data->settings,
                                SLEQP_SETTINGS_BOOL_ALWAYS_WARM_START_LP);

  if (warm_start)
  {
    switch (objective_type)
    {
    case SLEQP_CAUCHY_OBJTYPE_DEFAULT:
      // fallthrough
    case SLEQP_CAUCHY_OBJTYPE_FEAS:
      SLEQP_CALL(restore_basis(cauchy_data, objective_type));
      break;
    case SLEQP_CAUCHY_OBJTYPE_MIXED:
      if (cauchy_data->current_objective != SLEQP_CAUCHY_OBJTYPE_MIXED)
      {
        // restart from the default, this should be closer
        // to the initial mixed one
        SLEQP_CALL(restore_basis(cauchy_data, SLEQP_CAUCHY_OBJTYPE_DEFAULT));
      }
      break;
    default:
      break;
    }
  }

  cauchy_data->current_objective = objective_type;

  {
    SLEQP_CALL(sleqp_lpi_set_time_limit(cauchy_data->default_interface,
                                        cauchy_data->time_limit));

    SLEQP_CALL(sleqp_lpi_solve(cauchy_data->default_interface));

    SleqpTimer* lpi_timer
      = sleqp_lpi_solve_timer(cauchy_data->default_interface);

    const double elapsed_time = sleqp_timer_elapsed(lpi_timer);

    cauchy_data->time_limit
      = sleqp_remaining_time(elapsed_time, cauchy_data->time_limit);
  }

  cauchy_data->dirty_components      = ALL;
  cauchy_data->use_reduced_interface = false;

  SLEQP_LP_STATUS status = sleqp_lpi_status(cauchy_data->default_interface);

  if (status != SLEQP_LP_STATUS_OPTIMAL)
  {
    sleqp_raise(SLEQP_INTERNAL_ERROR, "Invalid LP status: %d", status);
  }

#if SLEQP_DEBUG
  {
    bool valid_basis = false;

    SLEQP_CALL(check_basis(cauchy_data, &valid_basis));

    assert(valid_basis);
  }

#endif

  if (warm_start)
  {
    SLEQP_CALL(
      sleqp_lpi_save_basis(cauchy_data->default_interface, objective_type));
  }

  cauchy_data->has_basis[objective_type] = true;

  const bool enable_resolves
    = sleqp_settings_bool_value(cauchy_data->settings,
                                SLEQP_SETTINGS_BOOL_LP_RESOLVES);

  if (enable_resolves && (objective_type == SLEQP_CAUCHY_OBJTYPE_DEFAULT))
  {
    bool resolve;

    SLEQP_CALL(needs_reduced_resolve(cauchy_data, &resolve));

    if (resolve)
    {
      sleqp_log_debug("Resolving reduced problem");
      SLEQP_CALL(switch_to_reduced_problem(cauchy_data));
    }
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
standard_cauchy_obj_val(double* objective_value, void* data)
{
  CauchyData* cauchy_data = (CauchyData*)data;

  SLEQP_CALL(sleqp_lpi_primal_sol(cauchy_data->default_interface,
                                  objective_value,
                                  NULL));

  (*objective_value) += sleqp_iterate_obj_val(cauchy_data->iterate);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
cauchy_working_set_add_vars(CauchyData* cauchy_data, SleqpIterate* iterate)
{
  SleqpProblem* problem = cauchy_data->problem;
  SleqpVec* x           = sleqp_iterate_primal(iterate);
  const SleqpVec* lb    = sleqp_problem_vars_lb(problem);
  const SleqpVec* ub    = sleqp_problem_vars_ub(problem);

  const int num_variables = sleqp_problem_num_vars(problem);

  SleqpWorkingSet* working_set = sleqp_iterate_working_set(iterate);

  const double trust_radius = cauchy_data->trust_radius;

  SLEQP_CALL(cauchy_fetch_components(cauchy_data, BASE_STATS));

  const double eps = sleqp_settings_real_value(cauchy_data->settings, SLEQP_SETTINGS_REAL_EPS);

  int k_x = 0, k_lb = 0, k_ub = 0;

  for (int i = 0; i < num_variables; ++i)
  {
    while (k_x < x->nnz && x->indices[k_x] < i)
    {
      ++k_x;
    }

    while (k_lb < lb->nnz && lb->indices[k_lb] < i)
    {
      ++k_lb;
    }

    while (k_ub < ub->nnz && ub->indices[k_ub] < i)
    {
      ++k_ub;
    }

    const bool valid_x  = (k_x < x->nnz) && (i == x->indices[k_x]);
    const bool valid_ub = (k_ub < ub->nnz) && (i == ub->indices[k_ub]);
    const bool valid_lb = (k_lb < lb->nnz) && (i == lb->indices[k_lb]);

    const double ubval = valid_ub ? ub->data[k_ub] : 0.;
    const double lbval = valid_lb ? lb->data[k_lb] : 0.;
    const double xval  = valid_x ? x->data[k_x] : 0.;

    assert(cauchy_data->var_stats[i] != SLEQP_BASESTAT_ZERO);

    sleqp_assert_is_leq(lbval, xval, eps);
    sleqp_assert_is_leq(xval, ubval, eps);

    if (sleqp_is_eq(lbval, ubval, eps))
    {
      SLEQP_CALL(sleqp_working_set_add_var(working_set, i, SLEQP_ACTIVE_BOTH));
    }
    else if ((cauchy_data->var_stats[i] == SLEQP_BASESTAT_LOWER)
             && ((xval - lbval) < trust_radius))
    {
      SLEQP_CALL(sleqp_working_set_add_var(working_set, i, SLEQP_ACTIVE_LOWER));
    }
    else if ((cauchy_data->var_stats[i] == SLEQP_BASESTAT_UPPER)
             && ((ubval - xval) < trust_radius))
    {
      SLEQP_CALL(sleqp_working_set_add_var(working_set, i, SLEQP_ACTIVE_UPPER));
    }
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
cauchy_working_set_add_cons(CauchyData* cauchy_data, SleqpIterate* iterate)
{
  SleqpProblem* problem = cauchy_data->problem;

  const int num_variables   = sleqp_problem_num_vars(problem);
  const int num_constraints = sleqp_problem_num_cons(problem);

  SleqpWorkingSet* working_set = sleqp_iterate_working_set(iterate);

  const double eps = sleqp_settings_real_value(cauchy_data->settings, SLEQP_SETTINGS_REAL_EPS);

  int k_lb = 0, k_ub = 0;

  SLEQP_CALL(cauchy_fetch_components(cauchy_data, BASE_STATS));

  SLEQP_BASESTAT* lower_slack_stats = cauchy_data->var_stats + num_variables;
  SLEQP_BASESTAT* upper_slack_stats = lower_slack_stats + num_constraints;

  const SleqpVec* lb = sleqp_problem_cons_lb(problem);
  const SleqpVec* ub = sleqp_problem_cons_ub(problem);

  for (int i = 0; i < num_constraints; ++i)
  {
    while (k_lb < lb->nnz && lb->indices[k_lb] < i)
    {
      ++k_lb;
    }

    while (k_ub < ub->nnz && ub->indices[k_ub] < i)
    {
      ++k_ub;
    }

    const bool valid_lb = (k_lb < lb->nnz) && (i == lb->indices[k_lb]);
    const bool valid_ub = (k_ub < ub->nnz) && (i == ub->indices[k_ub]);

    const double lb_val = valid_lb ? lb->data[k_lb] : 0.;
    const double ub_val = valid_ub ? ub->data[k_ub] : 0.;

    const bool equal_bounds = (sleqp_is_eq(lb_val, ub_val, eps));

    SLEQP_ACTIVE_STATE cons_stat = cauchy_cons_state(cauchy_data->cons_stats[i],
                                                     lower_slack_stats[i],
                                                     upper_slack_stats[i],
                                                     equal_bounds);

    if (cons_stat != SLEQP_INACTIVE)
    {
      SLEQP_CALL(sleqp_working_set_add_cons(working_set, i, cons_stat));
    }
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
get_reduced_working_set(CauchyData* cauchy_data, SleqpIterate* iterate)
{
  SLEQP_CALL(cauchy_fetch_components(cauchy_data, BASE_STATS));

  SleqpProblem* problem     = cauchy_data->problem;
  const int num_variables   = sleqp_problem_num_vars(problem);
  const int num_constraints = sleqp_problem_num_cons(problem);

  SLEQP_CALL(cauchy_working_set_add_vars(cauchy_data, iterate));

  SleqpWorkingSet* working_set = sleqp_iterate_working_set(iterate);

  const double* lower_slack_values = cauchy_data->primal_values + num_variables;
  const double* upper_slack_values = lower_slack_values + num_constraints;

  for (int i = 0; i < num_constraints; ++i)
  {
    const SLEQP_BASESTAT cons_stat = cauchy_data->cons_stats[i];

    if ((cons_stat == SLEQP_BASESTAT_BASIC)
        || (cons_stat == SLEQP_BASESTAT_ZERO))
    {
      continue;
    }

    if ((lower_slack_values[i] == 0) && (upper_slack_values[i] == 0.))
    {
      if (cons_stat == SLEQP_BASESTAT_LOWER)
      {
        SLEQP_CALL(
          sleqp_working_set_add_cons(working_set, i, SLEQP_ACTIVE_LOWER));
      }
      else if (cons_stat == SLEQP_BASESTAT_UPPER)
      {
        SLEQP_CALL(
          sleqp_working_set_add_cons(working_set, i, SLEQP_ACTIVE_UPPER));
      }
    }
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
get_simple_working_set(CauchyData* cauchy_data, SleqpIterate* iterate)
{
  SLEQP_CALL(cauchy_working_set_add_vars(cauchy_data, iterate));

  SLEQP_CALL(cauchy_working_set_add_cons(cauchy_data, iterate));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
standard_cauchy_working_set(SleqpIterate* iterate, void* data)
{
  CauchyData* cauchy_data = (CauchyData*)data;

  SleqpWorkingSet* working_set = sleqp_iterate_working_set(iterate);

  SLEQP_CALL(sleqp_working_set_reset(working_set));

  if (cauchy_data->use_reduced_interface)
  {
    SLEQP_CALL(get_reduced_working_set(cauchy_data, iterate));
  }
  else
  {
    SLEQP_CALL(get_simple_working_set(cauchy_data, iterate));
  }

  const int num_active_vars = sleqp_working_set_num_active_vars(working_set);
  const int num_active_cons = sleqp_working_set_num_active_cons(working_set);

  sleqp_log_debug("Created an active set with %d variables, %d constraints",
                  num_active_vars,
                  num_active_cons);

  assert(sleqp_working_set_valid(working_set));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
standard_cauchy_lp_step(SleqpVec* direction, void* data)
{
  CauchyData* cauchy_data = (CauchyData*)data;

  SLEQP_CALL(cauchy_fetch_components(cauchy_data, PRIMAL_VALS));

  const double zero_eps
    = sleqp_settings_real_value(cauchy_data->settings, SLEQP_SETTINGS_REAL_ZERO_EPS);

  const int num_variables = sleqp_problem_num_vars(cauchy_data->problem);

  SLEQP_CALL(sleqp_vec_set_from_raw(direction,
                                    cauchy_data->primal_values,
                                    num_variables,
                                    zero_eps));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
standard_cauchy_locally_infeasible(bool* locally_infeasible, void* data)
{
  CauchyData* cauchy_data = (CauchyData*)data;

  *locally_infeasible = false;

  SleqpIterate* iterate = cauchy_data->iterate;

  assert(iterate);
  assert(!cauchy_data->use_reduced_interface);

  SLEQP_CALL(cauchy_fetch_components(cauchy_data, BASE_STATS));

  const double trust_radius = cauchy_data->trust_radius;

  assert(trust_radius > 0.);

  SleqpProblem* problem = cauchy_data->problem;

  const int num_variables   = sleqp_problem_num_vars(problem);
  const int num_constraints = sleqp_problem_num_cons(problem);

  const double eps = sleqp_settings_real_value(cauchy_data->settings, SLEQP_SETTINGS_REAL_EPS);

  bool active_trust_region = false;

  // Check if trust region is active
  {
    SleqpVec* x        = sleqp_iterate_primal(iterate);
    const SleqpVec* lb = sleqp_problem_vars_lb(problem);
    const SleqpVec* ub = sleqp_problem_vars_ub(problem);

    int k_x = 0, k_lb = 0, k_ub = 0;

    for (int i = 0; i < num_variables; ++i)
    {
      while (k_x < x->nnz && x->indices[k_x] < i)
      {
        ++k_x;
      }

      while (k_lb < lb->nnz && lb->indices[k_lb] < i)
      {
        ++k_lb;
      }

      while (k_ub < ub->nnz && ub->indices[k_ub] < i)
      {
        ++k_ub;
      }

      const bool valid_x  = (k_x < x->nnz) && (i == x->indices[k_x]);
      const bool valid_ub = (k_ub < ub->nnz) && (i == ub->indices[k_ub]);
      const bool valid_lb = (k_lb < lb->nnz) && (i == lb->indices[k_lb]);

      const double ubval = valid_ub ? ub->data[k_ub] : 0.;
      const double lbval = valid_lb ? lb->data[k_lb] : 0.;
      const double xval  = valid_x ? x->data[k_x] : 0.;

      assert(cauchy_data->var_stats[i] != SLEQP_BASESTAT_ZERO);

      sleqp_assert_is_leq(lbval, xval, eps);
      sleqp_assert_is_leq(xval, ubval, eps);

      if (sleqp_is_eq(lbval, ubval, eps))
      {
        continue;
      }
      else if (cauchy_data->var_stats[i] == SLEQP_BASESTAT_LOWER
               && (xval - lbval) >= trust_radius)
      {
        active_trust_region = true;
        break;
      }
      else if (cauchy_data->var_stats[i] == SLEQP_BASESTAT_UPPER
               && ((ubval - xval) >= trust_radius))
      {
        active_trust_region = true;
        break;
      }
    }
  }

  bool feasible_direction = true;

  {
    int k_lb = 0, k_ub = 0;

    SLEQP_BASESTAT* lower_slack_stats = cauchy_data->var_stats + num_variables;
    SLEQP_BASESTAT* upper_slack_stats = lower_slack_stats + num_constraints;

    const SleqpVec* lb = sleqp_problem_cons_lb(problem);
    const SleqpVec* ub = sleqp_problem_cons_ub(problem);

    for (int i = 0; i < num_constraints; ++i)
    {
      const SLEQP_BASESTAT cons_stat = cauchy_data->cons_stats[i];

      if (cons_stat == SLEQP_BASESTAT_BASIC)
      {
        continue;
      }

      while (k_lb < lb->nnz && lb->indices[k_lb] < i)
      {
        ++k_lb;
      }

      while (k_ub < ub->nnz && ub->indices[k_ub] < i)
      {
        ++k_ub;
      }

      assert(lower_slack_stats[i] != SLEQP_BASESTAT_BASIC
             || upper_slack_stats[i] != SLEQP_BASESTAT_BASIC);

      bool zero_slack = lower_slack_stats[i] == SLEQP_BASESTAT_LOWER
                        && upper_slack_stats[i] == SLEQP_BASESTAT_LOWER;

      if (!zero_slack)
      {
        feasible_direction = false;
      }
    }
  }

  sleqp_log_debug("Trust region active: %s, feasible direction: %s",
                  sleqp_bool_string(active_trust_region),
                  sleqp_bool_string(feasible_direction));

  *locally_infeasible = !(feasible_direction || active_trust_region);

  return SLEQP_OKAY;
}

/*
 * Trims off dual values by
 * removing entries in the sparse vector
 * which are inactive or have the wrong sign
 */
static SLEQP_RETCODE
trim_duals_to_working_set(const SLEQP_ACTIVE_STATE* states,
                          SleqpVec* duals,
                          double zero_eps)
{
  int offset = 0;

  for (int k = 0; k < duals->nnz; ++k)
  {
    const int i = duals->indices[k];

    const SLEQP_ACTIVE_STATE state = states[i];

    if (state == SLEQP_INACTIVE)
    {
      ++offset;
    }
    else
    {
      duals->indices[k - offset] = duals->indices[k];
      duals->data[k - offset]    = duals->data[k];

      if (state == SLEQP_ACTIVE_UPPER)
      {
        // Should be zero up to dual feasibility tolerance
        sleqp_assert_is_geq(duals->data[k - offset], 0., zero_eps);

        if (duals->data[k - offset] < 0)
        {
          ++offset;
        }
      }
      else if (state == SLEQP_ACTIVE_LOWER)
      {
        // Should be zero up to dual feasibility tolerance
        sleqp_assert_is_leq(duals->data[k - offset], 0., zero_eps);

        if (duals->data[k - offset] > 0)
        {
          ++offset;
        }
      }
    }
  }

  duals->nnz -= offset;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
standard_cauchy_estimate_duals(const SleqpWorkingSet* working_set,
                               SleqpVec* cons_dual,
                               SleqpVec* vars_dual,
                               void* data)
{
  CauchyData* cauchy_data = (CauchyData*)data;

  const double zero_eps
    = sleqp_settings_real_value(cauchy_data->settings, SLEQP_SETTINGS_REAL_ZERO_EPS);

  SleqpProblem* problem = cauchy_data->problem;

  const int num_variables   = sleqp_problem_num_vars(problem);
  const int num_constraints = sleqp_problem_num_cons(problem);

  if (vars_dual)
  {
    assert(vars_dual->dim == num_variables);

    SLEQP_CALL(sleqp_lpi_dual_sol(cauchy_data->default_interface,
                                  cauchy_data->vars_dual,
                                  NULL));

    SLEQP_CALL(sleqp_vec_set_from_raw(vars_dual,
                                      cauchy_data->vars_dual,
                                      vars_dual->dim,
                                      zero_eps));

    // Note: We rescale here since sign conventions vary...
    SLEQP_CALL(sleqp_vec_scale(vars_dual, -1.));

    SLEQP_CALL(
      trim_duals_to_working_set(sleqp_working_set_var_states(working_set),
                                vars_dual,
                                zero_eps));
  }

  if (cons_dual)
  {
    assert(cons_dual->dim == num_constraints);

    SLEQP_CALL(sleqp_lpi_dual_sol(cauchy_data->default_interface,
                                  NULL,
                                  cauchy_data->vars_dual));

    SLEQP_CALL(sleqp_vec_set_from_raw(cons_dual,
                                      cauchy_data->vars_dual,
                                      cons_dual->dim,
                                      zero_eps));

    // Note: We rescale here since sign conventions vary...
    SLEQP_CALL(sleqp_vec_scale(cons_dual, -1.));

    SLEQP_CALL(
      trim_duals_to_working_set(sleqp_working_set_cons_states(working_set),
                                cons_dual,
                                zero_eps));
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
standard_cauchy_violation(double* violation, void* data)
{
  CauchyData* cauchy_data = (CauchyData*)data;

  SleqpProblem* problem = cauchy_data->problem;

  const int num_variables = sleqp_problem_num_vars(problem);

  SLEQP_CALL(sleqp_lpi_primal_sol(cauchy_data->default_interface,
                                  NULL,
                                  cauchy_data->primal_values));

  (*violation) = 0.;

  for (int j = num_variables; j < cauchy_data->num_lp_variables; ++j)
  {
    (*violation) += cauchy_data->primal_values[j];
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
standard_cauchy_set_time_limit(double time_limit, void* data)
{
  CauchyData* cauchy_data = (CauchyData*)data;

  cauchy_data->time_limit = time_limit;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
standard_cauchy_basis_condition(bool* exact, double* condition, void* data)
{
  CauchyData* cauchy_data = (CauchyData*)data;

  SLEQP_CALL(
    sleqp_lpi_basis_cond(cauchy_data->default_interface, exact, condition));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
standard_cauchy_print_stats(double total_elapsed, void* data)
{
  CauchyData* cauchy_data = (CauchyData*)data;

  SleqpTimer* default_timer
    = sleqp_lpi_solve_timer(cauchy_data->default_interface);

  SLEQP_CALL(sleqp_timer_display(default_timer, "Solved LPs", total_elapsed));

  if (cauchy_data->reduced_interface)
  {
    SleqpTimer* reduced_timer
      = sleqp_lpi_solve_timer(cauchy_data->reduced_interface);

    SLEQP_CALL(
      sleqp_timer_display(reduced_timer, "Solved reduced LPs", total_elapsed));
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
standard_cauchy_free(void* star)
{
  CauchyData* cauchy_data = (CauchyData*)star;

  sleqp_free(&cauchy_data->vars_dual);
  sleqp_free(&cauchy_data->primal_values);

  sleqp_free(&cauchy_data->vars_ub);
  sleqp_free(&cauchy_data->vars_lb);

  sleqp_free(&cauchy_data->cons_ub);
  sleqp_free(&cauchy_data->cons_lb);

  sleqp_free(&cauchy_data->objective);

  sleqp_free(&cauchy_data->reduced_cons_stats);

  SLEQP_CALL(sleqp_lpi_release(&cauchy_data->reduced_interface));
  SLEQP_CALL(sleqp_lpi_release(&cauchy_data->default_interface));

  sleqp_free(&cauchy_data->cons_stats);
  sleqp_free(&cauchy_data->var_stats);

  SLEQP_CALL(sleqp_iterate_release(&cauchy_data->iterate));

  SLEQP_CALL(sleqp_settings_release(&cauchy_data->settings));

  SLEQP_CALL(sleqp_problem_release(&cauchy_data->problem));

  sleqp_free(&star);

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_standard_cauchy_create(SleqpCauchy** star,
                             SleqpProblem* problem,
                             SleqpSettings* settings)
{
  CauchyData* cauchy_data;

  SLEQP_CALL(cauchy_data_create(&cauchy_data, problem, settings));

  SleqpCauchyCallbacks callbacks
    = {.set_iterate        = standard_cauchy_set_iterate,
       .set_trust_radius   = standard_cauchy_set_trust_radius,
       .solve              = standard_cauchy_solve,
       .obj_val            = standard_cauchy_obj_val,
       .working_set        = standard_cauchy_working_set,
       .lp_step            = standard_cauchy_lp_step,
       .locally_infeasible = standard_cauchy_locally_infeasible,
       .estimate_duals     = standard_cauchy_estimate_duals,
       .violation          = standard_cauchy_violation,
       .set_time_limit     = standard_cauchy_set_time_limit,
       .basis_condition    = standard_cauchy_basis_condition,
       .print_stats        = standard_cauchy_print_stats,
       .free               = standard_cauchy_free};

  SLEQP_CALL(sleqp_cauchy_create(star, &callbacks, (void*)cauchy_data));

  return SLEQP_OKAY;
}
