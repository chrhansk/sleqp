#include "sleqp_cauchy.h"

#include <math.h>

#include "sleqp_assert.h"
#include "sleqp_cmp.h"
#include "sleqp_mem.h"
#include "sleqp_merit.h"

struct SleqpCauchyData
{
  SleqpProblem* problem;
  SleqpParams* params;

  int num_lp_variables;
  int num_lp_constraints;

  double trust_radius;
  SleqpIterate* iterate;

  SLEQP_BASESTAT* var_stats;
  SLEQP_BASESTAT* cons_stats;

  SleqpLPi* lp_interface;

  double* objective;
  double* cons_lb;
  double* cons_ub;
  double* vars_lb;
  double* vars_ub;

  double* solution_values;
  double* dual_values;

  SleqpSparseVec* quadratic_gradient;
};

SLEQP_RETCODE sleqp_cauchy_data_create(SleqpCauchyData** star,
                                       SleqpProblem* problem,
                                       SleqpParams* params,
                                       SleqpLPi* lp_interface)
{
  SLEQP_CALL(sleqp_malloc(star));

  SleqpCauchyData* data = *star;

  *data = (SleqpCauchyData){0};

  data->problem = problem;
  data->params = params;

  data->num_lp_variables = problem->num_variables + 2 * problem->num_constraints;
  data->num_lp_constraints = problem->num_constraints;

  data->trust_radius = SLEQP_NONE;

  assert(data->num_lp_variables == sleqp_lpi_get_num_variables(lp_interface));
  assert(data->num_lp_constraints == sleqp_lpi_get_num_constraints(lp_interface));

  SLEQP_CALL(sleqp_calloc(&data->var_stats, data->num_lp_variables));
  SLEQP_CALL(sleqp_calloc(&data->cons_stats, data->num_lp_constraints));

  data->lp_interface = lp_interface;

  SLEQP_CALL(sleqp_calloc(&data->objective, data->num_lp_variables));

  SLEQP_CALL(sleqp_calloc(&data->cons_lb, data->num_lp_constraints));
  SLEQP_CALL(sleqp_calloc(&data->cons_ub, data->num_lp_constraints));

  SLEQP_CALL(sleqp_calloc(&data->vars_lb, data->num_lp_variables));
  SLEQP_CALL(sleqp_calloc(&data->vars_ub, data->num_lp_variables));

  SLEQP_CALL(sleqp_calloc(&data->solution_values, data->num_lp_variables));

  SLEQP_CALL(sleqp_calloc(&data->dual_values,
                          SLEQP_MAX(data->num_lp_constraints,
                                    data->num_lp_variables)));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&data->quadratic_gradient,
                                              problem->num_variables));

  const double inf = sleqp_infinity();

  for(int j = problem->num_variables; j < data->num_lp_variables; ++j)
  {
    data->vars_lb[j] = 0;
    data->vars_ub[j] = inf;
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE append_identities(SleqpSparseMatrix* cons_jac,
                                       int num_variables,
                                       int num_constraints)
{
  assert(num_constraints == sleqp_sparse_matrix_get_num_rows(cons_jac));
  assert(num_variables == sleqp_sparse_matrix_get_num_cols(cons_jac));

  const int nnz = sleqp_sparse_matrix_get_nnz(cons_jac);
  const int num_rows = sleqp_sparse_matrix_get_num_rows(cons_jac);
  const int num_cols = sleqp_sparse_matrix_get_num_cols(cons_jac);

  /*
   * Reserve a litte more so we can add the two
   * identity matrices afterwards
   */
  SLEQP_CALL(sleqp_sparse_matrix_reserve(cons_jac,
                                         nnz + 2*num_constraints));

  SLEQP_CALL(sleqp_sparse_matrix_resize(cons_jac,
                                        num_rows,
                                        num_cols + 2*num_constraints));

  // append a + id
  for(int i = 0; i < num_constraints; ++i)
  {
    const int con = num_variables + i;

    SLEQP_CALL(sleqp_sparse_matrix_push_column(cons_jac,
                                               con));

    SLEQP_CALL(sleqp_sparse_matrix_push(cons_jac,
                                        i,
                                        con,
                                        1.));
  }

  // append a - id
  for(int i = 0; i < num_constraints; ++i)
  {
    const int con = num_variables + num_constraints + i;

    SLEQP_CALL(sleqp_sparse_matrix_push_column(cons_jac,
                                               con));

    SLEQP_CALL(sleqp_sparse_matrix_push(cons_jac,
                                        i,
                                        con,
                                        -1.));
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE remove_identities(SleqpSparseMatrix* cons_jac,
                                       int num_variables,
                                       int num_constraints)
{
  const int num_rows = sleqp_sparse_matrix_get_num_rows(cons_jac);
  const int num_cols = sleqp_sparse_matrix_get_num_cols(cons_jac);

  assert(num_constraints == num_rows);
  assert(num_variables + 2*num_constraints == num_cols);

  SLEQP_CALL(sleqp_sparse_matrix_resize(cons_jac,
                                        num_constraints,
                                        num_variables));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE create_cons_bounds(SleqpCauchyData* cauchy_data,
                                        SleqpIterate* iterate,
                                        int num_variables,
                                        int num_constraints)
{
  int k_c = 0, k_lb = 0, k_ub = 0;

  SleqpSparseVec* lb = cauchy_data->problem->cons_lb;
  SleqpSparseVec* ub = cauchy_data->problem->cons_ub;

  SleqpSparseVec* val = sleqp_iterate_get_cons_val(iterate);

  const double eps = sleqp_params_get_eps(cauchy_data->params);
  const double inf = sleqp_infinity();

  for(int i = 0; i < num_constraints; ++i)
  {
    while(k_c < val->nnz && val->indices[k_c] < i)
    {
      ++k_c;
    }

    while(k_lb < lb->nnz && lb->indices[k_lb] < i)
    {
      ++k_lb;
    }

    while(k_ub < ub->nnz && ub->indices[k_ub] < i)
    {
      ++k_ub;
    }

    const bool valid_ub = (k_ub < ub->nnz && ub->indices[k_ub] == i);
    const bool valid_lb = (k_lb < lb->nnz && lb->indices[k_lb] == i);
    const bool valid_c = (k_c < val->nnz && val->indices[k_c] == i);

    const double ubval = valid_ub ? ub->data[k_ub] : 0;
    const double lbval = valid_lb ? lb->data[k_lb] : 0;
    const double cval = valid_c ? val->data[k_c] : 0;

    assert(!sleqp_is_inf(lbval));
    assert(!sleqp_is_inf(-ubval));
    assert(!sleqp_is_inf(SLEQP_ABS(cval)));

    if(sleqp_is_inf(ubval))
    {
      cauchy_data->cons_ub[i] = inf;
    }
    else
    {
      cauchy_data->cons_ub[i] = ubval - cval;
    }

    if(sleqp_is_inf(-lbval))
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

static SLEQP_RETCODE create_var_bounds(SleqpCauchyData* cauchy_data,
                                       SleqpIterate* iterate,
                                       int num_variables,
                                       int num_constraints)
{
  SleqpSparseVec* x = sleqp_iterate_get_primal(iterate);
  SleqpSparseVec* lb = cauchy_data->problem->var_lb;
  SleqpSparseVec* ub = cauchy_data->problem->var_ub;

  const double trust_radius = cauchy_data->trust_radius;

  assert(trust_radius > 0.);

  const double eps = sleqp_params_get_eps(cauchy_data->params);

  int k_x = 0, k_lb = 0, k_ub = 0;

  for(int j = 0; j < num_variables; ++j)
  {
    while(k_x < x->nnz && x->indices[k_x] < j)
    {
      ++k_x;
    }

    while(k_lb < lb->nnz && lb->indices[k_lb] < j)
    {
      ++k_lb;
    }

    while(k_ub < ub->nnz && ub->indices[k_ub] < j)
    {
      ++k_ub;
    }

    const bool valid_x = (k_x < x->nnz) && (j == x->indices[k_x]);
    const bool valid_ub = (k_ub < ub->nnz) && (j == ub->indices[k_ub]);
    const bool valid_lb = (k_lb < lb->nnz) && (j == lb->indices[k_lb]);

    const double ubval = valid_ub ? ub->data[k_ub] : 0.;
    const double lbval = valid_lb ? lb->data[k_lb] : 0.;
    const double xval = valid_x ? x->data[k_x] : 0.;

    assert(!sleqp_is_inf(lbval));
    assert(!sleqp_is_inf(-ubval));
    assert(!sleqp_is_inf(SLEQP_ABS(xval)));

    if(sleqp_is_inf(ubval))
    {
      cauchy_data->vars_ub[j] = trust_radius;
    }
    else
    {
      cauchy_data->vars_ub[j] = SLEQP_MIN(ubval - xval, trust_radius);
    }

    if(sleqp_is_inf(-lbval))
    {
      cauchy_data->vars_lb[j] = -trust_radius;
    }
    else
    {
      cauchy_data->vars_lb[j] = SLEQP_MAX(lbval - xval, -trust_radius);
    }

    assert(cauchy_data->vars_lb[j] <=cauchy_data->vars_ub[j]);
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE create_objective(SleqpCauchyData* cauchy_data,
                                      SleqpSparseVec* gradient,
                                      double penalty)
{
  SleqpProblem* problem = cauchy_data->problem;

  const int num_variables = problem->num_variables;
  const int num_constraints = problem->num_constraints;

  const int num_lp_variables = sleqp_lpi_get_num_variables(cauchy_data->lp_interface);

  assert(num_lp_variables == num_variables + 2*num_constraints);

  for(int j = 0; j < num_variables; ++j)
  {
    cauchy_data->objective[j] = 0.;
  }

  for(int j = num_variables; j < num_lp_variables; ++j)
  {
    cauchy_data->objective[j] = penalty;
  }

  if(gradient)
  {
    assert(gradient->dim == num_variables);

    for(int k = 0; k < gradient->nnz; ++k)
    {
      cauchy_data->objective[gradient->indices[k]] = gradient->data[k];
    }
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_cauchy_set_iterate(SleqpCauchyData* cauchy_data,
                                       SleqpIterate* iterate,
                                       double trust_radius)
{
  SleqpSparseMatrix* cons_jac = sleqp_iterate_get_cons_jac(iterate);

  assert(trust_radius > 0.);

  cauchy_data->trust_radius = trust_radius;

  SLEQP_CALL(sleqp_iterate_release(&cauchy_data->iterate));

  cauchy_data->iterate = iterate;

  SLEQP_CALL(sleqp_iterate_capture(cauchy_data->iterate));

  const int num_variables = sleqp_sparse_matrix_get_num_cols(cons_jac);
  const int num_constraints = sleqp_sparse_matrix_get_num_rows(cons_jac);

  SLEQP_CALL(append_identities(cons_jac,
                               num_variables,
                               num_constraints));

  assert(sleqp_sparse_matrix_valid(cons_jac));

  SLEQP_CALL(create_var_bounds(cauchy_data,
                               iterate,
                               num_variables,
                               num_constraints));

  SLEQP_CALL(create_cons_bounds(cauchy_data,
                                iterate,
                                num_variables,
                                num_constraints));

  SLEQP_CALL(sleqp_lpi_set_bounds(cauchy_data->lp_interface,
                                  cauchy_data->cons_lb,
                                  cauchy_data->cons_ub,
                                  cauchy_data->vars_lb,
                                  cauchy_data->vars_ub));

  SLEQP_CALL(sleqp_lpi_set_coefficients(cauchy_data->lp_interface,
                                        cons_jac));

  SLEQP_CALL(remove_identities(cons_jac,
                               num_variables,
                               num_constraints));

  assert(sleqp_sparse_matrix_valid(cons_jac));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE check_basis(SleqpCauchyData* cauchy_data,
                                 bool* valid_basis)
{
  (*valid_basis) = true;

  int basis_size = 0;

  SLEQP_CALL(sleqp_lpi_get_varstats(cauchy_data->lp_interface,
                                    cauchy_data->var_stats));

  for(int j = 0; j < cauchy_data->num_lp_variables; ++j)
  {
    SLEQP_BASESTAT var_stat = cauchy_data->var_stats[j];

    if(var_stat != SLEQP_BASESTAT_BASIC)
    {
      ++basis_size;
    }
  }

  SLEQP_CALL(sleqp_lpi_get_consstats(cauchy_data->lp_interface,
                                     cauchy_data->cons_stats));

  for(int i = 0; i < cauchy_data->num_lp_constraints; ++i)
  {
    SLEQP_BASESTAT cons_stat = cauchy_data->cons_stats[i];

    if(cons_stat != SLEQP_BASESTAT_BASIC)
    {
      ++basis_size;
    }
  }

  (*valid_basis) = (basis_size == cauchy_data->num_lp_variables);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE check_direction_bounds(SleqpCauchyData* cauchy_data,
                                            bool* valid_direction)
{
  SleqpProblem* problem = cauchy_data->problem;

  SleqpIterate* iterate = cauchy_data->iterate;

  assert(iterate);

  const double trust_radius = cauchy_data->trust_radius;

  assert(trust_radius > 0.);

  const double eps = sleqp_params_get_eps(cauchy_data->params);

  const int num_variables = problem->num_variables;

  SLEQP_CALL(sleqp_lpi_get_primal_sol(cauchy_data->lp_interface,
                                      NULL,
                                      cauchy_data->solution_values));

  (*valid_direction) = true;

  for(int j = 0; j < num_variables; ++j)
  {
    const double dval = cauchy_data->solution_values[j];

    if(sleqp_is_gt(dval, trust_radius, eps) || sleqp_is_lt(dval, -trust_radius, eps))
    {
      (*valid_direction) = false;
    }
  }

  {
    SleqpSparseVec* x = sleqp_iterate_get_primal(iterate);
    SleqpSparseVec* lb = problem->var_lb;
    SleqpSparseVec* ub = problem->var_ub;

    int k_x = 0, k_lb = 0, k_ub = 0;

    for(int j = 0; j < num_variables; ++j)
    {
      while(k_x < x->nnz && x->indices[k_x] < j)
      {
        ++k_x;
      }

      while(k_lb < lb->nnz && lb->indices[k_lb] < j)
      {
        ++k_lb;
      }

      while(k_ub < ub->nnz && ub->indices[k_ub] < j)
      {
        ++k_ub;
      }

      const bool valid_x = (k_x < x->nnz) && (j == x->indices[k_x]);
      const bool valid_ub = (k_ub < ub->nnz) && (j == ub->indices[k_ub]);
      const bool valid_lb = (k_lb < lb->nnz) && (j == lb->indices[k_lb]);

      const double ubval = valid_ub ? ub->data[k_ub] : 0.;
      const double lbval = valid_lb ? lb->data[k_lb] : 0.;
      const double xval = valid_x ? x->data[k_x] : 0.;
      const double dval = cauchy_data->solution_values[j];

      if(sleqp_is_gt(xval + dval, ubval, eps) || sleqp_is_lt(xval + dval, lbval, eps))
      {
        (*valid_direction) = false;
      }
    }
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_cauchy_solve(SleqpCauchyData* cauchy_data,
                                 SleqpSparseVec* gradient,
                                 double penalty)
{
  SLEQP_CALL(create_objective(cauchy_data,
                              gradient,
                              penalty));

  SLEQP_CALL(sleqp_lpi_set_objective(cauchy_data->lp_interface,
                                     cauchy_data->objective));

  SLEQP_CALL(sleqp_lpi_solve(cauchy_data->lp_interface));

#if !defined(NDEBUG)
  {
    bool valid_basis = false;

    SLEQP_CALL(check_basis(cauchy_data, &valid_basis));

    assert(valid_basis);
  }

  // TODO: Find out about tolerance guarantees given by the LP solvers
  /*
  {
    bool valid_direction = false;

    SLEQP_CALL(check_direction_bounds(cauchy_data, &valid_direction));

    assert(valid_direction);
  }
  */
#endif

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_cauchy_get_working_set(SleqpCauchyData* cauchy_data,
                                           SleqpIterate* iterate)
{
  SleqpWorkingSet* working_set = sleqp_iterate_get_working_set(iterate);

  SLEQP_CALL(sleqp_working_set_reset(working_set));

  assert(sleqp_working_set_valid(working_set));

  SLEQP_CALL(sleqp_lpi_get_varstats(cauchy_data->lp_interface,
                                    cauchy_data->var_stats));

  SLEQP_CALL(sleqp_lpi_get_consstats(cauchy_data->lp_interface,
                                     cauchy_data->cons_stats));

  const double trust_radius = cauchy_data->trust_radius;

  assert(trust_radius > 0.);

  SleqpProblem* problem = cauchy_data->problem;

  const int num_variables = problem->num_variables;
  const int num_constraints = problem->num_constraints;

  const double eps = sleqp_params_get_eps(cauchy_data->params);

  {
    SleqpSparseVec* x = sleqp_iterate_get_primal(iterate);
    SleqpSparseVec* lb = problem->var_lb;
    SleqpSparseVec* ub = problem->var_ub;

    int k_x = 0, k_lb = 0, k_ub = 0;

    for(int i = 0; i < num_variables; ++i)
    {
      while(k_x < x->nnz && x->indices[k_x] < i)
      {
        ++k_x;
      }

      while(k_lb < lb->nnz && lb->indices[k_lb] < i)
      {
        ++k_lb;
      }

      while(k_ub < ub->nnz && ub->indices[k_ub] < i)
      {
        ++k_ub;
      }

      const bool valid_x = (k_x < x->nnz) && (i == x->indices[k_x]);
      const bool valid_ub = (k_ub < ub->nnz) && (i == ub->indices[k_ub]);
      const bool valid_lb = (k_lb < lb->nnz) && (i == lb->indices[k_lb]);

      const double ubval = valid_ub ? ub->data[k_ub] : 0.;
      const double lbval = valid_lb ? lb->data[k_lb] : 0.;
      const double xval = valid_x ? x->data[k_x] : 0.;

      assert(cauchy_data->var_stats[i] != SLEQP_BASESTAT_ZERO);

      sleqp_assert_is_leq(lbval, xval, eps);
      sleqp_assert_is_leq(xval, ubval, eps);

      if(sleqp_is_eq(lbval, ubval, eps))
      {
        SLEQP_CALL(sleqp_working_set_add_variable(working_set,
                                                  i,
                                                  SLEQP_ACTIVE_BOTH));
      }
      else if((cauchy_data->var_stats[i] == SLEQP_BASESTAT_LOWER) &&
              ((xval - lbval) < trust_radius))
      {
        SLEQP_CALL(sleqp_working_set_add_variable(working_set,
                                                  i,
                                                  SLEQP_ACTIVE_LOWER));
      }
      else if((cauchy_data->var_stats[i] == SLEQP_BASESTAT_UPPER) &&
              ((ubval - xval) < trust_radius))
      {
        SLEQP_CALL(sleqp_working_set_add_variable(working_set,
                                                  i,
                                                  SLEQP_ACTIVE_UPPER));
      }
    }
  }

  {
    int k_lb = 0, k_ub = 0;

    SLEQP_BASESTAT* lower_slack_stats = cauchy_data->var_stats + num_variables;
    SLEQP_BASESTAT* upper_slack_stats = lower_slack_stats + num_constraints;

    SleqpSparseVec* lb = problem->cons_lb;
    SleqpSparseVec* ub = problem->cons_ub;

    for(int i = 0; i < num_constraints; ++i)
    {
      const SLEQP_BASESTAT cons_stat = cauchy_data->cons_stats[i];

      if(cons_stat == SLEQP_BASESTAT_BASIC)
      {
        continue;
      }

      while(k_lb < lb->nnz && lb->indices[k_lb] < i)
      {
        ++k_lb;
      }

      while(k_ub < ub->nnz && ub->indices[k_ub] < i)
      {
        ++k_ub;
      }

      const bool valid_ub = (k_ub < ub->nnz) && (i == ub->indices[k_ub]);
      const bool valid_lb = (k_lb < lb->nnz) && (i == lb->indices[k_lb]);

      const double ubval = valid_ub ? ub->data[k_ub] : 0.;
      const double lbval = valid_lb ? lb->data[k_lb] : 0.;

      assert(lower_slack_stats[i] != SLEQP_BASESTAT_BASIC ||
             upper_slack_stats[i] != SLEQP_BASESTAT_BASIC);

      const bool zero_slack = lower_slack_stats[i] == SLEQP_BASESTAT_LOWER &&
        upper_slack_stats[i] == SLEQP_BASESTAT_LOWER;

      if(cons_stat == SLEQP_BASESTAT_ZERO)
      {
        assert(sleqp_is_inf(ubval));
        assert(sleqp_is_inf(-lbval));

        continue;
      }


      if(sleqp_is_eq(lbval, ubval, eps))
      {
        if(zero_slack)
        {
          SLEQP_CALL(sleqp_working_set_add_constraint(working_set,
                                                      i,
                                                      SLEQP_ACTIVE_BOTH));
        }
      }
      else
      {
        if(cauchy_data->cons_stats[i] == SLEQP_BASESTAT_UPPER)
        {
          if(zero_slack)
          {

            // the constraint c(x) + grad(c(x))*d +I s_l -I s_u  <= u
            // is tight at i

            SLEQP_CALL(sleqp_working_set_add_constraint(working_set,
                                                        i,
                                                        SLEQP_ACTIVE_UPPER));
          }
        }
        else if(cauchy_data->cons_stats[i] == SLEQP_BASESTAT_LOWER)
        {
          if(zero_slack)
          {

            // the constraint l <= c(x) + grad(c(x))*d +I s_l -I s_u
            // is tight at i

            SLEQP_CALL(sleqp_working_set_add_constraint(working_set,
                                                        i,
                                                        SLEQP_ACTIVE_LOWER));
          }
        }
      }
    }
  }

  const int num_active_vars = sleqp_working_set_num_active_vars(working_set);
  const int num_active_cons = sleqp_working_set_num_active_cons(working_set);

  sleqp_log_debug("Created an active set with %d variables, %d constraints",
                  num_active_vars,
                  num_active_cons);

  assert(sleqp_working_set_valid(working_set));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_cauchy_get_direction(SleqpCauchyData* cauchy_data,
                                         SleqpSparseVec* direction)
{
  const double zero_eps = sleqp_params_get_zero_eps(cauchy_data->params);

  SLEQP_CALL(sleqp_lpi_get_primal_sol(cauchy_data->lp_interface,
                                      NULL,
                                      cauchy_data->solution_values));

  SLEQP_CALL(sleqp_sparse_vector_from_raw(direction,
                                          cauchy_data->solution_values,
                                          cauchy_data->problem->num_variables,
                                          zero_eps));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_cauchy_locally_infeasible(SleqpCauchyData* cauchy_data,
                                              bool* locally_infeasible)
{
  *locally_infeasible = false;

  SleqpIterate* iterate = cauchy_data->iterate;

  assert(iterate);

  SLEQP_CALL(sleqp_lpi_get_varstats(cauchy_data->lp_interface,
                                    cauchy_data->var_stats));

  SLEQP_CALL(sleqp_lpi_get_consstats(cauchy_data->lp_interface,
                                     cauchy_data->cons_stats));

  const double trust_radius = cauchy_data->trust_radius;

  assert(trust_radius > 0.);

  SleqpProblem* problem = cauchy_data->problem;

  const int num_variables = problem->num_variables;
  const int num_constraints = problem->num_constraints;

  const double eps = sleqp_params_get_eps(cauchy_data->params);

  bool active_trust_region = false;

  // Check if trust region is active
  {
    SleqpSparseVec* x = sleqp_iterate_get_primal(iterate);
    SleqpSparseVec* lb = problem->var_lb;
    SleqpSparseVec* ub = problem->var_ub;

    int k_x = 0, k_lb = 0, k_ub = 0;

    for(int i = 0; i < num_variables; ++i)
    {
      while(k_x < x->nnz && x->indices[k_x] < i)
      {
        ++k_x;
      }

      while(k_lb < lb->nnz && lb->indices[k_lb] < i)
      {
        ++k_lb;
      }

      while(k_ub < ub->nnz && ub->indices[k_ub] < i)
      {
        ++k_ub;
      }

      const bool valid_x = (k_x < x->nnz) && (i == x->indices[k_x]);
      const bool valid_ub = (k_ub < ub->nnz) && (i == ub->indices[k_ub]);
      const bool valid_lb = (k_lb < lb->nnz) && (i == lb->indices[k_lb]);

      const double ubval = valid_ub ? ub->data[k_ub] : 0.;
      const double lbval = valid_lb ? lb->data[k_lb] : 0.;
      const double xval = valid_x ? x->data[k_x] : 0.;

      assert(cauchy_data->var_stats[i] != SLEQP_BASESTAT_ZERO);

      sleqp_assert_is_leq(lbval, xval, eps);
      sleqp_assert_is_leq(xval, ubval, eps);

      if(sleqp_is_eq(lbval, ubval, eps))
      {
        continue;
      }
      else if(cauchy_data->var_stats[i] == SLEQP_BASESTAT_LOWER &&
              (xval - lbval) >= trust_radius)
      {
        active_trust_region = true;
        break;
      }
      else if(cauchy_data->var_stats[i] == SLEQP_BASESTAT_UPPER &&
              ((ubval - xval) >= trust_radius))
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

    SleqpSparseVec* lb = problem->cons_lb;
    SleqpSparseVec* ub = problem->cons_ub;

    for(int i = 0; i < num_constraints; ++i)
    {
      const SLEQP_BASESTAT cons_stat = cauchy_data->cons_stats[i];

      if(cons_stat == SLEQP_BASESTAT_BASIC)
      {
        continue;
      }

      while(k_lb < lb->nnz && lb->indices[k_lb] < i)
      {
        ++k_lb;
      }

      while(k_ub < ub->nnz && ub->indices[k_ub] < i)
      {
        ++k_ub;
      }

      const bool valid_ub = (k_ub < ub->nnz) && (i == ub->indices[k_ub]);
      const bool valid_lb = (k_lb < lb->nnz) && (i == lb->indices[k_lb]);

      const double ubval = valid_ub ? ub->data[k_ub] : 0.;
      const double lbval = valid_lb ? lb->data[k_lb] : 0.;

      assert(lower_slack_stats[i] != SLEQP_BASESTAT_BASIC ||
             upper_slack_stats[i] != SLEQP_BASESTAT_BASIC);

      bool zero_slack = lower_slack_stats[i] == SLEQP_BASESTAT_LOWER &&
        upper_slack_stats[i] == SLEQP_BASESTAT_LOWER;

      if(!zero_slack)
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

SLEQP_RETCODE sleqp_cauchy_get_dual_estimation(SleqpCauchyData* cauchy_data,
                                               SleqpIterate* iterate)
{
  SleqpSparseVec* vars_dual = sleqp_iterate_get_vars_dual(iterate);
  SleqpSparseVec* cons_dual = sleqp_iterate_get_cons_dual(iterate);

  const double zero_eps = sleqp_params_get_zero_eps(cauchy_data->params);

  SleqpProblem* problem = cauchy_data->problem;

  const int num_variables = problem->num_variables;
  const int num_constraints = problem->num_constraints;

  SleqpWorkingSet* working_set = sleqp_iterate_get_working_set(iterate);

  if(vars_dual)
  {
    assert(vars_dual->dim == problem->num_variables);

    SLEQP_CALL(sleqp_lpi_get_dual_sol(cauchy_data->lp_interface,
                                      cauchy_data->dual_values,
                                      NULL));

    SLEQP_CALL(sleqp_sparse_vector_from_raw(vars_dual,
                                            cauchy_data->dual_values,
                                            vars_dual->dim,
                                            zero_eps));

    // Note: We rescale here since sign conventions vary...
    SLEQP_CALL(sleqp_sparse_vector_scale(vars_dual, -1.));

    for(int k = 0; k < vars_dual->nnz; ++k)
    {
      const int i = vars_dual->indices[k];

      const SLEQP_ACTIVE_STATE var_state = sleqp_working_set_get_variable_state(working_set,
                                                                                i);

      if(var_state == SLEQP_INACTIVE)
      {
        vars_dual->data[k] = 0.;
      }
      else
      {
        if(var_state == SLEQP_ACTIVE_UPPER)
        {
          sleqp_assert_is_geq(vars_dual->data[k], 0., zero_eps);
        }
        else if(var_state == SLEQP_ACTIVE_LOWER)
        {
          sleqp_assert_is_leq(vars_dual->data[k], 0., zero_eps);
        }
      }
    }
  }

  if(cons_dual)
  {
    assert(cons_dual->dim == problem->num_constraints);

    SLEQP_CALL(sleqp_lpi_get_dual_sol(cauchy_data->lp_interface,
                                      NULL,
                                      cauchy_data->dual_values));

    SLEQP_CALL(sleqp_sparse_vector_from_raw(cons_dual,
                                            cauchy_data->dual_values,
                                            cons_dual->dim,
                                            zero_eps));

    // Note: We rescale here since sign conventions vary...
    SLEQP_CALL(sleqp_sparse_vector_scale(cons_dual, -1.));

    for(int k = 0; k < cons_dual->nnz; ++k)
    {
      const int i = cons_dual->indices[k];

      const SLEQP_ACTIVE_STATE var_state = sleqp_working_set_get_constraint_state(working_set,
                                                                                  i);

      if(var_state == SLEQP_INACTIVE)
      {
        cons_dual->data[k] = 0.;
      }
      else
      {
        if(var_state == SLEQP_ACTIVE_UPPER)
        {
          sleqp_assert_is_geq(cons_dual->data[k], 0., zero_eps);
        }
        else if(var_state == SLEQP_ACTIVE_LOWER)
        {
          sleqp_assert_is_leq(cons_dual->data[k], 0., zero_eps);
        }
      }
    }
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_cauchy_get_violation(SleqpCauchyData* cauchy_data,
                                         double* violation)
{
  SleqpProblem* problem = cauchy_data->problem;

  SLEQP_CALL(sleqp_lpi_get_primal_sol(cauchy_data->lp_interface,
                                      NULL,
                                      cauchy_data->solution_values));

  (*violation) = 0.;

  for(int j = problem->num_variables; j < cauchy_data->num_lp_variables; ++j)
  {
    (*violation) += cauchy_data->solution_values[j];
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_cauchy_data_free(SleqpCauchyData** star)
{
  SleqpCauchyData* data = *star;

  if(!data)
  {
    return SLEQP_OKAY;
  }

  //SLEQP_CALL(sleqp_merit_data_free(&data->merit_data));

  SLEQP_CALL(sleqp_sparse_vector_free(&data->quadratic_gradient));

  sleqp_free(&data->dual_values);
  sleqp_free(&data->solution_values);

  sleqp_free(&data->vars_ub);
  sleqp_free(&data->vars_lb);

  sleqp_free(&data->cons_ub);
  sleqp_free(&data->cons_lb);

  sleqp_free(&data->objective);

  sleqp_free(&data->cons_stats);
  sleqp_free(&data->var_stats);

  SLEQP_CALL(sleqp_iterate_release(&data->iterate));

  sleqp_free(star);

  return SLEQP_OKAY;
}
