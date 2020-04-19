#include "sleqp_cauchy.h"

#include <math.h>

#include "sleqp_cmp.h"
#include "sleqp_mem.h"
#include "sleqp_merit.h"

struct SleqpCauchyData
{
  SleqpProblem* problem;
  SleqpParams* params;

  int num_lp_variables;
  int num_lp_constraints;

  SLEQP_BASESTAT* var_stats;
  SLEQP_BASESTAT* cons_stats;

  SleqpLPi* lp_interface;

  double* objective;
  double* cons_lb;
  double* cons_ub;
  double* vars_lb;
  double* vars_ub;

  double* solution_values;

  SleqpSparseVec* quadratic_gradient;
  SleqpMeritData* merit_data;
};

SLEQP_RETCODE sleqp_cauchy_data_create(SleqpCauchyData** star,
                                       SleqpProblem* problem,
                                       SleqpParams* params,
                                       SleqpLPi* lp_interface)
{
  SLEQP_CALL(sleqp_malloc(star));

  SleqpCauchyData* data = *star;

  data->problem = problem;
  data->params = params;

  data->num_lp_variables = problem->num_variables + 2 * problem->num_constraints;
  data->num_lp_constraints = problem->num_constraints;

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

  SLEQP_CALL(sleqp_sparse_vector_create(&data->quadratic_gradient,
                                        problem->num_variables,
                                        0));

  SLEQP_CALL(sleqp_merit_data_create(&data->merit_data,
                                     problem,
                                     params));

  for(int i = problem->num_variables; i < data->num_lp_variables; ++i)
  {
    data->vars_lb[i] = 0;
    data->vars_ub[i] = sleqp_infinity();
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE append_identities(SleqpSparseMatrix* cons_jac,
                                       int num_variables,
                                       int num_constraints)
{
  assert(num_constraints == cons_jac->num_rows);
  assert(num_variables == cons_jac->num_cols);

  /*
   * Reserve a litte more so we can add the two
   * identity matrices afterwards
   */
  SLEQP_CALL(sleqp_sparse_matrix_reserve(cons_jac, cons_jac->nnz + 2*num_constraints));

  SLEQP_CALL(sleqp_sparse_matrix_resize(cons_jac,
                                        cons_jac->num_rows,
                                        cons_jac->num_cols + 2*num_constraints));

  // append a + id
  for(int i = 0; i < num_constraints; ++i)
  {
    int con = num_variables + i;

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
    int con = num_variables + num_constraints + i;

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
  assert(num_constraints == cons_jac->num_rows);
  assert(num_variables + 2*num_constraints == cons_jac->num_cols);

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

  SleqpSparseVec* val = iterate->cons_val;

  const double eps = sleqp_params_get_eps(cauchy_data->params);

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

    bool valid_ub = (k_ub < ub->nnz && ub->indices[k_ub] == i);
    bool valid_lb = (k_lb < lb->nnz && lb->indices[k_lb] == i);
    bool valid_c = (k_c < val->nnz && val->indices[k_c] == i);

    const double ubval = valid_ub ? ub->data[k_ub] : 0;
    const double lbval = valid_lb ? lb->data[k_lb] : 0;
    const double cval = valid_c ? val->data[k_c] : 0;

    cauchy_data->cons_ub[i] = ubval - cval;
    cauchy_data->cons_lb[i] = lbval - cval;

    assert(sleqp_le(cauchy_data->cons_lb[i],
                    cauchy_data->cons_ub[i],
                    eps));
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE create_var_bounds(SleqpCauchyData* cauchy_data,
                                       SleqpIterate* iterate,
                                       double trust_radius,
                                       int num_variables,
                                       int num_constraints)
{
  SleqpSparseVec* x = iterate->x;
  SleqpSparseVec* lb = cauchy_data->problem->var_lb;
  SleqpSparseVec* ub = cauchy_data->problem->var_ub;

  const double eps = sleqp_params_get_eps(cauchy_data->params);

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

    bool valid_x = (k_x < x->nnz) && (i == x->indices[k_x]);
    bool valid_ub = (k_ub < ub->nnz) && (i == ub->indices[k_ub]);
    bool valid_lb = (k_lb < lb->nnz) && (i == lb->indices[k_lb]);

    double ubval = valid_ub ? ub->data[k_ub] : 0.;
    double lbval = valid_lb ? lb->data[k_lb] : 0.;
    double xval = valid_x ? x->data[k_x] : 0.;

    cauchy_data->vars_ub[i] = SLEQP_MIN(ubval - xval, trust_radius);
    cauchy_data->vars_lb[i] = SLEQP_MAX(lbval - xval, -trust_radius);

    assert(sleqp_le(cauchy_data->vars_lb[i],
                    cauchy_data->vars_ub[i],
                    eps));
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE create_objective(SleqpCauchyData* cauchy_data,
                                      SleqpSparseVec* gradient,
                                      double penalty)
{
  SleqpProblem* problem = cauchy_data->problem;

  int num_variables = problem->num_variables;
  int num_constraints = problem->num_constraints;

  int num_lp_variables = sleqp_lpi_get_num_variables(cauchy_data->lp_interface);

  assert(num_lp_variables == num_variables + 2*num_constraints);

  for(int i = 0; i < num_lp_variables; ++i)
  {
    cauchy_data->objective[i] = (i >= num_variables) ? penalty : 0.;
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
  SleqpSparseMatrix* cons_jac = iterate->cons_jac;

  int num_variables = cons_jac->num_cols;
  int num_constraints = cons_jac->num_rows;

  SLEQP_CALL(append_identities(cons_jac,
                               num_variables,
                               num_constraints));

  assert(sleqp_sparse_matrix_valid(cons_jac));

  SLEQP_CALL(create_var_bounds(cauchy_data,
                               iterate,
                               trust_radius,
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

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_cauchy_get_active_set(SleqpCauchyData* cauchy_data,
                                          SleqpIterate* iterate,
                                          double trust_radius)
{
  SLEQP_CALL(sleqp_active_set_reset(iterate->active_set));

  SLEQP_CALL(sleqp_lpi_get_varstats(cauchy_data->lp_interface,
                                    cauchy_data->var_stats));

  SLEQP_CALL(sleqp_lpi_get_consstats(cauchy_data->lp_interface,
                                     cauchy_data->cons_stats));

  SleqpProblem* problem = cauchy_data->problem;

  const int num_variables = problem->num_variables;
  const int num_constraints = problem->num_constraints;

  const double eps = sleqp_params_get_eps(cauchy_data->params);

  {
    SleqpSparseVec* x = iterate->x;
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

      bool valid_x = (k_x < x->nnz) && (i == x->indices[k_x]);
      bool valid_ub = (k_ub < ub->nnz) && (i == ub->indices[k_ub]);
      bool valid_lb = (k_lb < lb->nnz) && (i == lb->indices[k_lb]);

      double ubval = valid_ub ? ub->data[k_ub] : 0.;
      double lbval = valid_lb ? lb->data[k_lb] : 0.;
      double xval = valid_x ? x->data[k_x] : 0.;

      assert(cauchy_data->var_stats[i] != SLEQP_BASESTAT_ZERO);

      assert(sleqp_le(lbval, xval, eps));
      assert(sleqp_le(xval, ubval, eps));

      if(sleqp_eq(lbval, ubval, eps))
      {
        SLEQP_CALL(sleqp_active_set_add_variable(iterate->active_set,
                                                 i,
                                                 SLEQP_ACTIVE_BOTH));
      }
      else if((cauchy_data->var_stats[i] == SLEQP_BASESTAT_LOWER) &&
         (-trust_radius < lbval - xval))
      {
        SLEQP_CALL(sleqp_active_set_add_variable(iterate->active_set,
                                                 i,
                                                 SLEQP_ACTIVE_LOWER));
      }
      else if((cauchy_data->var_stats[i] == SLEQP_BASESTAT_UPPER) &&
              (ubval - xval < trust_radius))
      {
        SLEQP_CALL(sleqp_active_set_add_variable(iterate->active_set,
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

      bool zero_slack = lower_slack_stats[i] == SLEQP_BASESTAT_LOWER &&
        upper_slack_stats[i] == SLEQP_BASESTAT_LOWER;

      if(cons_stat == SLEQP_BASESTAT_ZERO)
      {
        assert(sleqp_is_inf(ubval));
        assert(sleqp_is_inf(-lbval));

        continue;
      }


      if(sleqp_eq(lbval, ubval, eps))
      {
        if(zero_slack)
        {
          SLEQP_CALL(sleqp_active_set_add_constraint(iterate->active_set,
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

            SLEQP_CALL(sleqp_active_set_add_constraint(iterate->active_set,
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

            SLEQP_CALL(sleqp_active_set_add_constraint(iterate->active_set,
                                                       i,
                                                       SLEQP_ACTIVE_LOWER));
          }
        }
      }
    }
  }

  const int num_active_vars = sleqp_active_set_num_active_vars(iterate->active_set);
  const int num_active_conss = sleqp_active_set_num_active_conss(iterate->active_set);
  const int active_set_size = sleqp_active_set_size(iterate->active_set);

  assert(num_active_vars <= problem->num_variables);
  assert(num_active_conss <= problem->num_variables);

  assert(active_set_size <= problem->num_variables);

  sleqp_log_debug("Created an active set with %d variables, %d constraints",
                  num_active_vars,
                  num_active_conss);


  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_cauchy_get_direction(SleqpCauchyData* cauchy_data,
                                         SleqpIterate* iterate,
                                         SleqpSparseVec* direction)
{
  const double zero_eps = sleqp_params_get_zero_eps(cauchy_data->params);

  SLEQP_CALL(sleqp_lpi_get_solution(cauchy_data->lp_interface,
                                    NULL,
                                    cauchy_data->solution_values));

  SLEQP_CALL(sleqp_sparse_vector_from_raw(direction,
                                          cauchy_data->solution_values,
                                          cauchy_data->problem->num_variables,
                                          zero_eps));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_cauchy_get_violation(SleqpCauchyData* cauchy_data,
                                         double* violation)
{
  SleqpProblem* problem = cauchy_data->problem;

  SLEQP_CALL(sleqp_lpi_get_solution(cauchy_data->lp_interface,
                                    NULL,
                                    cauchy_data->solution_values));

  (*violation) = 0.;

  for(int i = problem->num_variables; i < cauchy_data->num_lp_variables; ++i)
  {
    (*violation) += cauchy_data->solution_values[i];
  }

  return SLEQP_OKAY;
}


SLEQP_RETCODE sleqp_cauchy_compute_step(SleqpCauchyData* cauchy_data,
                                        SleqpIterate* iterate,
                                        double penalty_parameter,
                                        double trust_radius,
                                        SleqpSparseVec* hessian_direction,
                                        SleqpSparseVec* direction,
                                        double* step_length)
{
  const double eps = sleqp_params_get_eps(cauchy_data->params);

  SleqpMeritData* merit_data = cauchy_data->merit_data;

  double exact_merit_value;

  SLEQP_CALL(sleqp_merit_func(merit_data,
                              iterate,
                              penalty_parameter,
                              &exact_merit_value));

  double hessian_product;

  SLEQP_CALL(sleqp_sparse_vector_dot(direction, hessian_direction, &hessian_product));

  double delta = 1.;

  {
    double direction_norm = sleqp_sparse_vector_norm(direction);

    double direction_factor = trust_radius / direction_norm;

    if(direction_factor < 1)
    {
      delta = direction_factor;

      hessian_product *= (delta*delta);

      SLEQP_CALL(sleqp_sparse_vector_scale(direction, delta));
    }
  }

  const double eta = sleqp_params_get_cauchy_eta(cauchy_data->params);
  const double tau = sleqp_params_get_cauchy_tau(cauchy_data->params);

  int max_it = 10000;
  int it = 0;

  for(it = 0; it < max_it; ++it)
  {
    // check

    {
      double linear_merit_value;

      SLEQP_CALL(sleqp_merit_linear(merit_data,
                                    iterate,
                                    direction,
                                    penalty_parameter,
                                    &linear_merit_value));

      double quadratic_merit_value = linear_merit_value + 0.5* hessian_product;

      if((exact_merit_value - quadratic_merit_value) >= eta*(exact_merit_value - linear_merit_value))
      {
        break;
      }
    }

    // update

    SLEQP_CALL(sleqp_sparse_vector_scale(direction, tau));

    hessian_product *= (tau*tau);

    delta *= tau;
  }

  assert(it != max_it);

  sleqp_log_debug("Cauchy line search converged after %d iterations (final value: %12e)",
                  it,
                  delta);

  if(step_length)
  {
    *step_length = delta;
  }

  assert(sleqp_le(sleqp_sparse_vector_norm(direction), trust_radius, eps));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_cauchy_data_free(SleqpCauchyData** star)
{
  SleqpCauchyData* data = *star;

  if(!data)
  {
    return SLEQP_OKAY;
  }

  SLEQP_CALL(sleqp_merit_data_free(&data->merit_data));

  SLEQP_CALL(sleqp_sparse_vector_free(&data->quadratic_gradient));

  sleqp_free(&data->solution_values);

  sleqp_free(&data->vars_ub);
  sleqp_free(&data->vars_lb);

  sleqp_free(&data->cons_ub);
  sleqp_free(&data->cons_lb);

  sleqp_free(&data->objective);

  sleqp_free(&data->cons_stats);
  sleqp_free(&data->var_stats);

  sleqp_free(star);

  return SLEQP_OKAY;
}
