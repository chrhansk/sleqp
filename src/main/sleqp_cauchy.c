#include "sleqp_cauchy.h"

#include "sleqp_cmp.h"
#include "sleqp_mem.h"
#include "sleqp_penalty.h"

struct SleqpCauchyData
{
  SleqpProblem* problem;

  int num_variables;
  int num_constraints;

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
  SleqpPenalty* penalty_data;
};

SLEQP_RETCODE sleqp_cauchy_data_create(SleqpCauchyData** star,
                                       SleqpProblem* problem,
                                       SleqpLPi* lp_interface)
{
  SLEQP_CALL(sleqp_malloc(star));

  SleqpCauchyData* data = *star;

  data->problem = problem;
  data->num_variables = problem->num_variables + 2 * problem->num_constraints;
  data->num_constraints = problem->num_constraints;

  SLEQP_CALL(sleqp_calloc(&data->var_stats, data->num_variables));
  SLEQP_CALL(sleqp_calloc(&data->cons_stats, data->num_constraints));

  data->lp_interface = lp_interface;

  SLEQP_CALL(sleqp_calloc(&data->objective, data->num_variables));

  SLEQP_CALL(sleqp_calloc(&data->cons_lb, data->num_constraints));
  SLEQP_CALL(sleqp_calloc(&data->cons_ub, data->num_constraints));

  SLEQP_CALL(sleqp_calloc(&data->vars_lb, data->num_variables));
  SLEQP_CALL(sleqp_calloc(&data->vars_ub, data->num_variables));

  SLEQP_CALL(sleqp_calloc(&data->solution_values, data->num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create(&data->quadratic_gradient,
                                        problem->num_variables,
                                        0));

  SLEQP_CALL(sleqp_penalty_create(&data->penalty_data,
                                  problem,
                                  problem->func));

  for(int i = problem->num_variables; i < data->num_variables; ++i)
  {
    data->vars_lb[i] = 0;
    data->vars_ub[i] = sleqp_infinity();
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_cauchy_data_free(SleqpCauchyData** star)
{
  SleqpCauchyData* data = *star;

  SLEQP_CALL(sleqp_penalty_free(&data->penalty_data));

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

    SLEQP_CALL(sleqp_sparse_matrix_add_column(cons_jac,
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

    SLEQP_CALL(sleqp_sparse_matrix_add_column(cons_jac,
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

static SLEQP_RETCODE append_penalties(SleqpSparseVec* func_grad,
                                      int num_variables,
                                      int num_constraints,
                                      double penalty)
{
  assert(num_variables == func_grad->dim);

  int size_increase = 2*num_constraints;

  /*
   * Reserve a litte more so we can add slack penalties
   * without reallocations
   */
  SLEQP_CALL(sleqp_sparse_vector_reserve(func_grad,
                                         func_grad->nnz + size_increase));

  func_grad->dim += size_increase;

  for(int i = 0; i < size_increase; ++i)
  {
    SLEQP_CALL(sleqp_sparse_vector_push(func_grad,
                                        num_variables + i,
                                        penalty));
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE create_cons_bounds(SleqpCauchyData* cauchy_data,
                                        SleqpIterate* iterate,
                                        int num_variables,
                                        int num_constraints)
{
  int k_y = 0, k_lb = 0, k_ub = 0;

  SleqpSparseVec* lb = cauchy_data->problem->cons_lb;
  SleqpSparseVec* ub = cauchy_data->problem->cons_ub;

  SleqpSparseVec* val = iterate->cons_val;

  for(int i = 0; i < num_constraints; ++i)
  {
    while(k_y < val->nnz && val->indices[k_y] < i)
    {
      ++k_y;
    }

    while(k_lb < lb->nnz && lb->indices[k_lb] < i)
    {
      ++k_lb;
    }

    while(k_ub < ub->nnz && ub->indices[k_ub] < i)
    {
      ++k_ub;
    }

    double ubval = (i == k_ub) ? ub->data[k_ub] : 0;
    double lbval = (i == k_lb) ? lb->data[k_lb] : 0;
    double yval = (i == k_y) ? val->data[k_y] : 0;


    cauchy_data->cons_ub[i] = ubval - yval;
    cauchy_data->cons_lb[i] = lbval - yval;

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

    double ubval = (i == k_ub) ? ub->data[k_ub] : 0;
    double lbval = (i == k_lb) ? lb->data[k_lb] : 0;
    double xval = (i == k_x) ? x->data[k_x] : 0;

    cauchy_data->vars_ub[i] = SLEQP_MIN(ubval - xval, trust_radius);
    cauchy_data->vars_lb[i] = SLEQP_MAX(lbval - xval, -trust_radius);
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE create_objective(SleqpCauchyData* cauchy_data,
                                      SleqpIterate* iterate,
                                      double trust_radius,
                                      int num_variables,
                                      int num_constraints)
{
  int index = 0;

  SleqpSparseVec* grad = iterate->func_grad;
  double* objective = cauchy_data->objective;

  for(int i = 0; i < num_variables; ++i)
  {
    if(index < grad->nnz && grad->indices[index] == i)
    {
      objective[i] = grad->data[index++];
    }
    else
    {
      objective[i] = 0.;
    }
  }

  for(int i = num_variables; i < cauchy_data->num_variables; ++i)
  {
    cauchy_data->objective[i] = trust_radius;
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_cauchy_compute_direction(SleqpCauchyData* cauchy_data,
                                             SleqpIterate* iterate,
                                             double penalty,
                                             double trust_radius)
{
  SleqpSparseMatrix* cons_jac = iterate->cons_jac;

  int num_variables = cons_jac->num_cols;
  int num_constraints = cons_jac->num_rows;

  SLEQP_CALL(append_penalties(iterate->func_grad,
                              num_variables,
                              num_constraints,
                              penalty));

  SLEQP_CALL(append_identities(cons_jac,
                               num_variables,
                               num_constraints));

  SLEQP_CALL(create_var_bounds(cauchy_data,
                               iterate,
                               trust_radius,
                               num_variables,
                               num_constraints));

  SLEQP_CALL(create_cons_bounds(cauchy_data,
                                iterate,
                                num_variables,
                                num_constraints));

  SLEQP_CALL(create_objective(cauchy_data,
                              iterate,
                              trust_radius,
                              num_variables,
                              num_constraints));

  SLEQP_CALL(sleqp_lpi_solve(cauchy_data->lp_interface,
                             cauchy_data->objective,
                             cons_jac,
                             cauchy_data->cons_lb,
                             cauchy_data->cons_ub,
                             cauchy_data->vars_lb,
                             cauchy_data->vars_ub));


  SLEQP_CALL(remove_identities(cons_jac,
                               num_variables,
                               num_constraints));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_cauchy_get_active_set(SleqpCauchyData* cauchy_data,
                                          SleqpIterate* iterate,
                                          double trust_radius)
{
  SleqpSparseVec* x = iterate->x;
  SleqpSparseVec* lb = cauchy_data->problem->var_lb;
  SleqpSparseVec* ub = cauchy_data->problem->var_ub;

  int num_variables = cauchy_data->problem->num_variables;
  int num_constraints = cauchy_data->problem->num_constraints;

  int k_x = 0, k_lb = 0, k_ub = 0;

  SLEQP_CALL(sleqp_lpi_get_varstats(cauchy_data->lp_interface,
                                    cauchy_data->num_variables,
                                    cauchy_data->var_stats));

  SLEQP_CALL(sleqp_lpi_get_consstats(cauchy_data->lp_interface,
                                     cauchy_data->num_constraints,
                                     cauchy_data->cons_stats));

  SLEQP_ACTIVE_STATE* var_states = sleqp_active_set_var_states(iterate->active_set);

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

    double ubval = (i == k_ub) ? ub->data[k_ub] : 0;
    double lbval = (i == k_lb) ? lb->data[k_lb] : 0;
    double xval = (i == k_x) ? x->data[k_x] : 0;

    if((cauchy_data->var_stats[i] == SLEQP_BASESTAT_UPPER) && sleqp_lt(ubval - xval, trust_radius))
    {
      var_states[i] = SLEQP_ACTIVE_UPPER;
    }
    else if((cauchy_data->var_stats[i] == SLEQP_BASESTAT_LOWER) && sleqp_lt(-trust_radius, lbval - xval))
    {
      var_states[i] = SLEQP_ACTIVE_LOWER;
    }
    else
    {
      var_states[i] = SLEQP_INACTIVE;
    }
  }

  SLEQP_ACTIVE_STATE* cons_states = sleqp_active_set_cons_states(iterate->active_set);

  SLEQP_BASESTAT* lower_slack_stats = cauchy_data->var_stats + num_variables;
  SLEQP_BASESTAT* upper_slack_stats = lower_slack_stats + num_variables;

  for(int i = 0; i < num_constraints; ++i)
  {
    cons_states[i] = SLEQP_INACTIVE;;

    if(cauchy_data->cons_stats[i] == SLEQP_BASESTAT_UPPER)
    {
      if(upper_slack_stats[i] == SLEQP_BASESTAT_LOWER)
      {

        // the constraint c(x) + grad(c(x))*d +I s_l -I s_u  <= u
        // is tight at i, and s_u is zero at i

        cons_states[i] = SLEQP_ACTIVE_UPPER;
      }
    }
    else if(cauchy_data->cons_stats[i] == SLEQP_BASESTAT_LOWER)
    {
      if(lower_slack_stats[i] == SLEQP_BASESTAT_LOWER)
      {

        // the constraint l <= c(x) + grad(c(x))*d +I s_l -I s_u
        // is tight at i, and s_l is zero at i

        cons_states[i] = SLEQP_ACTIVE_UPPER;
      }
    }
  }


  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_cauchy_get_direction(SleqpCauchyData* cauchy_data,
                                         SleqpIterate* iterate,
                                         SleqpSparseVec* direction)
{
  SLEQP_CALL(sleqp_lpi_get_solution(cauchy_data->lp_interface,
                                    cauchy_data->problem->num_variables,
                                    NULL,
                                    cauchy_data->solution_values));

  SLEQP_CALL(sleqp_sparse_vector_from_raw(direction,
                                          cauchy_data->solution_values,
                                          cauchy_data->problem->num_variables));

  return SLEQP_OKAY;
}


// TODO: Fix this to make it consistent with the paper
SLEQP_RETCODE sleqp_cauchy_compute_step(SleqpCauchyData* cauchy_data,
                                        SleqpIterate* iterate,
                                        double penalty_parameter,
                                        SleqpSparseVec* direction,
                                        double* predicted_reduction)
{
  double func_val = iterate->func_val;

  SleqpFunc* func = cauchy_data->problem->func;

  cauchy_data->quadratic_gradient->nnz = 0;

  double zero_value;

  SLEQP_CALL(sleqp_penalty_linear(cauchy_data->penalty_data,
                                  iterate,
                                  cauchy_data->quadratic_gradient,
                                  penalty_parameter,
                                  &zero_value));

  SLEQP_CALL(sleqp_penalty_quadratic_gradient(cauchy_data->penalty_data,
                                              iterate,
                                              penalty_parameter,
                                              cauchy_data->quadratic_gradient));

  double gradient_product;

  SLEQP_CALL(sleqp_sparse_vector_dot(direction,
                                     cauchy_data->quadratic_gradient,
                                     &gradient_product));

  double hessian_product;

  double func_dual = 1.;

  SLEQP_CALL(sleqp_hess_eval_bilinear(func,
                                      &func_dual,
                                      direction,
                                      iterate->cons_dual,
                                      &hessian_product));

  double alpha = 1.;
  double beta = 0.9;

  while(1)
  {
    double quadratic_penalty_value;

    SLEQP_CALL(sleqp_penalty_linear(cauchy_data->penalty_data,
                                    iterate,
                                    direction,
                                    penalty_parameter,
                                    &quadratic_penalty_value));

    quadratic_penalty_value += alpha*alpha*hessian_product;

    /*
    SLEQP_CALL(sleqp_penalty_quadratic(cauchy_data->penalty_data,
                                       iterate,
                                       direction,
                                       penalty_parameter,
                                       &quadratic_penalty_value));
    */

    *predicted_reduction = zero_value - quadratic_penalty_value;

    if(sleqp_le(quadratic_penalty_value, alpha * gradient_product))
    {
      break;
    }

    SLEQP_CALL(sleqp_sparse_vector_scale(direction, beta));

    alpha *= beta;
  }


  return SLEQP_OKAY;
}
