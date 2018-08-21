#include "sleqp_cauchy.h"

#include "sleqp_cmp.h"
#include "sleqp_mem.h"

struct SleqpCauchyData
{
  size_t num_variables;
  size_t num_constraints;

  SLEQP_BASESTAT* base_stats;

  SleqpLPi* lp_interface;

  double* objective;
  double* cons_lb;
  double* cons_ub;
  double* vars_lb;
  double* vars_ub;

  double* solution_values;
};

SLEQP_RETCODE sleqp_cauchy_data_create(SleqpCauchyData** star,
                                       SleqpProblem* problem,
                                       SleqpLPi* lp_interface)
{
  SLEQP_CALL(sleqp_malloc(star));

  SleqpCauchyData* data = *star;

  data->num_variables = problem->num_variables + 2 * problem->num_constraints;
  data->num_constraints = problem->num_constraints;

  SLEQP_CALL(sleqp_calloc(&data->base_stats, data->num_variables));

  data->lp_interface = lp_interface;

  SLEQP_CALL(sleqp_calloc(&data->objective, data->num_variables));

  SLEQP_CALL(sleqp_calloc(&data->cons_lb, data->num_constraints));
  SLEQP_CALL(sleqp_calloc(&data->cons_ub, data->num_constraints));

  SLEQP_CALL(sleqp_calloc(&data->vars_lb, data->num_variables));
  SLEQP_CALL(sleqp_calloc(&data->vars_ub, data->num_variables));

  SLEQP_CALL(sleqp_calloc(&data->solution_values, data->num_variables));

  for(size_t i = problem->num_variables; i < data->num_variables; ++i)
  {
    data->vars_lb[i] = 0;
    data->vars_ub[i] = sleqp_infinity();
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_cauchy_data_free(SleqpCauchyData** star)
{
  SleqpCauchyData* data = *star;

  sleqp_free(&data->solution_values);

  sleqp_free(&data->vars_ub);
  sleqp_free(&data->vars_lb);

  sleqp_free(&data->cons_ub);
  sleqp_free(&data->cons_lb);

  sleqp_free(&data->objective);

  sleqp_free(&data->base_stats);

  sleqp_free(star);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE append_identities(SleqpSparseMatrix* cons_jac,
                                       size_t num_variables,
                                       size_t num_constraints)
{
  assert(num_constraints == cons_jac->num_rows);
  assert(num_variables == cons_jac->num_cols);

  size_t nnz = cons_jac->nnz;

  assert(nnz + 2*num_constraints <= cons_jac->nnz);

  // append a + id
  for(size_t i = 0; i < num_constraints; ++i)
  {
    cons_jac->data[nnz + i] = +1.;
    cons_jac->rows[nnz + i] = i;
    cons_jac->cols[num_variables + i] = nnz + i;
  }

  nnz += num_constraints;

  // append a - id
  for(size_t i = 0; i < num_constraints; ++i)
  {
    cons_jac->data[nnz + i] = -1.;
    cons_jac->rows[nnz + i] = i;
    cons_jac->cols[num_variables + num_constraints + i] = nnz + i;
  }

  nnz += num_constraints;

  cons_jac->nnz = nnz;
  cons_jac->cols[num_variables + num_constraints] = nnz + 1;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE append_penalties(SleqpSparseVec* func_grad,
                                      size_t num_variables,
                                      size_t num_constraints,
                                      double penalty)
{
  assert(num_variables == func_grad->dim);

  size_t size_increase = 2*num_constraints;

  func_grad->dim += size_increase;

  for(size_t i = 0; i < size_increase; ++i)
  {
    SLEQP_CALL(sleqp_sparse_vector_push(func_grad,
                                        num_variables + i,
                                        penalty));
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE create_cons_bounds(SleqpProblem* problem,
                                        SleqpIterate* iterate,
                                        SleqpCauchyData* cauchy_data,
                                        size_t num_variables,
                                        size_t num_constraints)
{
  size_t k_y = 0, k_lb = 0, k_ub = 0;

  SleqpSparseVec* lb = problem->cons_lb;
  SleqpSparseVec* ub = problem->cons_ub;

  SleqpSparseVec* val = iterate->cons_val;

  for(size_t i = 0; i < num_constraints; ++i)
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

static SLEQP_RETCODE create_var_bounds(SleqpProblem* problem,
                                       SleqpIterate* iterate,
                                       SleqpCauchyData* cauchy_data,
                                       double trust_radius,
                                       size_t num_variables,
                                       size_t num_constraints)
{
  SleqpSparseVec* x = iterate->x;
  SleqpSparseVec* lb = problem->var_lb;
  SleqpSparseVec* ub = problem->var_ub;

  size_t k_x = 0, k_lb = 0, k_ub = 0;

  for(size_t i = 0; i < num_variables; ++i)
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

static SLEQP_RETCODE create_objective(SleqpProblem* problem,
                                      SleqpIterate* iterate,
                                      SleqpCauchyData* cauchy_data,
                                      double trust_radius,
                                      size_t num_variables,
                                      size_t num_constraints)
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

SLEQP_RETCODE sleqp_cauchy_compute_direction(SleqpProblem* problem,
                                             SleqpIterate* iterate,
                                             SleqpCauchyData* cauchy_data,
                                             double penalty,
                                             double trust_radius)
{
  SleqpSparseMatrix* cons_jac = iterate->cons_jac;

  size_t num_variables = cons_jac->num_cols;
  size_t num_constraints = cons_jac->num_rows;

  SLEQP_CALL(append_penalties(iterate->func_grad,
                              num_variables,
                              num_constraints,
                              penalty));

  SLEQP_CALL(append_identities(cons_jac,
                               num_variables,
                               num_constraints));

  SLEQP_CALL(create_var_bounds(problem,
                               iterate,
                               cauchy_data,
                               trust_radius,
                               num_variables,
                               num_constraints));

  SLEQP_CALL(create_cons_bounds(problem,
                                iterate,
                                cauchy_data,
                                num_variables,
                                num_constraints));

  SLEQP_CALL(create_objective(problem,
                              iterate,
                              cauchy_data,
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

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_cauchy_get_active_set(SleqpProblem* problem,
                                          SleqpIterate* iterate,
                                          SleqpCauchyData* cauchy_data,
                                          SleqpActiveSet* active_set,
                                          double trust_radius)
{
  SleqpSparseVec* x = iterate->x;
  SleqpSparseVec* lb = problem->var_lb;
  SleqpSparseVec* ub = problem->var_ub;

  size_t num_variables = problem->num_variables;
  size_t num_constraints = problem->num_constraints;

  size_t k_x = 0, k_lb = 0, k_ub = 0;

  SLEQP_CALL(sleqp_lpi_get_varstats(cauchy_data->lp_interface,
                                    cauchy_data->num_variables,
                                    cauchy_data->base_stats));

  SLEQP_ACTIVE_STATE* var_states = sleqp_active_set_var_states(active_set);

  for(size_t i = 0; i < num_variables; ++i)
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

    if((cauchy_data->base_stats[i] == SLEQP_BASESTAT_UPPER) && sleqp_lt(ubval - xval, trust_radius))
    {
      var_states[i] = SLEQP_ACTIVE_UPPER;
    }
    else if((cauchy_data->base_stats[i] == SLEQP_BASESTAT_LOWER) && sleqp_lt(-trust_radius, lbval - xval))
    {
      var_states[i] = SLEQP_ACTIVE_LOWER;
    }
    else
    {
      var_states[i] = SLEQP_INACTIVE;
    }
  }

  SLEQP_ACTIVE_STATE* cons_states = sleqp_active_set_cons_states(active_set);

  for(size_t i = 0; i < num_constraints; ++i)
  {
    if(cauchy_data->base_stats[num_variables + i] == SLEQP_BASESTAT_BASIC)
    {
      cons_states[i] = SLEQP_ACTIVE_UPPER;
    }
    else if(cauchy_data->base_stats[num_variables + num_constraints + i] == SLEQP_BASESTAT_BASIC)
    {
      cons_states[i] = SLEQP_ACTIVE_LOWER;
    }
  }


  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_cauchy_get_direction(SleqpProblem* problem,
                                         SleqpIterate* iterate,
                                         SleqpCauchyData* cauchy_data,
                                         SleqpSparseVec* direction)
{
  SLEQP_CALL(sleqp_lpi_get_solution(cauchy_data->lp_interface,
                                    problem->num_variables,
                                    NULL,
                                    cauchy_data->solution_values));

  SLEQP_CALL(sleqp_sparse_vector_from_raw(direction,
                                          cauchy_data->solution_values,
                                          problem->num_variables));

  return SLEQP_OKAY;
}
