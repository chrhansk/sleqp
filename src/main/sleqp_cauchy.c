#include "sleqp_cauchy.h"

#include "sleqp_cmp.h"
#include "sleqp_mem.h"

struct SleqpCauchyData
{
  size_t num_variables;
  size_t num_constraints;

  double* objective;
  double* cons_lb;
  double* cons_ub;
  double* vars_lb;
  double* vars_ub;
};

SLEQP_RETCODE sleqp_cauchy_data_create(SleqpCauchyData** star,
                                       SleqpProblem* problem)
{
  SLEQP_CALL(sleqp_malloc(star));

  SleqpCauchyData* data = *star;

  data->num_variables = problem->num_variables + 2 * problem->num_constraints;
  data->num_constraints = problem->num_constraints;

  SLEQP_CALL(sleqp_calloc(&data->objective, data->num_variables));

  SLEQP_CALL(sleqp_calloc(&data->cons_lb, data->num_constraints));
  SLEQP_CALL(sleqp_calloc(&data->cons_ub, data->num_constraints));

  SLEQP_CALL(sleqp_calloc(&data->vars_lb, data->num_variables));
  SLEQP_CALL(sleqp_calloc(&data->vars_ub, data->num_variables));

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

  sleqp_free(&data->vars_ub);
  sleqp_free(&data->vars_lb);

  sleqp_free(&data->cons_ub);
  sleqp_free(&data->cons_lb);

  sleqp_free(&data->objective);


  sleqp_free(star);

  return SLEQP_OKAY;
}

/*
  static SLEQP_RETCODE prepare_lp(SleqpSolver* solver)
  {
  const double radius = solver->lp_trust_radius;

  SleqpProblem* problem = solver->problem;

  // set variable lower bounds
  {
  size_t k_lb = 0;

  for(size_t j = 0; j < problem->num_variables; ++j)
  {
  while(k_lb < problem->var_lb->nnz && problem->var_lb->indices[k_lb] < j)
  {
  ++k_lb;
  }

  if(k_lb < problem->var_lb->nnz && problem->var_lb->indices[k_lb] == j)
  {
  solver->lp_vars_lb[j] = SLEQP_MAX(-radius, problem->var_lb->data[k_lb]);
  }
  else
  {
  solver->lp_vars_lb[j] = 0;
  }

  }
  }

  // set variable upper bounds
  {
  size_t k_ub = 0;

  for(size_t j = 0; j < problem->num_variables; ++j)
  {
  while(k_ub < problem->var_ub->nnz && problem->var_ub->indices[k_ub] < j)
  {
  ++k_ub;
  }

  if(k_ub < problem->var_ub->nnz && problem->var_ub->indices[k_ub] == j)
  {
  solver->lp_vars_ub[j] = SLEQP_MIN(radius, problem->var_ub->data[k_ub]);
  }
  else
  {
  solver->lp_vars_ub[j] = 0;
  }
  }
  }

  // append +/- id
  {

  }

  // set objective

  {
  size_t k = 0;

  for(size_t i = 0; i < problem->num_variables; ++i)
  {
  while(k < solver->func_grad->nnz && solver->func_grad->indices[k] < i)
  {
  ++k;
  }

  solver->lp_objective[i] = (k < solver->func_grad->nnz && solver->func_grad->indices[k] == i)
  ? solver->func_grad->data[k] : 0;
  }

  for(size_t i = problem->num_variables; i < solver->num_lp_variables; ++i)
  {
  solver->lp_objective[i] = solver->penalty;
  }
  }

  return SLEQP_OKAY;
  }
*/

static SLEQP_RETCODE append_identities(SleqpSparseMatrix* cons_jac,
                                       size_t num_variables,
                                       size_t num_constraints)
{
  assert(num_constraints == cons_jac->num_rows);
  assert(num_variables == cons_jac->num_cols);

  size_t nnz = cons_jac->nnz;

  assert(nnz + 2*num_constraints <= cons_jac->nnz);

  for(size_t i = 0; i < num_constraints; ++i)
  {
    cons_jac->data[nnz + i] = -1.;
    cons_jac->rows[nnz + i] = i;
    cons_jac->cols[num_variables + i] = nnz + i;
  }

  nnz += num_constraints;

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

SLEQP_RETCODE sleqp_cauchy_direction(SleqpProblem* problem,
                                     SleqpIterate* iterate,
                                     SleqpCauchyData* cauchy_data,
                                     SleqpLPi* lp_interface,
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

  SLEQP_CALL(sleqp_lpi_solve(lp_interface,
                             cauchy_data->objective,
                             cons_jac,
                             cauchy_data->cons_lb,
                             cauchy_data->cons_ub,
                             cauchy_data->vars_lb,
                             cauchy_data->vars_ub));

  return SLEQP_OKAY;
}
