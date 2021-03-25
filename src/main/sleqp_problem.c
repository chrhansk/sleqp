#include "sleqp_problem.h"

#include <assert.h>
#include <math.h>

#include "sleqp_assert.h"
#include "sleqp_cmp.h"
#include "sleqp_mem.h"

static SLEQP_RETCODE map_pos_inf(SleqpSparseVec* vec, double value)
{
  for(int k = 0; k < vec->nnz - 1; ++k)
  {
    if(isinf(vec->data[k]))
    {
      vec->data[k] = value;
    }
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE map_neg_inf(SleqpSparseVec* vec, double value)
{
  for(int k = 0; k < vec->nnz - 1; ++k)
  {
    if(isinf(-vec->data[k]))
    {
      vec->data[k] = value;
    }
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE check_bounds(const SleqpSparseVec* lb,
                                  const SleqpSparseVec* ub,
                                  bool cons_bounds)
{
  assert(lb->dim == ub->dim);

  int k_lb = 0, k_ub = 0;

  while(k_lb < lb->nnz || k_ub < ub->nnz)
  {
    double lb_val = 0;
    double ub_val = 0;

    bool valid_lb = (k_lb < lb->nnz);
    bool valid_ub = (k_ub < ub->nnz);

    int lb_i = valid_lb ? lb->indices[k_lb] : lb->dim + 1;
    int ub_i = valid_ub ? ub->indices[k_ub] : ub->dim + 1;

    int idx = SLEQP_MIN(lb_i, ub_i);

    valid_lb = valid_lb && (lb_i == idx);
    valid_ub = valid_ub && (ub_i == idx);

    if(valid_lb)
    {
      lb_val = lb->data[k_lb];
    }

    if(valid_ub)
    {
      ub_val = ub->data[k_ub];
    }

    if(lb_val > ub_val)
    {
      if(cons_bounds)
      {
        sleqp_log_error("Inconsistent constraint bound values at index %d: lower = %f > %f = upper",
                        idx,
                        lb_val,
                        ub_val);
      }
      else
      {
        sleqp_log_error("Inconsistent variable bound values at index %d: lower = %f > %f = upper",
                        idx,
                        lb_val,
                        ub_val);
      }

      return SLEQP_ILLEGAL_ARGUMENT;
    }

    if(valid_lb)
    {
      ++k_lb;
    }

    if(valid_ub)
    {
      ++k_ub;
    }
  }


  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_problem_create(SleqpProblem** star,
                                   SleqpFunc* func,
                                   SleqpSparseVec* var_lb,
                                   SleqpSparseVec* var_ub,
                                   SleqpSparseVec* cons_lb,
                                   SleqpSparseVec* cons_ub)
{
  assert(var_lb->dim == var_ub->dim);
  assert(cons_lb->dim == cons_ub->dim);

  sleqp_assert_msg(sleqp_sparse_vector_is_valid(var_lb), "Invalid variable bounds");
  sleqp_assert_msg(sleqp_sparse_vector_is_valid(var_ub), "Invalid variable bounds");

  sleqp_assert_msg(sleqp_sparse_vector_is_valid(cons_lb), "Invalid constraint bounds");
  sleqp_assert_msg(sleqp_sparse_vector_is_valid(cons_ub), "Invalid constraint bounds");

  const int num_constraints = sleqp_func_get_num_constraints(func);
  const int num_variables = sleqp_func_get_num_variables(func);

  SLEQP_CALL(sleqp_malloc(star));

  SleqpProblem* problem = *star;

  *problem = (SleqpProblem) {0};

  problem->num_variables = num_variables;
  problem->num_constraints = num_constraints;

  sleqp_assert_msg(var_lb->dim == num_variables, "Inconsisten variable dimensions");
  sleqp_assert_msg(var_ub->dim == num_variables, "Inconsisten variable dimensions");

  sleqp_assert_msg(cons_lb->dim == num_constraints, "Inconsisten constraint dimensions");
  sleqp_assert_msg(cons_ub->dim == num_constraints, "Inconsisten constraint dimensions");

  SLEQP_CALL(sleqp_func_capture(func));

  problem->func = func;

  SLEQP_CALL(sleqp_sparse_vector_create(&problem->var_lb,
                                        problem->num_variables,
                                        var_lb->nnz));

  SLEQP_CALL(sleqp_sparse_vector_create(&problem->var_ub,
                                        problem->num_variables,
                                        var_ub->nnz));

  SLEQP_CALL(sleqp_sparse_vector_create(&problem->cons_lb,
                                        problem->num_constraints,
                                        cons_lb->nnz));

  SLEQP_CALL(sleqp_sparse_vector_create(&problem->cons_ub,
                                        problem->num_constraints,
                                        cons_ub->nnz));

  SLEQP_CALL(sleqp_sparse_vector_copy(var_lb, problem->var_lb));
  SLEQP_CALL(sleqp_sparse_vector_copy(var_ub, problem->var_ub));

  SLEQP_CALL(sleqp_sparse_vector_copy(cons_lb, problem->cons_lb));
  SLEQP_CALL(sleqp_sparse_vector_copy(cons_ub, problem->cons_ub));

  const double inf = sleqp_infinity();

  SLEQP_CALL(map_pos_inf(problem->var_ub, inf));
  SLEQP_CALL(map_pos_inf(problem->cons_ub, inf));

  SLEQP_CALL(map_neg_inf(problem->var_lb, -inf));
  SLEQP_CALL(map_neg_inf(problem->cons_lb, -inf));

  SLEQP_CALL(check_bounds(problem->var_lb, problem->var_ub, false));
  SLEQP_CALL(check_bounds(problem->cons_lb, problem->cons_ub, true));

  sleqp_assert_msg(sleqp_sparse_vector_is_finite(problem->var_lb), "Infinite variable bounds");
  sleqp_assert_msg(sleqp_sparse_vector_is_finite(problem->var_ub), "Infinite variable bounds");

  sleqp_assert_msg(sleqp_sparse_vector_is_finite(problem->cons_lb), "Infinite constraint bounds");
  sleqp_assert_msg(sleqp_sparse_vector_is_finite(problem->cons_ub), "Infinite constraint bounds");

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_problem_free(SleqpProblem** star)
{
  SleqpProblem* problem = *star;

  if(!problem)
  {
    return SLEQP_OKAY;
  }

  SLEQP_CALL(sleqp_sparse_vector_free(&problem->var_lb));
  SLEQP_CALL(sleqp_sparse_vector_free(&problem->var_ub));

  SLEQP_CALL(sleqp_sparse_vector_free(&problem->cons_lb));
  SLEQP_CALL(sleqp_sparse_vector_free(&problem->cons_ub));

  SLEQP_CALL(sleqp_func_release(&problem->func));

  sleqp_free(star);

  return SLEQP_OKAY;
}
