#include "sleqp_problem_scaling.h"

#include "sleqp_cmp.h"
#include "sleqp_mem.h"

struct SleqpProblemScaling
{
  int refcount;

  SleqpScalingData* scaling_data;
  SleqpProblem* problem;
  SleqpParams* params;
  SleqpFunc* func;

  SleqpFunc* scaled_func;
  SleqpProblem* scaled_problem;

  SleqpSparseVec* unscaled_value;

  SleqpSparseVec* scaled_direction;
  SleqpSparseVec* scaled_cons_duals;
};

static SLEQP_RETCODE
scaled_func_set_value(SleqpSparseVec* scaled_value,
                      SLEQP_VALUE_REASON reason,
                      int num_variables,
                      int* func_grad_nnz,
                      int* cons_val_nnz,
                      int* cons_jac_nnz,
                      void* func_data)
{
  SleqpProblemScaling* problem_scaling = (SleqpProblemScaling*) func_data;
  SleqpScalingData* scaling_data = problem_scaling->scaling_data;

  SLEQP_CALL(sleqp_sparse_vector_copy(scaled_value,
                                      problem_scaling->unscaled_value));

  SLEQP_CALL(sleqp_unscale_point(scaling_data,
                                 problem_scaling->unscaled_value));

  SLEQP_CALL(sleqp_func_set_value(problem_scaling->func,
                                  problem_scaling->unscaled_value,
                                  reason,
                                  func_grad_nnz,
                                  cons_val_nnz,
                                  cons_jac_nnz));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
scaled_func_eval(int num_variables,
                 const SleqpSparseVec* cons_indices,
                 double* func_val,
                 SleqpSparseVec* func_grad,
                 SleqpSparseVec* cons_val,
                 SleqpSparseMatrix* cons_jac,
                 void* func_data)
{
  SleqpProblemScaling* problem_scaling = (SleqpProblemScaling*) func_data;
  SleqpScalingData* scaling_data = problem_scaling->scaling_data;

  SLEQP_CALL(sleqp_func_eval(problem_scaling->func,
                             cons_indices,
                             func_val,
                             func_grad,
                             cons_val,
                             cons_jac));

  if(func_val)
  {
    (*func_val) = sleqp_scale_func_val(scaling_data, (*func_val));
  }

  if(func_grad)
  {
    SLEQP_CALL(sleqp_scale_func_grad(scaling_data, func_grad));
  }

  if(cons_val)
  {
    SLEQP_CALL(sleqp_scale_cons_val(scaling_data, cons_val));
  }

  if(cons_jac)
  {
    SLEQP_CALL(sleqp_scale_cons_jac(scaling_data, cons_jac));
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
scaled_func_hess_prod(int num_variables,
                      const double* func_dual,
                      const SleqpSparseVec* direction,
                      const SleqpSparseVec* cons_duals,
                      SleqpSparseVec* product,
                      void* func_data)
{
  SleqpProblemScaling* problem_scaling = (SleqpProblemScaling*) func_data;
  SleqpScalingData* scaling_data = problem_scaling->scaling_data;

  SLEQP_CALL(sleqp_sparse_vector_copy(direction,
                                      problem_scaling->scaled_direction));

  SLEQP_CALL(sleqp_sparse_vector_copy(cons_duals,
                                      problem_scaling->scaled_cons_duals));

  SLEQP_CALL(sleqp_unscale_hessian_direction(scaling_data,
                                             problem_scaling->scaled_direction,
                                             problem_scaling->scaled_cons_duals));

  SLEQP_CALL(sleqp_func_hess_prod(problem_scaling->func,
                                  func_dual,
                                  problem_scaling->scaled_direction,
                                  problem_scaling->scaled_cons_duals,
                                  product));

  SLEQP_CALL(sleqp_scale_hessian_product(scaling_data,
                                         product));

  return SLEQP_OKAY;
}


SLEQP_RETCODE sleqp_problem_scaling_create(SleqpProblemScaling** star,
                                          SleqpScalingData* scaling_data,
                                          SleqpProblem* problem,
                                          SleqpParams* params)
{
  SLEQP_CALL(sleqp_malloc(star));

  SleqpProblemScaling* problem_scaling = *star;

  *problem_scaling = (SleqpProblemScaling) {0};

  const int num_variables = problem->num_variables;
  const int num_constraints = problem->num_constraints;

  if(num_variables != sleqp_scaling_get_num_variables(scaling_data))
  {
    sleqp_log_error("Invalid number of variables provided to scaled problem");
    return SLEQP_ILLEGAL_ARGUMENT;
  }

  if(num_constraints != sleqp_scaling_get_num_constraints(scaling_data))
  {
    sleqp_log_error("Invalid number of constraints provided to scaled problem");
    return SLEQP_ILLEGAL_ARGUMENT;
  }

  problem_scaling->refcount = 1;

  problem_scaling->problem = problem;

  problem_scaling->func = problem->func;

  problem_scaling->params = params;

  problem_scaling->scaling_data = scaling_data;

  SLEQP_CALL(sleqp_scaling_capture(problem_scaling->scaling_data));

  SleqpFuncCallbacks callbacks = {
    .set_value = scaled_func_set_value,
    .func_eval = scaled_func_eval,
    .hess_prod = scaled_func_hess_prod,
    .func_free = NULL
  };

  SLEQP_CALL(sleqp_func_create(&(problem_scaling->scaled_func),
                               &callbacks,
                               num_variables,
                               problem_scaling));

  SLEQP_CALL(sleqp_problem_create(&(problem_scaling->scaled_problem),
                                  problem_scaling->scaled_func,
                                  params,
                                  problem->var_lb,
                                  problem->var_ub,
                                  problem->cons_lb,
                                  problem->cons_ub));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&(problem_scaling->unscaled_value),
                                              num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&(problem_scaling->scaled_direction),
                                              num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&(problem_scaling->scaled_cons_duals),
                                              num_constraints));
}

SleqpProblem* sleqp_problem_scaling_get_problem(SleqpProblemScaling* problem_scaling)
{
  return problem_scaling->scaled_problem;
}

SLEQP_RETCODE sleqp_problem_scaling_set_func(SleqpProblemScaling* problem_scaling,
                                            SleqpFunc* func)
{
  problem_scaling->func = func;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_problem_scaling_flush(SleqpProblemScaling* problem_scaling)
{
  SleqpProblem* problem = problem_scaling->problem;
  SleqpScalingData* scaling_data = problem_scaling->scaling_data;
  SleqpProblem* scaled_problem = problem_scaling->scaled_problem;

  SLEQP_CALL(sleqp_sparse_vector_copy(problem->var_lb,
                                      scaled_problem->var_lb));

  SLEQP_CALL(sleqp_scale_point(scaling_data,
                               scaled_problem->var_lb));

  SLEQP_CALL(sleqp_sparse_vector_copy(problem->var_ub,
                                      scaled_problem->var_ub));

  SLEQP_CALL(sleqp_scale_point(scaling_data,
                               scaled_problem->var_ub));

  SLEQP_CALL(sleqp_sparse_vector_copy(problem->cons_lb,
                                      scaled_problem->cons_lb));

  SLEQP_CALL(sleqp_scale_cons_val(scaling_data,
                                  scaled_problem->cons_lb));

  SLEQP_CALL(sleqp_sparse_vector_copy(problem->cons_ub,
                                      scaled_problem->cons_ub));

  SLEQP_CALL(sleqp_scale_cons_val(scaling_data,
                                  scaled_problem->cons_ub));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE problem_scaling_free(SleqpProblemScaling** star)
{
  SleqpProblemScaling* problem_scaling = *star;

  if(!problem_scaling)
  {
    return SLEQP_OKAY;
  }

  SLEQP_CALL(sleqp_sparse_vector_free(&(problem_scaling->scaled_cons_duals)));

  SLEQP_CALL(sleqp_sparse_vector_free(&(problem_scaling->scaled_direction)));

  SLEQP_CALL(sleqp_sparse_vector_free(&(problem_scaling->unscaled_value)));

  SLEQP_CALL(sleqp_problem_free(&(problem_scaling->scaled_problem)));

  SLEQP_CALL(sleqp_func_release(&(problem_scaling->scaled_func)));

  SLEQP_CALL(sleqp_scaling_release(&problem_scaling->scaling_data));

  sleqp_free(star);

  *star = NULL;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_problem_scaling_capture(SleqpProblemScaling* scaling)
{
  ++scaling->refcount;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_problem_scaling_release(SleqpProblemScaling** star)
{
  SleqpProblemScaling* problem_scaling = *star;

  if(!problem_scaling)
  {
    return SLEQP_OKAY;
  }

  if(--problem_scaling->refcount == 0)
  {
    SLEQP_CALL(problem_scaling_free(star));
  }

  *star = NULL;

  return SLEQP_OKAY;
}
