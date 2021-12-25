#include "mex_solve.h"

#include <assert.h>

#include "mex_fields.h"
#include "mex_func.h"
#include "mex_output.h"

static int
array_or_cell_size(const mxArray* array)
{
  if (mxIsCell(array))
  {
    int size               = 0;
    const int num_elements = mxGetNumberOfElements(array);

    for (int i = 0; i < num_elements; ++i)
    {
      const mxArray* p               = mxGetCell(array, i);
      const int current_num_elements = mxGetNumberOfElements(p);

      size += current_num_elements;
    }

    return size;
  }
  else
  {
    return mxGetNumberOfElements(array);
  }
}

static int
num_vars_from_solution(const mxArray* x0)
{
  return array_or_cell_size(x0);
}

static int
num_cons_from_options(const mxArray* options)
{
  assert(mxIsStruct(options));

  const mxArray* cons_lb = mxGetField(options, 0, MEX_INPUT_CONS_LB);

  if (cons_lb)
  {
    return array_or_cell_size(cons_lb);
  }

  return 0;
}

static SLEQP_RETCODE
create_vec_from_array(SleqpSparseVec** star, const mxArray* array)
{
  const int dimension = array_or_cell_size(array);

  SLEQP_CALL(sleqp_sparse_vector_create_full(star, dimension));

  SleqpSparseVec* vec = *star;

  double* data = vec->data;

  if (mxIsCell(array))
  {
    const int num_elements = mxGetNumberOfElements(array);
    int offset             = 0;

    for (int i = 0; i < num_elements; ++i)
    {
      const mxArray* p           = mxGetCell(array, i);
      const int cur_num_elements = mxGetNumberOfElements(p);
      const double* ptr          = mxGetPr(p);

      for (int k = 0; k < cur_num_elements; ++k)
      {
        SLEQP_CALL(sleqp_sparse_vector_push(vec, k + offset, ptr[k]));
      }

      offset += cur_num_elements;
    }
  }
  else
  {
    const double* ptr = mxGetPr(array);

    for (int k = 0; k < dimension; ++k)
    {
      SLEQP_CALL(sleqp_sparse_vector_push(vec, k, ptr[k]));
    }
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
create_vec_from_array_default(SleqpSparseVec** star,
                              const mxArray* array,
                              const int dimension,
                              double value)
{
  if (array)
  {
    SLEQP_CALL(create_vec_from_array(star, array));
  }
  else
  {
    SLEQP_CALL(sleqp_sparse_vector_create_full(star, dimension));
    SLEQP_CALL(sleqp_sparse_vector_fill(*star, value));
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
create_bounds_from_options(SleqpSparseVec** lb_star,
                           SleqpSparseVec** ub_star,
                           const int dimension,
                           const mxArray* lb_array,
                           const mxArray* ub_array)
{
  const double inf = sleqp_infinity();

  SLEQP_CALL(create_vec_from_array_default(lb_star, lb_array, dimension, -inf));

  SLEQP_CALL(create_vec_from_array_default(ub_star, ub_array, dimension, inf));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
create_var_bounds_from_options(SleqpSparseVec** var_lb_star,
                               SleqpSparseVec** var_ub_star,
                               const int num_variables,
                               const mxArray* options)
{
  return create_bounds_from_options(var_lb_star,
                                    var_ub_star,
                                    num_variables,
                                    mxGetField(options, 0, MEX_INPUT_VAR_LB),
                                    mxGetField(options, 0, MEX_INPUT_VAR_UB));
  return SLEQP_OKAY;
}

static SLEQP_RETCODE
create_cons_bounds_from_options(SleqpSparseVec** cons_lb_star,
                                SleqpSparseVec** cons_ub_star,
                                const int num_constraints,
                                const mxArray* options)
{
  return create_bounds_from_options(cons_lb_star,
                                    cons_ub_star,
                                    num_constraints,
                                    mxGetField(options, 0, MEX_INPUT_CONS_LB),
                                    mxGetField(options, 0, MEX_INPUT_CONS_UB));
  return SLEQP_OKAY;
}

typedef struct
{
  SleqpSparseVec* var_lb;
  SleqpSparseVec* var_ub;

  SleqpSparseVec* cons_lb;
  SleqpSparseVec* cons_ub;

  SleqpSparseVec* initial;

  SleqpOptions* options;
  SleqpParams* params;

  SleqpFunc* func;
  SleqpProblem* problem;

  SleqpSolver* solver;

} Instance;

static SLEQP_RETCODE
create_instance(Instance* instance,
                const mxArray* mex_x0,
                const mxArray* mex_funcs,
                const mxArray* mex_options)
{
  *instance = (Instance){0};

  const int num_vars = num_vars_from_solution(mex_x0);
  const int num_cons = num_cons_from_options(mex_options);

  SLEQP_CALL(create_var_bounds_from_options(&instance->var_lb,
                                            &instance->var_ub,
                                            num_vars,
                                            mex_options));

  SLEQP_CALL(create_cons_bounds_from_options(&instance->cons_lb,
                                             &instance->cons_ub,
                                             num_cons,
                                             mex_options));

  SLEQP_CALL(create_vec_from_array(&instance->initial, mex_x0));

  SLEQP_CALL(sleqp_options_create(&instance->options));
  SLEQP_CALL(sleqp_params_create(&instance->params));

  SLEQP_CALL(sleqp_options_set_int_value(instance->options,
                                         SLEQP_OPTION_INT_HESS_EVAL,
                                         SLEQP_HESS_EVAL_DAMPED_BFGS));

  SLEQP_CALL(mex_func_create(&instance->func,
                             mex_funcs,
                             instance->params,
                             num_vars,
                             num_cons));

  SLEQP_CALL(sleqp_problem_create_simple(&instance->problem,
                                         instance->func,
                                         instance->params,
                                         instance->var_lb,
                                         instance->var_ub,
                                         instance->cons_lb,
                                         instance->cons_ub));

  SLEQP_CALL(sleqp_solver_create(&instance->solver,
                                 instance->problem,
                                 instance->params,
                                 instance->options,
                                 instance->initial,
                                 NULL));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
free_instance(Instance* instance)
{
  SLEQP_CALL(sleqp_solver_release(&instance->solver));

  SLEQP_CALL(sleqp_problem_release(&instance->problem));

  SLEQP_CALL(sleqp_func_release(&instance->func));

  SLEQP_CALL(sleqp_params_release(&instance->params));

  SLEQP_CALL(sleqp_options_release(&instance->options));

  SLEQP_CALL(sleqp_sparse_vector_free(&instance->initial));

  SLEQP_CALL(sleqp_sparse_vector_free(&instance->cons_ub));
  SLEQP_CALL(sleqp_sparse_vector_free(&instance->cons_lb));

  SLEQP_CALL(sleqp_sparse_vector_free(&instance->var_ub));
  SLEQP_CALL(sleqp_sparse_vector_free(&instance->var_lb));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
mex_command_solve(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
  assert(nlhs == 2);
  assert(nrhs == 3);

  const mxArray* mex_x0      = prhs[0];
  const mxArray* mex_funcs   = prhs[1];
  const mxArray* mex_options = prhs[2];

  Instance instance;

  SLEQP_CALL(create_instance(&instance, mex_x0, mex_funcs, mex_options));

  const int num_vars = sleqp_problem_num_vars(instance.problem);

  SleqpSolver* solver = instance.solver;

  SLEQP_CALL(sleqp_solver_solve(solver, SLEQP_NONE, SLEQP_NONE));

  SLEQP_CALL(
    create_mex_output(instance.problem, instance.solver, &plhs[0], &plhs[1]));

  SLEQP_CALL(free_instance(&instance));

  return SLEQP_OKAY;
}

