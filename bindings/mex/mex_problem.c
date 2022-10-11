#include "mex_problem.h"

#include <assert.h>

#include "mex_dyn.h"
#include "mex_error.h"
#include "mex_fields.h"
#include "mex_func.h"
#include "mex_lsq.h"

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
  MEX_EXPECT_STRUCT(options);

  const mxArray* cons_lb = mxGetField(options, 0, MEX_INPUT_CONS_LB);

  if (cons_lb)
  {
    return array_or_cell_size(cons_lb);
  }

  return 0;
}

SLEQP_RETCODE
mex_create_vec_from_array(SleqpVec** star, const mxArray* array)
{
  const int dimension = array_or_cell_size(array);

  SLEQP_CALL(sleqp_vec_create_full(star, dimension));

  SleqpVec* vec = *star;

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
        SLEQP_CALL(sleqp_vec_push(vec, k + offset, ptr[k]));
      }

      offset += cur_num_elements;
    }
  }
  else
  {
    const double* ptr = mxGetPr(array);

    for (int k = 0; k < dimension; ++k)
    {
      SLEQP_CALL(sleqp_vec_push(vec, k, ptr[k]));
    }
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
create_vec_from_array_default(SleqpVec** star,
                              const mxArray* array,
                              const int dimension,
                              double value)
{
  if (array)
  {
    SLEQP_CALL(mex_create_vec_from_array(star, array));
  }
  else
  {
    SLEQP_CALL(sleqp_vec_create_full(star, dimension));
    SLEQP_CALL(sleqp_vec_fill(*star, value));
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
create_bounds_from_options(SleqpVec** lb_star,
                           SleqpVec** ub_star,
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
create_var_bounds_from_options(SleqpVec** var_lb_star,
                               SleqpVec** var_ub_star,
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
create_cons_bounds_from_options(SleqpVec** cons_lb_star,
                                SleqpVec** cons_ub_star,
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

SLEQP_RETCODE
mex_problem_create(SleqpProblem** star,
                   SleqpParams* params,
                   SleqpOptions* options,
                   SLEQP_FUNC_TYPE func_type,
                   const mxArray* mex_x0,
                   const mxArray* mex_callbacks,
                   const mxArray* mex_options)
{
  const int num_vars = num_vars_from_solution(mex_x0);
  const int num_cons = num_cons_from_options(mex_options);

  SleqpVec* var_lb;
  SleqpVec* var_ub;

  SleqpVec* cons_lb;
  SleqpVec* cons_ub;

  SLEQP_CALL(
    create_var_bounds_from_options(&var_lb, &var_ub, num_vars, mex_options));

  SLEQP_CALL(
    create_cons_bounds_from_options(&cons_lb, &cons_ub, num_cons, mex_options));

  SleqpFunc* func;

  const SLEQP_HESS_EVAL hess_eval
    = sleqp_options_enum_value(options, SLEQP_OPTION_ENUM_HESS_EVAL);

  const bool with_hess = (hess_eval == SLEQP_HESS_EVAL_EXACT);

  switch (func_type)
  {
  case SLEQP_FUNC_TYPE_REGULAR:
    SLEQP_CALL(mex_func_create(&func,
                               mex_callbacks,
                               params,
                               num_vars,
                               num_cons,
                               with_hess));
    break;
  case SLEQP_FUNC_TYPE_LSQ:
    SLEQP_CALL(mex_lsq_func_create(&func,
                                   mex_x0,
                                   mex_callbacks,
                                   params,
                                   num_vars,
                                   num_cons));
    break;
  case SLEQP_FUNC_TYPE_DYNAMIC:
    SLEQP_CALL(mex_dyn_func_create(&func,
                                   mex_x0,
                                   mex_callbacks,
                                   params,
                                   num_vars,
                                   num_cons,
                                   with_hess));
    break;
  }

  SLEQP_CALL(sleqp_problem_create_simple(star,
                                         func,
                                         params,
                                         var_lb,
                                         var_ub,
                                         cons_lb,
                                         cons_ub));

  SLEQP_CALL(sleqp_vec_free(&cons_ub));
  SLEQP_CALL(sleqp_vec_free(&cons_lb));

  SLEQP_CALL(sleqp_vec_free(&var_ub));
  SLEQP_CALL(sleqp_vec_free(&var_lb));

  return SLEQP_OKAY;
}
