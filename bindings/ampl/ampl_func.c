#include "ampl_func.h"

#include <assert.h>

#include "ampl_util.h"

#define SLEQP_AMPL_IS_ERROR(errorptr) (errorptr && (*errorptr != 0))

#define SLEQP_AMPL_ERROR_CHECK(errorptr)                                       \
  do                                                                           \
  {                                                                            \
    if (SLEQP_AMPL_IS_ERROR(errorptr))                                         \
    {                                                                          \
      sleqp_raise(SLEQP_INTERNAL_ERROR,                                        \
                  "Error during evaluation. "                                  \
                  "Run with \"halt_on_ampl_error yes\" to see details.");      \
    }                                                                          \
  } while (false)

typedef struct AmplFuncData
{
  SleqpAmplData* ampl_data;
  double zero_eps;

  double* x;

  double obj_val;
  double* obj_grad;

  bool cached;

  enum
  {
    NONE     = 0,
    OBJ_VAL  = (1 << 0),
    CONS_VAL = (1 << 1),
  } evaluated;

  double* direction;
  double* multipliers;
  double* hessian_product;

  bool inverted_obj;
  double offset;

  fint error;
  fint* nerror;

} AmplFuncData;

static SLEQP_RETCODE
ampl_func_data_create(AmplFuncData** star,
                      SleqpAmplData* ampl_data,
                      double zero_eps,
                      bool halt_on_error)
{
  ASL* asl = ampl_data->asl;
  SLEQP_CALL(sleqp_malloc(star));

  AmplFuncData* data = *star;

  *data = (AmplFuncData){0};

  const int num_variables = ampl_data->num_variables;
  const int num_general   = ampl_data->num_constraints - ampl_data->num_linear;

  data->ampl_data = ampl_data;
  data->zero_eps  = zero_eps;

  data->inverted_obj = sleqp_ampl_max_problem(asl);

  data->cached    = false;
  data->evaluated = NONE;

  if (halt_on_error)
  {
    data->nerror = NULL;
  }
  else
  {
    data->nerror = &data->error;
  }

  if (n_obj > 0)
  {
    data->offset = objconst(0);
  }
  else
  {
    data->offset = 0.;
  }

  SLEQP_CALL(sleqp_alloc_array(&data->x, num_variables));
  SLEQP_CALL(sleqp_alloc_array(&data->multipliers, num_general));

  SLEQP_CALL(sleqp_alloc_array(&data->obj_grad, num_variables));

  SLEQP_CALL(sleqp_alloc_array(&data->direction, num_variables));
  SLEQP_CALL(sleqp_alloc_array(&data->hessian_product, num_variables));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
ampl_func_data_free(void* func_data)
{
  AmplFuncData* data  = (AmplFuncData*)func_data;
  AmplFuncData** star = &data;

  sleqp_free(&data->hessian_product);
  sleqp_free(&data->direction);

  sleqp_free(&data->obj_grad);
  sleqp_free(&data->multipliers);
  sleqp_free(&data->x);

  sleqp_free(star);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
eval_obj(AmplFuncData* data, bool* reject)
{
  ASL* asl = data->ampl_data->asl;

  data->obj_val = objval(0, data->x, data->nerror);

  *reject = SLEQP_AMPL_IS_ERROR(data->nerror);

  data->evaluated |= OBJ_VAL;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
eval_cons(AmplFuncData* data, bool* reject)
{
  SleqpAmplData* ampl_data = data->ampl_data;
  ASL* asl                 = data->ampl_data->asl;

  conval(data->x, ampl_data->cons_val, data->nerror);

  *reject = SLEQP_AMPL_IS_ERROR(data->nerror);

  data->evaluated |= CONS_VAL;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
set_eval(AmplFuncData* data, bool* reject)
{
  SleqpAmplData* ampl_data = data->ampl_data;

  *reject = false;

  SLEQP_CALL(eval_obj(data, reject));

  if (*reject)
  {
    return SLEQP_OKAY;
  }

  SLEQP_CALL(eval_cons(data, reject));

  if (*reject)
  {
    return SLEQP_OKAY;
  }

  if (!sleqp_is_finite(data->obj_val))
  {
    *reject = true;
    return SLEQP_OKAY;
  }

  const int num_constraints = data->ampl_data->num_constraints;

  for (int i = 0; i < num_constraints; ++i)
  {
    if (!sleqp_is_finite(ampl_data->cons_val[i]))
    {
      *reject = true;
      return SLEQP_OKAY;
    }
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
ampl_func_set(SleqpFunc* func,
              SleqpVec* x,
              SLEQP_VALUE_REASON reason,
              bool* reject,
              void* func_data)
{
  AmplFuncData* data = (AmplFuncData*)func_data;

  SLEQP_CALL(sleqp_vec_to_raw(x, data->x));

  data->evaluated = NONE;

  // pre-check for possible rejects, cache
  // objective and constraint values
  if ((reason == SLEQP_VALUE_REASON_TRYING_ITERATE)
      || (reason == SLEQP_VALUE_REASON_TRYING_SOC_ITERATE))
  {
    SLEQP_CALL(set_eval(data, reject));
    data->cached = true;
  }
  else
  {
    data->cached = false;
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE
ampl_func_nonzeros(SleqpFunc* func,
                   int* obj_grad_nnz,
                   int* cons_val_nnz,
                   int* cons_jac_nnz,
                   int* hess_prod_nnz,
                   void* func_data)
{
  AmplFuncData* data = (AmplFuncData*)func_data;

  *cons_jac_nnz = data->ampl_data->jac_nnz;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
ampl_obj_val(SleqpFunc* func, double* obj_val, void* func_data)
{
  AmplFuncData* data = (AmplFuncData*)func_data;
  ASL* asl           = data->ampl_data->asl;

  if (data->cached)
  {
    *obj_val = data->obj_val;
  }
  else
  {
    *obj_val = objval(0, data->x, data->nerror);

    data->evaluated |= OBJ_VAL;

    SLEQP_AMPL_ERROR_CHECK(data->nerror);
  }

  *obj_val += data->offset;

  if (data->inverted_obj)
  {
    *obj_val *= -1.;
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
ampl_obj_grad(SleqpFunc* func, SleqpVec* obj_grad, void* func_data)
{
  AmplFuncData* data = (AmplFuncData*)func_data;
  ASL* asl           = data->ampl_data->asl;

  objgrd(0, data->x, data->obj_grad, data->nerror);

  SLEQP_AMPL_ERROR_CHECK(data->nerror);

  SLEQP_CALL(sleqp_vec_set_from_raw(obj_grad,
                                    data->obj_grad,
                                    data->ampl_data->num_variables,
                                    data->zero_eps));

  if (data->inverted_obj)
  {
    SLEQP_CALL(sleqp_vec_scale(obj_grad, -1.));
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
ampl_cons_val(SleqpFunc* func, SleqpVec* cons_val, void* func_data)
{
  AmplFuncData* data       = (AmplFuncData*)func_data;
  SleqpAmplData* ampl_data = data->ampl_data;
  ASL* asl                 = data->ampl_data->asl;

  if (!data->cached)
  {
    conval(data->x, ampl_data->cons_val, data->nerror);

    data->evaluated |= CONS_VAL;

    SLEQP_AMPL_ERROR_CHECK(data->nerror);
  }

  const int num_general = ampl_data->num_constraints - ampl_data->num_linear;

  SLEQP_CALL(sleqp_vec_set_from_raw(cons_val,
                                    ampl_data->cons_val,
                                    num_general,
                                    data->zero_eps));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
ampl_cons_jac(SleqpFunc* func, SleqpMat* cons_jac, void* func_data)
{
  AmplFuncData* data       = (AmplFuncData*)func_data;
  ASL* asl                 = data->ampl_data->asl;
  SleqpAmplData* ampl_data = data->ampl_data;

  jacval(data->x, ampl_data->jac_vals, data->nerror);

  SLEQP_AMPL_ERROR_CHECK(data->nerror);

  int next_col = 0;

  SLEQP_CALL(sleqp_mat_reserve(cons_jac, ampl_data->jac_nnz));

  const int num_general = ampl_data->num_constraints - ampl_data->num_linear;

  for (int i = 0; i < ampl_data->jac_nnz; ++i)
  {
    const int row    = ampl_data->jac_rows[i];
    const int col    = ampl_data->jac_cols[i];
    const double val = ampl_data->jac_vals[i];

    while (col >= next_col)
    {
      SLEQP_CALL(sleqp_mat_push_col(cons_jac, next_col++));
    }

    if (row < num_general)
    {
      SLEQP_CALL(sleqp_mat_push(cons_jac, row, col, val));
    }
  }

  const int num_cols = data->ampl_data->num_variables;

  while (num_cols > next_col)
  {
    SLEQP_CALL(sleqp_mat_push_col(cons_jac, next_col++));
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
ensure_eval(AmplFuncData* data)
{
  if (!(data->evaluated & OBJ_VAL))
  {
    bool reject = false;

    SLEQP_CALL(eval_obj(data, &reject));

    assert(!reject);
  }

  if (!(data->evaluated & CONS_VAL))
  {
    bool reject = false;

    SLEQP_CALL(eval_cons(data, &reject));

    assert(!reject);
  }

  data->cached = true;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
ampl_func_hess_product(SleqpFunc* func,
                       const SleqpVec* direction,
                       const SleqpVec* cons_duals,
                       SleqpVec* product,
                       void* func_data)
{
  AmplFuncData* data = (AmplFuncData*)func_data;
  ASL* asl           = data->ampl_data->asl;

  // Need to ensure that objective and constraints
  // have been evaluated at the current primal point
  SLEQP_CALL(ensure_eval(data));

  SLEQP_CALL(sleqp_vec_to_raw(direction, data->direction));
  SLEQP_CALL(sleqp_vec_to_raw(cons_duals, data->multipliers));

  if (data->inverted_obj)
  {
    const int num_cons = sleqp_func_num_cons(func);

    for (int i = 0; i < num_cons; ++i)
    {
      data->multipliers[i] *= -1.;
    }
  }

  double one = 1.;

  hvcomp(data->hessian_product,
         data->direction,
         0,
         &one,
         data->multipliers);

  SLEQP_CALL(sleqp_vec_set_from_raw(product,
                                    data->hessian_product,
                                    data->ampl_data->num_variables,
                                    data->zero_eps));

  if (data->inverted_obj)
  {
    SLEQP_CALL(sleqp_vec_scale(product, -1.));
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_ampl_func_create(SleqpFunc** star,
                       SleqpAmplData* ampl_data,
                       SleqpSettings* settings,
                       bool halt_on_error)
{
  AmplFuncData* data;

  const int num_variables = ampl_data->num_variables;
  const int num_general   = ampl_data->num_constraints - ampl_data->num_linear;

  const double zero_eps = sleqp_settings_real_value(settings, SLEQP_SETTINGS_REAL_ZERO_EPS);

  SLEQP_CALL(ampl_func_data_create(&data, ampl_data, zero_eps, halt_on_error));

  SleqpFuncCallbacks callbacks
    = {.set_value = ampl_func_set,
       .nonzeros  = ampl_func_nonzeros,
       .obj_val   = ampl_obj_val,
       .obj_grad  = ampl_obj_grad,
       .cons_val  = ampl_data->is_constrained ? ampl_cons_val : NULL,
       .cons_jac  = ampl_data->is_constrained ? ampl_cons_jac : NULL,
       .hess_prod = ampl_func_hess_product,
       .func_free = ampl_func_data_free};

  SLEQP_CALL(
    sleqp_func_create(star, &callbacks, num_variables, num_general, data));

  return SLEQP_OKAY;
}
