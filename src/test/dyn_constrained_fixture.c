#include "dyn_constrained_fixture.h"

#include <stdlib.h>

#include "cmp.h"
#include "constrained_fixture.h"
#include "dyn.h"
#include "func.h"
#include "mem.h"

SleqpFunc* dyn_constrained_func = NULL;

struct
{
  double error_bound;
  double obj_weight;
  double* cons_weights;

  SleqpVec* orig_cons_val;
  SleqpVec* noise;

} func_data;

static SLEQP_RETCODE
dyn_func_set(SleqpFunc* func,
             SleqpVec* value,
             SLEQP_VALUE_REASON reason,
             bool* reject,
             void* data)
{
  SLEQP_CALL(sleqp_func_set_value(constrained_func, value, reason, reject));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
generate_noise(double* error)
{
  SLEQP_CALL(sleqp_vec_clear(func_data.noise));

  *error = 0.;

  for (int i = 0; i < constrained_num_constraints; ++i)
  {
    double noise = ((double)rand()) / ((double)RAND_MAX);

    noise = 2. * noise - 1.;

    double factor = func_data.error_bound;
    factor /= (func_data.cons_weights[i] * constrained_num_constraints);

    noise *= factor;

    SLEQP_CALL(sleqp_vec_push(func_data.noise, i, noise));

    *error += SLEQP_ABS(noise);
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
dyn_func_eval(SleqpFunc* func,
              double* obj_val,
              SleqpVec* cons_val,
              double* error,
              void* data)
{
  SLEQP_CALL(sleqp_func_obj_val(constrained_func, obj_val));

  SLEQP_CALL(sleqp_func_cons_val(constrained_func, func_data.orig_cons_val));

  generate_noise(error);

  SLEQP_CALL(
    sleqp_vec_add(func_data.orig_cons_val, func_data.noise, 0., cons_val));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
dyn_func_obj_grad(SleqpFunc* func, SleqpVec* obj_grad, void* data)
{
  SLEQP_CALL(sleqp_func_obj_grad(constrained_func, obj_grad));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
dyn_func_cons_jac(SleqpFunc* func, SleqpSparseMatrix* cons_jac, void* data)
{
  SLEQP_CALL(sleqp_func_cons_jac(constrained_func, cons_jac));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
dyn_func_hess_prod(SleqpFunc* func,
                   const double* obj_dual,
                   const SleqpVec* direction,
                   const SleqpVec* cons_duals,
                   SleqpVec* product,
                   void* data)
{
  SLEQP_CALL(sleqp_func_hess_prod(constrained_func,
                                  obj_dual,
                                  direction,
                                  cons_duals,
                                  product));
  return SLEQP_OKAY;
}

static SLEQP_RETCODE
dyn_func_set_error_bound(SleqpFunc* func, double error_bound, void* data)
{
  func_data.error_bound = error_bound;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
dyn_func_set_obj_weight(SleqpFunc* func, double obj_weight, void* data)
{
  func_data.obj_weight = obj_weight;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
dyn_func_set_cons_weights(SleqpFunc* func,
                          const double* cons_weights,
                          void* data)
{
  for (int i = 0; i < constrained_num_constraints; ++i)
  {
    func_data.cons_weights[i] = cons_weights[i];
  }

  return SLEQP_OKAY;
}

void
dyn_constrained_setup()
{
  constrained_setup();
  srand(42);

  SleqpDynFuncCallbacks callbacks
    = {.set_value        = dyn_func_set,
       .nonzeros         = NULL,
       .set_error_bound  = dyn_func_set_error_bound,
       .set_obj_weight   = dyn_func_set_obj_weight,
       .set_cons_weights = dyn_func_set_cons_weights,
       .eval             = dyn_func_eval,
       .obj_grad         = dyn_func_obj_grad,
       .cons_jac         = dyn_func_cons_jac,
       .hess_prod        = dyn_func_hess_prod,
       .func_free        = NULL};

  sleqp_alloc_array(&func_data.cons_weights, constrained_num_constraints);

  ASSERT_CALL(sleqp_vec_create_full(&func_data.orig_cons_val,
                                    constrained_num_constraints));

  ASSERT_CALL(
    sleqp_vec_create_full(&func_data.noise, constrained_num_constraints));

  ASSERT_CALL(sleqp_dyn_func_create(&dyn_constrained_func,
                                    &callbacks,
                                    constrained_num_variables,
                                    constrained_num_constraints,
                                    NULL));
}

void
dyn_constrained_teardown()
{
  ASSERT_CALL(sleqp_func_release(&dyn_constrained_func));

  sleqp_free(&func_data.cons_weights);

  ASSERT_CALL(sleqp_vec_free(&func_data.noise));
  ASSERT_CALL(sleqp_vec_free(&func_data.orig_cons_val));

  constrained_teardown();
}
