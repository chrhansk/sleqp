#include "dyn_rosenbrock_fixture.h"

#include <stdlib.h>

#include "dyn.h"

#include "func.h"
#include "mem.h"
#include "test_common.h"

SleqpFunc* dyn_rosenbrock_func;

typedef struct
{
  double error_bound;
  double obj_weight;

} FuncData;

static SLEQP_RETCODE
dyn_rosenbrock_set(SleqpFunc* func,
                   SleqpVec* x,
                   SLEQP_VALUE_REASON reason,
                   bool* reject,
                   void* func_data)
{
  SLEQP_CALL(sleqp_func_set_value(rosenbrock_func, x, reason, reject));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
dyn_rosenbrock_eval(SleqpFunc* func,
                    double* obj_val,
                    SleqpVec* cons_val,
                    double* error,
                    void* func_data)
{
  FuncData* data = (FuncData*)func_data;

  double orig_obj_val;

  SLEQP_CALL(sleqp_func_obj_val(rosenbrock_func, &orig_obj_val));

  // uniform in [0, 1]
  double noise = ((double)rand()) / ((double)RAND_MAX);

  // uniform in [-1, 1]
  noise = 2. * noise - 1;

  const double factor = data->error_bound / data->obj_weight;

  *obj_val = (orig_obj_val + factor * noise);

  *error = SLEQP_ABS(noise) * data->error_bound;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
dyn_rosenbrock_obj_grad(SleqpFunc* func, SleqpVec* obj_grad, void* func_data)
{
  SLEQP_CALL(sleqp_func_obj_grad(rosenbrock_func, obj_grad));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
dyn_rosenbrock_hess_prod(SleqpFunc* func,
                         const SleqpVec* direction,
                         const SleqpVec* cons_duals,
                         SleqpVec* product,
                         void* func_data)
{
  SLEQP_CALL(sleqp_func_hess_prod(rosenbrock_func,
                                  direction,
                                  cons_duals,
                                  product));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
dyn_rosenbrock_set_error_bound(SleqpFunc* func,
                               double error_bound,
                               void* func_data)
{
  FuncData* data = (FuncData*)func_data;

  assert(error_bound > 0.);

  data->error_bound = error_bound;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
dyn_rosenbrock_set_obj_weight(SleqpFunc* func,
                              double obj_weight,
                              void* func_data)
{
  FuncData* data = (FuncData*)func_data;

  assert(obj_weight > 0.);

  data->obj_weight = obj_weight;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
dyn_rosenbrock_free(void* func_data)
{
  FuncData* data = (FuncData*)func_data;

  sleqp_free(&data);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
dyn_rosenbrock_data_create(FuncData** star)
{
  SLEQP_CALL(sleqp_malloc(star));

  FuncData* func_data = *star;

  *func_data = (FuncData){0};

  return SLEQP_OKAY;
}

void
dyn_rosenbrock_setup()
{
  srand(42);

  rosenbrock_setup();

  SleqpDynFuncCallbacks callbacks
    = {.set_value        = dyn_rosenbrock_set,
       .set_error_bound  = dyn_rosenbrock_set_error_bound,
       .set_obj_weight   = dyn_rosenbrock_set_obj_weight,
       .set_cons_weights = NULL,
       .eval             = dyn_rosenbrock_eval,
       .obj_grad         = dyn_rosenbrock_obj_grad,
       .hess_prod        = dyn_rosenbrock_hess_prod,
       .func_free        = dyn_rosenbrock_free};

  FuncData* func_data;

  ASSERT_CALL(dyn_rosenbrock_data_create(&func_data));

  ASSERT_CALL(sleqp_dyn_func_create(&dyn_rosenbrock_func,
                                    &callbacks,
                                    rosenbrock_num_vars,
                                    rosenbrock_num_cons,
                                    func_data));
}

void
dyn_rosenbrock_teardown()
{
  ASSERT_CALL(sleqp_func_release(&dyn_rosenbrock_func));

  rosenbrock_teardown();
}
