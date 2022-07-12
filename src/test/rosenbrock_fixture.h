#ifndef ROSENBROCK_FIXTURE_H
#define ROSENBROCK_FIXTURE_H

#include "func.h"
#include "sparse/vec.h"

#include "test_common.h"

extern const int rosenbrock_num_vars;
extern const int rosenbrock_num_cons;

extern SleqpFunc* rosenbrock_func;

extern SleqpVec* rosenbrock_var_lb;
extern SleqpVec* rosenbrock_var_ub;
extern SleqpVec* rosenbrock_cons_lb;
extern SleqpVec* rosenbrock_cons_ub;
extern SleqpVec* rosenbrock_initial;
extern SleqpVec* rosenbrock_optimum;

SLEQP_RETCODE
rosenbrock_set(SleqpFunc* func,
               SleqpVec* x,
               SLEQP_VALUE_REASON reason,
               bool* reject,
               void* func_data);

SLEQP_RETCODE
rosenbrock_obj_val(SleqpFunc* func, double* obj_val, void* func_data);

SLEQP_RETCODE
rosenbrock_obj_grad(SleqpFunc* func, SleqpVec* obj_grad, void* func_data);

SLEQP_RETCODE
rosenbrock_hess_prod(SleqpFunc* func,
                     const double* obj_dual,
                     const SleqpVec* direction,
                     const SleqpVec* cons_duals,
                     SleqpVec* product,
                     void* func_data);

void
rosenbrock_create(SleqpFunc** fstar,
                  SleqpVec** var_lbstar,
                  SleqpVec** var_ubstar,
                  SleqpVec** cons_lbstar,
                  SleqpVec** cons_ubstar,
                  SleqpVec** init_star,
                  SleqpVec** opt_star);

void
rosenbrock_setup();

void
rosenbrock_teardown();

#endif /* ROSENBROCK_FIXTURE_H */
