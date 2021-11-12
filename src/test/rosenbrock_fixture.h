#ifndef ROSENBROCK_FIXTURE_H
#define ROSENBROCK_FIXTURE_H

#include "func.h"
#include "sparse/sparse_vec.h"

#include "test_common.h"

extern const int rosenbrock_num_variables;
extern const int rosenbrock_num_constraints;

extern SleqpFunc* rosenbrock_func;

extern SleqpSparseVec* rosenbrock_var_lb;
extern SleqpSparseVec* rosenbrock_var_ub;
extern SleqpSparseVec* rosenbrock_cons_lb;
extern SleqpSparseVec* rosenbrock_cons_ub;
extern SleqpSparseVec* rosenbrock_initial;
extern SleqpSparseVec* rosenbrock_optimal;

SLEQP_RETCODE
rosenbrock_set(SleqpFunc* func,
               SleqpSparseVec* x,
               SLEQP_VALUE_REASON reason,
               bool* reject,
               int* func_grad_nnz,
               int* cons_val_nnz,
               int* cons_jac_nnz,
               void* func_data);

SLEQP_RETCODE
rosenbrock_val(SleqpFunc* func, double* func_val, void* func_data);

SLEQP_RETCODE
rosenbrock_grad(SleqpFunc* func, SleqpSparseVec* func_grad, void* func_data);

SLEQP_RETCODE
rosenbrock_hess_prod(SleqpFunc* func,
                     const double* func_dual,
                     const SleqpSparseVec* direction,
                     const SleqpSparseVec* cons_duals,
                     SleqpSparseVec* product,
                     void* func_data);

void
rosenbrock_setup();

void
rosenbrock_teardown();

#endif /* ROSENBROCK_FIXTURE_H */
