#ifndef ROSENBROCK_FIXTURE_H
#define ROSENBROCK_FIXTURE_H

#include "cmp.h"
#include "func.h"
#include "mem.h"
#include "sparse/sparse_vec.h"

#include "test_common.h"

#ifdef __cplusplus
extern "C" {
#endif

  extern const int rosenbrock_num_variables;
  extern const int rosenbrock_num_constraints;

  extern SleqpFunc* rosenbrock_func;

  extern SleqpSparseVec* rosenbrock_var_lb;
  extern SleqpSparseVec* rosenbrock_var_ub;
  extern SleqpSparseVec* rosenbrock_cons_lb;
  extern SleqpSparseVec* rosenbrock_cons_ub;
  extern SleqpSparseVec* rosenbrock_initial;
  extern SleqpSparseVec* rosenbrock_optimal;

  void rosenbrock_setup();

  void rosenbrock_teardown();

#ifdef __cplusplus
}
#endif

#endif /* ROSENBROCK_FIXTURE_H */
