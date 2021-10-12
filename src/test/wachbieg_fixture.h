#ifndef WACHBIEG_FIXTURE_H
#define WACHBIEG_FIXTURE_H

#include "func.h"

#ifdef __cplusplus
extern "C" {
#endif

  extern const int wachbieg_num_variables;
  extern const int wachbieg_num_constraints;

  extern SleqpFunc* wachbieg_func;

  extern SleqpSparseVec* wachbieg_var_lb;
  extern SleqpSparseVec* wachbieg_var_ub;
  extern SleqpSparseVec* wachbieg_cons_lb;
  extern SleqpSparseVec* wachbieg_cons_ub;
  extern SleqpSparseVec* wachbieg_initial;
  extern SleqpSparseVec* wachbieg_optimal;

  extern SleqpFunc* wachbieg_func;

  void wachbieg_setup();

  void wachbieg_teardown();

#ifdef __cplusplus
}
#endif

#endif /* WACHBIEG_FIXTURE_H */
