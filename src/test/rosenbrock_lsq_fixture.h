#ifndef ROSENBROCK_LSQ_FIXTURE_H
#define ROSENBROCK_LSQ_FIXTURE_H

#include "cmp.h"
#include "func.h"
#include "mem.h"
#include "sparse/sparse_vec.h"

#include "test_common.h"

#ifdef __cplusplus
extern "C" {
#endif

  extern SleqpFunc* rosenbrock_lsq_func;

  void rosenbrock_lsq_setup();

  void rosenbrock_lsq_teardown();

#ifdef __cplusplus
}
#endif

#endif /* ROSENBROCK_LSQ_FIXTURE_H */
