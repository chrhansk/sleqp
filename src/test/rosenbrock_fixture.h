#ifndef ROSENBROCK_FIXTURE_H
#define ROSENBROCK_FIXTURE_H

#include "sleqp.h"
#include "sleqp_cmp.h"
#include "sleqp_mem.h"

#include "test_common.h"


extern SleqpFunc* rosenbrock_func;

extern SleqpSparseVec* rosenbrock_var_lb;
extern SleqpSparseVec* rosenbrock_var_ub;
extern SleqpSparseVec* rosenbrock_cons_lb;
extern SleqpSparseVec* rosenbrock_cons_ub;
extern SleqpSparseVec* rosenbrock_x;

void rosenbrock_setup();

void rosenbrock_teardown();

#endif /* ROSENBROCK_FIXTURE_H */
