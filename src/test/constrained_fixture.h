#ifndef CONSTRAINED_FIXTURE_H
#define CONSTRAINED_FIXTURE_H

#include "test_common.h"

#include "cmp.h"
#include "func.h"

extern int constrained_num_variables;
extern int constrained_num_constraints;

extern SleqpSparseVec* constrained_var_lb;
extern SleqpSparseVec* constrained_var_ub;
extern SleqpSparseVec* constrained_cons_lb;
extern SleqpSparseVec* constrained_cons_ub;
extern SleqpSparseVec* constrained_initial;
extern SleqpSparseVec* constrained_optimal;

extern SleqpFunc* constrained_func;

void constrained_setup();

void constrained_teardown();


#endif /* CONSTRAINED_FIXTURE_H */
