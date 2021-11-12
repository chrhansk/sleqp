#ifndef WACHBIEG_FIXTURE_H
#define WACHBIEG_FIXTURE_H

#include "func.h"

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

void
wachbieg_setup();

void
wachbieg_teardown();

#endif /* WACHBIEG_FIXTURE_H */
