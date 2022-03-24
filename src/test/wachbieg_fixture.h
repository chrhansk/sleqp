#ifndef WACHBIEG_FIXTURE_H
#define WACHBIEG_FIXTURE_H

#include "func.h"

extern const int wachbieg_num_variables;
extern const int wachbieg_num_constraints;

extern SleqpFunc* wachbieg_func;

extern SleqpVec* wachbieg_var_lb;
extern SleqpVec* wachbieg_var_ub;
extern SleqpVec* wachbieg_cons_lb;
extern SleqpVec* wachbieg_cons_ub;
extern SleqpVec* wachbieg_initial;
extern SleqpVec* wachbieg_optimal;

extern SleqpFunc* wachbieg_func;

void
wachbieg_setup();

void
wachbieg_teardown();

#endif /* WACHBIEG_FIXTURE_H */
