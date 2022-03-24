#ifndef QUADCONS_FIXTURE_H
#define QUADCONS_FIXTURE_H

#include "cmp.h"
#include "func.h"
#include "mem.h"
#include "sparse/vec.h"

#include "test_common.h"

extern SleqpFunc* quadconsfunc;

extern const int quadconsfunc_num_variables;
extern const int quadconsfunc_num_constraints;

extern SleqpVec* quadconsfunc_var_lb;
extern SleqpVec* quadconsfunc_var_ub;
extern SleqpVec* quadconsfunc_cons_lb;
extern SleqpVec* quadconsfunc_cons_ub;
extern SleqpVec* quadconsfunc_x;

void
quadconsfunc_setup();

void
quadconsfunc_teardown();

#endif /* QUADCONS_FIXTURE_H */
