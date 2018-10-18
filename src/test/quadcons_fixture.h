#ifndef QUADCONS_FIXTURE_H
#define QUADCONS_FIXTURE_H

#include "sleqp.h"
#include "sleqp_cmp.h"
#include "sleqp_mem.h"

#include "test_common.h"

extern SleqpFunc* quadconsfunc;

extern SleqpSparseVec* quadconsfunc_var_lb;
extern SleqpSparseVec* quadconsfunc_var_ub;
extern SleqpSparseVec* quadconsfunc_cons_lb;
extern SleqpSparseVec* quadconsfunc_cons_ub;
extern SleqpSparseVec* quadconsfunc_x;

void quadconsfunc_setup();

void quadconsfunc_teardown();

#endif /* QUADCONS_FIXTURE_H */
