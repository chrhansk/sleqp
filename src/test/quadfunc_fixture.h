#ifndef QUADFUNC_FIXTURE_H
#define QUADFUNC_FIXTURE_H

#include "cmp.h"
#include "func.h"
#include "mem.h"
#include "sparse/sparse_vec.h"

#include "test_common.h"

extern SleqpFunc* quadfunc;

extern SleqpSparseVec* quadfunc_var_lb;
extern SleqpSparseVec* quadfunc_var_ub;
extern SleqpSparseVec* quadfunc_cons_lb;
extern SleqpSparseVec* quadfunc_cons_ub;
extern SleqpSparseVec* quadfunc_x;

void
quadfunc_setup();

void
quadfunc_teardown();

#endif /* QUADFUNC_FIXTURE_H */
