#ifndef LOG_ROSENBROCK_FIXTURE_H
#define LOG_ROSENBROCK_FIXTURE_H

#include "cmp.h"
#include "func.h"
#include "mem.h"
#include "sparse/sparse_vec.h"

#include "test_common.h"

extern const int log_rosenbrock_num_variables;
extern const int log_rosenbrock_num_constraints;

extern SleqpFunc* log_rosenbrock_func;

extern SleqpSparseVec* log_rosenbrock_var_lb;
extern SleqpSparseVec* log_rosenbrock_var_ub;
extern SleqpSparseVec* log_rosenbrock_cons_lb;
extern SleqpSparseVec* log_rosenbrock_cons_ub;
extern SleqpSparseVec* log_rosenbrock_initial;
extern SleqpSparseVec* log_rosenbrock_optimal;

void log_rosenbrock_setup();

void log_rosenbrock_teardown();

#endif /* LOG_ROSENBROCK_FIXTURE_H */
