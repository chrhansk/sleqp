#ifndef LINEAR_LSQ_H
#define LINEAR_LSQ_H

#include "func.h"
#include "settings.h"
#include "sparse/vec.h"

#include "test_common.h"

extern const int linear_lsq_num_variables;
extern const int linear_lsq_num_constraints;
extern const int linear_lsq_num_residuals;

extern SleqpMat* linear_lsq_matrix;
extern SleqpVec* linear_lsq_rhs;

extern SleqpSettings* linear_lsq_settings;
extern SleqpFunc* linear_lsq_func;

extern SleqpVec* linear_lsq_var_lb;
extern SleqpVec* linear_lsq_var_ub;
extern SleqpVec* linear_lsq_cons_lb;
extern SleqpVec* linear_lsq_cons_ub;
extern SleqpVec* linear_lsq_initial;
extern SleqpVec* linear_lsq_optimal;

void
linear_lsq_setup();

void
linear_lsq_teardown();

#endif /* LINEAR_LSQ_H */
