#ifndef LINEAR_LSQ_H
#define LINEAR_LSQ_H

#include "func.h"
#include "params.h"
#include "sparse/sparse_vec.h"

#include "test_common.h"

extern const int linear_lsq_num_variables;
extern const int linear_lsq_num_constraints;
extern const int linear_lsq_num_residuals;

extern SleqpSparseMatrix* linear_lsq_matrix;
extern SleqpSparseVec* linear_lsq_rhs;

extern SleqpParams* linear_lsq_params;
extern SleqpFunc* linear_lsq_func;

extern SleqpSparseVec* linear_lsq_var_lb;
extern SleqpSparseVec* linear_lsq_var_ub;
extern SleqpSparseVec* linear_lsq_cons_lb;
extern SleqpSparseVec* linear_lsq_cons_ub;
extern SleqpSparseVec* linear_lsq_initial;
extern SleqpSparseVec* linear_lsq_optimal;

void linear_lsq_setup();

void linear_lsq_teardown();

#endif /* LINEAR_LSQ_H */
