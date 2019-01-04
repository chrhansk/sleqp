#!/usr/bin/python
#cython: language_level=3

from libc.stdio cimport FILE

cdef extern from "sleqp.h":

  cdef int sleqp_VERSION_MAJOR
  cdef int sleqp_VERSION_MINOR
  cdef int sleqp_VERSION_PATCH

  cdef char* sleqp_GIT_BRANCH
  cdef char* sleqp_GIT_COMMIT_HASH

  cdef enum SLEQP_Retcode:
    SLEQP_OKAY,
    SLEQP_NOMEM,
    SLEQP_ILLEGAL_ARGUMENT,
    SLEQP_INVALID_DERIV,
    SLEQP_INTERNAL_ERROR

  ctypedef SLEQP_Retcode SLEQP_RETCODE

  cdef enum SLEQP_Status:
    SLEQP_OPTIMAL,
    SLEQP_FEASIBLE,
    SLEQP_INFEASIBLE,
    SLEQP_INVALID

  ctypedef SLEQP_Status SLEQP_STATUS

  ctypedef struct SleqpSparseVec:
    double* data
    int* indices

    int dim
    int nnz
    int nnz_max

  ctypedef struct SleqpSparseMatrix:
    int num_rows
    int num_cols

    int nnz
    int nnz_max

    double* data
    int* cols
    int* rows

  ctypedef struct SleqpProblem:
    int num_variables
    int num_constraints

  ctypedef struct SleqpSolver:
    pass

  ctypedef struct SleqpParams:
    pass

  ctypedef struct SleqpActiveSet:
    pass

  ctypedef struct SleqpIterate:
    SleqpSparseVec* x
    double func_val
    SleqpSparseVec* func_grad
    SleqpSparseVec* cons_val
    SleqpSparseMatrix* cons_jac
    SleqpActiveSet* active_set
    SleqpSparseVec* cons_dual
    SleqpSparseVec* vars_dual

  cdef struct SleqpFunc:
    pass

  SLEQP_RETCODE sleqp_sparse_vector_create(SleqpSparseVec** vec,
                                           int dim,
                                           int nnz_max)

  SLEQP_RETCODE sleqp_sparse_vector_fprintf(SleqpSparseVec* vec,
                                            FILE* output)

  SLEQP_RETCODE sleqp_sparse_vector_free(SleqpSparseVec** vec)

  SLEQP_RETCODE sleqp_sparse_vector_clear(SleqpSparseVec* vec)

  SLEQP_RETCODE sleqp_sparse_vector_reserve(SleqpSparseVec* vec, int nnz)

  SLEQP_RETCODE sleqp_sparse_vector_resize(SleqpSparseVec* vec,
                                           int dim)

  SLEQP_RETCODE sleqp_sparse_vector_push(SleqpSparseVec* vec,
                                         int idx,
                                         double value)

  SLEQP_RETCODE sleqp_sparse_matrix_create(SleqpSparseMatrix** matrix,
                                           int num_rows,
                                           int num_cols,
                                           int nnz_max)

  SLEQP_RETCODE sleqp_sparse_matrix_reserve(SleqpSparseMatrix* matrix,
                                            int nnz)

  SLEQP_RETCODE sleqp_sparse_matrix_resize(SleqpSparseMatrix* matrix,
                                           int num_rows,
                                           int num_cols)

  SLEQP_RETCODE sleqp_sparse_matrix_push(SleqpSparseMatrix* matrix,
                                         int row,
                                         int col,
                                         double value)

  SLEQP_RETCODE sleqp_sparse_matrix_add_column(SleqpSparseMatrix* matrix,
                                               int col)

  SLEQP_RETCODE sleqp_sparse_matrix_free(SleqpSparseMatrix** matrix)

  ctypedef SLEQP_RETCODE (*SLEQP_FUNC_SET)(SleqpSparseVec* x,
                                           int num_variables,
                                           int* func_grad_nnz,
                                           int* cons_val_nnz,
                                           int* cons_jac_nnz,
                                           void* func_data)

  ctypedef SLEQP_RETCODE (*SLEQP_FUNC_EVAL)(int num_variables,
                                            SleqpSparseVec* cons_indices,
                                            double* func_val,
                                            SleqpSparseVec* func_grad,
                                            SleqpSparseVec* cons_val,
                                            SleqpSparseMatrix* cons_jac,
                                            void* func_data)

  ctypedef SLEQP_RETCODE (*SLEQP_HESS_PRODUCT)(int num_variables,
                                               double* func_dual,
                                               SleqpSparseVec* direction,
                                               SleqpSparseVec* cons_duals,
                                               SleqpSparseVec* product,
                                               void* func_data)

  SLEQP_RETCODE sleqp_func_create(SleqpFunc** fstar,
                                  SLEQP_FUNC_SET setx,
                                  SLEQP_FUNC_EVAL eval,
                                  SLEQP_HESS_PRODUCT eval_hess_prod,
                                  int num_variables,
                                  void* func_data)

  SLEQP_RETCODE sleqp_func_free(SleqpFunc** fstar)


  SLEQP_RETCODE sleqp_solver_create(SleqpSolver** star,
                                    SleqpProblem* problem,
                                    SleqpParams* params,
                                    SleqpSparseVec* x)

  SLEQP_RETCODE sleqp_solver_solve(SleqpSolver* solver,
                                   int max_num_iterations,
                                   double time_limit)

  SLEQP_STATUS sleqp_solver_get_status(SleqpSolver* solver)

  int sleqp_solver_get_iterations(SleqpSolver* solver)

  double sleqp_solver_get_elapsed_seconds(SleqpSolver* solver)

  SLEQP_RETCODE sleqp_solver_get_solution(SleqpSolver* solver,
                                          SleqpIterate** iterate)

  SLEQP_RETCODE sleqp_solver_free(SleqpSolver** star)

  SLEQP_RETCODE sleqp_problem_create(SleqpProblem** star,
                                     SleqpFunc* func,
                                     SleqpParams* params,
                                     SleqpSparseVec* var_lb,
                                     SleqpSparseVec* var_ub,
                                     SleqpSparseVec* cons_lb,
                                     SleqpSparseVec* cons_ub)

  SLEQP_RETCODE sleqp_problem_free(SleqpProblem** star)

  SLEQP_RETCODE sleqp_params_create(SleqpParams** star)

  SLEQP_RETCODE sleqp_params_free(SleqpParams** star)

  double sleqp_infinity()
