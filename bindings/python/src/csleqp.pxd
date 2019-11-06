#cython: language_level=3

cimport libc.stdio
cimport libc.time

cdef extern from "sleqp.h":

  cdef int SLEQP_VERSION_MAJOR
  cdef int SLEQP_VERSION_MINOR
  cdef int SLEQP_VERSION_PATCH

  cdef char* SLEQP_GIT_BRANCH
  cdef char* SLEQP_GIT_COMMIT_HASH

  ctypedef enum SLEQP_RETCODE:
    SLEQP_OKAY,
    SLEQP_NOMEM,
    SLEQP_ILLEGAL_ARGUMENT,
    SLEQP_INVALID_DERIV,
    SLEQP_INTERNAL_ERROR

  ctypedef enum SLEQP_STATUS:
    SLEQP_OPTIMAL,
    SLEQP_FEASIBLE,
    SLEQP_INFEASIBLE,
    SLEQP_INVALID

  ctypedef enum SLEQP_DERIV_CHECK:
    SLEQP_DERIV_CHECK_SKIP,
    SLEQP_DERIV_CHECK_FIRST,
    SLEQP_DERIV_CHECK_SEC,
    SLEQP_DERIV_CHECK_BOTH

  ctypedef enum SLEQP_HESSIAN_EVAL:
    SLEQP_HESSIAN_EVAL_EXACT,
    SLEQP_HESSIAN_EVAL_SR1,
    SLEQP_HESSIAN_EVAL_SIMPLE_BFGS,
    SLEQP_HESSIAN_EVAL_DAMPED_BFGS

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

  ctypedef struct SleqpScalingData:
    pass

  ctypedef struct SleqpParams:
    pass

  ctypedef struct SleqpOptions:
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

  ctypedef enum SLEQP_LOG_LEVEL:
    SLEQP_LOG_ERROR = 0,
    SLEQP_LOG_WARN = 1,
    SLEQP_LOG_INFO = 2,
    SLEQP_LOG_DEBUG = 3

  # Sparse vectors
  SLEQP_RETCODE sleqp_sparse_vector_create(SleqpSparseVec** vec,
                                           int dim,
                                           int nnz_max)

  SLEQP_RETCODE sleqp_sparse_vector_fprintf(SleqpSparseVec* vec,
                                            libc.stdio.FILE* output)

  SLEQP_RETCODE sleqp_sparse_vector_free(SleqpSparseVec** vec)

  SLEQP_RETCODE sleqp_sparse_vector_clear(SleqpSparseVec* vec)

  SLEQP_RETCODE sleqp_sparse_vector_reserve(SleqpSparseVec* vec, int nnz)

  SLEQP_RETCODE sleqp_sparse_vector_resize(SleqpSparseVec* vec,
                                           int dim)

  SLEQP_RETCODE sleqp_sparse_vector_push(SleqpSparseVec* vec,
                                         int idx,
                                         double value)

  # Sparse matrices
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

  SLEQP_RETCODE sleqp_sparse_matrix_push_column(SleqpSparseMatrix* matrix,
                                                int col)

  SLEQP_RETCODE sleqp_sparse_matrix_free(SleqpSparseMatrix** matrix)

  # Functions
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


  # Scaling
  SLEQP_RETCODE sleqp_scaling_create(SleqpScalingData** scaling,
                                     SleqpProblem* problem,
                                     SleqpParams* params);

  SLEQP_RETCODE sleqp_scaling_set_func_weight(SleqpScalingData* scaling,
                                              int weight);

  SLEQP_RETCODE sleqp_scaling_set_var_weight(SleqpScalingData* scaling,
                                             int index,
                                             int weight);

  SLEQP_RETCODE sleqp_scaling_set_cons_weight(SleqpScalingData* scaling,
                                              int index,
                                              int weight);

  SLEQP_RETCODE sleqp_func_scaling_from_gradient(SleqpScalingData* scaling,
                                                 SleqpSparseVec* gradient);

  SLEQP_RETCODE sleqp_scaling_from_cons_jac(SleqpScalingData* scaling,
                                            SleqpSparseMatrix* cons_jac);

  SLEQP_RETCODE sleqp_scaling_free(SleqpScalingData** scaling);

  # Solver
  SLEQP_RETCODE sleqp_solver_create(SleqpSolver** star,
                                    SleqpProblem* problem,
                                    SleqpParams* params,
                                    SleqpOptions* options,
                                    SleqpSparseVec* x,
                                    SleqpScalingData* scaling)

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

  # Parameters

  SLEQP_RETCODE sleqp_params_create(SleqpParams** star)

  double sleqp_params_get_zero_eps(SleqpParams* params)

  double sleqp_params_get_eps(SleqpParams* params)

  double sleqp_params_get_deriv_perturbation(SleqpParams* params)

  double sleqp_params_get_deriv_tolerance(SleqpParams* params)

  double sleqp_params_get_cauchy_tau(SleqpParams* params)
  double sleqp_params_get_cauchy_eta(SleqpParams* params)

  double sleqp_params_get_linesearch_tau(SleqpParams* params)
  double sleqp_params_get_linesearch_eta(SleqpParams* params)
  double sleqp_params_get_linesearch_cutoff(SleqpParams* params)

  double sleqp_params_get_optimality_tolerance(SleqpParams* params)

  double sleqp_params_get_accepted_reduction(SleqpParams* params)

  SLEQP_RETCODE sleqp_params_set_zero_eps(SleqpParams* params, double value)

  SLEQP_RETCODE sleqp_params_set_eps(SleqpParams* params, double value)

  SLEQP_RETCODE sleqp_params_set_deriv_perturbation(SleqpParams* params, double value)

  SLEQP_RETCODE sleqp_params_set_deriv_tolerance(SleqpParams* params, double value)

  SLEQP_RETCODE sleqp_params_set_cauchy_tau(SleqpParams* params, double value)
  SLEQP_RETCODE sleqp_params_set_cauchy_eta(SleqpParams* params, double value)

  SLEQP_RETCODE sleqp_params_set_linesearch_tau(SleqpParams* params, double value)
  SLEQP_RETCODE sleqp_params_set_linesearch_eta(SleqpParams* params, double value)
  SLEQP_RETCODE sleqp_params_set_linesearch_cutoff(SleqpParams* params, double value)

  SLEQP_RETCODE sleqp_params_set_optimality_tolerance(SleqpParams* params, double value)

  SLEQP_RETCODE sleqp_params_set_accepted_reduction(SleqpParams* params, double value)

  SLEQP_RETCODE sleqp_params_free(SleqpParams** star)

  # Options

  SLEQP_RETCODE sleqp_options_create(SleqpOptions** star)

  SLEQP_DERIV_CHECK sleqp_options_get_deriv_check(const SleqpOptions* options)

  SLEQP_HESSIAN_EVAL sleqp_options_get_hessian_eval(const SleqpOptions* options)

  SLEQP_RETCODE sleqp_options_set_deriv_check(SleqpOptions* options,
                                              SLEQP_DERIV_CHECK value)

  SLEQP_RETCODE sleqp_options_set_hessian_eval(SleqpOptions* options,
                                               SLEQP_HESSIAN_EVAL value)

  SLEQP_RETCODE sleqp_options_free(SleqpOptions** star)

  # Logging

  void sleqp_log_set_level(SLEQP_LOG_LEVEL value)

  ctypedef void (*SLEQP_LOG_HANDLER)(SLEQP_LOG_LEVEL level,
                                     libc.time.time_t time,
                                     const char* message)

  void sleqp_log_set_handler(SLEQP_LOG_HANDLER handler)

  SLEQP_LOG_LEVEL sleqp_log_level()

  # Numerics
  double sleqp_infinity()
