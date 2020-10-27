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

  ctypedef enum SLEQP_ACTIVE_STATE:
    SLEQP_INACTIVE,
    SLEQP_ACTIVE_LOWER,
    SLEQP_ACTIVE_UPPER,
    SLEQP_ACTIVE_BOTH

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

  ctypedef enum SLEQP_VALUE_REASON:
    SLEQP_VALUE_REASON_NONE,
    SLEQP_VALUE_REASON_INIT,
    SLEQP_VALUE_REASON_CHECKING_DERIV,
    SLEQP_VALUE_REASON_ACCEPTED_ITERATE,
    SLEQP_VALUE_REASON_TRYING_ITERATE,
    SLEQP_VALUE_REASON_TRYING_SOC_ITERATE,
    SLEQP_VALUE_REASON_REJECTED_ITERATE

  ctypedef enum SLEQP_DUAL_ESTIMATION_TYPE:
    SLEQP_DUAL_ESTIMATION_TYPE_LP,
    SLEQP_DUAL_ESTIMATION_TYPE_LSQ,

  ctypedef struct SleqpSparseVec:
    double* data
    int* indices

    int dim
    int nnz
    int nnz_max

  ctypedef struct SleqpSparseMatrix:
    pass

  ctypedef struct SleqpProblem:
    int num_variables
    int num_constraints

    SleqpSparseVec* var_lb
    SleqpSparseVec* var_ub

    SleqpSparseVec* cons_lb
    SleqpSparseVec* cons_ub

  ctypedef struct SleqpSolver:
    pass

  ctypedef struct SleqpScalingData:
    pass

  ctypedef struct SleqpParams:
    pass

  ctypedef struct SleqpOptions:
    pass

  ctypedef struct SleqpWorkingSet:
    pass

  ctypedef struct SleqpIterate:
    pass

  cdef struct SleqpFunc:
    pass

  cdef struct SleqpHessianStruct:
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

  SLEQP_RETCODE sleqp_sparse_vector_create_empty(SleqpSparseVec** vec,
                                                 int dim)

  SLEQP_RETCODE sleqp_sparse_vector_create_full(SleqpSparseVec** vec,
                                                int dim)

  SLEQP_RETCODE sleqp_sparse_vector_fprintf(const SleqpSparseVec* vec,
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

  int sleqp_sparse_matrix_get_num_cols(SleqpSparseMatrix* matrix)

  int sleqp_sparse_matrix_get_num_rows(SleqpSparseMatrix* matrix)

  int sleqp_sparse_matrix_get_nnz(SleqpSparseMatrix* matrix)

  int sleqp_sparse_matrix_get_nnz_max(SleqpSparseMatrix* matrix)

  double* sleqp_sparse_matrix_get_data(SleqpSparseMatrix* matrix)

  int* sleqp_sparse_matrix_get_cols(SleqpSparseMatrix* matrix)

  int* sleqp_sparse_matrix_get_rows(SleqpSparseMatrix* matrix)

  SLEQP_RETCODE sleqp_sparse_matrix_push(SleqpSparseMatrix* matrix,
                                         int row,
                                         int col,
                                         double value)

  SLEQP_RETCODE sleqp_sparse_matrix_push_column(SleqpSparseMatrix* matrix,
                                                int col)

  SLEQP_RETCODE sleqp_sparse_matrix_capture(SleqpSparseMatrix* matrix)

  SLEQP_RETCODE sleqp_sparse_matrix_release(SleqpSparseMatrix** matrix)

  # Functions
  ctypedef SLEQP_RETCODE (*SLEQP_FUNC_SET)(SleqpSparseVec* x,
                                           SLEQP_VALUE_REASON reason,
                                           int num_variables,
                                           int* func_grad_nnz,
                                           int* cons_val_nnz,
                                           int* cons_jac_nnz,
                                           void* func_data)

  ctypedef SLEQP_RETCODE (*SLEQP_FUNC_EVAL)(int num_variables,
                                            const SleqpSparseVec* cons_indices,
                                            double* func_val,
                                            SleqpSparseVec* func_grad,
                                            SleqpSparseVec* cons_val,
                                            SleqpSparseMatrix* cons_jac,
                                            void* func_data)

  ctypedef SLEQP_RETCODE (*SLEQP_HESS_PRODUCT)(int num_variables,
                                               const double* func_dual,
                                               const SleqpSparseVec* direction,
                                               const SleqpSparseVec* cons_duals,
                                               SleqpSparseVec* product,
                                               void* func_data)

  ctypedef SLEQP_RETCODE (*SLEQP_FUNC_FREE)(void* func_data)

  ctypedef struct SleqpFuncCallbacks:
    SLEQP_FUNC_SET set_value,
    SLEQP_FUNC_EVAL func_eval,
    SLEQP_HESS_PRODUCT hess_prod,
    SLEQP_FUNC_FREE func_free

  SLEQP_RETCODE sleqp_func_create(SleqpFunc** fstar,
                                  SleqpFuncCallbacks* callbacks,
                                  int num_variables,
                                  void* func_data)

  SLEQP_RETCODE sleqp_func_release(SleqpFunc** fstar)

  SleqpHessianStruct* sleqp_func_get_hess_struct(SleqpFunc* func)

  # Hessian struct

  int sleqp_hessian_struct_get_num_blocks(SleqpHessianStruct* hessian_struct)

  SLEQP_RETCODE sleqp_hessian_struct_get_block_range(SleqpHessianStruct* hessian_struct,
                                                     int block,
                                                     int* begin,
                                                     int* end)

  SLEQP_RETCODE sleqp_hessian_struct_push_block(SleqpHessianStruct* hessian_struct,
                                                int end)

  SLEQP_RETCODE sleqp_hessian_struct_clear(SleqpHessianStruct* hessian_struct)

  SLEQP_RETCODE sleqp_hessian_struct_get_linear_range(SleqpHessianStruct* hessian_struct,
                                                      int* begin,
                                                      int* end)

  SLEQP_RETCODE sleqp_hessian_struct_fprintf(SleqpHessianStruct* hessian_struct,
                                             libc.stdio.FILE* output)

  # Iterate

  SLEQP_RETCODE sleqp_iterate_capture(SleqpIterate* iterate)

  SleqpSparseVec* sleqp_iterate_get_primal(SleqpIterate* iterate)

  double sleqp_iterate_get_func_val(SleqpIterate* iterate)

  SLEQP_RETCODE sleqp_iterate_set_func_val(SleqpIterate* iterate,
                                           double value)

  SleqpSparseVec* sleqp_iterate_get_func_grad(SleqpIterate* iterate)

  SleqpSparseVec* sleqp_iterate_get_cons_val(SleqpIterate* iterate)

  SleqpSparseMatrix* sleqp_iterate_get_cons_jac(SleqpIterate* iterate)

  SleqpWorkingSet* sleqp_iterate_get_working_set(SleqpIterate* iterate)

  SleqpSparseVec* sleqp_iterate_get_cons_dual(SleqpIterate* iterate)

  SleqpSparseVec* sleqp_iterate_get_vars_dual(SleqpIterate* iterate)

  SLEQP_RETCODE sleqp_iterate_release(SleqpIterate** star)

  # Working set

  SLEQP_RETCODE sleqp_working_set_capture(SleqpWorkingSet* working_set)

  int sleqp_working_set_num_active_vars(const SleqpWorkingSet* working_set)

  int sleqp_working_set_num_active_cons(const SleqpWorkingSet* working_set)

  int sleqp_working_set_size(const SleqpWorkingSet* working_set)

  SLEQP_ACTIVE_STATE sleqp_working_set_get_variable_state(const SleqpWorkingSet* working_set,
                                                          int index)

  SLEQP_ACTIVE_STATE sleqp_working_set_get_constraint_state(const SleqpWorkingSet* working_set,
                                                            int index)

  SLEQP_RETCODE sleqp_working_set_release(SleqpWorkingSet** star)


  # LSQ

  ctypedef SLEQP_RETCODE (*SLEQP_LSQ_EVAL)(int num_variables,
                                           SleqpSparseVec* residual,
                                           void* func_data)

  ctypedef SLEQP_RETCODE (*SLEQP_LSQ_JAC_FORWARD)(int num_variables,
                                                  SleqpSparseVec* forward_direction,
                                                  SleqpSparseVec* product,
                                                  void* func_data)

  ctypedef SLEQP_RETCODE (*SLEQP_LSQ_JAC_ADJOINT)(int num_variables,
                                                  SleqpSparseVec* adjoint_direction,
                                                  SleqpSparseVec* product,
                                                  void* func_data)

  ctypedef struct SleqpLSQCallbacks:
    SLEQP_FUNC_SET set_value,
    SLEQP_LSQ_EVAL lsq_eval,
    SLEQP_LSQ_JAC_FORWARD lsq_jac_forward,
    SLEQP_LSQ_JAC_ADJOINT lsq_jac_adjoint,
    SLEQP_FUNC_EVAL eval_additional,
    SLEQP_HESS_PRODUCT hess_prod_additional,
    SLEQP_FUNC_FREE func_free

  SLEQP_RETCODE sleqp_lsq_func_create(SleqpFunc** fstar,
                                      SleqpLSQCallbacks* callbacks,
                                      int num_variables,
                                      int num_residuals,
                                      double levenberg_marquardt,
                                      SleqpParams* params,
                                      void* func_data)


  # Scaling

  SLEQP_RETCODE sleqp_scaling_create(SleqpScalingData** scaling,
                                     int num_variables,
                                     int num_constraints)

  int sleqp_scaling_get_num_variables(SleqpScalingData* scaling)
  int sleqp_scaling_get_num_constraints(SleqpScalingData* scaling)

  int sleqp_scaling_get_func_weight(SleqpScalingData* scaling)

  SLEQP_RETCODE sleqp_scaling_set_func_weight(SleqpScalingData* scaling,
                                              int weight)

  SLEQP_RETCODE sleqp_scaling_set_func_weight_from_nominal(SleqpScalingData* scaling,
                                                           double nominal_value)

  SLEQP_RETCODE sleqp_scaling_set_var_weight(SleqpScalingData* scaling,
                                             int index,
                                             int weight)

  SLEQP_RETCODE sleqp_scaling_set_var_weight_from_nominal(SleqpScalingData* scaling,
                                                          int index,
                                                          double nominal_value)

  SLEQP_RETCODE sleqp_scaling_set_var_weights_from_nominal(SleqpScalingData* scaling,
                                                           double* nominal_values)

  SLEQP_RETCODE sleqp_scaling_set_cons_weight(SleqpScalingData* scaling,
                                              int index,
                                              int weight)

  int* sleqp_scaling_get_var_weights(SleqpScalingData* scaling)

  int* sleqp_scaling_get_cons_weights(SleqpScalingData* scaling)

  SLEQP_RETCODE sleqp_scaling_set_cons_weights_from_nominal(SleqpScalingData* scaling,
                                                            double* nominal_values)

  SLEQP_RETCODE sleqp_scaling_set_cons_weight_from_nominal(SleqpScalingData* scaling,
                                                           int index,
                                                           double nominal_value);

  SLEQP_RETCODE sleqp_func_scaling_from_gradient(SleqpScalingData* scaling,
                                                 SleqpSparseVec* gradient,
                                                 double eps)

  SLEQP_RETCODE sleqp_scaling_from_cons_jac(SleqpScalingData* scaling,
                                            SleqpSparseMatrix* cons_jac,
                                            double eps)

  SLEQP_RETCODE sleqp_scaling_release(SleqpScalingData** scaling)

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

  SLEQP_RETCODE sleqp_solver_get_violated_constraints(SleqpSolver* solver,
                                                      SleqpIterate* iterate,
                                                      int* violated_constraints,
                                                      int* num_violated_constraints)

  SLEQP_RETCODE sleqp_solver_release(SleqpSolver** star)

  # Problem

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

  double sleqp_params_get_zero_eps(const SleqpParams* params)

  double sleqp_params_get_eps(const SleqpParams* params)

  double sleqp_params_get_deriv_perturbation(const SleqpParams* params)

  double sleqp_params_get_deriv_tolerance(const SleqpParams* params)

  double sleqp_params_get_cauchy_tau(const SleqpParams* params)
  double sleqp_params_get_cauchy_eta(const SleqpParams* params)

  double sleqp_params_get_linesearch_tau(const SleqpParams* params)
  double sleqp_params_get_linesearch_eta(const SleqpParams* params)
  double sleqp_params_get_linesearch_cutoff(const SleqpParams* params)

  double sleqp_params_get_optimality_tolerance(const SleqpParams* params)

  double sleqp_params_get_accepted_reduction(const SleqpParams* params)

  double sleqp_params_get_deadpoint_bound(const SleqpParams* params)

  double sleqp_params_get_newton_relative_tolerance(const SleqpParams* params)

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

  SLEQP_RETCODE sleqp_params_set_deadpoint_bound(SleqpParams* params, double value)

  SLEQP_RETCODE sleqp_params_set_newton_relative_tolerance(SleqpParams* params, double value)

  SLEQP_RETCODE sleqp_params_free(SleqpParams** star)

  # Options

  SLEQP_RETCODE sleqp_options_create(SleqpOptions** star)

  bint sleqp_options_get_perform_newton_step(const SleqpOptions* options)

  bint sleqp_options_get_perform_soc(const SleqpOptions* options)

  bint sleqp_options_get_use_quadratic_model(const SleqpOptions* options)

  SLEQP_DERIV_CHECK sleqp_options_get_deriv_check(const SleqpOptions* options)

  SLEQP_HESSIAN_EVAL sleqp_options_get_hessian_eval(const SleqpOptions* options)

  SLEQP_DUAL_ESTIMATION_TYPE sleqp_options_get_dual_estimation_type(const SleqpOptions* options)

  int sleqp_options_get_quasi_newton_num_iterates(const SleqpOptions* options)

  int sleqp_options_get_max_newton_iterations(const SleqpOptions* options)

  SLEQP_RETCODE sleqp_options_set_perform_newton_step(SleqpOptions* options, bint value)

  SLEQP_RETCODE sleqp_options_set_perform_soc(SleqpOptions* options, bint value)

  SLEQP_RETCODE sleqp_options_set_use_quadratic_model(SleqpOptions* options, bint value)

  SLEQP_RETCODE sleqp_options_set_deriv_check(SleqpOptions* options,
                                              SLEQP_DERIV_CHECK value)

  SLEQP_RETCODE sleqp_options_set_hessian_eval(SleqpOptions* options,
                                               SLEQP_HESSIAN_EVAL value)

  SLEQP_RETCODE sleqp_options_set_dual_estimation_type(SleqpOptions* options,
                                                       SLEQP_DUAL_ESTIMATION_TYPE dual_estimation_type)

  SLEQP_RETCODE sleqp_options_set_quasi_newton_num_iterates(SleqpOptions* options,
                                                            int size)

  SLEQP_RETCODE sleqp_options_set_max_newton_iterations(SleqpOptions* options, int iterations)

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
