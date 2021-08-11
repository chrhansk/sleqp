#cython: language_level=3

cimport libc.stdio
cimport libc.time

cdef extern from "sleqp.h":

  cdef int SLEQP_VERSION_MAJOR
  cdef int SLEQP_VERSION_MINOR
  cdef int SLEQP_VERSION_PATCH

  cdef char* SLEQP_GIT_BRANCH
  cdef char* SLEQP_GIT_COMMIT_HASH

  cdef enum:
    SLEQP_NONE

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
    SLEQP_DERIV_CHECK_SECOND_EXHAUSTIVE,
    SLEQP_DERIV_CHECK_SECOND_SIMPLE

  ctypedef enum SLEQP_HESSIAN_EVAL:
    SLEQP_HESSIAN_EVAL_EXACT,
    SLEQP_HESSIAN_EVAL_SR1,
    SLEQP_HESSIAN_EVAL_SIMPLE_BFGS,
    SLEQP_HESSIAN_EVAL_DAMPED_BFGS

  ctypedef enum SLEQP_BFGS_SIZING:
    SLEQP_BFGS_SIZING_NONE,
    SLEQP_BFGS_SIZING_CENTERED_OL

  ctypedef enum SLEQP_TR_SOLVER:
    SLEQP_TR_SOLVER_TRLIB
    SLEQP_TR_SOLVER_CG
    SLEQP_TR_SOLVER_LSQR
    SLEQP_TR_SOLVER_AUTO

  ctypedef enum SLEQP_PARAMETRIC_CAUCHY:
    SLEQP_PARAMETRIC_CAUCHY_DISABLED
    SLEQP_PARAMETRIC_CAUCHY_COARSE
    SLEQP_PARAMETRIC_CAUCHY_FINE

  ctypedef enum SLEQP_LINESEARCH:
    SLEQP_LINESEARCH_EXACT
    SLEQP_LINESEARCH_APPROX

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

  ctypedef enum SLEQP_PARAM:
    SLEQP_PARAM_ZERO_EPS,
    SLEQP_PARAM_EPS,
    SLEQP_PARAM_DERIV_PERTURBATION,
    SLEQP_PARAM_DERIV_TOL,
    SLEQP_PARAM_CAUCHY_TAU,
    SLEQP_PARAM_CAUCHY_ETA,
    SLEQP_PARAM_LINESEARCH_TAU,
    SLEQP_PARAM_LINESEARCH_ETA,
    SLEQP_PARAM_LINESEARCH_CUTOFF,
    SLEQP_PARAM_FEASIBILITY_TOL,
    SLEQP_PARAM_SLACKNESS_TOL,
    SLEQP_PARAM_STATIONARITY_TOL,
    SLEQP_PARAM_ACCEPTED_REDUCTION,
    SLEQP_PARAM_DEADPOINT_BOUND,
    SLEQP_PARAM_NEWTON_RELATIVE_TOL

  ctypedef enum SLEQP_OPTION_INT:
    SLEQP_OPTION_INT_DERIV_CHECK,
    SLEQP_OPTION_INT_HESSIAN_EVAL,
    SLEQP_OPTION_INT_DUAL_ESTIMATION_TYPE,
    SLEQP_OPTION_INT_NUM_QUASI_NEWTON_ITERATES,
    SLEQP_OPTION_INT_MAX_NEWTON_ITERATIONS,
    SLEQP_OPTION_INT_FLOAT_WARNING_FLAGS,
    SLEQP_OPTION_INT_FLOAT_ERROR_FLAGS,
    SLEQP_OPTION_INT_BFGS_SIZING,
    SLEQP_OPTION_INT_TR_SOLVER,
    SLEQP_OPTION_INT_LINESEARCH,
    SLEQP_OPTION_INT_PARAMETRIC_CAUCHY,
    SLEQP_OPTION_INT_NUM_THREADS

  ctypedef enum SLEQP_OPTION_BOOL:
    SLEQP_OPTION_BOOL_PERFORM_NEWTON_STEP,
    SLEQP_OPTION_BOOL_PERFORM_SOC,
    SLEQP_OPTION_BOOL_USE_QUADRATIC_MODEL,
    SLEQP_OPTION_BOOL_ENABLE_PREPROCESSOR

  ctypedef enum SLEQP_SOLVER_STATE_REAL:
    SLEQP_SOLVER_STATE_REAL_TRUST_RADIUS,
    SLEQP_SOLVER_STATE_REAL_LP_TRUST_RADIUS,
    SLEQP_SOLVER_STATE_REAL_SCALED_FUNC_VAL,
    SLEQP_SOLVER_STATE_REAL_SCALED_MERIT_VAL,
    SLEQP_SOLVER_STATE_REAL_SCALED_FEAS_RES,
    SLEQP_SOLVER_STATE_REAL_SCALED_STAT_RES,
    SLEQP_SOLVER_STATE_REAL_SCALED_SLACK_RES,
    SLEQP_SOLVER_STATE_REAL_PENALTY_PARAM,
    SLEQP_SOLVER_STATE_REAL_MIN_RAYLEIGH,
    SLEQP_SOLVER_STATE_REAL_MAX_RAYLEIGH

  ctypedef enum SLEQP_SOLVER_STATE_INT:
    SLEQP_SOLVER_STATE_INT_LAST_STEP_ON_BDRY
    SLEQP_SOLVER_STATE_INT_LAST_STEP_TYPE

  ctypedef enum SLEQP_SOLVER_STATE_VEC:
    SLEQP_SOLVER_STATE_VEC_SCALED_STAT_RESIDUALS,
    SLEQP_SOLVER_STATE_VEC_SCALED_FEAS_RESIDUALS,
    SLEQP_SOLVER_STATE_VEC_SCALED_CONS_SLACK_RESIDUALS,
    SLEQP_SOLVER_STATE_VEC_SCALED_VAR_SLACK_RESIDUALS

  ctypedef enum SLEQP_STEPTYPE:
    SLEQP_STEPTYPE_NONE,
    SLEQP_STEPTYPE_ACCEPTED,
    SLEQP_STEPTYPE_ACCEPTED_FULL,
    SLEQP_STEPTYPE_ACCEPTED_SOC,
    SLEQP_STEPTYPE_REJECTED

  ctypedef enum SLEQP_SOLVER_EVENT:
    SLEQP_SOLVER_EVENT_ACCEPTED_ITERATE,
    SLEQP_SOLVER_EVENT_PERFORMED_ITERATION,
    SLEQP_SOLVER_EVENT_FINISHED,

  ctypedef struct SleqpSparseVec:
    double* data
    int* indices

    int dim
    int nnz
    int nnz_max

  ctypedef struct SleqpSparseMatrix:
    pass

  ctypedef struct SleqpProblem:
    pass

  ctypedef struct SleqpSolver:
    pass

  ctypedef struct SleqpScaling:
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
  ctypedef SLEQP_RETCODE (*SLEQP_FUNC_SET)(SleqpFunc* func,
                                           SleqpSparseVec* x,
                                           SLEQP_VALUE_REASON reason,
                                           int* func_grad_nnz,
                                           int* cons_val_nnz,
                                           int* cons_jac_nnz,
                                           void* func_data)

  ctypedef SLEQP_RETCODE (*SLEQP_FUNC_VAL)(SleqpFunc* func,
                                           double* func_val,
                                           void* func_data)

  ctypedef SLEQP_RETCODE (*SLEQP_FUNC_GRAD)(SleqpFunc* func,
                                            SleqpSparseVec* func_grad,
                                            void* func_data)

  ctypedef SLEQP_RETCODE (*SLEQP_FUNC_CONS_VAL)(SleqpFunc* func,
                                                const SleqpSparseVec* cons_indices,
                                                SleqpSparseVec* cons_val,
                                                void* func_data)

  ctypedef SLEQP_RETCODE (*SLEQP_FUNC_CONS_JAC)(SleqpFunc* func,
                                                const SleqpSparseVec* cons_indices,
                                                SleqpSparseMatrix* cons_jac,
                                                void* func_data)

  ctypedef SLEQP_RETCODE (*SLEQP_HESS_PROD)(SleqpFunc* func,
                                            const double* func_dual,
                                            const SleqpSparseVec* direction,
                                            const SleqpSparseVec* cons_duals,
                                            SleqpSparseVec* product,
                                            void* func_data)

  ctypedef SLEQP_RETCODE (*SLEQP_FUNC_FREE)(void* func_data)

  ctypedef struct SleqpFuncCallbacks:
    SLEQP_FUNC_SET set_value
    SLEQP_FUNC_VAL func_val
    SLEQP_FUNC_GRAD func_grad
    SLEQP_FUNC_CONS_VAL cons_val
    SLEQP_FUNC_CONS_JAC cons_jac
    SLEQP_HESS_PROD hess_prod
    SLEQP_FUNC_FREE func_free

  SLEQP_RETCODE sleqp_func_create(SleqpFunc** fstar,
                                  SleqpFuncCallbacks* callbacks,
                                  int num_variables,
                                  int num_constraints,
                                  void* func_data)

  SLEQP_RETCODE sleqp_func_set_callbacks(SleqpFunc* func,
                                         SleqpFuncCallbacks* callbacks)

  int sleqp_func_get_num_variables(SleqpFunc* func)

  int sleqp_func_get_num_constraints(SleqpFunc* func)

  SLEQP_RETCODE sleqp_func_capture(SleqpFunc* func)

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

  SLEQP_RETCODE sleqp_hessian_struct_capture(SleqpHessianStruct* hessian_struct)

  SLEQP_RETCODE sleqp_hessian_struct_release(SleqpHessianStruct** star)

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

  SleqpProblem* sleqp_working_set_get_problem(const SleqpWorkingSet* working_set)

  SLEQP_ACTIVE_STATE sleqp_working_set_get_variable_state(const SleqpWorkingSet* working_set,
                                                          int index)

  SLEQP_ACTIVE_STATE sleqp_working_set_get_constraint_state(const SleqpWorkingSet* working_set,
                                                            int index)

  SLEQP_RETCODE sleqp_working_set_release(SleqpWorkingSet** star)


  # LSQ

  ctypedef SLEQP_RETCODE (*SLEQP_LSQ_RESIDUALS)(SleqpFunc* func,
                                                SleqpSparseVec* residual,
                                                void* func_data)

  ctypedef SLEQP_RETCODE (*SLEQP_LSQ_JAC_FORWARD)(SleqpFunc* func,
                                                  SleqpSparseVec* forward_direction,
                                                  SleqpSparseVec* product,
                                                  void* func_data)

  ctypedef SLEQP_RETCODE (*SLEQP_LSQ_JAC_ADJOINT)(SleqpFunc* func,
                                                  SleqpSparseVec* adjoint_direction,
                                                  SleqpSparseVec* product,
                                                  void* func_data)

  ctypedef struct SleqpLSQCallbacks:
    SLEQP_FUNC_SET        set_value,
    SLEQP_LSQ_RESIDUALS   lsq_residuals,
    SLEQP_LSQ_JAC_FORWARD lsq_jac_forward,
    SLEQP_LSQ_JAC_ADJOINT lsq_jac_adjoint,
    SLEQP_FUNC_CONS_VAL   cons_val,
    SLEQP_FUNC_CONS_JAC   cons_jac,
    SLEQP_FUNC_FREE       func_free

  SLEQP_RETCODE sleqp_lsq_func_create(SleqpFunc** fstar,
                                      SleqpLSQCallbacks* callbacks,
                                      int num_variables,
                                      int num_constraints,
                                      int num_residuals,
                                      double levenberg_marquardt,
                                      SleqpParams* params,
                                      void* func_data)

  SLEQP_RETCODE sleqp_lsq_func_set_callbacks(SleqpFunc* func,
                                             SleqpLSQCallbacks* callbacks)

  # Scaling

  SLEQP_RETCODE sleqp_scaling_create(SleqpScaling** scaling,
                                     int num_variables,
                                     int num_constraints)

  int sleqp_scaling_get_num_variables(SleqpScaling* scaling)
  int sleqp_scaling_get_num_constraints(SleqpScaling* scaling)

  int sleqp_scaling_get_func_weight(SleqpScaling* scaling)

  SLEQP_RETCODE sleqp_scaling_set_func_weight(SleqpScaling* scaling,
                                              int weight)

  SLEQP_RETCODE sleqp_scaling_set_func_weight_from_nominal(SleqpScaling* scaling,
                                                           double nominal_value)

  SLEQP_RETCODE sleqp_scaling_set_var_weight(SleqpScaling* scaling,
                                             int index,
                                             int weight)

  SLEQP_RETCODE sleqp_scaling_set_var_weight_from_nominal(SleqpScaling* scaling,
                                                          int index,
                                                          double nominal_value)

  SLEQP_RETCODE sleqp_scaling_set_var_weights_from_nominal(SleqpScaling* scaling,
                                                           double* nominal_values)

  SLEQP_RETCODE sleqp_scaling_set_cons_weight(SleqpScaling* scaling,
                                              int index,
                                              int weight)

  int* sleqp_scaling_get_var_weights(SleqpScaling* scaling)

  int* sleqp_scaling_get_cons_weights(SleqpScaling* scaling)

  SLEQP_RETCODE sleqp_scaling_set_cons_weights_from_nominal(SleqpScaling* scaling,
                                                            double* nominal_values)

  SLEQP_RETCODE sleqp_scaling_set_cons_weight_from_nominal(SleqpScaling* scaling,
                                                           int index,
                                                           double nominal_value);

  SLEQP_RETCODE sleqp_func_scaling_from_gradient(SleqpScaling* scaling,
                                                 SleqpSparseVec* gradient,
                                                 double eps)

  SLEQP_RETCODE sleqp_scaling_from_cons_jac(SleqpScaling* scaling,
                                            SleqpSparseMatrix* cons_jac,
                                            double eps)

  SLEQP_RETCODE sleqp_scaling_release(SleqpScaling** scaling)

  # Solver

  SLEQP_RETCODE sleqp_solver_create(SleqpSolver** star,
                                    SleqpProblem* problem,
                                    SleqpParams* params,
                                    SleqpOptions* options,
                                    SleqpSparseVec* x,
                                    SleqpScaling* scaling)

  const char* sleqp_solver_info(const SleqpSolver* solver)

  SLEQP_RETCODE sleqp_solver_solve(SleqpSolver* solver,
                                   int max_num_iterations,
                                   double time_limit) nogil

  SLEQP_RETCODE sleqp_solver_get_real_state(const SleqpSolver* solver,
                                            SLEQP_SOLVER_STATE_REAL state,
                                            double* value)

  SLEQP_RETCODE sleqp_solver_get_int_state(const SleqpSolver* solver,
                                           SLEQP_SOLVER_STATE_INT state,
                                           int* value)

  SLEQP_RETCODE sleqp_solver_get_vec_state(const SleqpSolver* solver,
                                           SLEQP_SOLVER_STATE_VEC value,
                                           SleqpSparseVec* result)

  SLEQP_STATUS sleqp_solver_get_status(SleqpSolver* solver)

  SLEQP_RETCODE sleqp_solver_abort(SleqpSolver* solver)

  int sleqp_solver_get_iterations(SleqpSolver* solver)

  double sleqp_solver_get_elapsed_seconds(SleqpSolver* solver)

  SLEQP_RETCODE sleqp_solver_add_callback(SleqpSolver* solver,
                                          SLEQP_SOLVER_EVENT solver_event,
                                          void* callback_func,
                                          void* callback_data)

  SLEQP_RETCODE sleqp_solver_remove_callback(SleqpSolver* solver,
                                             SLEQP_SOLVER_EVENT solver_event,
                                             void* callback_func,
                                             void* callback_data)

  SLEQP_RETCODE sleqp_solver_get_solution(SleqpSolver* solver,
                                          SleqpIterate** iterate)

  SLEQP_RETCODE sleqp_solver_get_violated_constraints(SleqpSolver* solver,
                                                      SleqpIterate* iterate,
                                                      int* violated_constraints,
                                                      int* num_violated_constraints)

  SLEQP_RETCODE sleqp_solver_release(SleqpSolver** star)

  # Problem

  SLEQP_RETCODE sleqp_problem_create_simple(SleqpProblem** star,
                                            SleqpFunc* func,
                                            SleqpParams* params,
                                            const SleqpSparseVec* var_lb,
                                            const SleqpSparseVec* var_ub,
                                            const SleqpSparseVec* cons_lb,
                                            const SleqpSparseVec* cons_ub)

  SLEQP_RETCODE sleqp_problem_create(SleqpProblem** star,
                                     SleqpFunc* func,
                                     SleqpParams* params,
                                     const SleqpSparseVec* var_lb,
                                     const SleqpSparseVec* var_ub,
                                     const SleqpSparseVec* genereal_lb,
                                     const SleqpSparseVec* genereal_ub,
                                     const SleqpSparseMatrix* linear_coeffs,
                                     const SleqpSparseVec* linear_lb,
                                     const SleqpSparseVec* linear_ub)

  int sleqp_problem_num_constraints(SleqpProblem* problem)

  int sleqp_problem_num_variables(SleqpProblem* problem)

  SleqpSparseVec* sleqp_problem_var_lb(SleqpProblem* problem)

  SleqpSparseVec* sleqp_problem_var_ub(SleqpProblem* problem)

  SleqpSparseVec* sleqp_problem_cons_lb(SleqpProblem* problem)

  SleqpSparseVec* sleqp_problem_cons_ub(SleqpProblem* problem)

  SLEQP_RETCODE sleqp_problem_release(SleqpProblem** star)

  # Parameters

  SLEQP_RETCODE sleqp_params_create(SleqpParams** star)

  double sleqp_params_get(const SleqpParams* params,
                          SLEQP_PARAM param)

  SLEQP_RETCODE sleqp_params_set(SleqpParams* params,
                                 SLEQP_PARAM param,
                                 double value)

  SLEQP_RETCODE sleqp_params_release(SleqpParams** star)

  # Options

  SLEQP_RETCODE sleqp_options_create(SleqpOptions** star)

  int sleqp_options_get_int(const SleqpOptions* options,
                            SLEQP_OPTION_INT option)

  SLEQP_RETCODE sleqp_options_set_int(SleqpOptions* options,
                                      SLEQP_OPTION_INT option,
                                      int value)


  bint sleqp_options_get_bool(const SleqpOptions* options,
                              SLEQP_OPTION_BOOL option)

  SLEQP_RETCODE sleqp_options_set_bool(SleqpOptions* options,
                                       SLEQP_OPTION_BOOL option,
                                       bint value)

  SLEQP_RETCODE sleqp_options_release(SleqpOptions** star)

  # Logging

  void sleqp_log_set_level(SLEQP_LOG_LEVEL value)

  ctypedef void (*SLEQP_LOG_HANDLER)(SLEQP_LOG_LEVEL level,
                                     libc.time.time_t time,
                                     const char* message)

  void sleqp_log_set_handler(SLEQP_LOG_HANDLER handler)

  SLEQP_LOG_LEVEL sleqp_log_level()

  # Numerics
  double sleqp_infinity()
