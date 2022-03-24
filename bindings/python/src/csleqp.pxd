#cython: language_level=3

cimport libc.stdio
cimport libc.time

cdef extern from "sleqp.h":

  ctypedef bint bool

  cdef int SLEQP_VERSION_MAJOR
  cdef int SLEQP_VERSION_MINOR
  cdef int SLEQP_VERSION_PATCH

  cdef char* SLEQP_GIT_BRANCH
  cdef char* SLEQP_GIT_COMMIT_HASH

  cdef enum:
    SLEQP_NONE

  ctypedef enum SLEQP_RETCODE:
    SLEQP_ERROR,
    SLEQP_OKAY

  ctypedef enum SLEQP_ERROR_TYPE:
    SLEQP_FAILED_ASSERTION,
    SLEQP_NOMEM,
    SLEQP_INTERNAL_ERROR,
    SLEQP_FUNC_EVAL_ERROR,
    SLEQP_CALLBACK_ERROR,
    SLEQP_MATH_ERROR,
    SLEQP_INVALID_DERIV,
    SLEQP_ILLEGAL_ARGUMENT

  ctypedef enum SLEQP_ACTIVE_STATE:
    SLEQP_INACTIVE,
    SLEQP_ACTIVE_LOWER,
    SLEQP_ACTIVE_UPPER,
    SLEQP_ACTIVE_BOTH

  ctypedef enum SLEQP_STATUS:
    SLEQP_STATUS_UNKNOWN,
    SLEQP_STATUS_RUNNING,
    SLEQP_STATUS_OPTIMAL,
    SLEQP_STATUS_INFEASIBLE,
    SLEQP_STATUS_UNBOUNDED,
    SLEQP_STATUS_ABORT_DEADPOINT,
    SLEQP_STATUS_ABORT_ITER,
    SLEQP_STATUS_ABORT_MANUAL,
    SLEQP_STATUS_ABORT_TIME

  ctypedef enum SLEQP_DERIV_CHECK:
    SLEQP_DERIV_CHECK_SKIP,
    SLEQP_DERIV_CHECK_FIRST_OBJ,
    SLEQP_DERIV_CHECK_FIRST_CONS,
    SLEQP_DERIV_CHECK_FIRST,
    SLEQP_DERIV_CHECK_SECOND_OBJ,
    SLEQP_DERIV_CHECK_SECOND_CONS,
    SLEQP_DERIV_CHECK_SECOND_EXHAUSTIVE,
    SLEQP_DERIV_CHECK_SECOND_SIMPLE

  ctypedef enum SLEQP_HESS_EVAL:
    SLEQP_HESS_EVAL_EXACT,
    SLEQP_HESS_EVAL_SR1,
    SLEQP_HESS_EVAL_SIMPLE_BFGS,
    SLEQP_HESS_EVAL_DAMPED_BFGS

  ctypedef enum SLEQP_BFGS_SIZING:
    SLEQP_BFGS_SIZING_NONE,
    SLEQP_BFGS_SIZING_CENTERED_OL

  ctypedef enum SLEQP_TR_SOLVER:
    SLEQP_TR_SOLVER_TRLIB
    SLEQP_TR_SOLVER_CG
    SLEQP_TR_SOLVER_LSQR
    SLEQP_TR_SOLVER_AUTO

  ctypedef enum SLEQP_POLISHING_TYPE:
    SLEQP_POLISHING_NONE
    SLEQP_POLISHING_ZERO_DUAL
    SLEQP_POLISHING_INACTIVE

  ctypedef enum SLEQP_STEP_RULE:
    SLEQP_STEP_RULE_DIRECT
    SLEQP_STEP_RULE_WINDOW
    SLEQP_STEP_RULE_MINSTEP

  ctypedef enum SLEQP_PARAMETRIC_CAUCHY:
    SLEQP_PARAMETRIC_CAUCHY_DISABLED
    SLEQP_PARAMETRIC_CAUCHY_COARSE
    SLEQP_PARAMETRIC_CAUCHY_FINE

  ctypedef enum SLEQP_INITIAL_TR_CHOICE:
    SLEQP_INITIAL_TR_CHOICE_NARROW
    SLEQP_INITIAL_TR_CHOICE_WIDE

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
    SLEQP_PARAM_FEAS_TOL,
    SLEQP_PARAM_SLACK_TOL,
    SLEQP_PARAM_STAT_TOL,
    SLEQP_PARAM_ACCEPTED_REDUCTION,
    SLEQP_PARAM_DEADPOINT_BOUND

  ctypedef enum SLEQP_OPTION_INT:
    SLEQP_OPTION_INT_NUM_QUASI_NEWTON_ITERATES,
    SLEQP_OPTION_INT_MAX_NEWTON_ITERATIONS,
    SLEQP_OPTION_INT_NUM_THREADS

  ctypedef enum SLEQP_OPTION_ENUM:
    SLEQP_OPTION_ENUM_DERIV_CHECK,
    SLEQP_OPTION_ENUM_HESS_EVAL,
    SLEQP_OPTION_ENUM_DUAL_ESTIMATION_TYPE,
    SLEQP_OPTION_ENUM_FLOAT_WARNING_FLAGS,
    SLEQP_OPTION_ENUM_FLOAT_ERROR_FLAGS,
    SLEQP_OPTION_ENUM_BFGS_SIZING,
    SLEQP_OPTION_ENUM_TR_SOLVER,
    SLEQP_OPTION_ENUM_POLISHING_TYPE,
    SLEQP_OPTION_ENUM_STEP_RULE,
    SLEQP_OPTION_ENUM_LINESEARCH,
    SLEQP_OPTION_ENUM_PARAMETRIC_CAUCHY,
    SLEQP_OPTION_ENUM_INITIAL_TR_CHOICE

  ctypedef enum SLEQP_OPTION_BOOL:
    SLEQP_OPTION_BOOL_PERFORM_NEWTON_STEP,
    SLEQP_OPTION_BOOL_GLOBAL_PENALTY_RESETS,
    SLEQP_OPTION_BOOL_PERFORM_SOC,
    SLEQP_OPTION_BOOL_USE_QUADRATIC_MODEL,
    SLEQP_OPTION_BOOL_ENABLE_PREPROCESSOR,
    SLEQP_OPTION_BOOL_LP_RESOLVES

  ctypedef enum SLEQP_SOLVER_STATE_REAL:
    SLEQP_SOLVER_STATE_REAL_TRUST_RADIUS,
    SLEQP_SOLVER_STATE_REAL_LP_TRUST_RADIUS,
    SLEQP_SOLVER_STATE_REAL_SCALED_OBJ_VAL,
    SLEQP_SOLVER_STATE_REAL_SCALED_MERIT_VAL,
    SLEQP_SOLVER_STATE_REAL_SCALED_FEAS_RES,
    SLEQP_SOLVER_STATE_REAL_SCALED_STAT_RES,
    SLEQP_SOLVER_STATE_REAL_SCALED_SLACK_RES,
    SLEQP_SOLVER_STATE_REAL_PENALTY_PARAM,
    SLEQP_SOLVER_STATE_REAL_MIN_RAYLEIGH,
    SLEQP_SOLVER_STATE_REAL_MAX_RAYLEIGH

  ctypedef enum SLEQP_SOLVER_STATE_INT:
    SLEQP_SOLVER_STATE_INT_LAST_STEP_ON_BDRY,
    SLEQP_SOLVER_STATE_INT_ITERATION,
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

  ctypedef struct SleqpVec:
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

  cdef struct SleqpHessStruct:
    pass


  ctypedef enum SLEQP_LOG_LEVEL:
    SLEQP_LOG_ERROR,
    SLEQP_LOG_WARN,
    SLEQP_LOG_INFO,
    SLEQP_LOG_DEBUG

  SLEQP_ERROR_TYPE sleqp_error_type()

  const char* sleqp_error_msg()

  # Sparse vectors
  SLEQP_RETCODE sleqp_vec_create(SleqpVec** vec,
                                           int dim,
                                           int nnz_max)

  SLEQP_RETCODE sleqp_vec_create_empty(SleqpVec** vec,
                                                 int dim)

  SLEQP_RETCODE sleqp_vec_create_full(SleqpVec** vec,
                                                int dim)

  SLEQP_RETCODE sleqp_vec_fprintf(const SleqpVec* vec,
                                            libc.stdio.FILE* output)

  SLEQP_RETCODE sleqp_vec_free(SleqpVec** vec)

  SLEQP_RETCODE sleqp_vec_clear(SleqpVec* vec)

  SLEQP_RETCODE sleqp_vec_reserve(SleqpVec* vec, int nnz)

  SLEQP_RETCODE sleqp_vec_resize(SleqpVec* vec,
                                           int dim)

  SLEQP_RETCODE sleqp_vec_push(SleqpVec* vec,
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

  int sleqp_sparse_matrix_num_cols(SleqpSparseMatrix* matrix)

  int sleqp_sparse_matrix_num_rows(SleqpSparseMatrix* matrix)

  int sleqp_sparse_matrix_nnz(SleqpSparseMatrix* matrix)

  int sleqp_sparse_matrix_nnz_max(SleqpSparseMatrix* matrix)

  double* sleqp_sparse_matrix_data(SleqpSparseMatrix* matrix)

  int* sleqp_sparse_matrix_cols(SleqpSparseMatrix* matrix)

  int* sleqp_sparse_matrix_rows(SleqpSparseMatrix* matrix)

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
                                           SleqpVec* x,
                                           SLEQP_VALUE_REASON reason,
                                           bool* reject,
                                           int* obj_grad_nnz,
                                           int* cons_val_nnz,
                                           int* cons_jac_nnz,
                                           void* func_data)

  ctypedef SLEQP_RETCODE (*SLEQP_FUNC_OBJ_VAL)(SleqpFunc* func,
                                               double* func_val,
                                               void* func_data)

  ctypedef SLEQP_RETCODE (*SLEQP_FUNC_OBJ_GRAD)(SleqpFunc* func,
                                                SleqpVec* func_grad,
                                                void* func_data)

  ctypedef SLEQP_RETCODE (*SLEQP_FUNC_CONS_VAL)(SleqpFunc* func,
                                                SleqpVec* cons_val,
                                                void* func_data)

  ctypedef SLEQP_RETCODE (*SLEQP_FUNC_CONS_JAC)(SleqpFunc* func,
                                                SleqpSparseMatrix* cons_jac,
                                                void* func_data)

  ctypedef SLEQP_RETCODE (*SLEQP_FUNC_HESS_PROD)(SleqpFunc* func,
                                                 const double* obj_dual,
                                                 const SleqpVec* direction,
                                                 const SleqpVec* cons_duals,
                                                 SleqpVec* product,
                                                 void* func_data)

  ctypedef SLEQP_RETCODE (*SLEQP_FUNC_FREE)(void* func_data)

  ctypedef struct SleqpFuncCallbacks:
    SLEQP_FUNC_SET       set_value
    SLEQP_FUNC_OBJ_VAL   obj_val
    SLEQP_FUNC_OBJ_GRAD  obj_grad
    SLEQP_FUNC_CONS_VAL  cons_val
    SLEQP_FUNC_CONS_JAC  cons_jac
    SLEQP_FUNC_HESS_PROD hess_prod
    SLEQP_FUNC_FREE      func_free

  SLEQP_RETCODE sleqp_func_create(SleqpFunc** fstar,
                                  SleqpFuncCallbacks* callbacks,
                                  int num_variables,
                                  int num_constraints,
                                  void* func_data)

  SLEQP_RETCODE sleqp_func_set_callbacks(SleqpFunc* func,
                                         SleqpFuncCallbacks* callbacks)

  int sleqp_func_num_vars(SleqpFunc* func)

  int sleqp_func_num_cons(SleqpFunc* func)

  SLEQP_RETCODE sleqp_func_capture(SleqpFunc* func)

  SLEQP_RETCODE sleqp_func_release(SleqpFunc** fstar)

  SleqpHessStruct* sleqp_func_hess_struct(SleqpFunc* func)

  # Hessian struct

  int sleqp_hess_struct_num_blocks(SleqpHessStruct* hess_struct)

  SLEQP_RETCODE sleqp_hess_struct_block_range(SleqpHessStruct* hess_struct,
                                              int block,
                                              int* begin,
                                              int* end)

  SLEQP_RETCODE sleqp_hess_struct_push_block(SleqpHessStruct* hess_struct,
                                             int end)

  SLEQP_RETCODE sleqp_hess_struct_clear(SleqpHessStruct* hess_struct)

  SLEQP_RETCODE sleqp_hess_struct_lin_range(SleqpHessStruct* hess_struct,
                                            int* begin,
                                            int* end)

  SLEQP_RETCODE sleqp_hess_struct_fprintf(SleqpHessStruct* hess_struct,
                                          libc.stdio.FILE* output)

  SLEQP_RETCODE sleqp_hess_struct_capture(SleqpHessStruct* hess_struct)

  SLEQP_RETCODE sleqp_hess_struct_release(SleqpHessStruct** star)

  # Iterate

  SLEQP_RETCODE sleqp_iterate_capture(SleqpIterate* iterate)

  SleqpVec* sleqp_iterate_primal(SleqpIterate* iterate)

  double sleqp_iterate_obj_val(SleqpIterate* iterate)

  SLEQP_RETCODE sleqp_iterate_set_obj_val(SleqpIterate* iterate,
                                          double value)

  SleqpVec* sleqp_iterate_obj_grad(SleqpIterate* iterate)

  SleqpVec* sleqp_iterate_cons_val(SleqpIterate* iterate)

  SleqpSparseMatrix* sleqp_iterate_cons_jac(SleqpIterate* iterate)

  SleqpWorkingSet* sleqp_iterate_working_set(SleqpIterate* iterate)

  SleqpVec* sleqp_iterate_cons_dual(SleqpIterate* iterate)

  SleqpVec* sleqp_iterate_vars_dual(SleqpIterate* iterate)

  SLEQP_RETCODE sleqp_iterate_release(SleqpIterate** star)

  # Working set

  SLEQP_RETCODE sleqp_working_set_capture(SleqpWorkingSet* working_set)

  int sleqp_working_set_num_active_vars(const SleqpWorkingSet* working_set)

  int sleqp_working_set_num_active_cons(const SleqpWorkingSet* working_set)

  int sleqp_working_set_size(const SleqpWorkingSet* working_set)

  SleqpProblem* sleqp_working_set_problem(const SleqpWorkingSet* working_set)

  SLEQP_ACTIVE_STATE sleqp_working_set_var_state(const SleqpWorkingSet* working_set,
                                                 int index)

  SLEQP_ACTIVE_STATE sleqp_working_set_cons_state(const SleqpWorkingSet* working_set,
                                                  int index)

  SLEQP_RETCODE sleqp_working_set_release(SleqpWorkingSet** star)


  # LSQ

  ctypedef SLEQP_RETCODE (*SLEQP_LSQ_RESIDUALS)(SleqpFunc* func,
                                                SleqpVec* residual,
                                                void* func_data)

  ctypedef SLEQP_RETCODE (*SLEQP_LSQ_JAC_FORWARD)(SleqpFunc* func,
                                                  SleqpVec* forward_direction,
                                                  SleqpVec* product,
                                                  void* func_data)

  ctypedef SLEQP_RETCODE (*SLEQP_LSQ_JAC_ADJOINT)(SleqpFunc* func,
                                                  SleqpVec* adjoint_direction,
                                                  SleqpVec* product,
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

  int sleqp_scaling_num_vars(SleqpScaling* scaling)
  int sleqp_scaling_num_cons(SleqpScaling* scaling)

  int sleqp_scaling_obj_weight(SleqpScaling* scaling)

  SLEQP_RETCODE sleqp_scaling_set_obj_weight(SleqpScaling* scaling,
                                             int weight)

  SLEQP_RETCODE sleqp_scaling_set_obj_weight_from_nominal(SleqpScaling* scaling,
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

  int* sleqp_scaling_var_weights(SleqpScaling* scaling)

  int* sleqp_scaling_cons_weights(SleqpScaling* scaling)

  SLEQP_RETCODE sleqp_scaling_set_cons_weights_from_nominal(SleqpScaling* scaling,
                                                            double* nominal_values)

  SLEQP_RETCODE sleqp_scaling_set_cons_weight_from_nominal(SleqpScaling* scaling,
                                                           int index,
                                                           double nominal_value);

  SLEQP_RETCODE sleqp_obj_scaling_from_grad(SleqpScaling* scaling,
                                            SleqpVec* gradient,
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
                                    SleqpVec* x,
                                    SleqpScaling* scaling)

  const char* sleqp_solver_info(const SleqpSolver* solver)

  SLEQP_RETCODE sleqp_solver_solve(SleqpSolver* solver,
                                   int max_num_iterations,
                                   double time_limit) nogil

  SLEQP_RETCODE sleqp_solver_real_state(const SleqpSolver* solver,
                                        SLEQP_SOLVER_STATE_REAL state,
                                        double* value)

  SLEQP_RETCODE sleqp_solver_int_state(const SleqpSolver* solver,
                                       SLEQP_SOLVER_STATE_INT state,
                                       int* value)

  SLEQP_RETCODE sleqp_solver_vec_state(const SleqpSolver* solver,
                                       SLEQP_SOLVER_STATE_VEC value,
                                       SleqpVec* result)

  SLEQP_STATUS sleqp_solver_status(SleqpSolver* solver)

  SLEQP_RETCODE sleqp_solver_abort(SleqpSolver* solver)

  int sleqp_solver_iterations(SleqpSolver* solver)

  double sleqp_solver_elapsed_seconds(SleqpSolver* solver)

  SLEQP_RETCODE sleqp_solver_add_callback(SleqpSolver* solver,
                                          SLEQP_SOLVER_EVENT solver_event,
                                          void* callback_func,
                                          void* callback_data)

  SLEQP_RETCODE sleqp_solver_remove_callback(SleqpSolver* solver,
                                             SLEQP_SOLVER_EVENT solver_event,
                                             void* callback_func,
                                             void* callback_data)

  SLEQP_RETCODE sleqp_solver_solution(SleqpSolver* solver,
                                      SleqpIterate** iterate)

  SLEQP_RETCODE sleqp_solver_violated_constraints(SleqpSolver* solver,
                                                  SleqpIterate* iterate,
                                                  int* violated_constraints,
                                                  int* num_violated_constraints)

  SLEQP_RETCODE sleqp_solver_release(SleqpSolver** star)

  # Problem

  SLEQP_RETCODE sleqp_problem_create_simple(SleqpProblem** star,
                                            SleqpFunc* func,
                                            SleqpParams* params,
                                            const SleqpVec* var_lb,
                                            const SleqpVec* var_ub,
                                            const SleqpVec* cons_lb,
                                            const SleqpVec* cons_ub)

  SLEQP_RETCODE sleqp_problem_create(SleqpProblem** star,
                                     SleqpFunc* func,
                                     SleqpParams* params,
                                     const SleqpVec* var_lb,
                                     const SleqpVec* var_ub,
                                     const SleqpVec* genereal_lb,
                                     const SleqpVec* genereal_ub,
                                     const SleqpSparseMatrix* linear_coeffs,
                                     const SleqpVec* linear_lb,
                                     const SleqpVec* linear_ub)

  int sleqp_problem_num_cons(SleqpProblem* problem)

  int sleqp_problem_num_vars(SleqpProblem* problem)

  SleqpVec* sleqp_problem_vars_lb(SleqpProblem* problem)

  SleqpVec* sleqp_problem_vars_ub(SleqpProblem* problem)

  SleqpVec* sleqp_problem_cons_lb(SleqpProblem* problem)

  SleqpVec* sleqp_problem_cons_ub(SleqpProblem* problem)

  SLEQP_RETCODE sleqp_problem_release(SleqpProblem** star)

  # Parameters

  SLEQP_RETCODE sleqp_params_create(SleqpParams** star)

  double sleqp_params_value(const SleqpParams* params,
                            SLEQP_PARAM param)

  SLEQP_RETCODE sleqp_params_set_value(SleqpParams* params,
                                       SLEQP_PARAM param,
                                       double value)

  SLEQP_RETCODE sleqp_params_release(SleqpParams** star)

  # Options

  SLEQP_RETCODE sleqp_options_create(SleqpOptions** star)

  int sleqp_options_enum_value(const SleqpOptions* options,
                               SLEQP_OPTION_ENUM option)

  SLEQP_RETCODE sleqp_options_set_enum_value(SleqpOptions* options,
                                             SLEQP_OPTION_ENUM option,
                                             int value)

  int sleqp_options_int_value(const SleqpOptions* options,
                              SLEQP_OPTION_INT option)

  SLEQP_RETCODE sleqp_options_set_int_value(SleqpOptions* options,
                                            SLEQP_OPTION_INT option,
                                            int value)


  bool sleqp_options_bool_value(const SleqpOptions* options,
                                SLEQP_OPTION_BOOL option)

  SLEQP_RETCODE sleqp_options_set_bool_value(SleqpOptions* options,
                                             SLEQP_OPTION_BOOL option,
                                             bool value)

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
