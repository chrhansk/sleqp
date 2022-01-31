#ifndef SLEQP_MEX_FIELDS
#define SLEQP_MEX_FIELDS

#define MEX_ERROR_PREFIX "SLEQP"

#define MEX_IDENTIFIER_FAILED_ASSERTION MEX_ERROR_PREFIX ":FailedAssertion"
#define MEX_IDENTIFIER_NOMEM MEX_ERROR_PREFIX ":OutOfMemory"
#define MEX_IDENTIFIER_INTERNAL_ERROR MEX_ERROR_PREFIX ":InternalError"
#define MEX_IDENTIFIER_FUNC_EVAL_ERROR MEX_ERROR_PREFIX ":FuncEvalError"
#define MEX_IDENTIFIER_CALLBACK_ERROR MEX_ERROR_PREFIX ":CallbackError"
#define MEX_IDENTIFIER_MATH_ERROR MEX_ERROR_PREFIX ":MathError"
#define MEX_IDENTIFIER_INVALID_DERIV MEX_ERROR_PREFIX ":InvalidDeriv"
#define MEX_IDENTIFIER_ILLEGAL_ARGUMENT MEX_ERROR_PREFIX ":IllegalArgument"

#define MEX_IDENTIFIER_DEFAULT MEX_ERROR_PREFIX ":Error"

#define MEX_COMMAND_INFO "info"
#define MEX_COMMAND_SOLVE "solve"
#define MEX_COMMAND_SOLVE_LSQ "solve_lsq"

#define MATLAB_FUNC_FEVAL "feval"
#define MATLAB_FUNC_DISP "disp"

#define MEX_INPUT_OBJ_VAL "obj_val"
#define MEX_INPUT_OBJ_GRAD "obj_grad"
#define MEX_INPUT_CONS_VAL "cons_val"
#define MEX_INPUT_CONS_JAC "cons_jac"
#define MEX_INPUT_HESS "hess"
#define MEX_INPUT_HESS_PROD "hess_prod"

#define MEX_INFO_VERSION "version"
#define MEX_INFO_VERSION_MAJOR "version_major"
#define MEX_INFO_VERSION_MINOR "version_minor"
#define MEX_INFO_VERSION_PATCH "version_patch"
#define MEX_INFO_GITHASH "git_hash"
#define MEX_INFO_FACT_NAME "fact_name"
#define MEX_INFO_FACT_VERSION "fact_version"
#define MEX_INFO_LPS_NAME "lps_name"
#define MEX_INFO_LPS_VERSION "lps_version"
#define MEX_INFO_TRLIB_VERSION "trlib_version"

#define MEX_INPUT_LSQ_RES "lsq_residuals"
#define MEX_INPUT_LSQ_JAC_FWD "lsq_jac_forward"
#define MEX_INPUT_LSQ_JAC_ADJ "lsq_jac_adjoint"

#define MEX_INPUT_CONS_LB "cons_lb"
#define MEX_INPUT_CONS_UB "cons_ub"

#define MEX_INPUT_VAR_LB "var_lb"
#define MEX_INPUT_VAR_UB "var_ub"

#define MEX_CALLBACK_ACCEPTED_ITERATE "accepted_iterate"

#define MEX_OUTPUT_PRIMAL "primal"
#define MEX_OUTPUT_CONS_DUAL "cons_dual"
#define MEX_OUTPUT_VARS_DUAL "vars_dal"
#define MEX_OUTPUT_ELAPSED "elapsed"
#define MEX_OUTPUT_ITER "iterations"
#define MEX_OUTPUT_STATUS "status"
#define MEX_OUTPUT_WORKING_VARS "working_vars"
#define MEX_OUTPUT_WORKING_CONS "working_cons"

#define MEX_PARAM_ZERO_EPS "zero_eps"
#define MEX_PARAM_EPS "eps"
#define MEX_PARAM_OBJ_LOWER "obj_lower"
#define MEX_PARAM_DERIV_PERTURBATION "deriv_pert"
#define MEX_PARAM_DERIV_TOL "deriv_tol"
#define MEX_PARAM_CAUCHY_TAU "cauchy_tau"
#define MEX_PARAM_CAUCHY_ETA "cauchy_eta"
#define MEX_PARAM_LINESEARCH_TAU "linesearch_tau"
#define MEX_PARAM_LINESEARCH_ETA "linesearch_eta"
#define MEX_PARAM_LINESEARCH_CUTOFF "linesearch_cutoff"
#define MEX_PARAM_FEASIBILITY_TOL "feas_tol"
#define MEX_PARAM_SLACKNESS_TOL "slack_tol"
#define MEX_PARAM_STATIONARITY_TOL "stat_tol"
#define MEX_PARAM_ACCEPTED_REDUCTION "accepted_reduction"
#define MEX_PARAM_DEADPOINT_BOUND "deadpoint_bound"

#define MEX_DERIV_CHECK "deriv_check"
#define MEX_HESS_EVAL "hess_eval"
#define MEX_DUAL_ESTIMATION_TYPE "dual_estimation_type"
#define MEX_BFGS_SIZING "bfgs_sizing"
#define MEX_TR_SOLVER "tr_solver"
#define MEX_POLISHING_TYPE "polishing_type"
#define MEX_STEP_RULE "step_rule"
#define MEX_LINESEARCH "linesearch"
#define MEX_PARAMETRIC_CAUCHY "parametric_cauchy"
#define MEX_INITIAL_TR_CHOICE "initial_tr_choice"

#define MEX_NUM_QUASI_NEWTON_ITERATES "num_quasi_newton_iterates"
#define MEX_MAX_NEWTON_ITERATIONS "max_newton_iterations"
#define MEX_NUM_THREADS "num_threads"

#define MEX_PERFORM_NEWTON_STEP "perform_newton_step"
#define MEX_PERFORM_SOC "perform_soc"
#define MEX_USE_QUADRATIC_MODEL "use_quadratic_model"
#define MEX_ALWAYS_WARM_START_LP "always_warm_start_lp"
#define MEX_ENABLE_RESTORATION_PHASE "enable_restoration_phase"
#define MEX_ENABLE_PREPROCESSOR "enable_preprocessor"

#endif /* SLEQP_MEX_FIELDS */
