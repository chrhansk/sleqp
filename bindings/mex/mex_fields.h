#ifndef SLEQP_MEX_FIELDS
#define SLEQP_MEX_FIELDS

#define MEX_MSG_IDENTIFIER "SLEQP:Error"

// #define MEX_COMMAND_INFO "info"
#define MEX_COMMAND_SOLVE "solve"
#define MEX_COMMAND_SOLVE_LSQ "solve_lsq"

#define MATLAB_FUNC_FEVAL "feval"
#define MATLAB_FUNC_DISP "disp"

#define MEX_INPUT_OBJ_VAL "obj_val"
#define MEX_INPUT_OBJ_GRAD "obj_grad"
#define MEX_INPUT_CONS_VAL "cons_val"
#define MEX_INPUT_CONS_JAC "cons_jac"
#define MEX_INPUT_HESS "hess"

#define MEX_INPUT_LSQ_RES "lsq_residuals"
#define MEX_INPUT_LSQ_JAC_FWD "lsq_jac_forward"
#define MEX_INPUT_LSQ_JAC_ADJ "lsq_jac_adjoint"

#define MEX_INPUT_CONS_LB "cons_lb"
#define MEX_INPUT_CONS_UB "cons_ub"

#define MEX_INPUT_VAR_LB "var_lb"
#define MEX_INPUT_VAR_UB "var_ub"

#define MEX_OUTPUT_PRIMAL "primal"
#define MEX_OUTPUT_CONS_DUAL "cons_dual"
#define MEX_OUTPUT_VARS_DUAL "vars_dal"
#define MEX_OUTPUT_ELAPSED "elapsed"
#define MEX_OUTPUT_ITER "iterations"
#define MEX_OUTPUT_STATUS "status"

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

#endif /* SLEQP_MEX_FIELDS */
