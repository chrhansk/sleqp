#include "problem_solver.h"

#include "error.h"
#include "feas.h"
#include "solver.h"

SLEQP_RETCODE
sleqp_problem_solver_get_real_state(const SleqpProblemSolver* solver,
                                    SLEQP_SOLVER_STATE_REAL state,
                                    double* value)
{
  SleqpTrialPointSolver* trial_point_solver = solver->trial_point_solver;

  double min_rayleigh, max_rayleigh;

  SLEQP_CALL(sleqp_trial_point_solver_rayleigh(trial_point_solver,
                                               &min_rayleigh,
                                               &max_rayleigh));

  switch (state)
  {
  case SLEQP_SOLVER_STATE_REAL_TRUST_RADIUS:
    (*value) = solver->trust_radius;
    break;
  case SLEQP_SOLVER_STATE_REAL_LP_TRUST_RADIUS:
    (*value) = solver->lp_trust_radius;
    break;
  case SLEQP_SOLVER_STATE_REAL_SCALED_OBJ_VAL:
    (*value) = sleqp_iterate_obj_val(solver->iterate);
    break;
  case SLEQP_SOLVER_STATE_REAL_SCALED_MERIT_VAL:
    (*value) = solver->current_merit_value;
    break;
  case SLEQP_SOLVER_STATE_REAL_SCALED_FEAS_RES:
    (*value) = solver->feas_res;
    break;
  case SLEQP_SOLVER_STATE_REAL_SCALED_STAT_RES:
    (*value) = solver->stat_res;
    break;
  case SLEQP_SOLVER_STATE_REAL_SCALED_SLACK_RES:
    (*value) = solver->slack_res;
    break;
  case SLEQP_SOLVER_STATE_REAL_PENALTY_PARAM:
    (*value) = solver->penalty_parameter;
    break;
  case SLEQP_SOLVER_STATE_REAL_MIN_RAYLEIGH:
    (*value) = min_rayleigh;
    break;
  case SLEQP_SOLVER_STATE_REAL_MAX_RAYLEIGH:
    (*value) = max_rayleigh;
    break;
  default:
    sleqp_raise(SLEQP_ILLEGAL_ARGUMENT, "Invalid state value (%d)", state);
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_problem_solver_get_int_state(const SleqpProblemSolver* solver,
                                   SLEQP_SOLVER_STATE_INT state,
                                   int* value)
{
  switch (state)
  {
  case SLEQP_SOLVER_STATE_INT_LAST_STEP_ON_BDRY:
    (*value) = solver->boundary_step;
    break;
  case SLEQP_SOLVER_STATE_INT_ITERATION:
    (*value) = solver->iteration;
    break;
  case SLEQP_SOLVER_STATE_INT_LAST_STEP_TYPE:
    (*value) = solver->last_step_type;
    break;
  default:
    sleqp_raise(SLEQP_ILLEGAL_ARGUMENT, "Invalid state value (%d)", state);
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_problem_solver_get_vec_state(const SleqpProblemSolver* solver,
                                   SLEQP_SOLVER_STATE_VEC value,
                                   SleqpVec* result)
{
  const double zero_eps
    = sleqp_params_value(solver->params, SLEQP_PARAM_ZERO_EPS);

  SleqpVec* cons_val = sleqp_iterate_cons_val(solver->iterate);

  switch (value)
  {
  case SLEQP_SOLVER_STATE_VEC_SCALED_STAT_RESIDUALS:
    SLEQP_CALL(sleqp_iterate_stationarity_residuals(solver->problem,
                                                    solver->iterate,
                                                    solver->dense_cache,
                                                    result,
                                                    zero_eps));
    break;
  case SLEQP_SOLVER_STATE_VEC_SCALED_FEAS_RESIDUALS:
    SLEQP_CALL(sleqp_violation_values(solver->problem, cons_val, result));
    break;
  case SLEQP_SOLVER_STATE_VEC_SCALED_CONS_SLACK_RESIDUALS:
    SLEQP_CALL(sleqp_iterate_cons_slackness_residuals(solver->problem,
                                                      solver->iterate,
                                                      result,
                                                      zero_eps));
    break;
  case SLEQP_SOLVER_STATE_VEC_SCALED_VAR_SLACK_RESIDUALS:
    SLEQP_CALL(sleqp_iterate_vars_slackness_residuals(solver->problem,
                                                      solver->iterate,
                                                      result,
                                                      zero_eps));
    break;
  default:
    sleqp_raise(SLEQP_ILLEGAL_ARGUMENT, "Invalid state value (%d)", value);
  }

  return SLEQP_OKAY;
}
