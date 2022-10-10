#include "measure.h"

#include <math.h>

#include "cmp.h"
#include "direction.h"
#include "fail.h"
#include "feas.h"
#include "mem.h"
#include "util.h"

#include "solver.h"
#include "sparse/sparse_matrix.h"

const double nonlinearity_cutoff = 1e-6;

struct SleqpMeasure
{
  int refcount;

  SleqpProblem* problem;
  SleqpParams* params;

  double penalty_parameter;

  double direction_norm_sq;

  double current_obj_val;
  double linear_obj_val;
  double hess_prod;
  double trial_obj_val;

  SleqpVec* expected_cons_val;

  double current_violation;
  double expected_violation;
  double trial_violation;

  double objective_nonlin;
  SleqpVec* cons_nonlin;
};

SLEQP_RETCODE
sleqp_measure_create(SleqpMeasure** star,
                     SleqpProblem* problem,
                     SleqpParams* params)
{
  SLEQP_CALL(sleqp_malloc(star));

  SleqpMeasure* measure = *star;

  *measure = (SleqpMeasure){0};

  measure->refcount = 1;

  SLEQP_CALL(sleqp_problem_capture(problem));
  measure->problem = problem;

  SLEQP_CALL(sleqp_params_capture(params));
  measure->params = params;

  const int num_cons = sleqp_problem_num_cons(problem);

  SLEQP_CALL(sleqp_vec_create_empty(&measure->expected_cons_val, num_cons));

  SLEQP_CALL(sleqp_vec_create_empty(&measure->cons_nonlin, num_cons));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
compute_obj_nonlin(SleqpMeasure* measure,
                   SleqpIterate* iterate,
                   SleqpIterate* trial_iterate,
                   SleqpDirection* direction)
{
  double objective_dot = *sleqp_direction_obj_grad(direction);

  measure->current_obj_val = sleqp_iterate_obj_val(iterate);

  measure->linear_obj_val = measure->current_obj_val + objective_dot;

  measure->trial_obj_val = sleqp_iterate_obj_val(trial_iterate);

  const double difference = (measure->linear_obj_val - measure->trial_obj_val);

  measure->objective_nonlin = difference * (2. / measure->direction_norm_sq);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
compute_lag_nonlin(SleqpMeasure* measure,
                   const SleqpVec* multipliers,
                   double* lag_nonlinearity)
{
  double lag_cons_nonlin;

  SLEQP_CALL(
    sleqp_vec_dot(measure->cons_nonlin, multipliers, &lag_cons_nonlin));

  *lag_nonlinearity = measure->objective_nonlin + lag_cons_nonlin;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
compute_cons_nonlin(SleqpMeasure* measure,
                    SleqpIterate* iterate,
                    SleqpIterate* trial_iterate,
                    SleqpDirection* direction)
{
  SleqpProblem* problem = measure->problem;

  SleqpVec* cons_nonlin = measure->cons_nonlin;

  const double zero_eps
    = sleqp_params_value(measure->params, SLEQP_PARAM_ZERO_EPS);

  SleqpVec* direction_cons_jac = sleqp_direction_cons_jac(direction);

  SLEQP_CALL(sleqp_vec_add(sleqp_iterate_cons_val(iterate),
                           direction_cons_jac,
                           zero_eps,
                           measure->expected_cons_val));

  SLEQP_CALL(sleqp_vec_add_scaled(measure->expected_cons_val,
                                  sleqp_iterate_cons_val(trial_iterate),
                                  1.,
                                  -1.,
                                  zero_eps,
                                  cons_nonlin));

  SLEQP_CALL(sleqp_vec_scale(cons_nonlin, (2. / measure->direction_norm_sq)));

  SLEQP_CALL(sleqp_total_violation(problem,
                                   sleqp_iterate_cons_val(iterate),
                                   &(measure->current_violation)));

  SLEQP_CALL(sleqp_total_violation(problem,
                                   measure->expected_cons_val,
                                   &(measure->expected_violation)));

  SLEQP_CALL(sleqp_total_violation(problem,
                                   sleqp_iterate_cons_val(trial_iterate),
                                   &(measure->trial_violation)));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_measure_set_iterates(SleqpMeasure* measure,
                           SleqpIterate* iterate,
                           SleqpIterate* trial_iterate,
                           SleqpDirection* direction)
{
  SleqpVec* direction_primal = sleqp_direction_primal(direction);

  measure->direction_norm_sq = sleqp_vec_norm_sq(direction_primal);

  SLEQP_CALL(sleqp_vec_dot(direction_primal,
                           sleqp_direction_hess(direction),
                           &measure->hess_prod));

  SLEQP_CALL(compute_obj_nonlin(measure, iterate, trial_iterate, direction));

  SLEQP_CALL(compute_cons_nonlin(measure, iterate, trial_iterate, direction));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_measure_set_penalty_parameter(SleqpMeasure* measure,
                                    double penalty_parameter)
{
  assert(penalty_parameter > 0.);

  measure->penalty_parameter = penalty_parameter;
  return SLEQP_OKAY;
}

SleqpVec*
sleqp_measure_cons_nonlin(SleqpMeasure* measure)
{
  return measure->cons_nonlin;
}

double
sleqp_measure_obj_nonlin(SleqpMeasure* measure)
{
  return measure->objective_nonlin;
}

double
sleqp_measure_step_norm(SleqpMeasure* measure)
{
  return sqrt(measure->direction_norm_sq);
}

SLEQP_RETCODE
sleqp_measure_lag_nonlin(SleqpMeasure* measure,
                         const SleqpVec* multipliers,
                         double* lag_nonlinearity)
{
  SLEQP_CALL(compute_lag_nonlin(measure, multipliers, lag_nonlinearity));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_measure_capture(SleqpMeasure* measure)
{
  ++measure->refcount;

  return SLEQP_OKAY;
}

static double
percent_reduction(double current, double trial)
{
  if (current == 0.)
  {
    return 0.;
  }

  double value = 100. * (current - trial) / current;

  if (current < 0.)
  {
    return -value;
  }

  return value;
}

static SLEQP_RETCODE
report_merit(SleqpMeasure* measure,
             double cur_merit,
             double exp_merit,
             double act_merit)
{
  const double model_red  = (cur_merit - exp_merit);
  const double actual_red = (cur_merit - act_merit);

  sleqp_log_debug("Merit: current: %14e, expected: %14e, actual: %14e, "
                  "prediced reduction: %9.4f%%, actual reduction: %9.4f%%, "
                  "reduction ratio: %e",
                  cur_merit,
                  exp_merit,
                  act_merit,
                  percent_reduction(cur_merit, exp_merit),
                  percent_reduction(cur_merit, act_merit),
                  sleqp_reduction_ratio(actual_red, model_red));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
report_lsq(SleqpMeasure* measure)
{
  const double cur_obj_val = measure->current_obj_val;
  const double exp_obj_val = measure->linear_obj_val + 0.5 * measure->hess_prod;
  const double act_obj_val = measure->trial_obj_val;

  sleqp_log_debug("Objective: current: %14e, expected: %14e, actual: %14e, "
                  "prediced reduction: %9.4f%%, actual reduction: %9.4f%%",
                  cur_obj_val,
                  exp_obj_val,
                  act_obj_val,
                  percent_reduction(cur_obj_val, exp_obj_val),
                  percent_reduction(cur_obj_val, act_obj_val));

  const double cur_violation = measure->current_violation;
  const double exp_violation = measure->expected_violation;
  const double act_violation = measure->trial_violation;

  sleqp_log_debug("Violation: current: %14e, expected: %14e, actual: %14e, "
                  "prediced reduction: %9.4f%%, actual reduction: %9.4f%%",
                  cur_violation,
                  exp_violation,
                  act_violation,
                  percent_reduction(cur_violation, exp_violation),
                  percent_reduction(cur_violation, act_violation));

  const double penalty   = measure->penalty_parameter;
  const double cur_merit = cur_obj_val + (penalty * cur_violation);
  const double exp_merit = exp_obj_val + (penalty * exp_violation);
  const double act_merit = act_obj_val + (penalty * act_violation);

  SLEQP_CALL(report_merit(measure, cur_merit, exp_merit, act_merit));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
report_nonlinearity(SleqpMeasure* measure, const SleqpVec* multipliers)
{
  double obj_nonlin = sleqp_measure_obj_nonlin(measure);

  SleqpVec* cons_nonlin = sleqp_measure_cons_nonlin(measure);

  double lag_nonlin;

  SLEQP_CALL(sleqp_measure_lag_nonlin(measure, multipliers, &lag_nonlin));

  const double max_cons_nonlin = sleqp_vec_inf_norm(cons_nonlin);

  const double step_norm = sleqp_measure_step_norm(measure);

  sleqp_log_debug("Objective nonlinearity: %g, "
                  "maximal constraint nonlinearity: %g, "
                  "Lagrangean nonlinearity: %g (step norm: %g)",
                  obj_nonlin,
                  max_cons_nonlin,
                  lag_nonlin,
                  step_norm);

  return SLEQP_OKAY;
}

/**
 * Other things to report (requiring additional Hessian products):
 *
 * - Quadratic vs actual objective
 * - Violation + (Hessian with penalty multipliers) vs actual merit
 * - Expected vs actual Lagrangean (without penalty multipliers)
 **/
static SLEQP_RETCODE
report_func(SleqpMeasure* measure)
{
  const double penalty = measure->penalty_parameter;

  const double cur_obj_val   = measure->current_obj_val;
  const double cur_violation = measure->current_violation;
  const double cur_merit     = cur_obj_val + (penalty * cur_violation);

  const double hess_prod     = measure->hess_prod;
  const double exp_obj_val   = measure->linear_obj_val;
  const double exp_violation = measure->expected_violation;

  double exp_merit = exp_obj_val + (penalty * exp_violation);
  exp_merit += 0.5 * hess_prod;

  const double act_obj_val   = measure->trial_obj_val;
  const double act_violation = measure->trial_violation;
  const double act_merit     = act_obj_val + (penalty * act_violation);

  SLEQP_CALL(report_merit(measure, cur_merit, exp_merit, act_merit));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_measure_report_trial_point(SleqpMeasure* measure,
                                 const SleqpVec* multipliers)
{
  SleqpProblem* problem = measure->problem;
  SleqpFunc* func       = sleqp_problem_func(problem);

  const double step_norm = sleqp_measure_step_norm(measure);

  if (step_norm > nonlinearity_cutoff)
  {
    SLEQP_CALL(report_nonlinearity(measure, multipliers));
  }
  else
  {
    sleqp_log_debug("Trial step norm: %e", sleqp_measure_step_norm(measure));
  }

  if (sleqp_func_get_type(func) == SLEQP_FUNC_TYPE_LSQ)
  {
    SLEQP_CALL(report_lsq(measure));
    return SLEQP_OKAY;
  }
  else
  {
    SLEQP_CALL(report_func(measure));
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_measure_report_soc_trial_point(SleqpMeasure* measure,
                                     SleqpIterate* soc_iterate)
{
  SleqpProblem* problem = measure->problem;

  const double cur_obj_val   = measure->current_obj_val;
  const double trial_obj_val = measure->trial_obj_val;
  const double soc_obj_val   = sleqp_iterate_obj_val(soc_iterate);

  sleqp_log_debug("Objective: current: %14e, trial: %14e, SOC: %14e, trial "
                  "reduction: %g%%, SOC reduction: %g%%",
                  cur_obj_val,
                  trial_obj_val,
                  soc_obj_val,
                  percent_reduction(cur_obj_val, trial_obj_val),
                  percent_reduction(cur_obj_val, soc_obj_val));

  const double cur_violation   = measure->current_violation;
  const double trial_violation = measure->trial_violation;

  double soc_violation;

  SLEQP_CALL(sleqp_total_violation(problem,
                                   sleqp_iterate_cons_val(soc_iterate),
                                   &(soc_violation)));

  sleqp_log_debug("Violation: current: %14e, trial: %14e, SOC: %14e, trial "
                  "reduction: %g%%, SOC reduction: %g%%",
                  cur_violation,
                  trial_violation,
                  soc_violation,
                  percent_reduction(cur_violation, trial_violation),
                  percent_reduction(cur_violation, soc_violation));

  const double penalty     = measure->penalty_parameter;
  const double cur_merit   = cur_obj_val + (penalty * cur_violation);
  const double trial_merit = trial_obj_val + (penalty * trial_violation);
  const double soc_merit   = soc_obj_val + (penalty * soc_violation);

  sleqp_log_debug("Merit: current: %14e, trial: %14e, SOC: %14e, trial "
                  "reduction: %g%%, SOC reduction: %g%%",
                  cur_merit,
                  trial_merit,
                  soc_merit,
                  percent_reduction(cur_merit, trial_merit),
                  percent_reduction(cur_merit, soc_merit));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
measure_free(SleqpMeasure** star)
{
  SleqpMeasure* measure = *star;

  SLEQP_CALL(sleqp_vec_free(&measure->cons_nonlin));

  SLEQP_CALL(sleqp_vec_free(&measure->expected_cons_val));

  SLEQP_CALL(sleqp_params_release(&measure->params));
  SLEQP_CALL(sleqp_problem_release(&measure->problem));

  sleqp_free(star);

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_measure_release(SleqpMeasure** star)
{
  SleqpMeasure* measure = *star;

  if (!measure)
  {
    return SLEQP_OKAY;
  }

  if (--measure->refcount == 0)
  {
    SLEQP_CALL(measure_free(star));
  }

  *star = NULL;

  return SLEQP_OKAY;
}
