#include "problem_solver.h"

#include <math.h>

#include "cmp.h"
#include "dyn.h"
#include "mem.h"

// penalty parameter as suggested:
const double penalty_parameter_default = 10.;
const double trust_region_factor       = .8;

SLEQP_RETCODE
sleqp_problem_solver_create(SleqpProblemSolver** star,
                            SLEQP_SOLVER_PHASE solver_phase,
                            SleqpProblem* problem,
                            SleqpSettings* settings)
{
  SLEQP_CALL(sleqp_malloc(star));

  SleqpProblemSolver* solver = *star;

  *solver = (SleqpProblemSolver){0};

  solver->refcount = 1;

  SLEQP_CALL(sleqp_problem_capture(problem));
  solver->problem = problem;

  SLEQP_CALL(sleqp_settings_capture(settings));
  solver->settings = settings;

  solver->solver_phase = solver_phase;

  SLEQP_CALL(sleqp_measure_create(&solver->measure, problem, settings));

  const int num_vars = sleqp_problem_num_vars(problem);
  const int num_cons = sleqp_problem_num_cons(problem);

  SLEQP_CALL(
    sleqp_alloc_array(&solver->dense_cache, SLEQP_MAX(num_vars, num_cons)));

  SLEQP_CALL(sleqp_vec_create_empty(&solver->primal_diff, num_vars));

  SLEQP_CALL(sleqp_vec_create_empty(&solver->cons_dual_diff, num_cons));

  SLEQP_CALL(sleqp_vec_create_empty(&solver->vars_dual_diff, num_vars));

  SleqpVec* var_lb = sleqp_problem_vars_lb(solver->problem);

  SLEQP_CALL(sleqp_iterate_create(&solver->iterate, solver->problem, var_lb));

  SLEQP_CALL(sleqp_iterate_create(&solver->trial_iterate,
                                  solver->problem,
                                  sleqp_iterate_primal(solver->iterate)));

  SLEQP_CALL(sleqp_timer_create(&solver->elapsed_timer));

  SLEQP_CALL(sleqp_trial_point_solver_create(&solver->trial_point_solver,
                                             problem,
                                             settings));

  SLEQP_CALL(
    sleqp_step_rule_create_default(&solver->step_rule, problem, settings));

  SLEQP_CALL(sleqp_deriv_checker_create(&solver->deriv_checker,
                                        solver->problem,
                                        settings));

  SLEQP_CALL(sleqp_merit_create(&solver->merit, solver->problem, settings));

  for (int i = 0; i < SLEQP_PROBLEM_SOLVER_NUM_EVENTS; ++i)
  {
    SLEQP_CALL(sleqp_callback_handler_create(solver->callback_handlers + i));
  }

  SLEQP_CALL(sleqp_problem_solver_reset(solver));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
compute_initial_trust_radius(SleqpProblemSolver* solver)
{
  const SLEQP_INITIAL_TR_CHOICE initial_tr_choice
    = sleqp_settings_enum_value(solver->settings,
                                SLEQP_SETTINGS_ENUM_INITIAL_TR_CHOICE);

  const int num_vars         = sleqp_problem_num_vars(solver->problem);
  const double sqrt_num_vars = sqrt((double)num_vars);

  switch (initial_tr_choice)
  {
  case SLEQP_INITIAL_TR_CHOICE_NARROW:
    // suggested in the original paper
    solver->trust_radius = 1.;
    solver->lp_trust_radius
      = trust_region_factor * solver->trust_radius / sqrt_num_vars;
    break;
  case SLEQP_INITIAL_TR_CHOICE_WIDE:
    // knitro default
    solver->trust_radius    = sqrt_num_vars;
    solver->lp_trust_radius = trust_region_factor;
    break;
  }

  assert(solver->lp_trust_radius < solver->trust_radius);

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_problem_solver_reset(SleqpProblemSolver* solver)
{
  SLEQP_CALL(compute_initial_trust_radius(solver));

  solver->penalty_parameter = penalty_parameter_default;

  solver->num_feasible_steps        = 0;
  solver->num_global_penalty_resets = 0;

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_problem_solver_abort(SleqpProblemSolver* solver)
{
  solver->abort_next = true;

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_problem_solver_set_iteration(SleqpProblemSolver* solver, int iteration)
{
  solver->iteration = iteration;

  return SLEQP_OKAY;
}

int
sleqp_problem_solver_elapsed_iterations(const SleqpProblemSolver* solver)
{
  return solver->elapsed_iterations;
}

double
sleqp_problem_solver_elapsed_seconds(const SleqpProblemSolver* solver)
{
  return sleqp_timer_get_ttl(solver->elapsed_timer);
}

SLEQP_RETCODE
sleqp_problem_solver_set_cons_weights(SleqpProblemSolver* solver)
{
  SleqpProblem* problem = solver->problem;
  SleqpFunc* func       = sleqp_problem_func(problem);

  SLEQP_CALL(sleqp_dyn_set_penalty_cons_weights(func,
                                                solver->penalty_parameter,
                                                solver->dense_cache));

  return SLEQP_OKAY;
}

SleqpIterate*
sleqp_problem_solver_iterate(const SleqpProblemSolver* solver)
{
  return solver->iterate;
}

SLEQP_PROBLEM_SOLVER_STATUS
sleqp_problem_solver_status(const SleqpProblemSolver* solver)
{
  return solver->status;
}

static SLEQP_RETCODE
problem_solver_free(SleqpProblemSolver** star)
{
  SleqpProblemSolver* solver = *star;

  for (int i = 0; i < SLEQP_PROBLEM_SOLVER_NUM_EVENTS; ++i)
  {
    SLEQP_CALL(sleqp_callback_handler_release(solver->callback_handlers + i));
  }

  SLEQP_CALL(sleqp_merit_release(&solver->merit));

  SLEQP_CALL(sleqp_deriv_checker_free(&solver->deriv_checker));

  SLEQP_CALL(sleqp_step_rule_release(&solver->step_rule));

  SLEQP_CALL(sleqp_trial_point_solver_release(&solver->trial_point_solver));

  SLEQP_CALL(sleqp_timer_free(&solver->elapsed_timer));

  SLEQP_CALL(sleqp_iterate_release(&solver->trial_iterate));
  SLEQP_CALL(sleqp_iterate_release(&solver->iterate));

  SLEQP_CALL(sleqp_vec_free(&solver->vars_dual_diff));

  SLEQP_CALL(sleqp_vec_free(&solver->cons_dual_diff));

  SLEQP_CALL(sleqp_vec_free(&solver->primal_diff));

  sleqp_free(&solver->dense_cache);

  SLEQP_CALL(sleqp_measure_release(&solver->measure));

  SLEQP_CALL(sleqp_settings_release(&solver->settings));

  SLEQP_CALL(sleqp_problem_release(&solver->problem));

  sleqp_free(&solver);

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_problem_solver_capture(SleqpProblemSolver* solver)
{
  ++solver->refcount;

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_problem_solver_release(SleqpProblemSolver** star)
{
  SleqpProblemSolver* solver = *star;

  if (!solver)
  {
    return SLEQP_OKAY;
  }

  if (--solver->refcount == 0)
  {
    SLEQP_CALL(problem_solver_free(star));
  }

  *star = NULL;

  return SLEQP_OKAY;
}
