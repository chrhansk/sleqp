#include "problem_solver.h"

#include <math.h>

#include "cmp.h"
#include "mem.h"

const double trust_region_factor = .8;

SLEQP_RETCODE
sleqp_problem_solver_create(SleqpProblemSolver** star,
                            SLEQP_SOLVER_PHASE solver_phase,
                            SleqpProblem* problem,
                            SleqpParams* params,
                            SleqpOptions* options,
                            SleqpVec* primal)
{
  SLEQP_CALL(sleqp_malloc(star));

  SleqpProblemSolver* solver = *star;

  *solver = (SleqpProblemSolver){0};

  solver->refcount = 1;

  SLEQP_CALL(sleqp_problem_capture(problem));
  solver->problem = problem;

  SLEQP_CALL(sleqp_params_capture(params));
  solver->params = params;

  SLEQP_CALL(sleqp_options_capture(options));
  solver->options = options;

  solver->solver_phase = solver_phase;

  SLEQP_CALL(sleqp_measure_create(&solver->measure, problem, params));

  const int num_variables   = sleqp_problem_num_vars(problem);
  const int num_constraints = sleqp_problem_num_cons(problem);

  SLEQP_CALL(sleqp_alloc_array(&solver->dense_cache,
                               SLEQP_MAX(num_variables, num_constraints)));

  SLEQP_CALL(sleqp_vec_create_empty(&solver->primal_diff, num_variables));

  SLEQP_CALL(sleqp_vec_create_empty(&solver->cons_dual_diff, num_constraints));

  SLEQP_CALL(sleqp_vec_create_empty(&solver->vars_dual_diff, num_variables));

  SLEQP_CALL(sleqp_iterate_create(&solver->iterate, solver->problem, primal));

  const double zero_eps
    = sleqp_params_value(solver->params, SLEQP_PARAM_ZERO_EPS);

  SLEQP_CALL(sleqp_vec_clip(primal,
                            sleqp_problem_vars_lb(solver->problem),
                            sleqp_problem_vars_ub(solver->problem),
                            zero_eps,
                            sleqp_iterate_primal(solver->iterate)));

  SLEQP_CALL(sleqp_iterate_create(&solver->trial_iterate,
                                  solver->problem,
                                  sleqp_iterate_primal(solver->iterate)));

  SLEQP_CALL(sleqp_timer_create(&solver->elapsed_timer));

  SLEQP_CALL(sleqp_trial_point_solver_create(&solver->trial_point_solver,
                                             problem,
                                             params,
                                             options));

  SLEQP_CALL(sleqp_step_rule_create_default(&solver->step_rule,
                                            problem,
                                            params,
                                            options));

  SLEQP_CALL(sleqp_deriv_checker_create(&solver->deriv_checker,
                                        solver->problem,
                                        params));

  SLEQP_CALL(sleqp_merit_create(&solver->merit, solver->problem, params));

  for (int i = 0; i < SLEQP_PROBLEM_SOLVER_NUM_EVENTS; ++i)
  {
    SLEQP_CALL(sleqp_callback_handler_create(solver->callback_handlers + i));
  }

  SLEQP_CALL(sleqp_problem_solver_reset(solver));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_problem_solver_set_primal(SleqpProblemSolver* solver,
                                const SleqpVec* primal)
{
  const int num_variables = sleqp_problem_num_vars(solver->problem);

  assert(primal->dim == num_variables);

  SLEQP_CALL(sleqp_vec_copy(primal, sleqp_iterate_primal(solver->iterate)));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
compute_initial_trust_radius(SleqpProblemSolver* solver)
{
  const SLEQP_INITIAL_TR_CHOICE initial_tr_choice
    = sleqp_options_enum_value(solver->options,
                               SLEQP_OPTION_ENUM_INITIAL_TR_CHOICE);

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

  // penalty parameter as suggested:
  solver->penalty_parameter = 10.;

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
sleqp_problem_solver_get_elapsed_seconds(const SleqpProblemSolver* solver)
{
  return sleqp_timer_get_ttl(solver->elapsed_timer);
}

SleqpIterate*
sleqp_problem_solver_get_iterate(const SleqpProblemSolver* solver)
{
  return solver->iterate;
}

SLEQP_PROBLEM_SOLVER_STATUS
sleqp_problem_solver_get_status(const SleqpProblemSolver* solver)
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

  SLEQP_CALL(sleqp_options_release(&solver->options));

  SLEQP_CALL(sleqp_params_release(&solver->params));

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
