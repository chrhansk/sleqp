#include "solver.h"

#include <math.h>
#include <string.h>

#include "cmp.h"
#include "defs.h"
#include "dyn.h"
#include "fail.h"
#include "feas.h"
#include "func.h"
#include "iterate.h"
#include "mem.h"
#include "problem.h"
#include "pub_mem.h"
#include "pub_settings.h"
#include "scale.h"

#include "timer.h"
#include "util.h"

#include "step/step_rule.h"

#define SOLVER_INFO_BUF_SIZE 400

_Thread_local char solver_info[SOLVER_INFO_BUF_SIZE];

static SLEQP_RETCODE
solver_convert_primal(SleqpSolver* solver,
                      const SleqpVec* source,
                      SleqpVec* target)
{
  assert(source->dim == sleqp_problem_num_vars(solver->original_problem));
  assert(target->dim == sleqp_problem_num_vars(solver->problem));

  SLEQP_CALL(sleqp_vec_copy(source, solver->scaled_primal));

  if (solver->scaling_data)
  {
    SLEQP_CALL(sleqp_scale_point(solver->scaling_data, solver->scaled_primal));
  }

  if (solver->preprocessor)
  {
    SLEQP_CALL(sleqp_preprocessor_transform_primal(solver->preprocessor,
                                                   solver->scaled_primal,
                                                   target));
  }
  else
  {
    SLEQP_CALL(sleqp_vec_copy(solver->scaled_primal, target));
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
do_restore_iterate(SleqpSolver* solver,
                   const SleqpIterate* source,
                   SleqpIterate* target)
{
  if (solver->preprocessor)
  {
    SLEQP_CALL(sleqp_preprocessor_restore_iterate(solver->preprocessor,
                                                  source,
                                                  solver->scaled_iterate));
  }
  else
  {
    SLEQP_CALL(sleqp_iterate_copy(source, solver->scaled_iterate));
  }

  SLEQP_CALL(sleqp_iterate_copy(solver->scaled_iterate, target));

  SleqpProblem* problem = solver->problem;
  SleqpFunc* func       = sleqp_problem_func(problem);
  const bool lsq        = sleqp_func_get_type(func) == SLEQP_FUNC_TYPE_LSQ;

  if (solver->scaling_data)
  {
    SLEQP_CALL(sleqp_unscale_iterate(solver->scaling_data, target, lsq));
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_solver_restore_original_iterate(SleqpSolver* solver)
{
  SLEQP_CALL(
    do_restore_iterate(solver,
                       sleqp_problem_solver_iterate(solver->problem_solver),
                       solver->original_iterate));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
set_obj_weight(SleqpProblem* problem)
{
  SleqpFunc* func = sleqp_problem_func(problem);

  if (sleqp_func_get_type(func) == SLEQP_FUNC_TYPE_DYNAMIC)
  {
    SLEQP_CALL(sleqp_dyn_func_set_obj_weight(func, 1.));
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
solver_create_problem(SleqpSolver* solver, SleqpProblem* problem)
{
  SleqpProblem* scaled_problem;

  SleqpSettings* settings = solver->settings;

  solver->original_problem = problem;
  SLEQP_CALL(sleqp_problem_capture(solver->original_problem));

  if (solver->scaling_data)
  {
    SLEQP_CALL(sleqp_problem_scaling_create(&solver->problem_scaling,
                                            solver->scaling_data,
                                            problem,
                                            settings));

    SLEQP_CALL(sleqp_problem_scaling_flush(solver->problem_scaling));

    scaled_problem = sleqp_problem_scaling_get_problem(solver->problem_scaling);
  }
  else
  {
    scaled_problem = problem;
  }

  SleqpFunc* func = sleqp_problem_func(scaled_problem);

  SLEQP_CALL(
    sleqp_quasi_newton_create_default(&solver->quasi_newton, func, settings));

  if (solver->quasi_newton)
  {
    func = sleqp_quasi_newton_get_func(solver->quasi_newton);
  }

  SLEQP_CALL(sleqp_problem_create(&solver->scaled_problem,
                                  func,
                                  sleqp_problem_vars_lb(scaled_problem),
                                  sleqp_problem_vars_ub(scaled_problem),
                                  sleqp_problem_general_lb(scaled_problem),
                                  sleqp_problem_general_ub(scaled_problem),
                                  sleqp_problem_linear_coeffs(scaled_problem),
                                  sleqp_problem_linear_lb(scaled_problem),
                                  sleqp_problem_linear_ub(scaled_problem),
                                  settings));

  const bool enable_preprocesor
    = sleqp_settings_bool_value(solver->settings,
                                SLEQP_SETTINGS_BOOL_ENABLE_PREPROCESSOR);

  if (enable_preprocesor)
  {
    SLEQP_CALL(sleqp_preprocessor_create(&solver->preprocessor,
                                         solver->scaled_problem,
                                         solver->settings));

    const SLEQP_PREPROCESSING_RESULT preprocessing_result
      = sleqp_preprocessor_result(solver->preprocessor);

    if (preprocessing_result == SLEQP_PREPROCESSING_RESULT_FAILURE)
    {
      solver->problem = solver->scaled_problem;
    }
    else
    {
      solver->problem
        = sleqp_preprocessor_transformed_problem(solver->preprocessor);

      if (preprocessing_result == SLEQP_PREPROCESSING_RESULT_INFEASIBLE)
      {
        solver->status = SLEQP_STATUS_INFEASIBLE;
      }
    }
  }
  else
  {
    solver->problem = solver->scaled_problem;
  }

  SLEQP_CALL(sleqp_problem_capture(solver->problem));

  // Only needs to be set once, does not change afterwards
  SLEQP_CALL(set_obj_weight(solver->problem));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
solver_create_iterates(SleqpSolver* solver, SleqpVec* primal)
{
  SleqpProblem* original_problem = solver->original_problem;

  {
    SleqpVec* var_lb = sleqp_problem_vars_lb(original_problem);
    SleqpVec* var_ub = sleqp_problem_vars_ub(original_problem);

    if (!sleqp_vec_is_boxed(primal, var_lb, var_ub))
    {
      sleqp_log_warn("Initial solution violates variable bounds");
    }
  }

  SLEQP_CALL(solver_convert_primal(solver, primal, solver->primal));

  SLEQP_CALL(sleqp_iterate_create(&solver->scaled_iterate,
                                  solver->scaled_problem,
                                  primal));

  SLEQP_CALL(sleqp_iterate_create(&solver->original_iterate,
                                  solver->original_problem,
                                  primal));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
on_problem_solver_performed_iteration(SleqpProblemSolver* problem_solver,
                                      void* callback_data)
{
  SleqpSolver* solver = (SleqpSolver*)callback_data;

  SLEQP_CALLBACK_EVENT(solver->callback_handlers,
                       SLEQP_SOLVER_EVENT_PERFORMED_ITERATION,
                       SLEQP_PERFORMED_ITERATION,
                       solver);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
on_problem_solver_accepted_iterate(SleqpProblemSolver* problem_solver,
                                   SleqpIterate* iterate,
                                   SleqpIterate* trial_iterate,
                                   void* callback_data)
{
  SleqpSolver* solver = (SleqpSolver*)callback_data;

  if (solver->quasi_newton)
  {
    SleqpVec* multipliers = sleqp_iterate_cons_dual(iterate);

    SLEQP_CALL(sleqp_quasi_newton_push(solver->quasi_newton,
                                       iterate,
                                       trial_iterate,
                                       multipliers));
  }

  SLEQP_SOLVER_EVENT event = SLEQP_SOLVER_EVENT_ACCEPTED_ITERATE;

  SleqpCallbackHandler* handler = solver->callback_handlers[event];

  if (sleqp_callback_handler_size(handler) != 0)
  {
    SLEQP_CALL(sleqp_solver_restore_original_iterate(solver));

    SLEQP_CALLBACK_EVENT(solver->callback_handlers,
                         SLEQP_SOLVER_EVENT_ACCEPTED_ITERATE,
                         SLEQP_ACCEPTED_ITERATE,
                         solver,
                         trial_iterate);
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_solver_create(SleqpSolver** star,
                    SleqpProblem* problem,
                    SleqpVec* primal,
                    SleqpScaling* scaling_data)
{
  assert(sleqp_vec_is_valid(primal));

  SLEQP_CALL(sleqp_malloc(star));

  SleqpSolver* solver = *star;

  *solver = (SleqpSolver){0};

  solver->refcount = 1;

  SleqpSettings* settings = sleqp_problem_settings(problem);

  SLEQP_CALL(sleqp_settings_capture(settings));
  solver->settings = settings;

  const int num_orig_vars = sleqp_problem_num_vars(problem);
  const int num_orig_cons = sleqp_problem_num_cons(problem);

  SLEQP_CALL(sleqp_alloc_array(&solver->dense_cache,
                               SLEQP_MAX(num_orig_vars, num_orig_cons)));

  SLEQP_CALL(sleqp_timer_create(&solver->elapsed_timer));

  if (scaling_data)
  {
    SLEQP_CALL(sleqp_scaling_capture(scaling_data));
    solver->scaling_data = scaling_data;
  }

  SLEQP_CALL(solver_create_problem(solver, problem));

  const int num_variables = sleqp_problem_num_vars(solver->problem);

  SLEQP_CALL(sleqp_vec_create_empty(&solver->scaled_primal, num_orig_vars));

  SLEQP_CALL(sleqp_vec_create_empty(&solver->primal, num_variables));

  SLEQP_CALL(solver_create_iterates(solver, primal));

  SLEQP_CALL(sleqp_problem_solver_create(&solver->problem_solver,
                                         SLEQP_SOLVER_PHASE_OPTIMIZATION,
                                         solver->problem,
                                         settings));

  SLEQP_CALL(sleqp_problem_solver_add_callback(
    solver->problem_solver,
    SLEQP_PROBLEM_SOLVER_EVENT_ACCEPTED_ITERATE,
    on_problem_solver_accepted_iterate,
    (void*)solver));

  SLEQP_CALL(sleqp_problem_solver_add_callback(
    solver->problem_solver,
    SLEQP_PROBLEM_SOLVER_EVENT_PERFORMED_ITERATION,
    on_problem_solver_performed_iteration,
    (void*)solver));

  SLEQP_CALL(
    sleqp_polishing_create(&solver->polishing, solver->problem, settings));

  for (int i = 0; i < SLEQP_SOLVER_NUM_EVENTS; ++i)
  {
    SLEQP_CALL(sleqp_callback_handler_create(solver->callback_handlers + i));
  }

  SLEQP_CALL(sleqp_solver_reset(solver));

  solver->time_limit = SLEQP_NONE;

  solver->abort_next = false;

  sleqp_log_debug("%s", sleqp_solver_info(solver));

  return SLEQP_OKAY;
}

const char*
sleqp_solver_info(const SleqpSolver* solver)
{
  snprintf(solver_info,
           SOLVER_INFO_BUF_SIZE,
           "SLEQP version %s [LP solver: %s] [Factorization: %s] [trlib: %s]",
           SLEQP_LONG_VERSION,
           SLEQP_LP_SOLVER_NAME " " SLEQP_LP_SOLVER_VERSION,
           SLEQP_FACT_NAME " " SLEQP_FACT_VERSION,
           SLEQP_TRLIB_VERSION);

  return solver_info;
}

SLEQP_RETCODE
sleqp_solver_solution(SleqpSolver* solver, SleqpIterate** iterate)
{
  SLEQP_CALL(sleqp_solver_restore_original_iterate(solver));

  (*iterate) = solver->original_iterate;

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_solver_violated_constraints(SleqpSolver* solver,
                                  SleqpIterate* iterate,
                                  int* violated_constraints,
                                  int* num_violated_constraints)
{
  const double feas_eps
    = sleqp_settings_real_value(solver->settings, SLEQP_SETTINGS_REAL_FEAS_TOL);

  SLEQP_CALL(sleqp_iterate_get_violated_constraints(solver->original_problem,
                                                    iterate,
                                                    violated_constraints,
                                                    num_violated_constraints,
                                                    feas_eps));

  return SLEQP_OKAY;
}

SLEQP_STATUS
sleqp_solver_status(const SleqpSolver* solver)
{
  return solver->status;
}

SLEQP_RETCODE
sleqp_solver_reset(SleqpSolver* solver)
{
  SLEQP_CALL(sleqp_problem_solver_reset(solver->problem_solver));

  if (solver->restoration_problem_solver)
  {
    SLEQP_CALL(sleqp_problem_solver_reset(solver->restoration_problem_solver));
  }

  if (solver->quasi_newton)
  {
    SLEQP_CALL(sleqp_quasi_newton_reset(solver->quasi_newton));
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_solver_abort(SleqpSolver* solver)
{
  solver->abort_next = true;

  SLEQP_CALL(sleqp_problem_solver_abort(solver->problem_solver));

  return SLEQP_OKAY;
}

int
sleqp_solver_iterations(const SleqpSolver* solver)
{
  return solver->iterations;
}

double
sleqp_solver_elapsed_seconds(const SleqpSolver* solver)
{
  return sleqp_timer_get_ttl(solver->elapsed_timer);
}

static SLEQP_RETCODE
solver_free(SleqpSolver** star)
{
  SleqpSolver* solver = *star;

  if (!solver)
  {
    return SLEQP_OKAY;
  }

  for (int i = 0; i < SLEQP_SOLVER_NUM_EVENTS; ++i)
  {
    SLEQP_CALL(sleqp_callback_handler_release(solver->callback_handlers + i));
  }

  SLEQP_CALL(sleqp_polishing_release(&solver->polishing));

  SLEQP_CALL(sleqp_timer_free(&solver->elapsed_timer));

  SLEQP_CALL(sleqp_iterate_release(&solver->scaled_iterate));

  SLEQP_CALL(sleqp_iterate_release(&solver->original_iterate));

  SLEQP_CALL(sleqp_problem_release(&solver->problem));

  SLEQP_CALL(sleqp_problem_scaling_release(&solver->problem_scaling));

  SLEQP_CALL(sleqp_preprocessor_release(&solver->preprocessor));

  SLEQP_CALL(sleqp_problem_solver_release(&solver->restoration_problem_solver));

  SLEQP_CALL(sleqp_problem_release(&solver->restoration_problem));

  SLEQP_CALL(sleqp_vec_free(&solver->restoration_primal));

  SLEQP_CALL(sleqp_problem_solver_release(&solver->problem_solver));

  SLEQP_CALL(sleqp_quasi_newton_release(&solver->quasi_newton));

  SLEQP_CALL(sleqp_problem_release(&solver->scaled_problem));

  SLEQP_CALL(sleqp_vec_free(&solver->primal));

  SLEQP_CALL(sleqp_vec_free(&solver->scaled_primal));

  SLEQP_CALL(sleqp_scaling_release(&solver->scaling_data));

  sleqp_free(&solver->dense_cache);

  SLEQP_CALL(sleqp_problem_release(&solver->original_problem));

  SLEQP_CALL(sleqp_settings_release(&solver->settings));

  sleqp_free(star);

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_solver_capture(SleqpSolver* solver)
{
  ++solver->refcount;

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_solver_release(SleqpSolver** star)
{
  SleqpSolver* solver = *star;

  if (!solver)
  {
    return SLEQP_OKAY;
  }

  if (--solver->refcount == 0)
  {
    SLEQP_CALL(solver_free(star));
  }

  *star = NULL;

  return SLEQP_OKAY;
}
