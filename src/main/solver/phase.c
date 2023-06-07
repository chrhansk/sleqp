#include "solver.h"

#include "feas.h"
#include "restoration.h"

static SLEQP_RETCODE
on_restoration_solver_accepted_iterate(SleqpProblemSolver* problem_solver,
                                       SleqpIterate* iterate,
                                       SleqpIterate* trial_iterate,
                                       void* callback_data)
{
  SleqpSolver* solver = (SleqpSolver*)callback_data;

  assert(problem_solver == solver->restoration_problem_solver);

  SleqpProblem* problem             = solver->problem;
  SleqpProblem* restoration_problem = solver->restoration_problem;
  SleqpFunc* restoration_func       = sleqp_problem_func(restoration_problem);

  const double feas_eps
    = sleqp_settings_real_value(solver->settings, SLEQP_SETTINGS_REAL_FEAS_TOL);

  SleqpVec* cons_val;

  SLEQP_CALL(sleqp_restoration_func_cons_val(restoration_func, &cons_val));

  double feas_res;

  SLEQP_CALL(sleqp_max_violation(problem, cons_val, &feas_res));

  if (feas_res <= feas_eps)
  {
    SLEQP_CALL(sleqp_problem_solver_abort(problem_solver));
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
create_restoration_primal(SleqpSolver* solver)
{
  SleqpIterate* iterate = sleqp_problem_solver_iterate(solver->problem_solver);

  SLEQP_CALL(
    sleqp_restoration_problem_transform(solver->problem,
                                        sleqp_iterate_primal(iterate),
                                        sleqp_iterate_cons_val(iterate),
                                        solver->restoration_primal));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
create_primal(SleqpSolver* solver)
{
  SleqpIterate* restoration_iterate
    = sleqp_problem_solver_iterate(solver->restoration_problem_solver);

  SLEQP_CALL(
    sleqp_restoration_problem_restore(solver->problem,
                                      sleqp_iterate_primal(restoration_iterate),
                                      solver->primal));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
create_restoration_solver(SleqpSolver* solver)
{
  if (solver->restoration_problem_solver)
  {
    return SLEQP_OKAY;
  }

  SLEQP_CALL(sleqp_restoration_problem_create(&solver->restoration_problem,
                                              solver->settings,
                                              solver->problem));

  const int num_variables = sleqp_problem_num_vars(solver->restoration_problem);

  SLEQP_CALL(
    sleqp_vec_create_empty(&solver->restoration_primal, num_variables));

  SLEQP_CALL(sleqp_problem_solver_create(&solver->restoration_problem_solver,
                                         SLEQP_SOLVER_PHASE_RESTORATION,
                                         solver->restoration_problem,
                                         solver->settings));

  SLEQP_CALL(sleqp_problem_solver_add_callback(
    solver->restoration_problem_solver,
    SLEQP_PROBLEM_SOLVER_EVENT_ACCEPTED_ITERATE,
    on_restoration_solver_accepted_iterate,
    (void*)solver));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
fill_optimization_iterate(SleqpSolver* solver)
{
  SLEQP_CALL(create_primal(solver));

  SleqpProblemSolver* problem_solver = solver->problem_solver;

  SleqpIterate* opt_iterate = sleqp_problem_solver_iterate(problem_solver);

  SleqpVec* primal = solver->primal;

  SLEQP_CALL(sleqp_vec_copy(primal, sleqp_iterate_primal(opt_iterate)));

  SleqpProblem* problem = solver->problem;
  SleqpFunc* func       = sleqp_problem_func(problem);

  SleqpProblem* rest_problem = solver->restoration_problem;
  SleqpFunc* rest_func       = sleqp_problem_func(rest_problem);

  {
    SleqpVec* cons_val;
    SleqpMat* cons_jac;

    SLEQP_CALL(sleqp_restoration_func_cons_val(rest_func, &cons_val));
    SLEQP_CALL(sleqp_restoration_func_cons_jac(rest_func, &cons_jac));

    SLEQP_CALL(sleqp_vec_copy(cons_val, sleqp_iterate_cons_val(opt_iterate)));

    SLEQP_CALL(sleqp_mat_copy(cons_jac, sleqp_iterate_cons_jac(opt_iterate)));
  }

  {
    double obj_val;

    SLEQP_CALL(sleqp_func_obj_val(func, &obj_val));

    SLEQP_CALL(sleqp_iterate_set_obj_val(opt_iterate, obj_val));
  }

  {
    SleqpVec* obj_grad = sleqp_iterate_obj_grad(opt_iterate);

    SLEQP_CALL(sleqp_func_obj_grad(func, obj_grad));
  }

  return SLEQP_OKAY;
}

// Assuming that we are switching to the restoration phase
// based on a *full* iterate (in particular, containing
// cons vals / cons Jac). Initialize iterate
// and restoration func based on these values
static SLEQP_RETCODE
fill_restoration_iterate(SleqpSolver* solver)
{
  SLEQP_CALL(create_restoration_primal(solver));

  SleqpProblemSolver* problem_solver = solver->problem_solver;
  SleqpProblemSolver* rest_solver    = solver->restoration_problem_solver;

  SleqpIterate* opt_iterate  = sleqp_problem_solver_iterate(problem_solver);
  SleqpIterate* rest_iterate = sleqp_problem_solver_iterate(rest_solver);

  SleqpVec* rest_primal = solver->restoration_primal;

  SLEQP_CALL(sleqp_vec_copy(rest_primal, sleqp_iterate_primal(rest_iterate)));

  SleqpProblem* rest_problem = solver->restoration_problem;
  SleqpFunc* rest_func       = sleqp_problem_func(rest_problem);

  SleqpVec* orig_cons_val = sleqp_iterate_cons_val(opt_iterate);
  SleqpMat* orig_cons_jac = sleqp_iterate_cons_jac(opt_iterate);

  SLEQP_CALL(sleqp_restoration_func_init(rest_func,
                                         rest_primal,
                                         orig_cons_val,
                                         orig_cons_jac));

  {
    double obj_val;

    SLEQP_CALL(sleqp_func_obj_val(rest_func, &obj_val));

    SLEQP_CALL(sleqp_iterate_set_obj_val(rest_iterate, obj_val));
  }

  {
    SleqpVec* obj_grad = sleqp_iterate_obj_grad(rest_iterate);

    SLEQP_CALL(sleqp_func_obj_grad(rest_func, obj_grad));
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_solver_toggle_phase(SleqpSolver* solver)
{
  SLEQP_SOLVER_PHASE current_phase = solver->solver_phase;

  SLEQP_CALL(create_restoration_solver(solver));

  if (current_phase == SLEQP_SOLVER_PHASE_OPTIMIZATION)
  {
    sleqp_log_info("Switching to restoration phase");

    SLEQP_CALL(create_restoration_primal(solver));

    SLEQP_CALL(fill_restoration_iterate(solver));

    solver->solver_phase = SLEQP_SOLVER_PHASE_RESTORATION;
  }
  else
  {
    sleqp_log_info("Switching to optimization phase");

    SLEQP_CALL(create_primal(solver));

    SLEQP_CALL(fill_optimization_iterate(solver));

    solver->solver_phase = SLEQP_SOLVER_PHASE_OPTIMIZATION;
  }

  return SLEQP_OKAY;
}
