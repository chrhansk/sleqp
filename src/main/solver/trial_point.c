#include "solver.h"

#include "fail.h"

SLEQP_RETCODE sleqp_solver_compute_trial_point_simple(SleqpSolver* solver,
                                                      double* cauchy_merit_value,
                                                      bool quadratic_model,
                                                      bool* full_step)
{
  SleqpProblem* problem = solver->problem;
  SleqpIterate* iterate = solver->iterate;

  const double eps = sleqp_params_get(solver->params,
                                      SLEQP_PARAM_EPS);

  SLEQP_NUM_ASSERT_PARAM(eps);

  const double zero_eps = sleqp_params_get(solver->params,
                                           SLEQP_PARAM_ZERO_EPS);

  SLEQP_CALL(sleqp_solver_compute_cauchy_step(solver,
                                              cauchy_merit_value,
                                              quadratic_model,
                                              full_step));

  const SleqpSparseVec* trial_step = solver->cauchy_step;

  // Compute quadratic merit value
  {
    SLEQP_CALL(sleqp_merit_linear(solver->merit_data,
                                  solver->iterate,
                                  trial_step,
                                  solver->penalty_parameter,
                                  cauchy_merit_value));
  }

  if(quadratic_model)
  {
    double hessian_prod;

    SLEQP_CALL(sleqp_sparse_vector_dot(trial_step,
                                       solver->cauchy_hessian_step,
                                       &hessian_prod));

    (*cauchy_merit_value) += .5 * hessian_prod;

#if !defined(NDEBUG)

    {
      double actual_quadratic_merit_value;

      double func_dual = 1.;

      SLEQP_CALL(sleqp_merit_quadratic(solver->merit_data,
                                       iterate,
                                       &func_dual,
                                       solver->cauchy_step,
                                       solver->multipliers,
                                       solver->penalty_parameter,
                                       &actual_quadratic_merit_value));

      sleqp_assert_is_eq(*cauchy_merit_value,
                         actual_quadratic_merit_value,
                         eps);
    }

#endif

  }

  SLEQP_CALL(sleqp_sparse_vector_add(sleqp_iterate_get_primal(iterate),
                                     solver->cauchy_step,
                                     zero_eps,
                                     solver->initial_trial_point));

  SLEQP_CALL(sleqp_sparse_vector_clip(solver->initial_trial_point,
                                      sleqp_problem_var_lb(problem),
                                      sleqp_problem_var_ub(problem),
                                      zero_eps,
                                      sleqp_iterate_get_primal(solver->trial_iterate)));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_solver_compute_trial_point_newton(SleqpSolver* solver,
                                                      double* trial_merit_value,
                                                      bool* full_step)
{
  SleqpProblem* problem = solver->problem;
  SleqpIterate* iterate = solver->iterate;

  const double eps = sleqp_params_get(solver->params, SLEQP_PARAM_EPS);

  SLEQP_NUM_ASSERT_PARAM(eps);

  const double zero_eps = sleqp_params_get(solver->params,
                                           SLEQP_PARAM_ZERO_EPS);

  const double one = 1.;

  double cauchy_merit_value;

  SLEQP_CALL(sleqp_solver_compute_cauchy_step(solver,
                                              &cauchy_merit_value,
                                              true,
                                              full_step));

  SLEQP_TR_SOLVER tr_solver = sleqp_options_get_int(solver->options,
                                                    SLEQP_OPTION_INT_TR_SOLVER);

  if(tr_solver == SLEQP_TR_SOLVER_LSQR)
  {
    SLEQP_CALL(sleqp_lsqr_solver_compute_step(solver->lsqr_solver,
                                              solver->newton_step));
  }
  else
  // compute Newton step
  {
    SLEQP_CALL(sleqp_newton_set_time_limit(solver->newton_data,
                                           sleqp_solver_remaining_time(solver)));

    SLEQP_CALL(sleqp_newton_compute_step(solver->newton_data,
                                         solver->multipliers,
                                         solver->newton_step));
  }

  SLEQP_CALL(sleqp_problem_hess_prod(problem,
                                     &one,
                                     solver->newton_step,
                                     solver->multipliers,
                                     solver->newton_hessian_step));

  {
    SLEQP_LINESEARCH lineserach = sleqp_options_get_int(solver->options,
                                                        SLEQP_OPTION_INT_LINESEARCH);

    double step_length;

    if(lineserach == SLEQP_LINESEARCH_EXACT)
    {
      SLEQP_CALL(sleqp_linesearch_trial_step_exact(solver->linesearch,
                                                   solver->cauchy_step,
                                                   solver->cauchy_hessian_step,
                                                   cauchy_merit_value,
                                                   solver->newton_step,
                                                   solver->newton_hessian_step,
                                                   solver->multipliers,
                                                   solver->trial_step,
                                                   &step_length,
                                                   trial_merit_value));
    }
    else
    {
      assert(lineserach == SLEQP_LINESEARCH_APPROX);

      SLEQP_CALL(sleqp_linesearch_trial_step(solver->linesearch,
                                             solver->cauchy_step,
                                             solver->cauchy_hessian_step,
                                             cauchy_merit_value,
                                             solver->newton_step,
                                             solver->newton_hessian_step,
                                             solver->multipliers,
                                             solver->trial_step,
                                             &step_length,
                                             trial_merit_value));
    }
  }

#if !defined(NDEBUG)

  {
    double actual_quadratic_merit_value;

    double func_dual = 1.;

    SLEQP_CALL(sleqp_merit_quadratic(solver->merit_data,
                                     iterate,
                                     &func_dual,
                                     solver->trial_step,
                                     solver->multipliers,
                                     solver->penalty_parameter,
                                     &actual_quadratic_merit_value));

    sleqp_assert_is_eq(*trial_merit_value,
                       actual_quadratic_merit_value,
                       eps);
  }

#endif

  SLEQP_CALL(sleqp_sparse_vector_add(sleqp_iterate_get_primal(iterate),
                                     solver->trial_step,
                                     zero_eps,
                                     solver->initial_trial_point));

  SLEQP_CALL(sleqp_sparse_vector_clip(solver->initial_trial_point,
                                      sleqp_problem_var_lb(problem),
                                      sleqp_problem_var_ub(problem),
                                      zero_eps,
                                      sleqp_iterate_get_primal(solver->trial_iterate)));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_solver_compute_trial_point_soc(SleqpSolver* solver)
{
  const double zero_eps = sleqp_params_get(solver->params,
                                           SLEQP_PARAM_ZERO_EPS);

  SleqpProblem* problem = solver->problem;

  SleqpIterate* iterate = solver->iterate;
  SleqpIterate* trial_iterate = solver->trial_iterate;

  SleqpSparseVec* trial_point = sleqp_iterate_get_primal(trial_iterate);

  SLEQP_CALL(sleqp_soc_compute_trial_point(solver->soc_data,
                                           solver->aug_jacobian,
                                           iterate,
                                           solver->trial_step,
                                           trial_iterate,
                                           solver->soc_trial_point,
                                           &solver->soc_step_norm));

  SLEQP_CALL(sleqp_sparse_vector_clip(solver->soc_trial_point,
                                      sleqp_problem_var_lb(problem),
                                      sleqp_problem_var_ub(problem),
                                      zero_eps,
                                      trial_point));

  return SLEQP_OKAY;
}
