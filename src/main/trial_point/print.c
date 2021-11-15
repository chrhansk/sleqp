#include "trial_point.h"

SLEQP_RETCODE
sleqp_trial_point_solver_print_stats(SleqpTrialPointSolver* solver,
                                     double elapsed_seconds)
{
  SLEQP_CALL(sleqp_timer_display(sleqp_aug_jac_creation_timer(solver->aug_jac),
                                 "Factorizations",
                                 elapsed_seconds));

  SLEQP_CALL(sleqp_timer_display(sleqp_aug_jac_solution_timer(solver->aug_jac),
                                 "Substitutions",
                                 elapsed_seconds));

  if (solver->lp_interface)
  {
    SLEQP_CALL(
      sleqp_timer_display(sleqp_lpi_get_solve_timer(solver->lp_interface),
                          "Solved LPs",
                          elapsed_seconds));
  }

  SLEQP_CALL(sleqp_timer_display(sleqp_eqp_solver_get_timer(solver->eqp_solver),
                                 "Solved EQPs",
                                 elapsed_seconds));

  SLEQP_CALL(sleqp_timer_display(sleqp_linesearch_get_timer(solver->linesearch),
                                 "Line searches",
                                 elapsed_seconds));

  return SLEQP_OKAY;
}
