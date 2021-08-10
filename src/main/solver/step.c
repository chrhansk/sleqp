#include "solver.h"

SLEQP_RETCODE sleqp_solver_set_func_value(SleqpSolver* solver,
                                          SleqpIterate* iterate,
                                          SLEQP_VALUE_REASON reason)
{
  SleqpProblem* problem = solver->problem;

  int func_grad_nnz = 0;
  int cons_val_nnz = 0;
  int cons_jac_nnz = 0;

  SLEQP_CALL(sleqp_problem_set_value(problem,
                                     sleqp_iterate_get_primal(iterate),
                                     reason,
                                     &func_grad_nnz,
                                     &cons_val_nnz,
                                     &cons_jac_nnz));

  SLEQP_CALL(sleqp_sparse_vector_reserve(sleqp_iterate_get_func_grad(iterate),
                                         func_grad_nnz));

  SLEQP_CALL(sleqp_sparse_vector_reserve(sleqp_iterate_get_cons_val(iterate),
                                         cons_val_nnz));

  SLEQP_CALL(sleqp_sparse_matrix_reserve(sleqp_iterate_get_cons_jac(iterate),
                                         cons_jac_nnz));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_solver_reject_step(SleqpSolver* solver)
{
  SleqpIterate* iterate = solver->iterate;

  SLEQP_CALL(sleqp_solver_set_func_value(solver,
                                         iterate,
                                         SLEQP_VALUE_REASON_REJECTED_ITERATE));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_solver_accept_step(SleqpSolver* solver)
{
  SleqpProblem* problem = solver->problem;
  SleqpIterate* iterate = solver->iterate;
  SleqpIterate* trial_iterate = solver->trial_iterate;

  SLEQP_CALL(sleqp_solver_set_func_value(solver,
                                         trial_iterate,
                                         SLEQP_VALUE_REASON_ACCEPTED_ITERATE));

  // get the remaining data to fill the iterate

  SLEQP_CALL(sleqp_problem_eval(problem,
                                NULL,
                                NULL,
                                sleqp_iterate_get_func_grad(trial_iterate),
                                sleqp_iterate_get_cons_val(trial_iterate),
                                sleqp_iterate_get_cons_jac(trial_iterate)));

  if(solver->bfgs_data)
  {
    SLEQP_CALL(sleqp_bfgs_push(solver->bfgs_data,
                               solver->iterate,
                               solver->trial_iterate,
                               solver->multipliers));
  }

  if(solver->sr1_data)
  {
    SLEQP_CALL(sleqp_sr1_push(solver->sr1_data,
                              solver->iterate,
                              solver->trial_iterate,
                              solver->multipliers));
  }

  // perform simple swaps
  solver->trial_iterate = iterate;
  solver->iterate = trial_iterate;

  // ensure that the unscaled iterate is kept up to date
  if(solver->scaling_data || solver->preprocessor)
  {
    solver->restore_original_iterate = true;
  }
  else
  {
    solver->original_iterate = solver->iterate;
  }

  SleqpCallbackHandler* handler = solver->callback_handlers[SLEQP_SOLVER_EVENT_ACCEPTED_ITERATE];

  if (sleqp_callback_handler_size(handler) != 0) {
    SLEQP_CALL(sleqp_solver_restore_original_iterate(solver));

    SLEQP_CALLBACK_HANDLER_EXECUTE(handler, SLEQP_ACCEPTED_ITERATE, solver,
                                   solver->original_iterate);
  }

  return SLEQP_OKAY;
}
