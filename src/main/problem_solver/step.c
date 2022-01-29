#include "problem_solver.h"

#include "error.h"

SLEQP_RETCODE
sleqp_problem_solver_set_func_value(SleqpProblemSolver* solver,
                                    SleqpIterate* iterate,
                                    SLEQP_VALUE_REASON reason,
                                    bool* reject)
{
  SleqpProblem* problem = solver->problem;

  int obj_grad_nnz = 0;
  int cons_val_nnz = 0;
  int cons_jac_nnz = 0;

  bool manual_reject = false;

  SLEQP_CALL(sleqp_problem_set_value(problem,
                                     sleqp_iterate_primal(iterate),
                                     reason,
                                     &manual_reject,
                                     &obj_grad_nnz,
                                     &cons_val_nnz,
                                     &cons_jac_nnz));

  if (reject)
  {
    *reject = manual_reject;
  }

  if (manual_reject)
  {
    if (reason != SLEQP_VALUE_REASON_TRYING_ITERATE
        && reason != SLEQP_VALUE_REASON_TRYING_SOC_ITERATE)
    {
      sleqp_raise(SLEQP_FUNC_EVAL_ERROR,
                  "Function can only reject trial steps");
    }

    return SLEQP_OKAY;
  }

  SLEQP_CALL(
    sleqp_sparse_vector_reserve(sleqp_iterate_obj_grad(iterate), obj_grad_nnz));

  SLEQP_CALL(
    sleqp_sparse_vector_reserve(sleqp_iterate_cons_val(iterate), cons_val_nnz));

  SLEQP_CALL(
    sleqp_sparse_matrix_reserve(sleqp_iterate_cons_jac(iterate), cons_jac_nnz));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_problem_solver_reject_step(SleqpProblemSolver* solver)
{
  SleqpIterate* iterate = solver->iterate;

  SLEQP_CALL(
    sleqp_problem_solver_set_func_value(solver,
                                        iterate,
                                        SLEQP_VALUE_REASON_REJECTED_ITERATE,
                                        NULL));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_problem_solver_accept_step(SleqpProblemSolver* solver)
{
  SleqpProblem* problem       = solver->problem;
  SleqpIterate* iterate       = solver->iterate;
  SleqpIterate* trial_iterate = solver->trial_iterate;

  // SleqpTrialPointSolver* trial_point_solver = solver->trial_point_solver;

  // SleqpSparseVec* multipliers =
  // sleqp_trial_point_solver_get_multipliers(trial_point_solver);

  SLEQP_CALL(
    sleqp_problem_solver_set_func_value(solver,
                                        trial_iterate,
                                        SLEQP_VALUE_REASON_ACCEPTED_ITERATE,
                                        NULL));

  // get the remaining data to fill the iterate

  SLEQP_CALL(sleqp_problem_eval(problem,
                                NULL,
                                sleqp_iterate_obj_grad(trial_iterate),
                                sleqp_iterate_cons_val(trial_iterate),
                                sleqp_iterate_cons_jac(trial_iterate)));

  SleqpCallbackHandler* handler
    = solver->callback_handlers[SLEQP_PROBLEM_SOLVER_EVENT_ACCEPTED_ITERATE];

  SLEQP_CALLBACK_HANDLER_EXECUTE(handler,
                                 SLEQP_PROBLEM_SOLVER_ACCEPTED_ITERATE,
                                 solver,
                                 iterate,
                                 trial_iterate);

  /*
  if(solver->quasi_newton)
  {
    SLEQP_CALL(sleqp_quasi_newton_push(solver->quasi_newton,
                                       solver->iterate,
                                       solver->trial_iterate,
                                       multipliers));
  }
  */

  // perform simple swaps
  solver->trial_iterate = iterate;
  solver->iterate       = trial_iterate;

  // ensure that the unscaled iterate is kept up to date

  /*
  if(solver->scaling_data || solver->preprocessor)
  {
    solver->restore_original_iterate = true;
  }
  else
  {
    solver->original_iterate = solver->iterate;
  }
  */

  /*
  if (sleqp_callback_handler_size(handler) != 0)
  {
    SLEQP_CALL(sleqp_solver_restore_original_iterate(solver));

    SLEQP_CALLBACK_HANDLER_EXECUTE(handler,
                                   SLEQP_ACCEPTED_ITERATE,
                                   solver,
                                   solver->original_iterate);
  }
  */

  return SLEQP_OKAY;
}
