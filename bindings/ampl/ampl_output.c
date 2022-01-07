#include "ampl_output.h"
#include "ampl_mem.h"

SLEQP_RETCODE
sleqp_ampl_report(SleqpProblem* problem,
                  SleqpSolver* solver,
                  ASL* asl,
                  Option_Info* option_info)
{
  const char* message = "";

  const int num_vars = sleqp_problem_num_vars(problem);
  const int num_cons = sleqp_problem_num_cons(problem);

  SleqpIterate* iterate;

  SLEQP_CALL(sleqp_solver_solution(solver, &iterate));

  double* primal;
  double* cons_dual;

  SLEQP_CALL(sleqp_ampl_alloc_array(&primal, num_vars));

  SLEQP_CALL(sleqp_sparse_vector_to_raw(sleqp_iterate_primal(iterate), primal));

  SLEQP_CALL(sleqp_ampl_alloc_array(&cons_dual, num_cons));

  SLEQP_CALL(
    sleqp_sparse_vector_to_raw(sleqp_iterate_cons_dual(iterate), cons_dual));

  write_sol(message, primal, cons_dual, option_info);

  sleqp_ampl_free(&cons_dual);
  sleqp_ampl_free(&primal);

  return SLEQP_OKAY;
}
