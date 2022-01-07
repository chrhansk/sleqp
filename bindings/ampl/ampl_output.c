#include "ampl_output.h"
#include "ampl_mem.h"

#include <asl.h>

// solve_result values as suggested
// in AMPL book, section 14.2

enum AMPL_CODE
{
  AMPL_OPTIMAL    = 0,
  AMPL_INFEASIBLE = 200,
  AMPL_UNBOUNDED  = 300,
  AMPL_TIME_LIMIT = 400,
  AMPL_ITER_LIMIT = 410,
  AMPL_DEADPOINT  = 500,
  AMPL_UNKNOWN    = 501
};

static SLEQP_RETCODE
report_with_status_message(SleqpSolver* solver,
                           ASL* asl,
                           Option_Info* option_info,
                           double* primal,
                           double* cons_dual)
{
  const char* message = NULL;

  SLEQP_STATUS status = sleqp_solver_status(solver);

  switch (status)
  {
  case SLEQP_STATUS_OPTIMAL:
    message          = "Optimal solution found";
    solve_result_num = AMPL_OPTIMAL;
    break;
  case SLEQP_STATUS_INFEASIBLE:
    message          = "Failed to find feasible solution";
    solve_result_num = AMPL_INFEASIBLE;
    break;
  case SLEQP_STATUS_UNBOUNDED:
    message          = "Problem appears unbounded";
    solve_result_num = AMPL_UNBOUNDED;
    break;
  case SLEQP_STATUS_ABORT_DEADPOINT:
    message          = "Reached dead point";
    solve_result_num = AMPL_DEADPOINT;
    break;
  case SLEQP_STATUS_ABORT_ITER:
    message          = "Reached iteration limit";
    solve_result_num = AMPL_ITER_LIMIT;
    break;
  case SLEQP_STATUS_ABORT_TIME:
    message          = "Reached time limit";
    solve_result_num = AMPL_TIME_LIMIT;
    break;
  default:
    message          = "Unknown status";
    solve_result_num = AMPL_UNKNOWN;
    break;
  }

  write_sol(message, primal, cons_dual, option_info);

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_ampl_report(SleqpProblem* problem,
                  SleqpSolver* solver,
                  ASL* asl,
                  Option_Info* option_info)
{

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

  SLEQP_CALL(
    report_with_status_message(solver, asl, option_info, primal, cons_dual));

  sleqp_ampl_free(&cons_dual);
  sleqp_ampl_free(&primal);

  return SLEQP_OKAY;
}
