#include "ampl_output.h"
#include "ampl_suffix.h"
#include "ampl_util.h"

#include <assert.h>

#define BUF_SIZE 512

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

// AMPL variable / constraint states,
// as defined in the sstatus_table
typedef enum
{
  AMPL_STATE_NONE       = 0,
  AMPL_STATE_BASIC      = 1,
  AMPL_STATE_SUPERBASIC = 2,
  AMPL_STATE_LOWER      = 3,
  AMPL_STATE_UPPER      = 4,
  AMPL_STATE_EQUAL      = 5,
  AMPL_STATE_BETWEEN    = 6,
} AMPL_STATE;

static AMPL_STATE
map_state(SLEQP_ACTIVE_STATE state)
{
  switch (state)
  {
  case SLEQP_INACTIVE:
    return AMPL_STATE_BETWEEN;
  case SLEQP_ACTIVE_LOWER:
    return AMPL_STATE_LOWER;
  case SLEQP_ACTIVE_UPPER:
    return AMPL_STATE_UPPER;
  case SLEQP_ACTIVE_BOTH:
    return AMPL_STATE_EQUAL;
  }

  assert(0);
  return AMPL_STATE_NONE;
}

static SLEQP_RETCODE
set_suffixes(SleqpProblem* problem, SleqpSolver* solver, ASL* asl)
{
  SufDecl* ampl_suffixes = sleqp_ampl_suffixes();

  int* var_stats;
  int* cons_stats;

  const int num_vars = sleqp_problem_num_vars(problem);
  const int num_cons = sleqp_problem_num_cons(problem);

  SLEQP_CALL(sleqp_alloc_array(&var_stats, num_vars));
  SLEQP_CALL(sleqp_alloc_array(&cons_stats, num_cons));

  SleqpIterate* iterate;

  SLEQP_CALL(sleqp_solver_solution(solver, &iterate));

  SleqpWorkingSet* working_set = sleqp_iterate_working_set(iterate);

  for (int j = 0; j < num_vars; ++j)
  {
    var_stats[j] = map_state(sleqp_working_set_var_state(working_set, j));
  }

  suf_iput(ampl_suffixes[AMPL_SUFFIX_VARSTAT].name,
           ampl_suffixes[AMPL_SUFFIX_VARSTAT].kind,
           var_stats);

  for (int i = 0; i < num_cons; ++i)
  {
    cons_stats[i] = map_state(sleqp_working_set_cons_state(working_set, i));
  }

  suf_iput(ampl_suffixes[AMPL_SUFFIX_CONSSTAT].name,
           ampl_suffixes[AMPL_SUFFIX_CONSSTAT].kind,
           cons_stats);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
report_with_status_message(SleqpSolver* solver,
                           ASL* asl,
                           Option_Info* option_info,
                           double* primal,
                           double* cons_dual)
{
  char message[BUF_SIZE];
  char sol_val[BUF_SIZE];
  const char* status_message = NULL;

  SLEQP_STATUS status = sleqp_solver_status(solver);

  switch (status)
  {
  case SLEQP_STATUS_OPTIMAL:
    status_message   = "Optimal solution found";
    solve_result_num = AMPL_OPTIMAL;
    break;
  case SLEQP_STATUS_INFEASIBLE:
    status_message   = "Failed to find feasible solution";
    solve_result_num = AMPL_INFEASIBLE;
    break;
  case SLEQP_STATUS_UNBOUNDED:
    status_message   = "Problem appears unbounded";
    solve_result_num = AMPL_UNBOUNDED;
    break;
  case SLEQP_STATUS_ABORT_DEADPOINT:
    status_message   = "Reached dead point";
    solve_result_num = AMPL_DEADPOINT;
    break;
  case SLEQP_STATUS_ABORT_ITER:
    status_message   = "Reached iteration limit";
    solve_result_num = AMPL_ITER_LIMIT;
    break;
  case SLEQP_STATUS_ABORT_TIME:
    status_message   = "Reached time limit";
    solve_result_num = AMPL_TIME_LIMIT;
    break;
  default:
    status_message   = "Unknown status";
    solve_result_num = AMPL_UNKNOWN;
    break;
  }

  if (solve_result_num == AMPL_UNKNOWN)
  {
    snprintf(sol_val, BUF_SIZE, "");
  }
  else
  {
    SleqpIterate* iterate;

    SLEQP_CALL(sleqp_solver_solution(solver, &iterate));

    const int iterations = sleqp_solver_iterations(solver);

    snprintf(sol_val,
             BUF_SIZE,
             ", objective: %.*g, iterations: %d",
             obj_prec(),
             sleqp_iterate_obj_val(iterate),
             iterations);
  }

  snprintf(message,
           BUF_SIZE,
           "%s: %s%s",
           "SLEQP " SLEQP_LONG_VERSION,
           status_message,
           sol_val);

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

  SLEQP_CALL(sleqp_alloc_array(&primal, num_vars));

  SLEQP_CALL(sleqp_sparse_vector_to_raw(sleqp_iterate_primal(iterate), primal));

  SLEQP_CALL(sleqp_alloc_array(&cons_dual, num_cons));

  SLEQP_CALL(
    sleqp_sparse_vector_to_raw(sleqp_iterate_cons_dual(iterate), cons_dual));

  if (!sleqp_ampl_max_problem(asl))
  {
    for (int i = 0; i < num_cons; ++i)
    {
      cons_dual[i] *= -1.;
    }
  }

  SLEQP_CALL(set_suffixes(problem, solver, asl));

  SLEQP_CALL(
    report_with_status_message(solver, asl, option_info, primal, cons_dual));

  sleqp_free(&cons_dual);
  sleqp_free(&primal);

  return SLEQP_OKAY;
}
