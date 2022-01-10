#include <getstub.h>

#include "sleqp.h"

#include "ampl_data.h"
#include "ampl_keywords.h"
#include "ampl_output.h"
#include "ampl_problem.h"
#include "ampl_suffix.h"

int
main(int argc, char* argv[])
{
  SleqpOptions* options;
  SleqpParams* params;

  SLEQP_CALL(sleqp_options_create(&options));
  SLEQP_CALL(sleqp_params_create(&params));

  ASL* asl;
  // allocate ASL
  asl = ASL_alloc(ASL_read_pfgh);

  // set ampl solver options: require primal starting, point, objective,
  // derivatives
  want_xpi0   = 1;
  obj_no      = 0;
  want_derivs = 1;

  SleqpAmplKeywords* ampl_keywords;

  SLEQP_CALL(sleqp_ampl_keywords_create(&ampl_keywords, options, params));

  keyword* keywords;
  int num_keywords;

  SLEQP_CALL(sleqp_ampl_keywords_get(ampl_keywords, &keywords, &num_keywords));

  Option_Info Oinfo = {"sleqp",
                       "SLEQP",
                       "sleqp_options",
                       keywords,
                       num_keywords,
                       1,
                       "SLEQP " SLEQP_LONG_VERSION,
                       NULL,
                       NULL,
                       NULL,
                       keywords,
                       num_keywords};

  // get stub, options
  char* stub = getstops(argv, &Oinfo);

  if (stub == NULL)
  {
    sleqp_log_error("Failed to open nl stub.");
    return -1;
  }
  // get problem dimensions from stub
  FILE* nl = jac0dim(stub, strlen(stub));
  if (nl == NULL)
  {
    sleqp_log_error("Failed to read nl stub.");
    return EXIT_FAILURE;
  }

  sleqp_log_info(
    "Problem has %d variables, %d constraints, %d of which are linear",
    n_var,
    n_con,
    n_con - nlc);

  suf_declare(sleqp_ampl_suffixes(), AMPL_NUM_SUFFIX);

  // allocate problem data
  SleqpAmplData* data;
  SLEQP_CALL(sleqp_ampl_data_create(&data, asl));

  const double zero_eps = sleqp_params_value(params, SLEQP_PARAM_ZERO_EPS);

  SleqpProblem* problem;
  SleqpSolver* solver;

  SLEQP_CALL(sleqp_ampl_problem_create(&problem, data, nl, params));

  SleqpSparseVec* x;
  SLEQP_CALL(sleqp_sparse_vector_create(&x, n_var, 0));
  SLEQP_CALL(sleqp_sparse_vector_from_raw(x, data->x, n_var, zero_eps));

  SLEQP_CALL(sleqp_solver_create(&solver, problem, params, options, x, NULL));

  const int iter_limit    = sleqp_ampl_keywords_iter_limit(ampl_keywords);
  const double time_limit = sleqp_ampl_keywords_time_limit(ampl_keywords);

  bool success = true;

  SLEQP_RETCODE retcode = sleqp_solver_solve(solver, iter_limit, time_limit);

  success = (retcode == SLEQP_OKAY);

  SLEQP_CALL(sleqp_ampl_report(problem, solver, asl, &Oinfo));

  // free data
  SLEQP_CALL(sleqp_solver_release(&solver));
  SLEQP_CALL(sleqp_problem_release(&problem));
  SLEQP_CALL(sleqp_sparse_vector_free(&x));
  ASL_free(&asl);
  SLEQP_CALL(sleqp_ampl_data_free(&data));

  SLEQP_CALL(sleqp_ampl_keywords_free(&ampl_keywords));

  SLEQP_CALL(sleqp_params_release(&params));
  SLEQP_CALL(sleqp_options_release(&options));

  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
