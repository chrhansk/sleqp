#include "sleqp.h"

#include "sleqp_cutest_defs.h"

#include "sleqp_cutest_constrained.h"
#include "sleqp_cutest_unconstrained.h"

int main(int argc, char *argv[])
{
  //char *fname = "OUTSDIF.d"; /* CUTEst data file */

  char *fname = PROBLEM_OUTSDIF_D; /* CUTEst data file */

  integer funit = 42;        /* FORTRAN unit number for OUTSDIF.d */
  integer iout = 6;          /* FORTRAN unit number for error output */
  integer ierr;              /* Exit flag from OPEN and CLOSE */
  integer io_buffer = 11;    /* FORTRAN unit internal input/output */
  integer cutest_status;            /* Exit flag from CUTEst tools */


  integer CUTEst_nvar;        /* number of variables */
  integer CUTEst_ncons;       /* number of constraints */

  bool CUTest_constrained = false;

  ierr = 0;
  FORTRAN_open(&funit, fname, &ierr);

  if(ierr)
  {
    sleqp_log_error("Failed to open %s, aborting.", fname);
    return 1;
  }

  CUTEST_cdimen( &cutest_status, &funit, &CUTEst_nvar, &CUTEst_ncons);

  sleqp_log_info("Problem has %d variables, %d constraints",
                 CUTEst_nvar,
                 CUTEst_ncons);

  if(CUTEst_ncons)
  {
    CUTest_constrained = true;
  }

  double* x_dense, *var_ub_dense, *var_lb_dense;
  double* cons_ub_dense = NULL, *cons_lb_dense  = NULL;

  double *v = NULL;

  logical *equatn = NULL, *linear = NULL;

  // we don't care about the ordering (or indeed the variable /
  // constraint types) at this point
  integer e_order = 0, l_order = 0, v_order = 0;

  SLEQP_CALL(sleqp_calloc(&x_dense, CUTEst_nvar));
  SLEQP_CALL(sleqp_calloc(&var_lb_dense, CUTEst_nvar));
  SLEQP_CALL(sleqp_calloc(&var_ub_dense, CUTEst_nvar));

  if(CUTest_constrained)
  {
    SLEQP_CALL(sleqp_calloc(&cons_lb_dense, CUTEst_ncons + 1));
    SLEQP_CALL(sleqp_calloc(&cons_ub_dense, CUTEst_ncons + 1));

    SLEQP_CALL(sleqp_calloc(&equatn, CUTEst_ncons + 1));
    SLEQP_CALL(sleqp_calloc(&linear, CUTEst_ncons + 1));

    SLEQP_CALL(sleqp_calloc(&v, CUTEst_ncons + CUTEst_nvar + 1));

    CUTEST_csetup(&cutest_status, &funit, &iout, &io_buffer,
                  &CUTEst_nvar, &CUTEst_ncons,
                  x_dense,
                  var_lb_dense,
                  var_ub_dense,
                  v,
                  cons_lb_dense,
                  cons_ub_dense,
                  equatn,
                  linear,
                  &e_order, &l_order, &v_order);
  }
  else
  {
    CUTEST_usetup(&cutest_status, &funit, &iout, &io_buffer,
                  &CUTEst_nvar, x_dense, var_lb_dense, var_ub_dense);
  }

  {
    const double inf = sleqp_infinity();

    for(int i = 0; i < CUTEst_nvar; i++)
    {
      if(var_lb_dense[i] == -CUTE_INF)
      {
        var_lb_dense[i] = -inf;
      }
      if(var_ub_dense[i] ==  CUTE_INF)
      {
        var_ub_dense[i] = inf;
      }
    }

    for (int i = 0; i < CUTEst_ncons; i++)
    {
      if(cons_lb_dense[i] == -CUTE_INF)
      {
        cons_lb_dense[i] = -inf;
      }
      if(cons_ub_dense[i] == CUTE_INF)
      {
        cons_ub_dense[i] = inf;
      }
    }
  }

  SleqpParams* params;

  SLEQP_CALL(sleqp_params_create(&params));

  const double eps = sleqp_params_get_eps(params);

  SleqpSparseVec* var_lb;
  SleqpSparseVec* var_ub;
  SleqpSparseVec* x;

  SleqpSparseVec* cons_lb;
  SleqpSparseVec* cons_ub;
  SleqpFunc* func;

  SleqpProblem* problem;
  SleqpSolver* solver;

  SLEQP_CALL(sleqp_sparse_vector_create(&var_lb, CUTEst_nvar, 0));
  SLEQP_CALL(sleqp_sparse_vector_create(&var_ub, CUTEst_nvar, 0));
  SLEQP_CALL(sleqp_sparse_vector_create(&x, CUTEst_nvar, 0));

  SLEQP_CALL(sleqp_sparse_vector_create(&cons_lb, CUTEst_ncons, 0));
  SLEQP_CALL(sleqp_sparse_vector_create(&cons_ub, CUTEst_ncons, 0));

  SLEQP_CALL(sleqp_sparse_vector_from_raw(var_lb, var_lb_dense, CUTEst_nvar, eps));
  SLEQP_CALL(sleqp_sparse_vector_from_raw(var_ub, var_ub_dense, CUTEst_nvar, eps));

  SLEQP_CALL(sleqp_sparse_vector_from_raw(x, x_dense, CUTEst_nvar, eps));

  SLEQP_CALL(sleqp_sparse_vector_from_raw(cons_lb, cons_lb_dense, CUTEst_ncons, eps));
  SLEQP_CALL(sleqp_sparse_vector_from_raw(cons_ub, cons_ub_dense, CUTEst_ncons, eps));

  if(CUTest_constrained)
  {
    SLEQP_CALL(sleqp_cutest_cons_func_create(&func, CUTEst_nvar, CUTEst_ncons, eps));
  }
  else
  {
    SLEQP_CALL(sleqp_cutest_uncons_func_create(&func, CUTEst_nvar, eps));
  }

  SLEQP_CALL(sleqp_problem_create(&problem,
                                  func,
                                  params,
                                  var_lb,
                                  var_ub,
                                  cons_lb,
                                  cons_ub));

  SLEQP_CALL(sleqp_solver_create(&solver,
                                 problem,
                                 params,
                                 x));

  /**/

  int max_num_iterations = 1000;

  bool success = true;

  SLEQP_RETCODE status = sleqp_solver_solve(solver, max_num_iterations);

  if(status != SLEQP_OKAY)
  {
    success = false;
  }

  if(success && (sleqp_solver_get_status(solver) != SLEQP_OPTIMAL))
  {
    success = false;
  }

  /**/

  SLEQP_CALL(sleqp_solver_free(&solver));

  SLEQP_CALL(sleqp_problem_free(&problem));

  if(CUTest_constrained)
  {
    SLEQP_CALL(sleqp_cutest_cons_func_free(&func));
  }
  else
  {
    SLEQP_CALL(sleqp_cutest_uncons_func_free(&func));
  }

  SLEQP_CALL(sleqp_sparse_vector_free(&cons_ub));
  SLEQP_CALL(sleqp_sparse_vector_free(&cons_lb));

  SLEQP_CALL(sleqp_sparse_vector_free(&x));
  SLEQP_CALL(sleqp_sparse_vector_free(&var_ub));
  SLEQP_CALL(sleqp_sparse_vector_free(&var_lb));

  SLEQP_CALL(sleqp_params_free(&params));

  sleqp_free(&v);
  sleqp_free(&linear);
  sleqp_free(&equatn);

  sleqp_free(&cons_ub_dense);
  sleqp_free(&cons_lb_dense);

  sleqp_free(&cons_ub_dense);
  sleqp_free(&cons_lb_dense);

  sleqp_free(&var_ub_dense);
  sleqp_free(&var_lb_dense);
  sleqp_free(&x_dense);


  FORTRAN_close(&funit, &ierr);

  if(ierr)
  {
    sleqp_log_error("Error closing %s on unit %d.\n",
                    fname,
                    funit);
  }

  return !(success);
}
