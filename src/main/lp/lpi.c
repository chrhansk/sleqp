#include "lpi.h"

#include <assert.h>
#include <string.h>

#include "log.h"
#include "mem.h"

struct SleqpLPi
{
  int refcount;
  // data
  void* lp_data;

  char* name;
  char* version;

  SleqpTimer* timer;

  int num_variables, num_constraints;

  double time_limit;

  // callbacks
  SleqpLPiCallbacks callbacks;
};

SLEQP_RETCODE
sleqp_lpi_create(SleqpLPi** lp_star,
                 const char* name,
                 const char* version,
                 int num_variables,
                 int num_constraints,
                 SleqpParams* params,
                 SleqpOptions* options,
                 SleqpLPiCallbacks* callbacks)
{
  SLEQP_CALL(sleqp_malloc(lp_star));

  SleqpLPi* lp_interface = *lp_star;

  *lp_interface = (SleqpLPi){0};

  lp_interface->refcount = 1;

  lp_interface->name    = strdup(name);
  lp_interface->version = strdup(version);

  SLEQP_CALL(sleqp_timer_create(&lp_interface->timer));

  lp_interface->callbacks = *callbacks;

  lp_interface->num_variables   = num_variables;
  lp_interface->num_constraints = num_constraints;

  lp_interface->time_limit = SLEQP_NONE;

  SLEQP_CALL(lp_interface->callbacks.create_problem(&lp_interface->lp_data,
                                                    num_variables,
                                                    num_constraints,
                                                    params,
                                                    options));

  return SLEQP_OKAY;
}

const char*
sleqp_lpi_name(SleqpLPi* lp_interface)
{
  return lp_interface->name;
}

const char*
sleqp_lpi_version(SleqpLPi* lp_interface)
{
  return lp_interface->version;
}

int
sleqp_lpi_num_vars(SleqpLPi* lp_interface)
{
  return lp_interface->num_variables;
}

int
sleqp_lpi_num_cons(SleqpLPi* lp_interface)
{
  return lp_interface->num_constraints;
}

SLEQP_RETCODE
sleqp_lpi_solve(SleqpLPi* lp_interface)
{
  SLEQP_CALL(sleqp_timer_start(lp_interface->timer));

  SLEQP_RETCODE retcode
    = lp_interface->callbacks.solve(lp_interface->lp_data,
                                    lp_interface->num_variables,
                                    lp_interface->num_constraints,
                                    lp_interface->time_limit);

  SLEQP_CALL(sleqp_timer_stop(lp_interface->timer));

  return retcode;
}

SLEQP_LP_STATUS
sleqp_lpi_status(SleqpLPi* lp_interface)
{
  return lp_interface->callbacks.status(lp_interface->lp_data);
}

SLEQP_RETCODE
sleqp_lpi_set_bounds(SleqpLPi* lp_interface,
                     double* cons_lb,
                     double* cons_ub,
                     double* vars_lb,
                     double* vars_ub)
{
  return lp_interface->callbacks.set_bounds(lp_interface->lp_data,
                                            lp_interface->num_variables,
                                            lp_interface->num_constraints,
                                            cons_lb,
                                            cons_ub,
                                            vars_lb,
                                            vars_ub);
}

SLEQP_RETCODE
sleqp_lpi_set_coeffs(SleqpLPi* lp_interface, SleqpSparseMatrix* coeff_matrix)
{
  return lp_interface->callbacks.set_coefficients(lp_interface->lp_data,
                                                  lp_interface->num_variables,
                                                  lp_interface->num_constraints,
                                                  coeff_matrix);
}

SLEQP_RETCODE
sleqp_lpi_set_objective(SleqpLPi* lp_interface, double* objective)
{
  return lp_interface->callbacks.set_objective(lp_interface->lp_data,
                                               lp_interface->num_variables,
                                               lp_interface->num_constraints,
                                               objective);
}

SLEQP_RETCODE
sleqp_lpi_set_time_limit(SleqpLPi* lp_interface, double time_limit)
{
  lp_interface->time_limit = time_limit;

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_lpi_save_basis(SleqpLPi* lp_interface, int index)
{
  return lp_interface->callbacks.save_basis(lp_interface->lp_data, index);
}

SLEQP_RETCODE
sleqp_lpi_restore_basis(SleqpLPi* lp_interface, int index)
{
  return lp_interface->callbacks.restore_basis(lp_interface->lp_data, index);
}

SLEQP_RETCODE
sleqp_lpi_primal_sol(SleqpLPi* lp_interface,
                     double* objective_value,
                     double* solution_values)
{
  return lp_interface->callbacks.primal_sol(lp_interface->lp_data,
                                            lp_interface->num_variables,
                                            lp_interface->num_constraints,
                                            objective_value,
                                            solution_values);
}

SLEQP_RETCODE
sleqp_lpi_dual_sol(SleqpLPi* lp_interface, double* vars_dual, double* cons_dual)
{
  return lp_interface->callbacks.dual_sol(lp_interface->lp_data,
                                          lp_interface->num_variables,
                                          lp_interface->num_constraints,
                                          vars_dual,
                                          cons_dual);
}

SLEQP_RETCODE
sleqp_lpi_vars_stats(SleqpLPi* lp_interface, SLEQP_BASESTAT* variable_stats)
{
  return lp_interface->callbacks.vars_stats(lp_interface->lp_data,
                                            lp_interface->num_variables,
                                            lp_interface->num_constraints,
                                            variable_stats);
}

SLEQP_RETCODE
sleqp_lpi_cons_stats(SleqpLPi* lp_interface, SLEQP_BASESTAT* constraint_stats)
{
  return lp_interface->callbacks.cons_stats(lp_interface->lp_data,
                                            lp_interface->num_variables,
                                            lp_interface->num_constraints,
                                            constraint_stats);
}

SleqpTimer*
sleqp_lpi_solve_timer(SleqpLPi* lp_interface)
{
  return lp_interface->timer;
}

SLEQP_RETCODE
sleqp_lpi_basis_condition_estimate(SleqpLPi* lp_interface,
                                   bool* exact,
                                   double* condition)
{
  return lp_interface->callbacks.basis_condition_estimate(lp_interface->lp_data,
                                                          exact,
                                                          condition);
}

static SLEQP_RETCODE
lpi_free(SleqpLPi** lp_star)
{
  SleqpLPi* lp_interface = *lp_star;

  lp_interface->callbacks.free_problem(&lp_interface->lp_data);

  SLEQP_CALL(sleqp_timer_free(&lp_interface->timer));

  sleqp_free(&lp_interface->version);
  sleqp_free(&lp_interface->name);

  sleqp_free(lp_star);

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_lpi_capture(SleqpLPi* lp_interface)
{
  ++lp_interface->refcount;

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_lpi_release(SleqpLPi** star)
{
  SleqpLPi* lp_interface = *star;

  if (!lp_interface)
  {
    return SLEQP_OKAY;
  }

  if (--lp_interface->refcount == 0)
  {
    SLEQP_CALL(lpi_free(star));
  }

  *star = NULL;

  return SLEQP_OKAY;
}
