#include "sleqp_lpi.h"

#include <assert.h>

#include "sleqp_mem.h"

struct SleqpLPi
{
  // data
  void* lp_data;

  SleqpTimer* timer;

  int num_variables, num_constraints;

  double time_limit;

  // callbacks
  SleqpLPiCallbacks callbacks;
};

SLEQP_RETCODE sleqp_lpi_create_interface(SleqpLPi** lp_star,
                                         int num_variables,
                                         int num_constraints,
                                         SleqpParams* params,
                                         SleqpLPiCallbacks* callbacks)
{
  SLEQP_CALL(sleqp_malloc(lp_star));

  SleqpLPi* lp_interface = *lp_star;

  *lp_interface = (SleqpLPi) {0};

  SLEQP_CALL(sleqp_timer_create(&lp_interface->timer));

  lp_interface->callbacks = *callbacks;

  lp_interface->num_variables = num_variables;
  lp_interface->num_constraints = num_constraints;

  lp_interface->time_limit = -1;

  SLEQP_CALL(lp_interface->callbacks.create_problem(&lp_interface->lp_data,
                                                    num_variables,
                                                    num_constraints,
                                                    params));

  return SLEQP_OKAY;
}

int sleqp_lpi_get_num_variables(SleqpLPi* lp_interface)
{
  return lp_interface->num_variables;
}

int sleqp_lpi_get_num_constraints(SleqpLPi* lp_interface)
{
  return lp_interface->num_constraints;
}

SLEQP_RETCODE sleqp_lpi_solve(SleqpLPi* lp_interface)
{
  SLEQP_CALL(sleqp_timer_start(lp_interface->timer));

  SLEQP_CALL(lp_interface->callbacks.solve(lp_interface->lp_data,
                                           lp_interface->num_variables,
                                           lp_interface->num_constraints,
                                           lp_interface->time_limit));

  SLEQP_CALL(sleqp_timer_stop(lp_interface->timer));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_lpi_set_bounds(SleqpLPi* lp_interface,
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

SLEQP_RETCODE sleqp_lpi_set_coefficients(SleqpLPi* lp_interface,
                                         SleqpSparseMatrix* coeff_matrix)
{
  return lp_interface->callbacks.set_coefficients(lp_interface->lp_data,
                                                  lp_interface->num_variables,
                                                  lp_interface->num_constraints,
                                                  coeff_matrix);
}

SLEQP_RETCODE sleqp_lpi_set_objective(SleqpLPi* lp_interface,
                                      double* objective)
{
  return lp_interface->callbacks.set_objective(lp_interface->lp_data,
                                               lp_interface->num_variables,
                                               lp_interface->num_constraints,
                                               objective);
}

SLEQP_RETCODE sleqp_lpi_set_time_limit(SleqpLPi* lp_interface,
                                       double time_limit)
{
  lp_interface->time_limit = time_limit;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_lpi_get_solution(SleqpLPi* lp_interface,
                                     double* objective_value,
                                     double* solution_values)
{
  return lp_interface->callbacks.get_solution(lp_interface->lp_data,
                                              lp_interface->num_variables,
                                              lp_interface->num_constraints,
                                              objective_value,
                                              solution_values);
}

SLEQP_RETCODE sleqp_lpi_get_varstats(SleqpLPi* lp_interface,
                                     SLEQP_BASESTAT* variable_stats)
{
  return lp_interface->callbacks.get_varstats(lp_interface->lp_data,
                                              lp_interface->num_variables,
                                              lp_interface->num_constraints,
                                              variable_stats);
}

SLEQP_RETCODE sleqp_lpi_get_consstats(SleqpLPi* lp_interface,
                                      SLEQP_BASESTAT* constraint_stats)
{
  return lp_interface->callbacks.get_consstats(lp_interface->lp_data,
                                               lp_interface->num_variables,
                                               lp_interface->num_constraints,
                                               constraint_stats);
}

SleqpTimer* sleqp_lpi_get_solve_timer(SleqpLPi* lp_interface)
{
  return lp_interface->timer;
}

SLEQP_RETCODE sleqp_lpi_get_basis_condition(SleqpLPi* lp_interface,
                                            bool* exact,
                                            double* condition)
{
  return lp_interface->callbacks.get_basis_condition(lp_interface->lp_data,
                                                     exact,
                                                     condition);
}

SLEQP_RETCODE sleqp_lpi_free(SleqpLPi** lp_star)
{
  SleqpLPi* lp_interface = *lp_star;

  if(!lp_interface)
  {
    return SLEQP_OKAY;
  }

  lp_interface->callbacks.free_problem(&lp_interface->lp_data);

  SLEQP_CALL(sleqp_timer_free(&lp_interface->timer));

  sleqp_free(lp_star);

  return SLEQP_OKAY;
}
