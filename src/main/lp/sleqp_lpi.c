#include "sleqp_lpi.h"

#include <assert.h>

#include "sleqp_mem.h"

struct SleqpLPi
{
  // data
  void* lp_data;

  // callbacks
  SLEQP_LPI_CREATE create_problem;
  SLEQP_LPI_SOLVE solve;
  SLEQP_LPI_GET_SOLUTION get_solution;
  SLEQP_LPI_GET_VARSTATS get_varstats;
  SLEQP_LPI_FREE free_problem;
};

SLEQP_RETCODE sleqp_lpi_create_interface(SleqpLPi** lp_star,
                                         SLEQP_LPI_CREATE create_problem,
                                         SLEQP_LPI_SOLVE solve,
                                         SLEQP_LPI_GET_SOLUTION get_solution,
                                         SLEQP_LPI_GET_VARSTATS get_varstats,
                                         SLEQP_LPI_FREE free_problem)
{
  sleqp_malloc(lp_star);

  SleqpLPi* lp_interface = *lp_star;

  lp_interface->lp_data = NULL;

  lp_interface->create_problem = create_problem;
  lp_interface->solve = solve;
  lp_interface->get_solution = get_solution;
  lp_interface->get_varstats = get_varstats;
  lp_interface->free_problem = free_problem;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_lpi_create_problem(SleqpLPi* lp_interface,
                                       int num_variables,
                                       int num_constraints)
{
  assert(lp_interface);

  return lp_interface->create_problem(&lp_interface->lp_data,
                                      num_variables,
                                      num_constraints);
}

SLEQP_RETCODE sleqp_lpi_solve(SleqpLPi* lp_interface,
                              double* objective,
                              SleqpSparseMatrix* cons_matrix,
                              double* cons_lb,
                              double* cons_ub,
                              double* vars_lb,
                              double* vars_ub)
{
  assert(lp_interface);

  return lp_interface->solve(lp_interface->lp_data,
                             objective,
                             cons_matrix,
                             cons_lb,
                             cons_ub,
                             vars_lb,
                             vars_ub);
}

SLEQP_RETCODE sleqp_lpi_get_solution(SleqpLPi* lp_interface,
                                     int num_variables,
                                     double* objective_value,
                                     double* solution_values)
{
  assert(lp_interface);

  return lp_interface->get_solution(lp_interface->lp_data,
                                    num_variables,
                                    objective_value,
                                    solution_values);
}

SLEQP_RETCODE sleqp_lpi_get_varstats(SleqpLPi* lp_interface,
                                     int num_variables,
                                     SLEQP_BASESTAT* variable_stats)
{
  assert(lp_interface);

  return lp_interface->get_varstats(lp_interface->lp_data,
                                    num_variables,
                                    variable_stats);
}

SLEQP_RETCODE sleqp_lpi_free(SleqpLPi** lp_star)
{
  SleqpLPi* lp_interface = *lp_star;

  assert(lp_interface);

  lp_interface->free_problem(&lp_interface->lp_data);

  sleqp_free(lp_star);

  return SLEQP_OKAY;
}
