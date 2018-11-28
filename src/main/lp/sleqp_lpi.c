#include "sleqp_lpi.h"

#include <assert.h>

#include "sleqp_mem.h"

struct SleqpLPi
{
  // data
  void* lp_data;

  int num_variables, num_constraints;

  // callbacks
  SLEQP_LPI_CREATE create_problem;
  SLEQP_LPI_SOLVE solve;
  SLEQP_LPI_SET_BOUNDS set_bounds;
  SLEQP_LPI_SET_COEFFICIENTS set_coefficients;
  SLEQP_LPI_SET_OBJECTIVE set_objective;
  SLEQP_LPI_GET_SOLUTION get_solution;
  SLEQP_LPI_GET_VARSTATS get_varstats;
  SLEQP_LPI_GET_CONSSTATS get_consstats;
  SLEQP_LPI_FREE free_problem;
};

SLEQP_RETCODE sleqp_lpi_create_interface(SleqpLPi** lp_star,
                                         int num_variables,
                                         int num_constraints,
                                         SleqpParams* params,
                                         SLEQP_LPI_CREATE create_problem,
                                         SLEQP_LPI_SOLVE solve,
                                         SLEQP_LPI_SET_BOUNDS set_bounds,
                                         SLEQP_LPI_SET_COEFFICIENTS set_coefficients,
                                         SLEQP_LPI_SET_OBJECTIVE set_objective,
                                         SLEQP_LPI_GET_SOLUTION get_solution,
                                         SLEQP_LPI_GET_VARSTATS get_varstats,
                                         SLEQP_LPI_GET_CONSSTATS get_consstats,
                                         SLEQP_LPI_FREE free_problem)
{
  SLEQP_CALL(sleqp_malloc(lp_star));

  SleqpLPi* lp_interface = *lp_star;

  lp_interface->lp_data = NULL;

  lp_interface->create_problem = create_problem;

  lp_interface->solve = solve;

  lp_interface->set_bounds = set_bounds;
  lp_interface->set_coefficients = set_coefficients;
  lp_interface->set_objective = set_objective;

  lp_interface->get_solution = get_solution;
  lp_interface->get_varstats = get_varstats;
  lp_interface->get_consstats = get_consstats;
  lp_interface->free_problem = free_problem;

  lp_interface->num_variables = num_variables;
  lp_interface->num_constraints = num_constraints;

  SLEQP_CALL(create_problem(&lp_interface->lp_data,
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
  return lp_interface->solve(lp_interface->lp_data,
                             lp_interface->num_variables,
                             lp_interface->num_constraints);
}

SLEQP_RETCODE sleqp_lpi_set_bounds(SleqpLPi* lp_interface,
                                   double* cons_lb,
                                   double* cons_ub,
                                   double* vars_lb,
                                   double* vars_ub)
{
  return lp_interface->set_bounds(lp_interface->lp_data,
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
  return lp_interface->set_coefficients(lp_interface->lp_data,
                                        lp_interface->num_variables,
                                        lp_interface->num_constraints,
                                        coeff_matrix);
}

SLEQP_RETCODE sleqp_lpi_set_objective(SleqpLPi* lp_interface,
                                      double* objective)
{
  return lp_interface->set_objective(lp_interface->lp_data,
                                     lp_interface->num_variables,
                                     lp_interface->num_constraints,
                                     objective);
}

SLEQP_RETCODE sleqp_lpi_get_solution(SleqpLPi* lp_interface,
                                     double* objective_value,
                                     double* solution_values)
{
  return lp_interface->get_solution(lp_interface->lp_data,
                                    lp_interface->num_variables,
                                    lp_interface->num_constraints,
                                    objective_value,
                                    solution_values);
}

SLEQP_RETCODE sleqp_lpi_get_varstats(SleqpLPi* lp_interface,
                                     SLEQP_BASESTAT* variable_stats)
{
  return lp_interface->get_varstats(lp_interface->lp_data,
                                    lp_interface->num_variables,
                                    lp_interface->num_constraints,
                                    variable_stats);
}

SLEQP_RETCODE sleqp_lpi_get_consstats(SleqpLPi* lp_interface,
                                      SLEQP_BASESTAT* constraint_stats)
{
  return lp_interface->get_consstats(lp_interface->lp_data,
                                     lp_interface->num_variables,
                                     lp_interface->num_constraints,
                                     constraint_stats);
}

SLEQP_RETCODE sleqp_lpi_free(SleqpLPi** lp_star)
{
  SleqpLPi* lp_interface = *lp_star;

  lp_interface->free_problem(&lp_interface->lp_data);

  sleqp_free(lp_star);

  return SLEQP_OKAY;
}
