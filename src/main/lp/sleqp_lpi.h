#ifndef SLEQP_LPI_H
#define SLEQP_LPI_H

/**
 * @file sleqp_lpi.h
 * @brief Definition of the LP interface.
 **/

#include "sparse/sleqp_sparse_matrix.h"
#include "sparse/sleqp_sparse_vec.h"

#include "sleqp_types.h"
#include "sleqp_lpi_types.h"

#ifdef __cplusplus
extern "C" {
#endif

  SLEQP_RETCODE sleqp_lpi_create_interface(SleqpLPi** lp_interface,
                                           int num_variables,
                                           int num_constraints,
                                           SLEQP_LPI_CREATE create_problem,
                                           SLEQP_LPI_SOLVE solve,
                                           SLEQP_LPI_SET_BOUNDS set_bounds,
                                           SLEQP_LPI_SET_COEFFICIENTS set_coefficients,
                                           SLEQP_LPI_SET_OBJECTIVE set_objective,
                                           SLEQP_LPI_GET_SOLUTION get_solution,
                                           SLEQP_LPI_GET_VARSTATS get_varstats,
                                           SLEQP_LPI_GET_CONSSTATS get_consstats,
                                           SLEQP_LPI_FREE free_problem);

  int sleqp_lpi_get_num_variables(SleqpLPi* lp_interface);

  int sleqp_lpi_get_num_constraints(SleqpLPi* lp_interface);

  SLEQP_RETCODE sleqp_lpi_solve(SleqpLPi* lp_interface);

  SLEQP_RETCODE sleqp_lpi_set_bounds(SleqpLPi* lp_interface,
                                     double* cons_lb,
                                     double* cons_ub,
                                     double* vars_lb,
                                     double* vars_ub);

  SLEQP_RETCODE sleqp_lpi_set_coefficients(SleqpLPi* lp_interface,
                                           SleqpSparseMatrix* coeff_matrix);

  SLEQP_RETCODE sleqp_lpi_set_objective(SleqpLPi* lp_interface,
                                        double* objective);

  SLEQP_RETCODE sleqp_lpi_get_solution(SleqpLPi* lp_interface,
                                       double* objective_value,
                                       double* solution_values);

  SLEQP_RETCODE sleqp_lpi_get_varstats(SleqpLPi* lp_interface,
                                       SLEQP_BASESTAT* variable_stats);

  SLEQP_RETCODE sleqp_lpi_get_consstats(SleqpLPi* lp_interface,
                                        SLEQP_BASESTAT* constraint_stats);

  SLEQP_RETCODE sleqp_lpi_free(SleqpLPi** lp_interface);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_LPI_H */
