#ifndef SLEQP_LPI_H
#define SLEQP_LPI_H

#include "sparse/sleqp_sparse_matrix.h"
#include "sparse/sleqp_sparse_vec.h"

#include "sleqp_types.h"
#include "sleqp_lpi_types.h"

#ifdef __cplusplus
extern "C" {
#endif

  SLEQP_RETCODE sleqp_lpi_create_interface(SleqpLPi** lp_interface,
                                           SLEQP_LPI_CREATE create_problem,
                                           SLEQP_LPI_SOLVE solve,
                                           SLEQP_LPI_GET_SOLUTION get_solution,
                                           SLEQP_LPI_GET_VARSTATS get_varstats,
                                           SLEQP_LPI_GET_CONSSTATS get_consstats,
                                           SLEQP_LPI_FREE free_problem);

  SLEQP_RETCODE sleqp_lpi_create_problem(SleqpLPi* lp_interface,
                                         int num_variables,
                                         int num_constraints);

  SLEQP_RETCODE sleqp_lpi_solve(SleqpLPi* lp_interface,
                                double* objective,
                                SleqpSparseMatrix* cons_matrix,
                                double* cons_lb,
                                double* cons_ub,
                                double* vars_lb,
                                double* vars_ub);

  SLEQP_RETCODE sleqp_lpi_get_solution(SleqpLPi* lp_interface,
                                       int num_variables,
                                       double* objective_value,
                                       double* solution_values);

  SLEQP_RETCODE sleqp_lpi_get_varstats(SleqpLPi* lp_interface,
                                       int num_variables,
                                       SLEQP_BASESTAT* variable_stats);

  SLEQP_RETCODE sleqp_lpi_get_consstats(SleqpLPi* lp_interface,
                                        int num_constraints,
                                        SLEQP_BASESTAT* constraint_stats);

  SLEQP_RETCODE sleqp_lpi_free(SleqpLPi** lp_interface);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_LPI_H */
