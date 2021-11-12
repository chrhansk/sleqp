#ifndef SLEQP_LPI_H
#define SLEQP_LPI_H

/**
 * @file lpi.h
 * @brief Definition of the LP interface.
 *
 * The LP interface solves problem of the form
 * \f$ \min \langle c, x \rangle, \textrm{ s.t. } l_x \leq x \leq u_x, l \leq A x \leq u \f$
 *
 **/

#include "timer.h"
#include "types.h"
#include "lpi_types.h"

typedef struct SleqpLPi SleqpLPi;

SLEQP_NODISCARD
SLEQP_RETCODE sleqp_lpi_create_interface(SleqpLPi** lp_interface,
                                         const char* name,
                                         const char* version,
                                         int num_variables,
                                         int num_constraints,
                                         SleqpParams* params,
                                         SleqpOptions* options,
                                         SleqpLPiCallbacks* callbacks);

const char* sleqp_lpi_get_name(SleqpLPi* lp_interface);

const char* sleqp_lpi_get_version(SleqpLPi* lp_interface);

int sleqp_lpi_get_num_variables(SleqpLPi* lp_interface);

int sleqp_lpi_get_num_constraints(SleqpLPi* lp_interface);

SLEQP_NODISCARD
SLEQP_RETCODE sleqp_lpi_solve(SleqpLPi* lp_interface);

SLEQP_NODISCARD
SLEQP_LPI_STATUS sleqp_get_status(SleqpLPi* lp_interface);

SLEQP_NODISCARD
SLEQP_RETCODE sleqp_lpi_set_bounds(SleqpLPi* lp_interface,
                                   double* cons_lb,
                                   double* cons_ub,
                                   double* vars_lb,
                                   double* vars_ub);

SLEQP_NODISCARD
SLEQP_RETCODE sleqp_lpi_set_coefficients(SleqpLPi* lp_interface,
                                         SleqpSparseMatrix* coeff_matrix);

SLEQP_NODISCARD
SLEQP_RETCODE sleqp_lpi_set_objective(SleqpLPi* lp_interface,
                                      double* objective);

SLEQP_NODISCARD
SLEQP_RETCODE sleqp_lpi_set_time_limit(SleqpLPi* lp_interface,
                                       double time_limit);

SLEQP_NODISCARD
SLEQP_RETCODE sleqp_lpi_save_basis(SleqpLPi* lp_interface,
                                   int index);

SLEQP_NODISCARD
SLEQP_RETCODE sleqp_lpi_restore_basis(SleqpLPi* lp_interface,
                                      int index);

SLEQP_NODISCARD
SLEQP_RETCODE sleqp_lpi_get_primal_sol(SleqpLPi* lp_interface,
                                       double* objective_value,
                                       double* solution_values);

SLEQP_NODISCARD
SLEQP_RETCODE sleqp_lpi_get_dual_sol(SleqpLPi* lp_interface,
                                     double* vars_dual,
                                     double* cons_dual);

SLEQP_NODISCARD
SLEQP_RETCODE sleqp_lpi_get_varstats(SleqpLPi* lp_interface,
                                     SLEQP_BASESTAT* variable_stats);

SLEQP_NODISCARD
SLEQP_RETCODE sleqp_lpi_get_consstats(SleqpLPi* lp_interface,
                                      SLEQP_BASESTAT* constraint_stats);

SleqpTimer* sleqp_lpi_get_solve_timer(SleqpLPi* lp_interface);

SLEQP_NODISCARD
SLEQP_RETCODE sleqp_lpi_get_basis_condition(SleqpLPi* lp_interface,
                                            bool* exact,
                                            double* condition);

SLEQP_NODISCARD
SLEQP_RETCODE sleqp_lpi_capture(SleqpLPi* lp_interface);

SLEQP_NODISCARD
SLEQP_RETCODE sleqp_lpi_release(SleqpLPi** star);

SLEQP_NODISCARD
SLEQP_RETCODE sleqp_lpi_create_default_interface(SleqpLPi** lp_interface,
                                                 int num_variables,
                                                 int num_constraints,
                                                 SleqpParams* params,
                                                 SleqpOptions* options);


#endif /* SLEQP_LPI_H */
