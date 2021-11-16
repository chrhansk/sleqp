#ifndef SLEQP_LPI_GUROBI_H
#define SLEQP_LPI_GUROBI_H

/**
 * @file lpi_gurobi.h
 * @brief Definition of the Gurobi LP interface.
 **/

#include "lpi.h"

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_lpi_gurobi_create(SleqpLPi** lp_star,
                        int num_variables,
                        int num_constraints,
                        SleqpParams* params,
                        SleqpOptions* options);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_lpi_create_default(SleqpLPi** lp_interface,
                         int num_variables,
                         int num_constraints,
                         SleqpParams* params,
                         SleqpOptions* options);

#endif /* SLEQP_LPI_GUROBI_H */
