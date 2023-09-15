#ifndef SLEQP_LPI_GUROBI_H
#define SLEQP_LPI_GUROBI_H

/**
 * @file lpi_gurobi.h
 * @brief Definition of the Gurobi LP interface.
 **/

#include "lpi.h"

SLEQP_WARNUNUSED
SLEQP_RETCODE
sleqp_lpi_gurobi_create(SleqpLPi** lp_star,
                        int num_variables,
                        int num_constraints,
                        SleqpSettings* settings);

SLEQP_WARNUNUSED
SLEQP_RETCODE
sleqp_lpi_create_default(SleqpLPi** lp_interface,
                         int num_variables,
                         int num_constraints,
                         SleqpSettings* settings);

#endif /* SLEQP_LPI_GUROBI_H */
