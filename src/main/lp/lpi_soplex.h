#ifndef SLEQP_LPI_SOPLEX_H
#define SLEQP_LPI_SOPLEX_H

/**
 * @file lpi_soplex.h
 * @brief Definition of the SoPlex LP interface.
 **/

#include "lpi.h"

SLEQP_WARNUNUSED
SLEQP_RETCODE
sleqp_lpi_soplex_create(SleqpLPi** lp_star,
                        int num_variables,
                        int num_constraints,
                        SleqpSettings* settings);

SLEQP_WARNUNUSED
SLEQP_RETCODE
sleqp_lpi_create_default(SleqpLPi** lp_interface,
                         int num_variables,
                         int num_constraints,
                         SleqpSettings* settings);

#endif /* SLEQP_LPI_SOPLEX_H */
