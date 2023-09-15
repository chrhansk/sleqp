#ifndef SLEQP_LPI_HIGHS_H
#define SLEQP_LPI_HIGHS_H

/**
 * @file lpi_highs.h
 * @brief Definition of the HiGHS LP interface.
 **/

#include "lpi.h"

SLEQP_WARNUNUSED
SLEQP_RETCODE
sleqp_lpi_highs_create(SleqpLPi** lp_star,
                       int num_variables,
                       int num_constraints,
                       SleqpSettings* settings);

SLEQP_WARNUNUSED
SLEQP_RETCODE
sleqp_lpi_create_default(SleqpLPi** lp_interface,
                         int num_variables,
                         int num_constraints,
                         SleqpSettings* settings);

#endif /* SLEQP_LPI_HIGHS_H */
