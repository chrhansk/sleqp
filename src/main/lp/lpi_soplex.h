#ifndef SLEQP_LPI_SOPLEX_H
#define SLEQP_LPI_SOPLEX_H

/**
 * @file lpi_soplex.h
 * @brief Definition of the SoPlex LP interface.
 **/

#include "lpi.h"

#ifdef __cplusplus
extern "C" {
#endif

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_lpi_soplex_create_interface(SleqpLPi** lp_star,
                                                  int num_variables,
                                                  int num_constraints,
                                                  SleqpParams* params);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_lpi_create_default_interface(SleqpLPi** lp_interface,
                                                   int num_variables,
                                                   int num_constraints,
                                                   SleqpParams* params);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_LPI_SOPLEX_H */
