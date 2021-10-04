#ifndef SLEQP_LPI_HIGHS_H
#define SLEQP_LPI_HIGHS_H

/**
 * @file lpi_highs.h
 * @brief Definition of the HiGHS LP interface.
 **/

#include "lpi.h"

#ifdef __cplusplus
extern "C" {
#endif

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_lpi_highs_create_interface(SleqpLPi** lp_star,
                                                 int num_variables,
                                                 int num_constraints,
                                                 SleqpParams* params,
                                                 SleqpOptions* options);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_lpi_create_default_interface(SleqpLPi** lp_interface,
                                                   int num_variables,
                                                   int num_constraints,
                                                   SleqpParams* params,
                                                   SleqpOptions* options);


#ifdef __cplusplus
}
#endif

#endif /* SLEQP_LPI_HIGHS_H */
