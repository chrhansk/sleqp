#ifndef SLEQP_LPI_GUROBI_H
#define SLEQP_LPI_GUROBI_H

/**
 * @file sleqp_lpi_gurobi.h
 * @brief Definition of the Gurobi LP interface.
 **/

#include "sleqp_lpi.h"

#ifdef __cplusplus
extern "C" {
#endif

  SLEQP_RETCODE sleqp_lpi_gurobi_create_interface(SleqpLPi** lp_star,
                                                  int num_variables,
                                                  int num_constraints,
                                                  SleqpParams* params);

  SLEQP_RETCODE sleqp_lpi_create_default_interface(SleqpLPi** lp_interface,
                                                   int num_variables,
                                                   int num_constraints,
                                                   SleqpParams* params);


#ifdef __cplusplus
}
#endif

#endif /* SLEQP_LPI_GUROBI_H */
