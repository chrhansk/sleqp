#ifndef SLEQP_CUTEST_DATA_H
#define SLEQP_CUTEST_DATA_H

#include "pub_types.h"

#include "sleqp_cutest_defs.h"

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct
  {
    int num_variables;
    int num_constraints;
    int num_linear;

    double* var_lb;
    double* var_ub;

    double* cons_lb;
    double* cons_ub;

    logical* equatn;
    logical* linear;
    double* v;
    double* x;

  } SleqpCutestData;


  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_cutest_data_create(SleqpCutestData** star,
                                         integer funit,
                                         int num_variables,
                                         int num_constraints);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_cutest_data_free(SleqpCutestData** star);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_CUTEST_DATA_H */
