#ifndef SLEQP_AMPL
#define SLEQP_AMPL

#include "sleqp.h"
#include "asl.h"

#ifdef __cplusplus
extern "C"
{
#endif

  typedef struct
  {
    ASL *asl;
    bool is_constrained;
    int num_variables;
    int num_constraints;
    int num_linear;

    double *var_lb;
    double *var_ub;

    double *cons_lb;
    double *cons_ub;

    double *x;

  } SleqpAmplData;

  SLEQP_NODISCARD
  SLEQP_RETCODE map_ampl_inf(double *values, int num_values);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_ampl_data_create(SleqpAmplData **star,
                                       ASL *asl);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_ampl_data_free(SleqpAmplData **star);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_AMPL_DATA_H */
