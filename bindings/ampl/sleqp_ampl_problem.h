#ifndef SLEQP_AMPL_PROBLEM_H
#define SLEQP_AMPL_PROBLEM_H

#include "sleqp.h"

#include "sleqp_ampl_data.h"

#ifdef __cplusplus
extern "C"
{
#endif

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_ampl_problem_create(SleqpProblem **star,
                                          SleqpAmplData *data,
                                          FILE *nl,
                                          SleqpParams *params);

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_AMPL_PROBLEM_H */
