#ifndef SLEQP_AMPL_PROBLEM_H
#define SLEQP_AMPL_PROBLEM_H

#include "sleqp.h"

#include "ampl_data.h"

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_ampl_problem_create(SleqpProblem** star,
                          SleqpAmplData* data,
                          FILE* nl,
                          SleqpParams* params);

#endif /* SLEQP_AMPL_PROBLEM_H */
