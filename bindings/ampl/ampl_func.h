#ifndef SLEQP_AMPL_FUNC_H
#define SLEQP_AMPL_FUNC_H

#include "sleqp.h"

#include "ampl_data.h"

SLEQP_RETCODE
sleqp_ampl_func_create(SleqpFunc** star,
                       SleqpAmplData* ampl_data,
                       SleqpParams* params,
                       bool halt_on_error);

#endif /* SLEQP_AMPL_FUNC_H */
