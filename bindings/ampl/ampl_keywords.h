#ifndef SLEQP_AMPL_KEYWORDS_H
#define SLEQP_AMPL_KEYWORDS_H

#include <asl.h>
#include <getstub.h>

#include "sleqp.h"

typedef struct SleqpAmplKeywords SleqpAmplKeywords;

SLEQP_RETCODE
sleqp_ampl_keywords_create(SleqpAmplKeywords** star,
                           SleqpOptions* options,
                           SleqpParams* params);

SLEQP_RETCODE
sleqp_ampl_keywords_get(SleqpAmplKeywords* ampl_keywords,
                        keyword** star,
                        int* num_keywords);

double
sleqp_ampl_keywords_iter_limit(SleqpAmplKeywords* ampl_keywords);

double
sleqp_ampl_keywords_time_limit(SleqpAmplKeywords* ampl_keywords);

SLEQP_RETCODE
sleqp_ampl_keywords_free(SleqpAmplKeywords** star);

#endif /* SLEQP_AMPL_KEYWORDS_H */
