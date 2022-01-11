#ifndef SLEQP_AMPL_SUFFIX_H
#define SLEQP_AMPL_SUFFIX_H

#include <asl.h>

typedef enum
{
  AMPL_SUFFIX_VARSTAT = 0,
  AMPL_SUFFIX_CONSSTAT,
  AMPL_NUM_SUFFIX
} AMPL_SUFFIX;

SufDecl*
sleqp_ampl_suffixes();

#endif /* SLEQP_AMPL_SUFFIX_H */
