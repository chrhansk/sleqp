#include "ampl_suffix.h"

#define AMPL_SUFFIX_STAT_NAME "sstatus"

static SufDecl suffixes[AMPL_NUM_SUFFIX]
  = {[AMPL_SUFFIX_VARSTAT]  = {.name   = AMPL_SUFFIX_STAT_NAME,
                               .table  = NULL,
                               .kind   = ASL_Sufkind_var,
                               .nextra = 0},
     [AMPL_SUFFIX_CONSSTAT] = {.name   = AMPL_SUFFIX_STAT_NAME,
                               .table  = NULL,
                               .kind   = ASL_Sufkind_con,
                               .nextra = 0}};

SufDecl*
sleqp_ampl_suffixes()
{
  return suffixes;
}
