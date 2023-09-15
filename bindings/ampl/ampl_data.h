#ifndef SLEQP_AMPL_DATA_H
#define SLEQP_AMPL_DATA_H

#include <asl.h>

#include "sleqp.h"

typedef struct
{
  ASL* asl;
  bool is_constrained;
  int num_variables;
  int num_constraints;
  int num_linear;

  double* var_lb;
  double* var_ub;

  double* cons_val;

  double* cons_lb;
  double* cons_ub;

  int jac_nnz;
  int* jac_rows;
  int* jac_cols;
  double* jac_vals;

  double* x;

} SleqpAmplData;

SLEQP_WARNUNUSED
SLEQP_RETCODE
map_ampl_inf(double* values, int num_values);

SLEQP_WARNUNUSED
SLEQP_RETCODE
sleqp_ampl_data_create(SleqpAmplData** star, ASL* asl, FILE* nl);

SLEQP_WARNUNUSED
SLEQP_RETCODE
sleqp_ampl_data_free(SleqpAmplData** star);

#endif /* SLEQP_AMPL_DATA_H */
