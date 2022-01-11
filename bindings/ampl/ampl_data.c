#include "ampl_data.h"
#include "ampl_mem.h"

#include <assert.h>

SLEQP_RETCODE
map_ampl_inf(double* values, int num_values)
{
  const double inf = sleqp_infinity();

  for (int i = 0; i < num_values; i++)
  {
    if (values[i] == negInfinity)
    {
      values[i] = -inf;
    }
    else if (values[i] == Infinity)
    {
      values[i] = inf;
    }
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_ampl_data_create(SleqpAmplData** star, ASL* asl)
{

  SLEQP_CALL(sleqp_ampl_malloc(star));

  SleqpAmplData* data = *star;

  *data = (SleqpAmplData){0};

  data->asl = asl;

  int num_variables   = n_var;
  int num_constraints = n_con;
  int num_general     = nlc;
  int num_linear      = n_con - nlc;

  data->num_variables   = num_variables;
  data->num_constraints = num_constraints;
  data->num_linear      = num_linear;
  data->is_constrained  = (num_constraints != 0);

  SLEQP_CALL(sleqp_ampl_alloc_array(&data->var_lb, num_variables));
  SLEQP_CALL(sleqp_ampl_alloc_array(&data->var_ub, num_variables));

  SLEQP_CALL(sleqp_ampl_alloc_array(&data->x, num_variables));

  if (data->is_constrained)
  {
    SLEQP_CALL(sleqp_ampl_alloc_array(&data->cons_lb, num_constraints));
    SLEQP_CALL(sleqp_ampl_alloc_array(&data->cons_ub, num_constraints));
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_ampl_data_free(SleqpAmplData** star)
{
  SleqpAmplData* data = *star;

  sleqp_ampl_free(&data->x);

  sleqp_ampl_free(&data->cons_ub);
  sleqp_ampl_free(&data->cons_lb);

  sleqp_ampl_free(&data->var_ub);
  sleqp_ampl_free(&data->var_lb);

  sleqp_ampl_free(star);

  return SLEQP_OKAY;
}
