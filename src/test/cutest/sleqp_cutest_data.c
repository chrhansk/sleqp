#include "sleqp_cutest_data.h"

#include <assert.h>

#include "cmp.h"
#include "log.h"
#include "mem.h"


static
SLEQP_RETCODE map_cutest_inf(double* values, int num_values)
{
  const double inf = sleqp_infinity();

  for(int i = 0; i < num_values; i++)
  {
    if(values[i] == -CUTE_INF)
    {
      values[i] = -inf;
    }
    else if(values[i] == CUTE_INF)
    {
      values[i] = inf;
    }
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_cutest_data_create(SleqpCutestData** star,
                                       integer funit,
                                       int num_variables,
                                       int num_constraints)
{
  integer cutest_status;

  // l_order = 2: general linear constraints should follow the
  //              general nonlinear ones
  integer e_order = 0, l_order = 2, v_order = 0;

  const bool is_constrained = (num_constraints != 0);

  SLEQP_CALL(sleqp_malloc(star));

  SleqpCutestData* data = *star;

  *data = (SleqpCutestData) {0};

  data->num_variables = num_variables;
  data->num_constraints = num_constraints;

  SLEQP_CALL(sleqp_alloc_array(&data->var_lb, num_variables));
  SLEQP_CALL(sleqp_alloc_array(&data->var_ub, num_variables));

  SLEQP_CALL(sleqp_alloc_array(&data->x, num_variables));

  if(is_constrained)
  {
    SLEQP_CALL(sleqp_alloc_array(&data->cons_lb, num_constraints));
    SLEQP_CALL(sleqp_alloc_array(&data->cons_ub, num_constraints));

    SLEQP_CALL(sleqp_alloc_array(&data->equatn, num_constraints));
    SLEQP_CALL(sleqp_alloc_array(&data->linear, num_constraints));

    SLEQP_CALL(sleqp_alloc_array(&data->v, num_constraints));

    CUTEST_csetup(&cutest_status, &funit, &cutest_iout, &cutest_io_buffer,
                  &num_variables, &num_constraints,
                  data->x,
                  data->var_lb,
                  data->var_ub,
                  data->v,
                  data->cons_lb,
                  data->cons_ub,
                  data->equatn,
                  data->linear,
                  &e_order, &l_order, &v_order);
  }
  else
  {
    CUTEST_usetup(&cutest_status, &funit, &cutest_iout, &cutest_io_buffer,
                  &num_variables, data->x, data->var_lb, data->var_ub);
  }

  int num_general = 0, i = 0;

  for(i = 0; i < num_constraints; ++i, ++num_general)
  {
    if(data->linear[i])
    {
      break;
    }
  }

  assert(num_general <= num_constraints);

  data->num_linear = num_constraints - num_general;

  // ensure that l_order = 2 is satisfied
  for(; i < num_constraints; ++i)
  {
    assert(data->linear[i]);
  }

  SLEQP_CALL(map_cutest_inf(data->var_lb, num_variables));
  SLEQP_CALL(map_cutest_inf(data->var_ub, num_variables));

  SLEQP_CALL(map_cutest_inf(data->cons_lb, num_constraints));
  SLEQP_CALL(map_cutest_inf(data->cons_ub, num_constraints));

  return SLEQP_OKAY;
}


SLEQP_RETCODE sleqp_cutest_data_free(SleqpCutestData** star)
{
  SleqpCutestData* data = *star;

  sleqp_free(&data->x);
  sleqp_free(&data->v);

  sleqp_free(&data->linear);
  sleqp_free(&data->equatn);

  sleqp_free(&data->cons_ub);
  sleqp_free(&data->cons_lb);

  sleqp_free(&data->var_ub);
  sleqp_free(&data->var_lb);

  sleqp_free(star);
  
  return SLEQP_OKAY;
}
