#include "ampl_data.h"

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
sleqp_ampl_data_create(SleqpAmplData** star, ASL* asl, FILE* nl)
{
  SLEQP_CALL(sleqp_malloc(star));

  SleqpAmplData* data = *star;

  *data = (SleqpAmplData){0};

  data->asl = asl;

  int num_variables   = n_var;
  int num_constraints = n_con;
  int num_linear      = n_con - nlc;

  data->num_variables   = num_variables;
  data->num_constraints = num_constraints;
  data->num_linear      = num_linear;
  data->is_constrained  = (num_constraints != 0);

  SLEQP_CALL(sleqp_alloc_array(&data->var_lb, num_variables));
  SLEQP_CALL(sleqp_alloc_array(&data->var_ub, num_variables));

  if (data->is_constrained)
  {
    SLEQP_CALL(sleqp_alloc_array(&data->cons_lb, num_constraints));
    SLEQP_CALL(sleqp_alloc_array(&data->cons_ub, num_constraints));
  }

  SLEQP_CALL(sleqp_alloc_array(&data->x, num_variables));

  // set ASL pointer to allocated data
  X0    = data->x;
  LUv   = data->var_lb;
  Uvx   = data->var_ub;
  LUrhs = data->cons_lb;
  Urhsx = data->cons_ub;

  // read remaining data from stub and set functions
  int retcode = pfgh_read(nl, ASL_return_read_err);

  if (retcode != ASL_readerr_none)
  {
    sleqp_raise(SLEQP_INTERNAL_ERROR, "Error %d in reading nl file", retcode);
  }

  if (data->is_constrained)
  {
    data->jac_nnz = nzc;

    SLEQP_CALL(sleqp_alloc_array(&data->cons_val, num_constraints));

    SLEQP_CALL(sleqp_alloc_array(&data->jac_rows, data->jac_nnz));
    SLEQP_CALL(sleqp_alloc_array(&data->jac_cols, data->jac_nnz));
    SLEQP_CALL(sleqp_alloc_array(&data->jac_vals, data->jac_nnz));

    for (int i = 0; i < num_constraints; ++i)
    {
      for (cgrad* cg = Cgrad[i]; cg; cg = cg->next)
      {
        assert(cg->goff >= 0);
        assert(cg->goff < data->jac_nnz);

        data->jac_rows[cg->goff] = i;
        data->jac_cols[cg->goff] = cg->varno;
      }
    }

    for (int i = 0; i + 1 < data->jac_nnz; ++i)
    {
      assert(data->jac_cols[i] >= 0);
      assert(data->jac_cols[i] < num_variables);

      assert(data->jac_rows[i] >= 0);
      assert(data->jac_rows[i] < num_constraints);

      assert(data->jac_cols[i] <= data->jac_cols[i + 1]);

      if (data->jac_cols[i] == data->jac_cols[i + 1])
      {
        assert(data->jac_rows[i] < data->jac_rows[i + 1]);
      }
    }
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_ampl_data_free(SleqpAmplData** star)
{
  SleqpAmplData* data = *star;

  sleqp_free(&data->x);

  sleqp_free(&data->jac_vals);
  sleqp_free(&data->jac_cols);
  sleqp_free(&data->jac_rows);

  sleqp_free(&data->cons_ub);
  sleqp_free(&data->cons_lb);
  sleqp_free(&data->cons_val);

  sleqp_free(&data->var_ub);
  sleqp_free(&data->var_lb);

  sleqp_free(star);

  return SLEQP_OKAY;
}
