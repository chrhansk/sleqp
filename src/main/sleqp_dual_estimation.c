#include "sleqp_dual_estimation.h"

#include "sleqp_mem.h"

struct SleqpDualEstimationData
{
  size_t num_variables;
  size_t num_constraints;

  int* active_vars;
  int* active_cons;
};


SLEQP_RETCODE sleqp_dual_estimation_data_create(SleqpDualEstimationData** star,
                                                SleqpProblem* problem)
{
  SLEQP_CALL(sleqp_malloc(star));

  SleqpDualEstimationData* estimation_data = *star;

  estimation_data->num_variables = problem->num_variables;
  estimation_data->num_constraints = problem->num_constraints;

  SLEQP_CALL(sleqp_calloc(&estimation_data->active_vars, estimation_data->num_variables));
  SLEQP_CALL(sleqp_calloc(&estimation_data->active_cons, estimation_data->num_constraints));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_dual_estimation_compute(SleqpDualEstimationData* data,
                                            SleqpIterate* iterate,
                                            SleqpActiveSet* active_set)
{
  for(size_t i = 0; i < data->num_variables; ++i)
  {

  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_dual_estimation_data_free(SleqpDualEstimationData** star)
{
  SleqpDualEstimationData* estimation_data = *star;

  sleqp_free(&estimation_data->active_vars);
  sleqp_free(&estimation_data->active_cons);

  sleqp_free(star);

  return SLEQP_OKAY;
}
