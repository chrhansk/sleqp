#include "dual_estimation_mixed.h"

#include "cmp.h"
#include "dual_estimation_lp.h"
#include "dual_estimation_lsq.h"
#include "fail.h"
#include "mem.h"

typedef struct
{
  SleqpDualEstimation* estimation_lsq;
  SleqpDualEstimation* estimation_lp;

  SleqpVec* cons_dual_lsq;
  SleqpVec* vars_dual_lsq;

  SleqpVec* cons_dual_lp;
  SleqpVec* vars_dual_lp;

} EstimationData;

static SLEQP_RETCODE
merge(const SleqpVec* dual_lsq,
      int num_active,
      const SleqpVec* dual_lp,
      SleqpVec* dual)
{
  SLEQP_CALL(sleqp_vec_clear(dual));
  SLEQP_CALL(sleqp_vec_reserve(dual, num_active));

  const int dim = dual_lsq->dim;

  assert(dim == dual_lp->dim);
  assert(dim == dual->dim);

  int k_lsq = 0, k_lp = 0;

  while (true)
  {
    const bool valid_lsq = (k_lsq < dual_lsq->nnz);
    const bool valid_lp  = (k_lp < dual_lp->nnz);

    if (!(valid_lp || valid_lsq))
    {
      break;
    }

    const int i_lsq = (valid_lsq) ? dual_lsq->indices[k_lsq] : dim;
    const int i_lp  = (valid_lp) ? dual_lp->indices[k_lp] : dim;

    const int i = SLEQP_MIN(i_lsq, i_lp);

    if (i_lsq <= i_lp)
    {
      assert(dual_lsq->data[k_lsq] != 0.);
      SLEQP_CALL(sleqp_vec_push(dual, i_lsq, dual_lsq->data[k_lsq]));
    }
    else
    {
      assert(dual_lp->data[k_lp] != 0.);
      SLEQP_CALL(sleqp_vec_push(dual, i_lp, dual_lp->data[k_lp]));
    }

    if (i == i_lsq)
    {
      ++k_lsq;
    }
    if (i == i_lp)
    {
      ++k_lp;
    }
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
copy_or_merge(bool clipped,
              int num_active,
              SleqpVec* dual_lsq,
              SleqpVec* dual_lp,
              SleqpVec* dual)
{
  if (!clipped)
  {
    return sleqp_vec_copy(dual_lsq, dual);
  }

  return merge(dual_lsq, num_active, dual_lp, dual);
}

static SLEQP_RETCODE
estimate_duals(const SleqpIterate* iterate,
               SleqpVec* cons_dual,
               SleqpVec* vars_dual,
               void* data)
{
  EstimationData* estimation_data = (EstimationData*)data;

  SleqpWorkingSet* working_set = sleqp_iterate_working_set(iterate);

  int num_clipped_vars;
  int num_clipped_cons;

  SLEQP_CALL(sleqp_estimate_duals_lsq(estimation_data->estimation_lsq,
                                      iterate,
                                      estimation_data->cons_dual_lsq,
                                      estimation_data->vars_dual_lsq,
                                      &num_clipped_vars,
                                      &num_clipped_cons));

  const bool clipped = (num_clipped_cons > 0) || (num_clipped_vars > 0);

  if (clipped)
  {
    SLEQP_CALL(sleqp_estimate_duals(estimation_data->estimation_lp,
                                    iterate,
                                    estimation_data->cons_dual_lp,
                                    estimation_data->vars_dual_lp));
  }

  const int num_active_vars = sleqp_working_set_num_active_vars(working_set);

  SLEQP_CALL(copy_or_merge(num_clipped_vars > 0,
                           num_active_vars,
                           estimation_data->vars_dual_lsq,
                           estimation_data->vars_dual_lp,
                           vars_dual));

  const int num_active_cons = sleqp_working_set_num_active_cons(working_set);

  SLEQP_CALL(copy_or_merge(num_clipped_cons > 0,
                           num_active_cons,
                           estimation_data->cons_dual_lsq,
                           estimation_data->cons_dual_lp,
                           cons_dual));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
estimation_free(void* data)
{
  EstimationData* estimation_data = (EstimationData*)data;

  SLEQP_CALL(sleqp_vec_free(&estimation_data->vars_dual_lp));
  SLEQP_CALL(sleqp_vec_free(&estimation_data->cons_dual_lp));

  SLEQP_CALL(sleqp_vec_free(&estimation_data->vars_dual_lsq));
  SLEQP_CALL(sleqp_vec_free(&estimation_data->cons_dual_lsq));

  SLEQP_CALL(sleqp_dual_estimation_release(&(estimation_data->estimation_lp)));

  SLEQP_CALL(sleqp_dual_estimation_release(&(estimation_data->estimation_lsq)));

  sleqp_free(&estimation_data);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
estimation_data_create(EstimationData** star,
                       SleqpProblem* problem,
                       SleqpCauchy* cauchy,
                       SleqpAugJac* aug_jac)
{
  SLEQP_CALL(sleqp_malloc(star));

  const int num_vars = sleqp_problem_num_vars(problem);
  const int num_cons = sleqp_problem_num_cons(problem);

  EstimationData* estimation_data = *star;

  *estimation_data = (EstimationData){0};

  SLEQP_CALL(sleqp_dual_estimation_lsq_create(&estimation_data->estimation_lsq,
                                              problem,
                                              aug_jac));

  SLEQP_CALL(
    sleqp_dual_estimation_lp_create(&estimation_data->estimation_lp, cauchy));

  SLEQP_CALL(sleqp_vec_create_empty(&estimation_data->cons_dual_lsq, num_cons));

  SLEQP_CALL(sleqp_vec_create_empty(&estimation_data->vars_dual_lsq, num_vars));

  SLEQP_CALL(sleqp_vec_create_empty(&estimation_data->cons_dual_lp, num_cons));

  SLEQP_CALL(sleqp_vec_create_empty(&estimation_data->vars_dual_lp, num_vars));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_dual_estimation_mixed_create(SleqpDualEstimation** star,
                                   SleqpProblem* problem,
                                   SleqpCauchy* cauchy,
                                   SleqpAugJac* aug_jac)
{
  SleqpDualEstimationCallbacks callbacks = {
    .estimate_duals  = estimate_duals,
    .estimation_free = estimation_free,
  };

  EstimationData* estimation_data;

  SLEQP_CALL(
    estimation_data_create(&estimation_data, problem, cauchy, aug_jac));

  SLEQP_CALL(
    sleqp_dual_estimation_create(star, &callbacks, (void*)estimation_data));

  return SLEQP_OKAY;
}
