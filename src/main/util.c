#include "util.h"

#include <float.h>

#include "cmp.h"
#include "error.h"
#include "fail.h"
#include "log.h"

static const double eps_factor = 10.;

SLEQP_RETCODE
sleqp_set_and_evaluate(SleqpProblem* problem,
                       SleqpIterate* iterate,
                       SLEQP_VALUE_REASON reason,
                       bool* reject)
{
  bool manual_reject = false;

  SLEQP_CALL(sleqp_problem_set_value(problem,
                                     sleqp_iterate_primal(iterate),
                                     reason,
                                     &manual_reject));

  if (reject)
  {
    *reject = manual_reject;
  }
  else if (manual_reject)
  {
    sleqp_raise(SLEQP_ILLEGAL_ARGUMENT, "Function is not allowed to raise");
  }

  SLEQP_CALL(sleqp_iterate_reserve(iterate, problem));

  SleqpVec* obj_grad = sleqp_iterate_obj_grad(iterate);
  SleqpMat* cons_jac = sleqp_iterate_cons_jac(iterate);
  SleqpVec* cons_val = sleqp_iterate_cons_val(iterate);

  double obj_val;

  SLEQP_CALL(
    sleqp_problem_eval(problem, &obj_val, obj_grad, cons_val, cons_jac));

  SLEQP_CALL(sleqp_iterate_set_obj_val(iterate, obj_val));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_direction_in_working_set(SleqpProblem* problem,
                               const SleqpIterate* iterate,
                               const SleqpVec* direction,
                               double* cache,
                               double eps,
                               bool* in_working_set)
{
  (*in_working_set) = true;

  SleqpMat* cons_jac = sleqp_iterate_cons_jac(iterate);

  SLEQP_CALL(sleqp_mat_mult_vec(cons_jac, direction, cache));

  const SleqpVec* lb = sleqp_problem_cons_lb(problem);
  const SleqpVec* ub = sleqp_problem_cons_ub(problem);

  SleqpWorkingSet* working_set = sleqp_iterate_working_set(iterate);

  SleqpVec* c = sleqp_iterate_cons_val(iterate);

  const int dim = lb->dim;

  int k_lb = 0, k_c = 0, k_ub = 0;

  while (k_lb < lb->nnz || k_c < c->nnz || k_ub < ub->nnz)
  {
    double lb_val = 0., c_val = 0, ub_val = 0.;

    bool valid_lb = (k_lb < lb->nnz);
    bool valid_c  = (k_c < c->nnz);
    bool valid_ub = (k_ub < ub->nnz);

    int idx = valid_c ? c->indices[k_c] : dim + 1;
    idx     = SLEQP_MIN(idx, valid_lb ? lb->indices[k_lb] : dim + 1);
    idx     = SLEQP_MIN(idx, valid_ub ? ub->indices[k_ub] : dim + 1);

    if (valid_lb && idx == lb->indices[k_lb])
    {
      lb_val = lb->data[k_lb++];
    }

    if (valid_c && idx == c->indices[k_c])
    {
      c_val = c->data[k_c++];
    }

    if (valid_ub && idx == ub->indices[k_ub])
    {
      ub_val = ub->data[k_ub++];
    }

    const SLEQP_ACTIVE_STATE state
      = sleqp_working_set_cons_state(working_set, idx);

    const double prod_val = c_val + cache[idx];

    if (state == SLEQP_INACTIVE)
    {
      continue;
    }
    else if (state == SLEQP_ACTIVE_UPPER && !sleqp_is_eq(prod_val, ub_val, eps))
    {
      (*in_working_set) = false;
      return SLEQP_OKAY;
    }
    else if (state == SLEQP_ACTIVE_LOWER && !sleqp_is_eq(prod_val, lb_val, eps))
    {
      (*in_working_set) = false;
      return SLEQP_OKAY;
    }
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_max_step_length(const SleqpVec* x,
                      const SleqpVec* d,
                      const SleqpVec* l,
                      const SleqpVec* u,
                      double* max_step_length)
{
  const int dim = x->dim;

  assert(l->dim == dim);
  assert(u->dim == dim);
  assert(d->dim == dim);

  assert((*max_step_length) > 0.);

  // upper bound

  {
    int k_x = 0, k_d = 0, k_u = 0;

    while (k_x < x->nnz || k_d < d->nnz || k_u < u->nnz)
    {
      bool valid_x = (k_x < x->nnz);
      bool valid_d = (k_d < d->nnz);
      bool valid_u = (k_u < u->nnz);

      int i_combined = valid_x ? x->indices[k_x] : dim + 1;
      i_combined = SLEQP_MIN(i_combined, valid_d ? d->indices[k_d] : dim + 1);
      i_combined = SLEQP_MIN(i_combined, valid_u ? u->indices[k_u] : dim + 1);

      valid_x = valid_x && x->indices[k_x] == i_combined;
      valid_d = valid_d && d->indices[k_d] == i_combined;
      valid_u = valid_u && u->indices[k_u] == i_combined;

      double x_value = valid_x ? x->data[k_x] : 0.;
      double d_value = valid_d ? d->data[k_d] : 0.;
      double u_value = valid_u ? u->data[k_u] : 0.;

      double diff = u_value - x_value;

      if ((d_value > 0. && diff > 0.) || (d_value < 0. && diff < 0.))
      {
        double current_bound = diff / d_value;

        (*max_step_length) = SLEQP_MIN(*max_step_length, current_bound);
      }

      if (valid_x)
      {
        ++k_x;
      }

      if (valid_d)
      {
        ++k_d;
      }

      if (valid_u)
      {
        ++k_u;
      }
    }
  }

  // lower bound

  {
    int k_x = 0, k_d = 0, k_l = 0;

    while (k_x < x->nnz || k_d < d->nnz || k_l < l->nnz)
    {
      bool valid_x = (k_x < x->nnz);
      bool valid_d = (k_d < d->nnz);
      bool valid_l = (k_l < l->nnz);

      int i_combined = valid_x ? x->indices[k_x] : dim + 1;
      i_combined = SLEQP_MIN(i_combined, valid_d ? d->indices[k_d] : dim + 1);
      i_combined = SLEQP_MIN(i_combined, valid_l ? l->indices[k_l] : dim + 1);

      valid_x = valid_x && x->indices[k_x] == i_combined;
      valid_d = valid_d && d->indices[k_d] == i_combined;
      valid_l = valid_l && l->indices[k_l] == i_combined;

      double x_value = valid_x ? x->data[k_x] : 0.;
      double d_value = valid_d ? d->data[k_d] : 0.;
      double l_value = valid_l ? l->data[k_l] : 0.;

      double diff = l_value - x_value;

      if ((d_value < 0. && diff < 0.) || (d_value > 0. && diff > 0.))
      {
        double current_bound = diff / d_value;

        (*max_step_length) = SLEQP_MIN(*max_step_length, current_bound);
      }

      if (valid_x)
      {
        ++k_x;
      }

      if (valid_d)
      {
        ++k_d;
      }

      if (valid_l)
      {
        ++k_l;
      }
    }
  }

  assert((*max_step_length) >= 0.);

  return SLEQP_OKAY;
}

double
sleqp_reduction_ratio(const double exact_reduction,
                      const double model_reduction)
{
  const double eps = eps_factor * DBL_EPSILON;

  const double corr_model_reduction = model_reduction - eps;
  const double corr_exact_reduction = exact_reduction - eps;

  // Safeguard against roundoff errors
  if (SLEQP_ABS(corr_model_reduction) <= eps
      && SLEQP_ABS(corr_exact_reduction) <= eps)
  {
    return 1.;
  }

  return corr_exact_reduction / corr_model_reduction;
}
