#include "sleqp_sr1.h"

#include <math.h>

#include "sleqp_mem.h"

/*
 * Our SR1 implementation works in the following way:
 * The Hessian approximation \f$ B_k \f$ in the k-th iteration
 * is given as
 *
 * \f[ B_k := B_k^{0} +  \sum_{i=k-m}^{k - 1} \frac{a_i a_i^{T}}{a_i^{T} s_i}.  \f]
 *
 * where the inner products \f$ a_i \f$ are given as \f$a _i := y_i - B_i s_i \f$,
 * where the \f$ s_i \f$ are the step differences (\f$ s_i := x_{i + 1} - x_{i} \f$)
 * and the \f$  y_i \f$ are the differences of the gradients of the Lagrangian
 * with multipliers according to \f$ x_i \f$ of the accepted steps.
 *
 * We recompute the values of \f$ a_i \f$ and the
 * products \f$ a_i^{T} s_i \f$ whenever a new pair \f$ a_i, y_i \f$
 * is pushed. Note that some pairs are actually discarded when computing \f$ B_k \f$
 * (we don't want to divide by zero). The criterion is that we only use
 * pairs where
 *
 * \f[ |a_i^{T} s_i| >= r \|s_i\| \| a_i \|  \f]
 *
 * with a safeguard factor \f$ r \in (0, 1) \f$.
 */

static const double safeguard_factor = 1e-8;

struct SleqpSR1Data
{
  int num_variables;
  SleqpParams* params;

  SleqpSparseVec** step_diffs;
  SleqpSparseVec** grad_diffs;

  SleqpSparseVec** inner_prods;

  SleqpSparseVec* inner_cache;
  SleqpSparseVec* outer_cache;

  double* inner_dots;

  double initial_scale;

  int num;
  int len;
  int curr;

  SleqpFunc* sr1_func;
  SleqpFunc* func;
};

static SLEQP_RETCODE
sr1_func_set_value(SleqpSparseVec* x,
                   int num_variables,
                   int* func_grad_nnz,
                   int* cons_val_nnz,
                   int* cons_jac_nnz,
                   void* func_data)
{
  SleqpSR1Data* sr1_data = (SleqpSR1Data*) func_data;

  SLEQP_CALL(sleqp_func_set_value(sr1_data->func,
                                  x,
                                  func_grad_nnz,
                                  cons_val_nnz,
                                  cons_jac_nnz));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
sr1_func_eval(int num_variables,
              SleqpSparseVec* cons_indices,
              double* func_val,
              SleqpSparseVec* func_grad,
              SleqpSparseVec* cons_val,
              SleqpSparseMatrix* cons_jac,
              void* func_data)
{
  SleqpSR1Data* sr1_data = (SleqpSR1Data*) func_data;

  SLEQP_CALL(sleqp_func_eval(sr1_data->func,
                             cons_indices,
                             func_val,
                             func_grad,
                             cons_val,
                             cons_jac));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
sr1_func_hess_product(int num_variables,
                      double* func_dual,
                      SleqpSparseVec* direction,
                      SleqpSparseVec* cons_duals,
                      SleqpSparseVec* product,
                      void* func_data)
{
  SleqpSR1Data* sr1_data = (SleqpSR1Data*) func_data;

  SLEQP_CALL(sleqp_sr1_data_hess_prod(sr1_data,
                                      direction,
                                      product));

  return SLEQP_OKAY;
}


SLEQP_RETCODE sr1_func_create(SleqpFunc** fstar,
                              SleqpFunc* func,
                              SleqpSR1Data* sr1_data)
{

  const int num_variables = sleqp_func_get_num_variables(func);

  SLEQP_CALL(sleqp_func_create(fstar,
                               sr1_func_set_value,
                               sr1_func_eval,
                               sr1_func_hess_product,
                               num_variables,
                               sr1_data));
  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_sr1_data_create(SleqpSR1Data** star,
                                    SleqpFunc* func,
                                    SleqpParams* params,
                                    int num)
{
  assert(num > 0);

  SLEQP_CALL(sleqp_malloc(star));

  SleqpSR1Data* data = *star;

  *data = (SleqpSR1Data) {0};

  const int num_variables = sleqp_func_get_num_variables(func);

  data->num_variables = num_variables;
  data->params = params;

  data->num = num;

  data->len = 0;
  data->curr = -1;

  SLEQP_CALL(sleqp_calloc(&(data->step_diffs), num));

  SLEQP_CALL(sleqp_calloc(&(data->grad_diffs), num));

  SLEQP_CALL(sleqp_calloc(&(data->inner_prods), num));

  for(int i = 0; i < num; ++i)
  {
    SLEQP_CALL(sleqp_sparse_vector_create(data->step_diffs + i,
                                          num_variables,
                                          0));

    SLEQP_CALL(sleqp_sparse_vector_create(data->grad_diffs + i,
                                          num_variables,
                                          0));

    SLEQP_CALL(sleqp_sparse_vector_create(data->inner_prods + i,
                                          num_variables,
                                          0));
  }

  SLEQP_CALL(sleqp_sparse_vector_create(&(data->inner_cache),
                                        num_variables,
                                        0));

  SLEQP_CALL(sleqp_sparse_vector_create(&(data->outer_cache),
                                        num_variables,
                                        0));

  SLEQP_CALL(sleqp_calloc(&(data->inner_dots), num));

  SLEQP_CALL(sr1_func_create(&(data->sr1_func),
                             func,
                             data));

  data->func = func;

  return SLEQP_OKAY;
}

static int data_index(SleqpSR1Data* data, int index)
{
  if(data->len == 0)
  {
    return 0;
  }

  int current_index = index % (data->len);

  return (current_index < 0) ? (current_index + data->len) : current_index;
}

static
SLEQP_RETCODE sr1_initial_scale(SleqpSparseVec* step_diff,
                                SleqpSparseVec* grad_diff,
                                double* initial_scale)
{
  double step_dot;

  SLEQP_CALL(sleqp_sparse_vector_dot(grad_diff,
                                     step_diff,
                                     &step_dot));

  const double step_diff_normsq = sleqp_sparse_vector_normsq(step_diff);

  (*initial_scale) = step_dot / step_diff_normsq;

  return SLEQP_OKAY;
}

static
SLEQP_RETCODE sr1_compute_inner_products(SleqpSR1Data* data)
{
  assert(data->len > 0);

  const double eps = sleqp_params_get_eps(data->params);

  {
    SLEQP_CALL(sr1_initial_scale(data->step_diffs[data->curr],
                                 data->grad_diffs[data->curr],
                                 &data->initial_scale));
  }

  for(int val = data->curr - data->len + 1; val <= data->curr; ++val)
  {
    int i = data_index(data, val);

    SleqpSparseVec* current_step_diff = data->step_diffs[i];
    SleqpSparseVec* current_grad_diff = data->grad_diffs[i];

    SleqpSparseVec* current_inner_prod = data->inner_prods[i];

    {
      SLEQP_CALL(sleqp_sparse_vector_copy(data->step_diffs[i],
                                          data->outer_cache));

      SLEQP_CALL(sleqp_sparse_vector_scale(data->outer_cache,
                                           data->initial_scale));
    }

    for(int prev = data->curr - data->len + 1; prev < val; ++prev)
    {
      int j = data_index(data, prev);

      SleqpSparseVec* inner_prod = data->inner_prods[j];

      double inner_step_dot;

      SLEQP_CALL(sleqp_sparse_vector_dot(current_step_diff,
                                         inner_prod,
                                         &inner_step_dot));

      if(data->inner_dots[j] == 0.)
      {
        // this marks that the safeguard tripped
        continue;
      }

      double combined_factor = inner_step_dot / (data->inner_dots[j]);

      SLEQP_CALL(sleqp_sparse_vector_add_scaled(data->outer_cache,
                                                inner_prod,
                                                1.,
                                                combined_factor,
                                                eps,
                                                data->inner_cache));

      SLEQP_CALL(sleqp_sparse_vector_copy(data->inner_cache, data->outer_cache));

    }

    SLEQP_CALL(sleqp_sparse_vector_add_scaled(data->outer_cache,
                                              current_grad_diff,
                                              -1.,
                                              1.,
                                              eps,
                                              current_inner_prod));

    // Check the safeguard
    {
      double inner_dot;

      SLEQP_CALL(sleqp_sparse_vector_dot(current_step_diff,
                                         current_inner_prod,
                                         &inner_dot));

      double current_inner_norm = sleqp_sparse_vector_norm(current_inner_prod);
      double current_step_norm = sleqp_sparse_vector_norm(current_step_diff);

      assert(current_inner_norm >= 0.);
      assert(current_step_norm >= 0.);

      if(fabs(inner_dot) < safeguard_factor * current_inner_norm * current_step_norm)
      {
        data->inner_dots[i] = 0.;
      }
      else
      {
        assert(inner_dot != 0.);

        data->inner_dots[i] = inner_dot;
      }
    }

  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_sr1_data_push(SleqpSR1Data* data,
                                  SleqpIterate* old_iterate,
                                  SleqpIterate* new_iterate)
{
  SleqpSparseVec* cons_dual = old_iterate->cons_dual;

  const double eps = sleqp_params_get_eps(data->params);

  const int next = data_index(data, data->curr + 1);


  {
    SLEQP_CALL(sleqp_sparse_matrix_trans_vector_product(old_iterate->cons_jac,
                                                        cons_dual,
                                                        eps,
                                                        data->outer_cache));

    SLEQP_CALL(sleqp_sparse_vector_add(data->outer_cache,
                                       old_iterate->func_grad,
                                       eps,
                                       data->inner_cache));
  }

  {
    SLEQP_CALL(sleqp_sparse_matrix_trans_vector_product(new_iterate->cons_jac,
                                                        cons_dual,
                                                        eps,
                                                        data->outer_cache));

    SLEQP_CALL(sleqp_sparse_vector_add(data->outer_cache,
                                       new_iterate->func_grad,
                                       eps,
                                       data->step_diffs[next]));
  }

  SLEQP_CALL(sleqp_sparse_vector_add_scaled(data->inner_cache,
                                            data->step_diffs[next],
                                            -1.,
                                            1.,
                                            eps,
                                            data->grad_diffs[next]));

  SLEQP_CALL(sleqp_sparse_vector_add_scaled(old_iterate->x,
                                            new_iterate->x,
                                            -1.,
                                            1.,
                                            eps,
                                            data->step_diffs[next]));

  if(data->len < data->num)
  {
    ++data->len;
  }

  data->curr = next;

  SLEQP_CALL(sr1_compute_inner_products(data));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_sr1_data_hess_prod(SleqpSR1Data* data,
                                       SleqpSparseVec* direction,
                                       SleqpSparseVec* product)
{
  SLEQP_CALL(sleqp_sparse_vector_copy(direction, product));

  if(data->len == 0)
  {
    /* If we have not acquired any information yet,
     * then we simply use the identity without any scaling.
     */

    return SLEQP_OKAY;
  }

  const double eps = sleqp_params_get_eps(data->params);

  // Note: We have already computed the required inner products / dots

  SLEQP_CALL(sleqp_sparse_vector_scale(product, data->initial_scale));

  for(int val = data->curr - data->len + 1; val <= data->curr; ++val)
  {
    int i = data_index(data, val);

    if(data->inner_dots[i] == 0.)
    {
      // this marks that the safeguard tripped
      continue;
    }

    SleqpSparseVec* current_inner_prod = data->inner_prods[i];

    double direction_dot;

    SLEQP_CALL(sleqp_sparse_vector_dot(current_inner_prod, direction, &direction_dot));

    double combined_factor = direction_dot / (data->inner_dots[i]);

    SLEQP_CALL(sleqp_sparse_vector_add_scaled(product,
                                              current_inner_prod,
                                              1.,
                                              combined_factor,
                                              eps,
                                              data->inner_cache));

    SLEQP_CALL(sleqp_sparse_vector_copy(data->inner_cache, product));
  }

  return SLEQP_OKAY;
}

SleqpFunc* sleqp_sr1_get_func(SleqpSR1Data* data)
{
  return data->sr1_func;
}

SLEQP_RETCODE sleqp_sr1_data_free(SleqpSR1Data** star)
{
  SleqpSR1Data* data = *star;

  if(!data)
  {
    return SLEQP_OKAY;
  }

  SLEQP_CALL(sleqp_func_free(&(data->sr1_func)));

  sleqp_free(&(data->inner_dots));

  const int num = data->num;

  SLEQP_CALL(sleqp_sparse_vector_free(&(data->outer_cache)));
  SLEQP_CALL(sleqp_sparse_vector_free(&(data->inner_cache)));

  for(int i = 0; i < num; ++i)
  {
    SLEQP_CALL(sleqp_sparse_vector_free(data->inner_prods + i));

    SLEQP_CALL(sleqp_sparse_vector_free(data->grad_diffs + i));
    SLEQP_CALL(sleqp_sparse_vector_free(data->step_diffs + i));
  }

  sleqp_free(data->grad_diffs);
  sleqp_free(data->step_diffs);

  sleqp_free(star);

  return SLEQP_OKAY;
}
