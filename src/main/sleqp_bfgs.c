#include "sleqp_bfgs.h"

#include <assert.h>
#include <math.h>

#include "sleqp_cmp.h"
#include "sleqp_mem.h"

static const double damping_factor = 0.2;

struct SleqpBFGSData
{
  int num_variables;
  SleqpParams* params;

  SleqpSparseVec** step_diffs;
  SleqpSparseVec** grad_diffs;

  SleqpSparseVec** inner_prods;
  SleqpSparseVec** outer_prods;

  SleqpSparseVec* inner_cache;
  SleqpSparseVec* outer_cache;

  SleqpSparseVec* step_cache;

  SleqpSparseVec* prod_cache;
  SleqpSparseVec* combined_cache;

  double initial_scale;

  // max size
  int num;
  // curr size
  int len;
  // curr index
  int curr;

  bool damped;

  SleqpFunc* bfgs_func;
  SleqpFunc* func;
};

static SLEQP_RETCODE
bfgs_func_set_value(SleqpSparseVec* x,
                    int num_variables,
                    int* func_grad_nnz,
                    int* cons_val_nnz,
                    int* cons_jac_nnz,
                    void* func_data)
{
  SleqpBFGSData* bfgs_data = (SleqpBFGSData*) func_data;

  SLEQP_CALL(sleqp_func_set_value(bfgs_data->func,
                                  x,
                                  func_grad_nnz,
                                  cons_val_nnz,
                                  cons_jac_nnz));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
bfgs_func_eval(int num_variables,
               SleqpSparseVec* cons_indices,
               double* func_val,
               SleqpSparseVec* func_grad,
               SleqpSparseVec* cons_val,
               SleqpSparseMatrix* cons_jac,
               void* func_data)
{
  SleqpBFGSData* bfgs_data = (SleqpBFGSData*) func_data;

  SLEQP_CALL(sleqp_func_eval(bfgs_data->func,
                             cons_indices,
                             func_val,
                             func_grad,
                             cons_val,
                             cons_jac));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
bfgs_func_hess_product(int num_variables,
                       double* func_dual,
                       SleqpSparseVec* direction,
                       SleqpSparseVec* cons_duals,
                       SleqpSparseVec* product,
                       void* func_data)
{
  SleqpBFGSData* bfgs_data = (SleqpBFGSData*) func_data;

  SLEQP_CALL(sleqp_bfgs_data_hess_prod(bfgs_data,
                                       direction,
                                       product));

  return SLEQP_OKAY;
}

SLEQP_RETCODE bfgs_func_create(SleqpFunc** fstar,
                               SleqpFunc* func,
                               SleqpBFGSData* bfgs_data)
{

  const int num_variables = sleqp_func_get_num_variables(func);

  SLEQP_CALL(sleqp_func_create(fstar,
                               bfgs_func_set_value,
                               bfgs_func_eval,
                               bfgs_func_hess_product,
                               num_variables,
                               bfgs_data));
  return SLEQP_OKAY;
}

SleqpFunc* sleqp_bfgs_get_func(SleqpBFGSData* data)
{
  return data->bfgs_func;
}

SLEQP_RETCODE sleqp_bfgs_data_create(SleqpBFGSData** star,
                                     SleqpFunc* func,
                                     SleqpParams* params,
                                     int num,
                                     bool damped)
{
  assert(num > 0);

  SLEQP_CALL(sleqp_malloc(star));

  SleqpBFGSData* data = *star;

  *data = (SleqpBFGSData) {0};

  const int num_variables = sleqp_func_get_num_variables(func);

  data->num_variables = num_variables;
  data->params = params;

  data->num = num;
  data->len = 0;
  data->curr = -1;
  data->damped = damped;

  SLEQP_CALL(sleqp_calloc(&(data->step_diffs), num));

  SLEQP_CALL(sleqp_calloc(&(data->grad_diffs), num));

  SLEQP_CALL(sleqp_calloc(&(data->inner_prods), num));

  SLEQP_CALL(sleqp_calloc(&(data->outer_prods), num));

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

    SLEQP_CALL(sleqp_sparse_vector_create(data->outer_prods + i,
                                          num_variables,
                                          0));
  }

  SLEQP_CALL(sleqp_sparse_vector_create(&(data->inner_cache),
                                        num_variables,
                                        0));

  SLEQP_CALL(sleqp_sparse_vector_create(&(data->outer_cache),
                                        num_variables,
                                        0));

  SLEQP_CALL(sleqp_sparse_vector_create(&(data->step_cache),
                                        num_variables,
                                        0));

  SLEQP_CALL(sleqp_sparse_vector_create(&(data->prod_cache),
                                        num_variables,
                                        0));

  SLEQP_CALL(sleqp_sparse_vector_create(&(data->combined_cache),
                                        num_variables,
                                        0));

  SLEQP_CALL(bfgs_func_create(&data->bfgs_func,
                              func,
                              data));

  data->func = func;

  return SLEQP_OKAY;
}

static int data_index(SleqpBFGSData* data, int index)
{
  if(data->len == 0)
  {
    return 0;
  }

  int current_index = index % (data->num);

  return (current_index < 0) ? (current_index + data->num) : current_index;
}

static
SLEQP_RETCODE bfgs_hess_prod_range(SleqpBFGSData* data,
                                   SleqpSparseVec* direction,
                                   SleqpSparseVec* product,
                                   int final)
{
  const double eps = sleqp_params_get_eps(data->params);
  const int begin = data->curr - data->len + 1;

  SLEQP_CALL(sleqp_sparse_vector_copy(direction, product));

  SLEQP_CALL(sleqp_sparse_vector_scale(product, data->initial_scale));

  for(int prev = begin; prev <= final; ++prev)
  {
    int j = data_index(data, prev);

    SleqpSparseVec* inner_prod = data->inner_prods[j];
    SleqpSparseVec* outer_prod = data->outer_prods[j];

    double inner_dot, outer_dot;

    SLEQP_CALL(sleqp_sparse_vector_dot(inner_prod,
                                       direction,
                                       &inner_dot));

    SLEQP_CALL(sleqp_sparse_vector_dot(outer_prod,
                                       direction,
                                       &outer_dot));

    SLEQP_CALL(sleqp_sparse_vector_add_scaled(product,
                                              inner_prod,
                                              1.,
                                              -1. * inner_dot,
                                              eps,
                                              data->inner_cache));

    SLEQP_CALL(sleqp_sparse_vector_add_scaled(data->inner_cache,
                                              outer_prod,
                                              1.,
                                              1. * outer_dot,
                                              eps,
                                              product));
  }

  return SLEQP_OKAY;
}

static
SLEQP_RETCODE bfgs_initial_scale(SleqpBFGSData* data,
                                 SleqpSparseVec* step_diff,
                                 SleqpSparseVec* grad_diff,
                                 double* initial_scale)
{
  double step_dot;

  SLEQP_CALL(sleqp_sparse_vector_dot(grad_diff,
                                     step_diff,
                                     &step_dot));

  const double grad_diff_normsq = sleqp_sparse_vector_normsq(grad_diff);

  assert(grad_diff_normsq > 0.);

  (*initial_scale) = grad_diff_normsq / step_dot;

  if(data->damped)
  {
    // TODO: Find out if there is smoe better way
    // of applying the damping to the initial approximation
    (*initial_scale) = SLEQP_MAX((*initial_scale), 1.);
  }

  return SLEQP_OKAY;
}

static
SLEQP_RETCODE bfgs_compute_products(SleqpBFGSData* data)
{
  const double eps = sleqp_params_get_eps(data->params);

  assert(data->len > 0);

  {
    SLEQP_CALL(bfgs_initial_scale(data,
                                  data->step_diffs[data->curr],
                                  data->grad_diffs[data->curr],
                                  &data->initial_scale));

    assert(sleqp_pos(data->initial_scale, eps));
  }

  const int begin = data->curr - data->len + 1;

  for(int val = begin; val <= data->curr; ++val)
  {
    int i = data_index(data, val);

    SleqpSparseVec* current_step_diff = data->step_diffs[i];
    SleqpSparseVec* current_grad_diff = data->grad_diffs[i];

    SleqpSparseVec* direction = current_step_diff;
    SleqpSparseVec* product = data->prod_cache;

    SLEQP_CALL(bfgs_hess_prod_range(data, direction, product, val - 1));

    double bidir_product;

    SLEQP_CALL(sleqp_sparse_vector_dot(product,
                                       direction,
                                       &bidir_product));

    assert(bidir_product > 0);

    double dot_product;

    SLEQP_CALL(sleqp_sparse_vector_dot(current_step_diff,
                                       current_grad_diff,
                                       &dot_product));

    SleqpSparseVec* current_inner_prod = data->inner_prods[i];
    SleqpSparseVec* current_outer_prod = data->outer_prods[i];

    if(data->damped && (dot_product < damping_factor*bidir_product))
    {
      double combination_factor = (1. - damping_factor) *
        bidir_product / (bidir_product - dot_product);

      assert(sleqp_pos(combination_factor, eps));
      assert(sleqp_lt(combination_factor, 1., eps));

      SLEQP_CALL(sleqp_sparse_vector_add_scaled(current_grad_diff,
                                                product,
                                                combination_factor,
                                                1. - combination_factor,
                                                eps,
                                                current_outer_prod));

      SLEQP_CALL(sleqp_sparse_vector_dot(current_outer_prod,
                                         current_step_diff,
                                         &dot_product));
    }
    else
    {
      SLEQP_CALL(sleqp_sparse_vector_copy(current_grad_diff, current_outer_prod));
    }

    // set outer product
    {
      assert(dot_product > 0);

      SLEQP_CALL(sleqp_sparse_vector_scale(current_outer_prod, 1./sqrt(dot_product)));
    }

    // set inner product
    {
      SLEQP_CALL(sleqp_sparse_vector_copy(product, current_inner_prod));

      SLEQP_CALL(sleqp_sparse_vector_scale(current_inner_prod, 1./sqrt(bidir_product)));
    }
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_bfgs_data_push(SleqpBFGSData* data,
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

  SleqpSparseVec* next_grad_diff = data->grad_diffs[next];
  SleqpSparseVec* next_step_diff = data->step_diffs[next];

  {
  SLEQP_CALL(sleqp_sparse_matrix_trans_vector_product(new_iterate->cons_jac,
                                                      cons_dual,
                                                      eps,
                                                      data->outer_cache));

  SLEQP_CALL(sleqp_sparse_vector_add(data->outer_cache,
                                     new_iterate->func_grad,
                                     eps,
                                     data->step_cache));
  }

  SLEQP_CALL(sleqp_sparse_vector_add_scaled(data->inner_cache,
                                            data->step_cache,
                                            -1.,
                                            1.,
                                            eps,
                                            next_grad_diff));

  SLEQP_CALL(sleqp_sparse_vector_add_scaled(old_iterate->primal,
                                            new_iterate->primal,
                                            -1.,
                                            1.,
                                            eps,
                                            next_step_diff));

  if(data->len < data->num)
  {
    ++data->len;
  }

  data->curr = next;

  SLEQP_CALL(bfgs_compute_products(data));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_bfgs_data_hess_prod(SleqpBFGSData* data,
                                        SleqpSparseVec* direction,
                                        SleqpSparseVec* product)
{
  if(data->len == 0)
  {
    /* If we have not acquired any information yet,
     * then we simply use the identity without any scaling.
     */

    SLEQP_CALL(sleqp_sparse_vector_copy(direction, product));

    return SLEQP_OKAY;
  }

  assert(data->len > 0);

  SLEQP_CALL(bfgs_hess_prod_range(data, direction, product, data->curr));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_bfgs_data_free(SleqpBFGSData** star)
{
  SleqpBFGSData* data = *star;

  if(!data)
  {
    return SLEQP_OKAY;
  }

  const int num = data->num;

  SLEQP_CALL(sleqp_func_free(&(data->bfgs_func)));

  SLEQP_CALL(sleqp_sparse_vector_free(&(data->combined_cache)));

  SLEQP_CALL(sleqp_sparse_vector_free(&(data->prod_cache)));

  SLEQP_CALL(sleqp_sparse_vector_free(&(data->step_cache)));

  SLEQP_CALL(sleqp_sparse_vector_free(&(data->outer_cache)));
  SLEQP_CALL(sleqp_sparse_vector_free(&(data->inner_cache)));

  for(int i = 0; i < num; ++i)
  {
    SLEQP_CALL(sleqp_sparse_vector_free(data->outer_prods + i));
    SLEQP_CALL(sleqp_sparse_vector_free(data->inner_prods + i));

    SLEQP_CALL(sleqp_sparse_vector_free(data->grad_diffs + i));
    SLEQP_CALL(sleqp_sparse_vector_free(data->step_diffs + i));
  }

  sleqp_free(data->grad_diffs);
  sleqp_free(data->step_diffs);

  sleqp_free(star);

  return SLEQP_OKAY;
}
