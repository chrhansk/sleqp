#include "sleqp_bfgs.h"

#include <assert.h>
#include <math.h>

#include "sleqp_mem.h"

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

  double* grad_step_dots;

  int num;
  int len;
  int curr;

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
                                     int num)
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

  SLEQP_CALL(sleqp_calloc(&(data->step_diffs), num));

  SLEQP_CALL(sleqp_calloc(&(data->grad_diffs), num));

  SLEQP_CALL(sleqp_calloc(&(data->inner_prods), num));

  SLEQP_CALL(sleqp_calloc(&(data->outer_prods), num));

  SLEQP_CALL(sleqp_calloc(&(data->grad_step_dots), num));

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

  int current_index = index % (data->len);

  return (current_index < 0) ? (current_index + data->len) : current_index;
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

  SLEQP_CALL(sleqp_sparse_vector_dot(data->grad_diffs[next],
                                     data->step_diffs[next],
                                     data->grad_step_dots + next));

  if(data->len < data->num)
  {
    ++data->len;
  }

  data->curr = next;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_bfgs_data_hess_prod(SleqpBFGSData* data,
                                        SleqpSparseVec* direction,
                                        SleqpSparseVec* product)
{
  if(data->len == 0)
  {
    SLEQP_CALL(sleqp_sparse_vector_copy(direction, product));

    return SLEQP_OKAY;
  }

  assert(data->len > 0);

  const double eps = sleqp_params_get_eps(data->params);

  double initial_scale;

  {
    double step_dot;

    SleqpSparseVec* current_step_diff = data->step_diffs[data->curr];
    SleqpSparseVec* current_grad_diff = data->grad_diffs[data->curr];

    SLEQP_CALL(sleqp_sparse_vector_dot(current_grad_diff,
                                       current_step_diff,
                                       &step_dot));

    const double grad_diff_normsq = sleqp_sparse_vector_normsq(current_grad_diff);

    assert(grad_diff_normsq > 0);

    initial_scale = step_dot / grad_diff_normsq;
  }

  SLEQP_CALL(sleqp_sparse_vector_copy(direction, product));

  SLEQP_CALL(sleqp_sparse_vector_scale(product, initial_scale));

  for(int val = data->curr - data->len + 1; val <= data->curr; ++val)
  {
    int i = data_index(data, val);

    SleqpSparseVec* current_step_diff = data->step_diffs[i];
    SleqpSparseVec* current_grad_diff = data->grad_diffs[i];

    {
      SLEQP_CALL(sleqp_sparse_vector_copy(data->grad_diffs[i],
                                          data->outer_prods[i]));

      SLEQP_CALL(sleqp_sparse_vector_scale(data->outer_prods[i],
                                           1./sqrt(data->grad_step_dots[i])));
    }

    {
      SLEQP_CALL(sleqp_sparse_vector_copy(data->step_diffs[i],
                                          data->inner_prods[i]));

      SLEQP_CALL(sleqp_sparse_vector_scale(data->inner_prods[i],
                                           initial_scale));
    }

    SLEQP_CALL(sleqp_sparse_vector_clear(data->outer_cache));

    for(int prev = data->curr - data->len + 1; prev < val; ++prev)
    {
      int j = data_index(data, prev);

      SleqpSparseVec* old_inner_prod = data->inner_prods[j];
      SleqpSparseVec* old_outer_prod = data->outer_prods[j];

      double old_inner_dot, old_outer_dot;

      SLEQP_CALL(sleqp_sparse_vector_dot(current_step_diff,
                                         old_inner_prod,
                                         &old_inner_dot));

      SLEQP_CALL(sleqp_sparse_vector_dot(current_step_diff,
                                         old_outer_prod,
                                         &old_outer_dot));

      SLEQP_CALL(sleqp_sparse_vector_add_scaled(old_inner_prod,
                                                old_outer_prod,
                                                -1.*old_inner_dot,
                                                old_outer_dot,
                                                eps,
                                                data->inner_cache));

      SLEQP_CALL(sleqp_sparse_vector_copy(data->inner_prods[i], data->outer_cache));


      SLEQP_CALL(sleqp_sparse_vector_add(data->inner_cache,
                                         data->outer_cache,
                                         eps,
                                         data->inner_prods[i]));

    }

    {
      double inner_dot;

      SLEQP_CALL(sleqp_sparse_vector_dot(current_step_diff,
                                         data->inner_prods[i],
                                         &inner_dot));

      SLEQP_CALL(sleqp_sparse_vector_scale(data->inner_prods[i],
                                           1./sqrt(inner_dot)));
    }

    {
      double inner_dir_dot, outer_dir_dot;

      SLEQP_CALL(sleqp_sparse_vector_dot(direction,
                                         data->inner_prods[i],
                                         &inner_dir_dot));

      SLEQP_CALL(sleqp_sparse_vector_dot(direction,
                                         data->outer_prods[i],
                                         &outer_dir_dot));

      SLEQP_CALL(sleqp_sparse_vector_add_scaled(product,
                                                data->inner_prods[i],
                                                1.,
                                                -1.*inner_dir_dot,
                                                eps,
                                                data->outer_cache));

      SLEQP_CALL(sleqp_sparse_vector_add_scaled(data->outer_cache,
                                                data->outer_prods[i],
                                                1.,
                                                outer_dir_dot,
                                                eps,
                                                product));

    }

  }

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

  SLEQP_CALL(sleqp_sparse_vector_free(&(data->outer_cache)));
  SLEQP_CALL(sleqp_sparse_vector_free(&(data->inner_cache)));

  for(int i = 0; i < num; ++i)
  {
    SLEQP_CALL(sleqp_sparse_vector_free(data->outer_prods + i));
    SLEQP_CALL(sleqp_sparse_vector_free(data->inner_prods + i));

    SLEQP_CALL(sleqp_sparse_vector_free(data->grad_diffs + i));
    SLEQP_CALL(sleqp_sparse_vector_free(data->step_diffs + i));
  }

  sleqp_free(&(data->grad_step_dots));

  sleqp_free(data->grad_diffs);
  sleqp_free(data->step_diffs);

  sleqp_free(star);

  return SLEQP_OKAY;
}
