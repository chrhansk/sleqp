#include "sleqp_bfgs.h"

#include <assert.h>
#include <math.h>

#include "sleqp_cmp.h"
#include "sleqp_mem.h"

static const double damping_factor = 0.2;

typedef struct BFGSBlock
{
  int dimension;

  SleqpSparseVec** step_diffs;
  SleqpSparseVec** grad_diffs;

  SleqpSparseVec** inner_prods;
  SleqpSparseVec** outer_prods;

  SleqpSparseVec* inner_cache;
  SleqpSparseVec* outer_cache;

  SleqpSparseVec* prod_cache;

  double initial_scale;

  // max size
  int num;
  // curr size
  int len;
  // curr index
  int curr;

  bool damped;

} BFGSBlock;

struct SleqpBFGSData
{
  int refcount;

  int num_variables;

  const SleqpParams* params;

  int num_blocks;
  BFGSBlock* blocks;

  SleqpSparseVec* grad_diff;
  SleqpSparseVec* step_diff;

  SleqpSparseVec* previous_grad;
  SleqpSparseVec* current_grad;

  SleqpSparseVec* block_grad_diff;
  SleqpSparseVec* block_step_diff;

  SleqpSparseVec* prod_cache;

  SleqpSparseVec* block_direction;
  SleqpSparseVec* block_prod;

  SleqpFunc* bfgs_func;
  SleqpFunc* func;
};

static SLEQP_RETCODE
bfgs_func_set_value(SleqpSparseVec* x,
                    SLEQP_VALUE_REASON reason,
                    int num_variables,
                    int* func_grad_nnz,
                    int* cons_val_nnz,
                    int* cons_jac_nnz,
                    void* func_data)
{
  SleqpBFGSData* bfgs_data = (SleqpBFGSData*) func_data;

  SLEQP_CALL(sleqp_func_set_value(bfgs_data->func,
                                  x,
                                  reason,
                                  func_grad_nnz,
                                  cons_val_nnz,
                                  cons_jac_nnz));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
bfgs_func_eval(int num_variables,
               const SleqpSparseVec* cons_indices,
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
                       const double* func_dual,
                       const SleqpSparseVec* direction,
                       const SleqpSparseVec* cons_duals,
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

  SleqpFuncCallbacks callbacks = {
    .set_value = bfgs_func_set_value,
    .func_eval = bfgs_func_eval,
    .hess_prod = bfgs_func_hess_product,
    .func_free = NULL
  };

  SLEQP_CALL(sleqp_func_create(fstar,
                               &callbacks,
                               num_variables,
                               bfgs_data));
  return SLEQP_OKAY;
}

SleqpFunc* sleqp_bfgs_get_func(SleqpBFGSData* data)
{
  return data->bfgs_func;
}

static SLEQP_RETCODE bfgs_block_create_at(BFGSBlock* block,
                                          int dimension,
                                          const SleqpParams* params,
                                          int num,
                                          bool damped)
{
  assert(dimension > 0);
  assert(num > 0);

  *block = (BFGSBlock) {0};

  block->dimension = dimension;

  block->num = num;
  block->len = 0;
  block->curr = -1;
  block->damped = damped;

  SLEQP_CALL(sleqp_calloc(&(block->step_diffs), num));

  SLEQP_CALL(sleqp_calloc(&(block->grad_diffs), num));

  SLEQP_CALL(sleqp_calloc(&(block->inner_prods), num));

  SLEQP_CALL(sleqp_calloc(&(block->outer_prods), num));

  for(int i = 0; i < num; ++i)
  {
    SLEQP_CALL(sleqp_sparse_vector_create(block->step_diffs + i,
                                          dimension,
                                          0));

    SLEQP_CALL(sleqp_sparse_vector_create(block->grad_diffs + i,
                                          dimension,
                                          0));

    SLEQP_CALL(sleqp_sparse_vector_create(block->inner_prods + i,
                                          dimension,
                                          0));

    SLEQP_CALL(sleqp_sparse_vector_create(block->outer_prods + i,
                                          dimension,
                                          0));
  }

  SLEQP_CALL(sleqp_sparse_vector_create(&(block->inner_cache),
                                        dimension,
                                        0));

  SLEQP_CALL(sleqp_sparse_vector_create(&(block->outer_cache),
                                        dimension,
                                        0));

  SLEQP_CALL(sleqp_sparse_vector_create(&(block->prod_cache),
                                        dimension,
                                        0));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE bfgs_block_free_at(BFGSBlock* block)
{
  const int num = block->num;

  SLEQP_CALL(sleqp_sparse_vector_free(&(block->prod_cache)));

  SLEQP_CALL(sleqp_sparse_vector_free(&(block->outer_cache)));
  SLEQP_CALL(sleqp_sparse_vector_free(&(block->inner_cache)));

  for(int i = 0; i < num; ++i)
  {
    SLEQP_CALL(sleqp_sparse_vector_free(block->outer_prods + i));
    SLEQP_CALL(sleqp_sparse_vector_free(block->inner_prods + i));

    SLEQP_CALL(sleqp_sparse_vector_free(block->grad_diffs + i));
    SLEQP_CALL(sleqp_sparse_vector_free(block->step_diffs + i));
  }

  sleqp_free(&block->outer_prods);
  sleqp_free(&block->inner_prods);

  sleqp_free(&block->grad_diffs);
  sleqp_free(&block->step_diffs);

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_bfgs_data_create(SleqpBFGSData** star,
                                     SleqpFunc* func,
                                     const SleqpParams* params,
                                     int num,
                                     bool damped)
{
  assert(num > 0);

  SLEQP_CALL(sleqp_malloc(star));

  SleqpBFGSData* data = *star;

  *data = (SleqpBFGSData) {0};

  data->refcount = 1;

  const int num_variables = sleqp_func_get_num_variables(func);

  SleqpHessianStruct* hessian_struct = sleqp_func_get_hess_struct(func);

  const int num_blocks = sleqp_hessian_struct_get_num_blocks(hessian_struct);

  data->num_variables = num_variables;
  data->num_blocks = num_blocks;
  data->func = func;
  data->params = params;

  sleqp_calloc(&data->blocks, num_blocks);

  for(int block = 0; block < num_blocks; ++block)
  {
    int begin, end;

    SLEQP_CALL(sleqp_hessian_struct_get_block_range(hessian_struct,
                                                    block,
                                                    &begin,
                                                    &end));

    int block_dimension = end - begin;

    SLEQP_CALL(bfgs_block_create_at(data->blocks + block,
                                    block_dimension,
                                    params,
                                    num,
                                    damped));
  }

  SLEQP_CALL(sleqp_sparse_vector_create(&data->grad_diff,
                                        num_variables,
                                        num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create(&data->step_diff,
                                        num_variables,
                                        num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create(&data->previous_grad,
                                        num_variables,
                                        num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create(&data->current_grad,
                                        num_variables,
                                        num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create(&data->block_grad_diff,
                                        num_variables,
                                        num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create(&data->block_step_diff,
                                        num_variables,
                                        num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create(&data->prod_cache,
                                        num_variables,
                                        num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create(&data->block_direction,
                                        num_variables,
                                        num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create(&data->block_prod,
                                        num_variables,
                                        num_variables));

  SLEQP_CALL(bfgs_func_create(&data->bfgs_func,
                              data->func,
                              data));

  return SLEQP_OKAY;
}

static int data_index(BFGSBlock* block, int index)
{
  if(block->len == 0)
  {
    return 0;
  }

  int current_index = index % (block->num);

  return (current_index < 0) ? (current_index + block->num) : current_index;
}

static
SLEQP_RETCODE bfgs_hess_prod_range(BFGSBlock* block,
                                   const SleqpParams* params,
                                   const SleqpSparseVec* direction,
                                   SleqpSparseVec* product,
                                   int final)
{
  const double eps = sleqp_params_get_eps(params);
  const int begin = block->curr - block->len + 1;

  SLEQP_CALL(sleqp_sparse_vector_copy(direction, product));

  SLEQP_CALL(sleqp_sparse_vector_scale(product, block->initial_scale));

  for(int prev = begin; prev <= final; ++prev)
  {
    int j = data_index(block, prev);

    SleqpSparseVec* inner_prod = block->inner_prods[j];
    SleqpSparseVec* outer_prod = block->outer_prods[j];

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
                                              block->inner_cache));

    SLEQP_CALL(sleqp_sparse_vector_add_scaled(block->inner_cache,
                                              outer_prod,
                                              1.,
                                              1. * outer_dot,
                                              eps,
                                              product));
  }

  return SLEQP_OKAY;
}

static
SLEQP_RETCODE bfgs_initial_scale(BFGSBlock* block,
                                 SleqpSparseVec* step_diff,
                                 SleqpSparseVec* grad_diff,
                                 double* initial_scale)
{
  double grad_step_dot;

  SLEQP_CALL(sleqp_sparse_vector_dot(grad_diff,
                                     step_diff,
                                     &grad_step_dot));

  const double step_diff_normsq = sleqp_sparse_vector_norm_sq(step_diff);

  assert(step_diff_normsq >= 0.);

  if(grad_step_dot == 0.)
  {
    (*initial_scale) = 1.;
    return SLEQP_OKAY;
  }

  (*initial_scale) = step_diff_normsq / grad_step_dot;

  if(block->damped)
  {
    // TODO: Find out if there is smoe better way
    // of applying the damping to the initial approximation
    (*initial_scale) = SLEQP_MAX((*initial_scale), 1.);
  }

  return SLEQP_OKAY;
}

static
SLEQP_RETCODE bfgs_compute_products(BFGSBlock* block,
                                    const SleqpParams* params)
{
  const double eps = sleqp_params_get_eps(params);

  assert(block->len > 0);

  {
    SLEQP_CALL(bfgs_initial_scale(block,
                                  block->step_diffs[block->curr],
                                  block->grad_diffs[block->curr],
                                  &block->initial_scale));

    assert(sleqp_pos(block->initial_scale, eps));
  }

  const int begin = block->curr - block->len + 1;

  for(int val = begin; val <= block->curr; ++val)
  {
    int i = data_index(block, val);

    SleqpSparseVec* current_step_diff = block->step_diffs[i];
    SleqpSparseVec* current_grad_diff = block->grad_diffs[i];

    SleqpSparseVec* direction = current_step_diff;
    SleqpSparseVec* product = block->prod_cache;

    SLEQP_CALL(bfgs_hess_prod_range(block,
                                    params,
                                    direction,
                                    product,
                                    val - 1));

    double bidir_product;

    SLEQP_CALL(sleqp_sparse_vector_dot(product,
                                       direction,
                                       &bidir_product));

    assert(bidir_product > 0);

    double dot_product;

    SLEQP_CALL(sleqp_sparse_vector_dot(current_step_diff,
                                       current_grad_diff,
                                       &dot_product));

    SleqpSparseVec* current_inner_prod = block->inner_prods[i];
    SleqpSparseVec* current_outer_prod = block->outer_prods[i];

    if(block->damped && (dot_product < damping_factor*bidir_product))
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

static
SLEQP_RETCODE bfgs_block_push(BFGSBlock* block,
                              const SleqpParams* params,
                              const SleqpSparseVec* step_diff,
                              const SleqpSparseVec* grad_diff)
{
  const int next = data_index(block, block->curr + 1);

  SLEQP_CALL(sleqp_sparse_vector_copy(step_diff,
                                      block->step_diffs[next]));

  SLEQP_CALL(sleqp_sparse_vector_copy(grad_diff,
                                      block->grad_diffs[next]));

  if(block->len < block->num)
  {
    ++block->len;
  }

  block->curr = next;

  SLEQP_CALL(bfgs_compute_products(block, params));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_bfgs_data_push(SleqpBFGSData* data,
                                   SleqpIterate* previous_iterate,
                                   SleqpIterate* current_iterate)
{
  SleqpSparseVec* cons_dual = sleqp_iterate_get_cons_dual(previous_iterate);

  const double eps = sleqp_params_get_eps(data->params);

  const int num_blocks = data->num_blocks;

  // Compute gradient difference
  {
    SLEQP_CALL(sleqp_sparse_matrix_trans_vector_product(sleqp_iterate_get_cons_jac(previous_iterate),
                                                        cons_dual,
                                                        eps,
                                                        data->prod_cache));

    SLEQP_CALL(sleqp_sparse_vector_add(data->prod_cache,
                                       sleqp_iterate_get_func_grad(previous_iterate),
                                       eps,
                                       data->previous_grad));
  }

  {
    SLEQP_CALL(sleqp_sparse_matrix_trans_vector_product(sleqp_iterate_get_cons_jac(current_iterate),
                                                        cons_dual,
                                                        eps,
                                                        data->prod_cache));

    SLEQP_CALL(sleqp_sparse_vector_add(data->prod_cache,
                                       sleqp_iterate_get_func_grad(current_iterate),
                                       eps,
                                       data->current_grad));
  }

  {
    SLEQP_CALL(sleqp_sparse_vector_add_scaled(data->previous_grad,
                                              data->current_grad,
                                              -1.,
                                              1.,
                                              eps,
                                              data->grad_diff));
  }

  // Compute primal difference
  SLEQP_CALL(sleqp_sparse_vector_add_scaled(sleqp_iterate_get_primal(previous_iterate),
                                            sleqp_iterate_get_primal(current_iterate),
                                            -1.,
                                            1.,
                                            eps,
                                            data->step_diff));

  int k_step = 0;
  int k_grad = 0;

  int offset = 0;

  for(int i = 0; i < num_blocks; ++i)
  {
    BFGSBlock* block = data->blocks + i;

    int next_offset = offset + block->dimension;

    SLEQP_CALL(sleqp_sparse_vector_clear(data->block_grad_diff));
    SLEQP_CALL(sleqp_sparse_vector_clear(data->block_step_diff));

    data->block_grad_diff->dim = block->dimension;
    data->block_step_diff->dim = block->dimension;

    while(k_step < data->step_diff->nnz &&
          data->step_diff->indices[k_step] < next_offset)
    {
      SLEQP_CALL(sleqp_sparse_vector_push(data->block_step_diff,
                                          data->step_diff->indices[k_step] - offset,
                                          data->step_diff->data[k_step]));

      ++k_step;
    }

    while(k_grad < data->grad_diff->nnz &&
          data->grad_diff->indices[k_grad] < next_offset)
    {
      SLEQP_CALL(sleqp_sparse_vector_push(data->block_grad_diff,
                                          data->grad_diff->indices[k_grad] - offset,
                                          data->grad_diff->data[k_grad]));

      ++k_grad;
    }

    const double step_normsq = sleqp_sparse_vector_norm_sq(data->block_step_diff);

    if(!sleqp_zero(step_normsq, eps))
    {
      SLEQP_CALL(bfgs_block_push(block,
                                 data->params,
                                 data->block_step_diff,
                                 data->block_grad_diff));
    }


    offset = next_offset;
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_bfgs_data_hess_prod(SleqpBFGSData* data,
                                        const SleqpSparseVec* direction,
                                        SleqpSparseVec* product)
{
  const int num_blocks = data->num_blocks;

  int offset = 0;
  int k_direction = 0, k_product = 0;

  SLEQP_CALL(sleqp_sparse_vector_reserve(product, data->num_variables));

  for(int i = 0; i < num_blocks; ++i)
  {
    BFGSBlock* block = data->blocks + i;

    int next_offset = offset + block->dimension;

    if(block->len == 0)
    {
      while(k_direction < direction->nnz &&
            direction->indices[k_direction] < next_offset)
      {
        SLEQP_CALL(sleqp_sparse_vector_push(product,
                                            direction->indices[k_direction],
                                            direction->data[k_direction]));

        ++k_direction;
        ++k_product;
      }
    }
    else
    {
      data->block_direction->dim = block->dimension;
      data->block_prod->dim = block->dimension;

      SLEQP_CALL(sleqp_sparse_vector_clear(data->block_direction));
      SLEQP_CALL(sleqp_sparse_vector_clear(data->block_prod));

      while(k_direction < direction->nnz &&
            direction->indices[k_direction] < next_offset)
      {
        SLEQP_CALL(sleqp_sparse_vector_push(data->block_direction,
                                            direction->indices[k_direction] - offset,
                                            direction->data[k_direction]));

        ++k_direction;
      }

      SLEQP_CALL(bfgs_hess_prod_range(block,
                                      data->params,
                                      data->block_direction,
                                      data->block_prod,
                                      block->curr));

      for(int k_block_prod = 0;
          k_block_prod < data->block_prod->nnz;
          ++k_block_prod)
      {
        SLEQP_CALL(sleqp_sparse_vector_push(product,
                                            data->block_prod->indices[k_block_prod] + offset,
                                            data->block_prod->data[k_block_prod]));

        ++k_product;
      }
    }

    offset = next_offset;
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE bfgs_data_free(SleqpBFGSData** star)
{
  SleqpBFGSData* data = *star;

  if(!data)
  {
    return SLEQP_OKAY;
  }

  SLEQP_CALL(sleqp_func_release(&data->bfgs_func));

  SLEQP_CALL(sleqp_sparse_vector_free(&data->block_prod));
  SLEQP_CALL(sleqp_sparse_vector_free(&data->block_direction));

  SLEQP_CALL(sleqp_sparse_vector_free(&data->prod_cache));

  SLEQP_CALL(sleqp_sparse_vector_free(&data->block_prod));
  SLEQP_CALL(sleqp_sparse_vector_free(&data->block_prod));

  SLEQP_CALL(sleqp_sparse_vector_free(&data->block_step_diff));
  SLEQP_CALL(sleqp_sparse_vector_free(&data->block_grad_diff));

  SLEQP_CALL(sleqp_sparse_vector_free(&data->current_grad));
  SLEQP_CALL(sleqp_sparse_vector_free(&data->previous_grad));

  SLEQP_CALL(sleqp_sparse_vector_free(&data->step_diff));
  SLEQP_CALL(sleqp_sparse_vector_free(&data->grad_diff));

  int num_blocks = data->num_blocks;

  for(int block = 0; block < num_blocks; ++block)
  {
    SLEQP_CALL(bfgs_block_free_at(data->blocks + block));
  }

  sleqp_free(&data->blocks);

  sleqp_free(star);

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_bfgs_data_capture(SleqpBFGSData* data)
{
  ++data->refcount;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_bfgs_data_release(SleqpBFGSData** star)
{
  SleqpBFGSData* bfgs_data = *star;

  if(!bfgs_data)
  {
    return SLEQP_OKAY;
  }

  if(--bfgs_data->refcount == 0)
  {
    SLEQP_CALL(bfgs_data_free(star));
  }

  *star = NULL;

  return SLEQP_OKAY;
}
