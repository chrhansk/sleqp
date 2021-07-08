#include "sr1.h"

#include <math.h>

#include "cmp.h"
#include "log.h"
#include "mem.h"

#include "sparse/sparse_matrix.h"

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
 * products \f$ a_i^{T} s_i \f$ whenever a new pair \f$ s_i, y_i \f$
 * is pushed. Note that some pairs are actually discarded when computing \f$ B_k \f$
 * (we don't want to divide by zero). The criterion is that we only use
 * pairs where
 *
 * \f[ |a_i^{T} s_i| >= r \|s_i\| \| a_i \|  \f]
 *
 * with a safeguard factor \f$ r \in (0, 1) \f$.
 *
 * As initial approximation, we choose the identity scaled by
 * a factor of \f$ (s_k^{T} s_k) / (y_k^{T} s_k) \f$, as suggested in
 *
 * "An SR1/BFGS SQP algorithm for nonconvex nonlinear
 *  programs with block-diagonal Hessian matrix"
 *
 */

static const double safeguard_factor = 1e-8;

typedef struct SR1Block
{
  int dimension;

  SleqpSparseVec** step_diffs;
  SleqpSparseVec** grad_diffs;

  SleqpSparseVec** inner_prods;

  double* inner_dots;

  double initial_scale;

  int num;
  int len;
  int curr;

} SR1Block;

struct SleqpSR1
{
  int refcount;

  int num_variables;

  SleqpParams* params;
  SleqpOptions* options;

  SleqpSparseVec* grad_diff;
  SleqpSparseVec* step_diff;

  SleqpSparseVec* previous_grad;
  SleqpSparseVec* current_grad;

  SleqpSparseVec* block_grad_diff;
  SleqpSparseVec* block_step_diff;

  SleqpSparseVec* prod_cache;

  SleqpSparseVec* inner_cache;
  SleqpSparseVec* outer_cache;

  SleqpSparseVec* block_direction;
  SleqpSparseVec* block_prod;

  SR1Block* blocks;
  int num_blocks;

  SleqpFunc* sr1_func;
  SleqpFunc* func;

  SleqpTimer* update_timer;
};

static SLEQP_RETCODE
sr1_func_set_value(SleqpFunc* func,
                   SleqpSparseVec* x,
                   SLEQP_VALUE_REASON reason,
                   int* func_grad_nnz,
                   int* cons_val_nnz,
                   int* cons_jac_nnz,
                   void* func_data)
{
  SleqpSR1* sr1 = (SleqpSR1*) func_data;

  SLEQP_CALL(sleqp_func_set_value(sr1->func,
                                  x,
                                  reason,
                                  func_grad_nnz,
                                  cons_val_nnz,
                                  cons_jac_nnz));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
sr1_func_val(SleqpFunc* func,
             double* func_val,
             void* func_data)
{
  SleqpSR1* sr1 = (SleqpSR1*) func_data;

  SLEQP_CALL(sleqp_func_val(sr1->func,
                            func_val));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
sr1_func_grad(SleqpFunc* func,
              SleqpSparseVec* func_grad,
              void* func_data)
{
  SleqpSR1* sr1 = (SleqpSR1*) func_data;

  SLEQP_CALL(sleqp_func_grad(sr1->func,
                             func_grad));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
sr1_func_cons_val(SleqpFunc* func,
                  const SleqpSparseVec* cons_indices,
                  SleqpSparseVec* cons_val,
                  void* func_data)
{
  SleqpSR1* sr1 = (SleqpSR1*) func_data;

  SLEQP_CALL(sleqp_func_cons_val(sr1->func,
                                 cons_indices,
                                 cons_val));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
sr1_func_cons_jac(SleqpFunc* func,
                  const SleqpSparseVec* cons_indices,
                  SleqpSparseMatrix* cons_jac,
                  void* func_data)
{
  SleqpSR1* sr1 = (SleqpSR1*) func_data;

  SLEQP_CALL(sleqp_func_cons_jac(sr1->func,
                                 cons_indices,
                                 cons_jac));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
sr1_func_hess_prod(SleqpFunc* func,
                   const double* func_dual,
                   const SleqpSparseVec* direction,
                   const SleqpSparseVec* cons_duals,
                   SleqpSparseVec* product,
                   void* func_data)
{
  SleqpSR1* sr1 = (SleqpSR1*) func_data;

  SLEQP_CALL(sleqp_sr1_hess_prod(sr1,
                                 direction,
                                 product));

  return SLEQP_OKAY;
}


static SLEQP_RETCODE sr1_func_create(SleqpFunc** fstar,
                                     SleqpFunc* func,
                                     SleqpSR1* sr1)
{

  const int num_variables = sleqp_func_get_num_variables(func);
  const int num_constraints = sleqp_func_get_num_constraints(func);

  SleqpFuncCallbacks callbacks = {
    .set_value = sr1_func_set_value,
    .func_val = sr1_func_val,
    .func_grad = sr1_func_grad,
    .cons_val = sr1_func_cons_val,
    .cons_jac = sr1_func_cons_jac,
    .hess_prod = sr1_func_hess_prod,
    .func_free = NULL
  };

  SLEQP_CALL(sleqp_func_create(fstar,
                               &callbacks,
                               num_variables,
                               num_constraints,
                               sr1));

  SleqpFunc* sr1_func = *fstar;

  SLEQP_CALL(sleqp_hessian_struct_copy(sleqp_func_get_hess_struct(func),
                                       sleqp_func_get_hess_struct(sr1_func)));

  return SLEQP_OKAY;
}

static
SLEQP_RETCODE sr1_block_create_at(SR1Block* block,
                                  int dimension,
                                  int num)
{
  assert(dimension > 0);
  assert(num > 0);

  *block= (SR1Block) {0};

  block->num = num;
  block->dimension = dimension;

  block->len = 0;
  block->curr = -1;

  SLEQP_CALL(sleqp_alloc_array(&(block->step_diffs), num));

  SLEQP_CALL(sleqp_alloc_array(&(block->grad_diffs), num));

  SLEQP_CALL(sleqp_alloc_array(&(block->inner_prods), num));

  SLEQP_CALL(sleqp_alloc_array(&(block->inner_dots), num));

  for(int i = 0; i < num; ++i)
  {
    SLEQP_CALL(sleqp_sparse_vector_create_empty(block->step_diffs + i,
                                                dimension));

    SLEQP_CALL(sleqp_sparse_vector_create_empty(block->grad_diffs + i,
                                                dimension));

    SLEQP_CALL(sleqp_sparse_vector_create_empty(block->inner_prods + i,
                                                dimension));
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_sr1_create(SleqpSR1** star,
                               SleqpFunc* func,
                               SleqpParams* params,
                               SleqpOptions* options)
{
  SLEQP_CALL(sleqp_malloc(star));

  SleqpSR1* data = *star;

  *data = (SleqpSR1) {0};

  data->refcount = 1;

  SLEQP_CALL(sleqp_params_capture(params));
  data->params = params;

  SLEQP_CALL(sleqp_options_capture(options));
  data->options = options;

  const int num_iter = sleqp_options_get_int(options,
                                             SLEQP_OPTION_INT_NUM_QUASI_NEWTON_ITERATES);

  assert(num_iter > 0);

  const int num_variables = sleqp_func_get_num_variables(func);

  SleqpHessianStruct* hessian_struct = sleqp_func_get_hess_struct(func);

  const int num_blocks = sleqp_hessian_struct_get_num_blocks(hessian_struct);

  data->num_blocks = num_blocks;
  data->num_variables = num_variables;

  SLEQP_CALL(sleqp_alloc_array(&data->blocks, num_blocks));

  for(int block = 0; block < num_blocks; ++block)
  {
    int begin, end;

    SLEQP_CALL(sleqp_hessian_struct_get_block_range(hessian_struct,
                                                    block,
                                                    &begin,
                                                    &end));

    int block_dimension = end - begin;

    SLEQP_CALL(sr1_block_create_at(data->blocks + block,
                                   block_dimension,
                                   num_iter));
  }

  SLEQP_CALL(sleqp_sparse_vector_create_full(&(data->grad_diff),
                                             num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_full(&(data->step_diff),
                                             num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_full(&(data->previous_grad),
                                             num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_full(&(data->current_grad),
                                             num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_full(&(data->block_grad_diff),
                                             num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_full(&(data->block_step_diff),
                                             num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_full(&(data->prod_cache),
                                             num_variables));


  SLEQP_CALL(sleqp_sparse_vector_create_full(&(data->inner_cache),
                                             num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_full(&(data->outer_cache),
                                             num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_full(&(data->block_direction),
                                             num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_full(&(data->block_prod),
                                             num_variables));

  SLEQP_CALL(sr1_func_create(&(data->sr1_func),
                             func,
                             data));

  data->func = func;

  SLEQP_CALL(sleqp_timer_create(&(data->update_timer)));

  return SLEQP_OKAY;
}

static int data_index(SR1Block* block, int index)
{
  if(block->len == 0)
  {
    return 0;
  }

  int current_index = index % (block->num);

  return (current_index < 0) ? (current_index + block->num) : current_index;
}

static
SLEQP_RETCODE sr1_initial_scale(SleqpSparseVec* step_diff,
                                SleqpSparseVec* grad_diff,
                                double* initial_scale)
{
  double grad_step_dot;

  SLEQP_CALL(sleqp_sparse_vector_dot(grad_diff,
                                     step_diff,
                                     &grad_step_dot));

  const double grad_diff_normsq = sleqp_sparse_vector_norm_sq(grad_diff);

  assert(grad_diff_normsq >= 0.);

  if(grad_step_dot > 0)
  {
    (*initial_scale) = grad_diff_normsq / grad_step_dot;
  }
  else
  {
    (*initial_scale) = 1e-4;
  }

  return SLEQP_OKAY;
}

static
SLEQP_RETCODE sr1_compute_inner_products(SleqpSR1* sr1,
                                         SR1Block* block)
{
  assert(block->len > 0);

  const int first = block->curr - block->len + 1;

  const double zero_eps = sleqp_params_get(sr1->params, SLEQP_PARAM_ZERO_EPS);

  {
    int i = data_index(block, block->curr);

    SLEQP_CALL(sr1_initial_scale(block->step_diffs[i],
                                 block->grad_diffs[i],
                                 &block->initial_scale));
  }

  sr1->inner_cache->dim = block->dimension;
  sr1->outer_cache->dim = block->dimension;

  for(int val = first; val <= block->curr; ++val)
  {
    int i = data_index(block, val);

    SleqpSparseVec* current_step_diff = block->step_diffs[i];
    SleqpSparseVec* current_grad_diff = block->grad_diffs[i];

    SleqpSparseVec* current_inner_prod = block->inner_prods[i];

    {
      SLEQP_CALL(sleqp_sparse_vector_copy(block->step_diffs[i],
                                          sr1->outer_cache));

      SLEQP_CALL(sleqp_sparse_vector_scale(sr1->outer_cache,
                                           block->initial_scale));
    }

    for(int prev = first; prev < val; ++prev)
    {
      int j = data_index(block, prev);

      if(block->inner_dots[j] == 0.)
      {
        // this marks that the safeguard tripped
        continue;
      }

      SleqpSparseVec* inner_prod = block->inner_prods[j];

      double inner_step_dot;

      SLEQP_CALL(sleqp_sparse_vector_dot(current_step_diff,
                                         inner_prod,
                                         &inner_step_dot));

      double combined_factor = inner_step_dot / (block->inner_dots[j]);

      SLEQP_CALL(sleqp_sparse_vector_add_scaled(sr1->outer_cache,
                                                inner_prod,
                                                1.,
                                                combined_factor,
                                                zero_eps,
                                                sr1->inner_cache));

      SLEQP_CALL(sleqp_sparse_vector_copy(sr1->inner_cache, sr1->outer_cache));

    }

    SLEQP_CALL(sleqp_sparse_vector_add_scaled(sr1->outer_cache,
                                              current_grad_diff,
                                              -1.,
                                              1.,
                                              zero_eps,
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
        block->inner_dots[i] = 0.;
      }
      else
      {
        block->inner_dots[i] = inner_dot;
      }
    }

  }

  return SLEQP_OKAY;
}

static
SLEQP_RETCODE sr1_block_push(SleqpSR1* sr1,
                             SR1Block* block,
                             SleqpSparseVec* step_diff,
                             SleqpSparseVec* grad_diff)
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

  SLEQP_CALL(sr1_compute_inner_products(sr1, block));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_sr1_push(SleqpSR1* sr1,
                             SleqpIterate* previous_iterate,
                             SleqpIterate* current_iterate,
                             SleqpSparseVec* multipliers)
{
  const double eps = sleqp_params_get(sr1->params, SLEQP_PARAM_EPS);

  const double zero_eps = sleqp_params_get(sr1->params, SLEQP_PARAM_ZERO_EPS);

  const int num_blocks = sr1->num_blocks;

  SLEQP_CALL(sleqp_timer_start(sr1->update_timer));

  // Compute gradient difference
  {
    SLEQP_CALL(sleqp_sparse_matrix_trans_vector_product(sleqp_iterate_get_cons_jac(previous_iterate),
                                                        multipliers,
                                                        zero_eps,
                                                        sr1->prod_cache));

    SLEQP_CALL(sleqp_sparse_vector_add(sr1->prod_cache,
                                       sleqp_iterate_get_func_grad(previous_iterate),
                                       zero_eps,
                                       sr1->previous_grad));
  }

  {
    SLEQP_CALL(sleqp_sparse_matrix_trans_vector_product(sleqp_iterate_get_cons_jac(current_iterate),
                                                        multipliers,
                                                        zero_eps,
                                                        sr1->prod_cache));

    SLEQP_CALL(sleqp_sparse_vector_add(sr1->prod_cache,
                                       sleqp_iterate_get_func_grad(current_iterate),
                                       zero_eps,
                                       sr1->current_grad));
  }

  {
    SLEQP_CALL(sleqp_sparse_vector_add_scaled(sr1->previous_grad,
                                              sr1->current_grad,
                                              -1.,
                                              1.,
                                              zero_eps,
                                              sr1->grad_diff));
  }

  // Compute primal difference
  SLEQP_CALL(sleqp_sparse_vector_add_scaled(sleqp_iterate_get_primal(previous_iterate),
                                            sleqp_iterate_get_primal(current_iterate),
                                            -1.,
                                            1.,
                                            zero_eps,
                                            sr1->step_diff));

  int k_step = 0;
  int k_grad = 0;

  int offset = 0;

  for(int i = 0; i < num_blocks; ++i)
  {
    SR1Block* block = sr1->blocks + i;

    int next_offset = offset + block->dimension;

    SLEQP_CALL(sleqp_sparse_vector_clear(sr1->block_grad_diff));
    SLEQP_CALL(sleqp_sparse_vector_clear(sr1->block_step_diff));

    sr1->block_grad_diff->dim = block->dimension;
    sr1->block_step_diff->dim = block->dimension;

    while(k_step < sr1->step_diff->nnz &&
          sr1->step_diff->indices[k_step] < next_offset)
    {
      SLEQP_CALL(sleqp_sparse_vector_push(sr1->block_step_diff,
                                          sr1->step_diff->indices[k_step] - offset,
                                          sr1->step_diff->data[k_step]));

      ++k_step;
    }

    while(k_grad < sr1->grad_diff->nnz &&
          sr1->grad_diff->indices[k_grad] < next_offset)
    {
      SLEQP_CALL(sleqp_sparse_vector_push(sr1->block_grad_diff,
                                          sr1->grad_diff->indices[k_grad] - offset,
                                          sr1->grad_diff->data[k_grad]));

      ++k_grad;
    }

    const double step_normsq = sleqp_sparse_vector_norm_sq(sr1->block_step_diff);

    assert(sleqp_sparse_vector_is_finite(sr1->block_step_diff));
    assert(sleqp_sparse_vector_is_finite(sr1->block_grad_diff));

    if(!sleqp_is_zero(step_normsq, eps))
    {
      SLEQP_CALL(sr1_block_push(sr1,
                                block,
                                sr1->block_step_diff,
                                sr1->block_grad_diff));
    }


    offset = next_offset;
  }

  SLEQP_CALL(sleqp_timer_stop(sr1->update_timer));

  return SLEQP_OKAY;
}

static
SLEQP_RETCODE sr1_block_hess_prod(SleqpSR1* sr1,
                                  SR1Block* block,
                                  const SleqpSparseVec* direction,
                                  SleqpSparseVec* product)
{
  SLEQP_CALL(sleqp_sparse_vector_copy(direction, product));

  if(block->len == 0)
  {
    // If we have not acquired any information yet,
    // then we simply use the identity without any scaling.

    return SLEQP_OKAY;
  }

  const double zero_eps = sleqp_params_get(sr1->params, SLEQP_PARAM_ZERO_EPS);

  // Note: We have already computed the required inner products / dots

  SLEQP_CALL(sleqp_sparse_vector_scale(product, block->initial_scale));

  for(int val = block->curr - block->len + 1; val <= block->curr; ++val)
  {
    int i = data_index(block, val);

    if(block->inner_dots[i] == 0.)
    {
      // this marks that the safeguard tripped
      continue;
    }

    SleqpSparseVec* current_inner_prod = block->inner_prods[i];

    double direction_dot;

    SLEQP_CALL(sleqp_sparse_vector_dot(current_inner_prod, direction, &direction_dot));

    double combined_factor = direction_dot / (block->inner_dots[i]);

    SLEQP_CALL(sleqp_sparse_vector_add_scaled(product,
                                              current_inner_prod,
                                              1.,
                                              combined_factor,
                                              zero_eps,
                                              sr1->inner_cache));

    SLEQP_CALL(sleqp_sparse_vector_copy(sr1->inner_cache, product));
  }

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_sr1_hess_prod(SleqpSR1* sr1,
                                  const SleqpSparseVec* direction,
                                  SleqpSparseVec* product)
{
  const int num_blocks = sr1->num_blocks;

  int offset = 0;
  int k_direction = 0, k_product = 0;

  SLEQP_CALL(sleqp_sparse_vector_clear(product));
  SLEQP_CALL(sleqp_sparse_vector_reserve(product, sr1->num_variables));

  for(int i = 0; i < num_blocks; ++i)
  {
    SR1Block* block = sr1->blocks + i;

    int next_offset = offset + block->dimension;

    if(block->len == 0)
    {
      while(k_direction < direction->nnz &&
            direction->indices[k_direction] < next_offset)
      {
        assert(direction->indices[k_direction] >= offset);

        SLEQP_CALL(sleqp_sparse_vector_push(product,
                                            direction->indices[k_direction],
                                            direction->data[k_direction]));

        ++k_direction;
        ++k_product;
      }
    }
    else
    {
      sr1->block_direction->dim = block->dimension;
      sr1->block_prod->dim = block->dimension;

      SLEQP_CALL(sleqp_sparse_vector_clear(sr1->block_direction));
      SLEQP_CALL(sleqp_sparse_vector_clear(sr1->block_prod));

      while(k_direction < direction->nnz &&
            direction->indices[k_direction] < next_offset)
      {
        SLEQP_CALL(sleqp_sparse_vector_push(sr1->block_direction,
                                            direction->indices[k_direction] - offset,
                                            direction->data[k_direction]));

        ++k_direction;
      }

      assert(sleqp_sparse_vector_is_finite(sr1->block_direction));

      SLEQP_CALL(sr1_block_hess_prod(sr1,
                                     block,
                                     sr1->block_direction,
                                     sr1->block_prod));

      assert(sleqp_sparse_vector_is_finite(sr1->block_prod));

      for(int k_block_prod = 0;
          k_block_prod < sr1->block_prod->nnz;
          ++k_block_prod)
      {
        assert(sr1->block_prod->indices[k_block_prod] + offset < next_offset);

        SLEQP_CALL(sleqp_sparse_vector_push(product,
                                            sr1->block_prod->indices[k_block_prod] + offset,
                                            sr1->block_prod->data[k_block_prod]));

        ++k_product;
      }
    }

    offset = next_offset;
  }

  return SLEQP_OKAY;
}

SleqpFunc* sleqp_sr1_get_func(SleqpSR1* sr1)
{
  return sr1->sr1_func;
}

static SLEQP_RETCODE sr1_block_free_at(SR1Block* block)
{
  const int num = block->num;

  for(int i = 0; i < num; ++i)
  {
    SLEQP_CALL(sleqp_sparse_vector_free(block->inner_prods + i));

    SLEQP_CALL(sleqp_sparse_vector_free(block->grad_diffs + i));
    SLEQP_CALL(sleqp_sparse_vector_free(block->step_diffs + i));
  }

  sleqp_free(&(block->inner_dots));

  sleqp_free(&block->inner_prods);

  sleqp_free(&block->grad_diffs);
  sleqp_free(&block->step_diffs);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE sr1_free(SleqpSR1** star)
{
  SleqpSR1* sr1 = *star;

  if(!sr1)
  {
    return SLEQP_OKAY;
  }

  SLEQP_CALL(sleqp_timer_free(&(sr1->update_timer)));

  SLEQP_CALL(sleqp_func_release(&(sr1->sr1_func)));

  int num_blocks = sr1->num_blocks;

  for(int block = 0; block < num_blocks; ++block)
  {
    SLEQP_CALL(sr1_block_free_at(sr1->blocks + block));
  }

  sleqp_free(&sr1->blocks);

  SLEQP_CALL(sleqp_sparse_vector_free(&(sr1->block_prod)));
  SLEQP_CALL(sleqp_sparse_vector_free(&(sr1->block_direction)));

  SLEQP_CALL(sleqp_sparse_vector_free(&(sr1->outer_cache)));
  SLEQP_CALL(sleqp_sparse_vector_free(&(sr1->inner_cache)));

  SLEQP_CALL(sleqp_sparse_vector_free(&(sr1->prod_cache)));

  SLEQP_CALL(sleqp_sparse_vector_free(&(sr1->block_step_diff)));
  SLEQP_CALL(sleqp_sparse_vector_free(&(sr1->block_grad_diff)));

  SLEQP_CALL(sleqp_sparse_vector_free(&(sr1->current_grad)));
  SLEQP_CALL(sleqp_sparse_vector_free(&(sr1->previous_grad)));

  SLEQP_CALL(sleqp_sparse_vector_free(&(sr1->step_diff)));
  SLEQP_CALL(sleqp_sparse_vector_free(&(sr1->grad_diff)));

  SLEQP_CALL(sleqp_options_release(&sr1->options));
  SLEQP_CALL(sleqp_params_release(&sr1->params));

  sleqp_free(star);

  return SLEQP_OKAY;
}

SleqpTimer* sleqp_sr1_update_timer(SleqpSR1* data)
{
  return data->update_timer;
}

SLEQP_RETCODE sleqp_sr1_capture(SleqpSR1* data)
{
  ++data->refcount;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_sr1_release(SleqpSR1** star)
{
  SleqpSR1* data = *star;

  if(!data)
  {
    return SLEQP_OKAY;
  }

  if(--data->refcount == 0)
  {
    SLEQP_CALL(sr1_free(star));
  }

  *star = NULL;

  return SLEQP_OKAY;
}
