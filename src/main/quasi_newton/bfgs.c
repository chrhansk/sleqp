#include "bfgs.h"

#include <assert.h>
#include <math.h>

#include "fail.h"
#include "cmp.h"
#include "mem.h"
#include "sparse/sparse_matrix.h"

#include "quasi_newton.h"

static const double damping_factor = 0.2;
static const double sizing_cutoff = 0.1;

static const double initial_scale_min = 1e-6;
static const double damped_initial_scale_max = 1.;

/**
 * We follow the notation in "Numerical Optimization"
 **/

typedef struct
{
  int dimension;

  // The vectors \f$ s_k \f$
  SleqpSparseVec** point_diffs;
  // The vectors \f$ y_k \f$
  SleqpSparseVec** grad_diffs;

  // The products \f$ B_k s_k \f$
  SleqpSparseVec** point_products;

  // The vectors \f$ r_k \f$
  SleqpSparseVec** damped_grad_diffs;

  // The products \f$ s_k^{T} r_k \f$
  double* damped_grad_point_diff_dots;

  // The products \f$ s_k^{T} y_k \f$
  double* grad_point_diff_dots;

  // The products \f$ s_k^{T} s_k \f$
  double* point_diff_inner_dots;

  // The products \f$ s_k^{T} B_k s_k \f$
  double* bidir_products;

  SleqpSparseVec* inner_cache;

  SleqpSparseVec* prod_cache;

  double initial_scale;

  double* sizing_factors;

  // max size
  int num;
  // curr size
  int len;
  // curr index
  int curr;

  bool damped;

  SLEQP_BFGS_SIZING sizing;

} BFGSBlock;

typedef struct
{
  int refcount;

  int num_variables;

  SleqpParams* params;
  SleqpOptions* options;

  int num_blocks;
  BFGSBlock* blocks;

  SleqpSparseVec* grad_diff;
  SleqpSparseVec* point_diff;

  SleqpSparseVec* previous_grad;
  SleqpSparseVec* current_grad;

  SleqpSparseVec* block_grad_diff;
  SleqpSparseVec* block_point_diff;

  SleqpSparseVec* prod_cache;

  SleqpSparseVec* block_direction;
  SleqpSparseVec* block_prod;
} BFGS;

static SLEQP_RETCODE bfgs_block_create_at(BFGSBlock* block,
                                          int dimension,
                                          const SleqpParams* params,
                                          int num,
                                          bool damped,
                                          SLEQP_BFGS_SIZING sizing)
{
  assert(dimension > 0);
  assert(num > 0);

  *block = (BFGSBlock) {0};

  block->dimension = dimension;

  block->num = num;
  block->len = 0;
  block->curr = -1;
  block->damped = damped;
  block->sizing = sizing;

  SLEQP_CALL(sleqp_alloc_array(&(block->point_diffs), num));

  SLEQP_CALL(sleqp_alloc_array(&(block->grad_diffs), num));

  SLEQP_CALL(sleqp_alloc_array(&(block->point_products), num));

  SLEQP_CALL(sleqp_alloc_array(&(block->damped_grad_diffs), num));

  SLEQP_CALL(sleqp_alloc_array(&(block->damped_grad_point_diff_dots), num));

  if(sizing != SLEQP_BFGS_SIZING_NONE)
  {
    SLEQP_CALL(sleqp_alloc_array(&(block->grad_point_diff_dots), num));

    SLEQP_CALL(sleqp_alloc_array(&(block->point_diff_inner_dots), num));
  }

  SLEQP_CALL(sleqp_alloc_array(&(block->bidir_products), num));

  SLEQP_CALL(sleqp_alloc_array(&(block->sizing_factors), num));

  for(int i = 0; i < num; ++i)
  {
    SLEQP_CALL(sleqp_sparse_vector_create_empty(block->point_diffs + i,
                                                dimension));

    SLEQP_CALL(sleqp_sparse_vector_create_empty(block->grad_diffs + i,
                                                dimension));

    SLEQP_CALL(sleqp_sparse_vector_create_empty(block->point_products + i,
                                                dimension));

    SLEQP_CALL(sleqp_sparse_vector_create_empty(block->damped_grad_diffs + i,
                                                dimension));

    block->damped_grad_point_diff_dots[i] = 0.;
    block->sizing_factors[i] = 1.;
  }

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&(block->inner_cache),
                                              dimension));

  SLEQP_CALL(sleqp_sparse_vector_create_empty(&(block->prod_cache),
                                              dimension));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE bfgs_block_free_at(BFGSBlock* block)
{
  const int num = block->num;

  SLEQP_CALL(sleqp_sparse_vector_free(&(block->prod_cache)));

  SLEQP_CALL(sleqp_sparse_vector_free(&(block->inner_cache)));

  for(int i = 0; i < num; ++i)
  {
    SLEQP_CALL(sleqp_sparse_vector_free(block->damped_grad_diffs + i));
    SLEQP_CALL(sleqp_sparse_vector_free(block->point_products + i));

    SLEQP_CALL(sleqp_sparse_vector_free(block->grad_diffs + i));
    SLEQP_CALL(sleqp_sparse_vector_free(block->point_diffs + i));
  }

  sleqp_free(&block->damped_grad_diffs);
  sleqp_free(&block->point_products);

  sleqp_free(&block->point_diff_inner_dots);

  sleqp_free(&block->grad_point_diff_dots);

  sleqp_free(&block->damped_grad_point_diff_dots);

  sleqp_free(&block->sizing_factors);

  sleqp_free(&block->bidir_products);

  sleqp_free(&block->grad_diffs);
  sleqp_free(&block->point_diffs);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
bfgs_create(BFGS** star,
            SleqpFunc* func,
            SleqpParams* params,
            SleqpOptions* options)
{
  SLEQP_CALL(sleqp_malloc(star));

  BFGS* data = *star;

  *data = (BFGS) {0};

  data->refcount = 1;

  SLEQP_CALL(sleqp_params_capture(params));
  data->params = params;

  SLEQP_CALL(sleqp_options_capture(options));
  data->options = options;

  const SLEQP_HESSIAN_EVAL hessian_eval = sleqp_options_get_int(options, SLEQP_OPTION_INT_HESSIAN_EVAL);

  assert(hessian_eval == SLEQP_HESSIAN_EVAL_SIMPLE_BFGS ||
         hessian_eval == SLEQP_HESSIAN_EVAL_DAMPED_BFGS);

  const bool damped = (hessian_eval == SLEQP_HESSIAN_EVAL_DAMPED_BFGS);

  const int num = sleqp_options_get_int(options, SLEQP_OPTION_INT_NUM_QUASI_NEWTON_ITERATES);

  const SLEQP_BFGS_SIZING sizing = sleqp_options_get_int(options, SLEQP_OPTION_INT_BFGS_SIZING);

  assert(num > 0);

  const int num_variables = sleqp_func_get_num_variables(func);

  SleqpHessianStruct* hessian_struct = sleqp_func_get_hess_struct(func);

  const int num_blocks = sleqp_hessian_struct_get_num_blocks(hessian_struct);

  data->num_variables = num_variables;
  data->num_blocks = num_blocks;

  SLEQP_CALL(sleqp_alloc_array(&data->blocks, num_blocks));

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
                                    damped,
                                    sizing));
  }

  SLEQP_CALL(sleqp_sparse_vector_create_full(&data->grad_diff,
                                             num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_full(&data->point_diff,
                                             num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_full(&data->previous_grad,
                                             num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_full(&data->current_grad,
                                             num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_full(&data->block_grad_diff,
                                             num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_full(&data->block_point_diff,
                                             num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_full(&data->prod_cache,
                                             num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_full(&data->block_direction,
                                             num_variables));

  SLEQP_CALL(sleqp_sparse_vector_create_full(&data->block_prod,
                                             num_variables));

  return SLEQP_OKAY;
}

static int
data_index(BFGSBlock* block, int index)
{
  if(block->len == 0)
  {
    return 0;
  }

  int current_index = index % (block->num);

  return (current_index < 0) ? (current_index + block->num) : current_index;
}

static SLEQP_RETCODE
bfgs_hess_prod_range(BFGSBlock* block,
                     const SleqpParams* params,
                     const SleqpSparseVec* direction,
                     SleqpSparseVec* product,
                     int final)
{
  const double zero_eps = sleqp_params_get(params, SLEQP_PARAM_ZERO_EPS);
  const int begin = block->curr - block->len + 1;

  // Initially apply scaled identity
  {
    SLEQP_CALL(sleqp_sparse_vector_copy(direction, product));

    SLEQP_CALL(sleqp_sparse_vector_scale(product, block->initial_scale));
  }

  // invariant: after iteration k, "product" contains
  // product of given "direction" with the approximate
  // Hessian consisting of the first k terms
  for(int prev = begin; prev <= final; ++prev)
  {
    int j = data_index(block, prev);

    SleqpSparseVec* current_point_product = block->point_products[j];
    SleqpSparseVec* current_damped_grad_diff = block->damped_grad_diffs[j];
    const double sizing_factor = block->sizing_factors[j];
    const double bidir_product = block->bidir_products[j];

    double direction_product_dot, direction_damped_grad_dot;

    SLEQP_CALL(sleqp_sparse_vector_dot(current_point_product,
                                       direction,
                                       &direction_product_dot));

    SLEQP_CALL(sleqp_sparse_vector_dot(current_damped_grad_diff,
                                       direction,
                                       &direction_damped_grad_dot));

    SLEQP_CALL(sleqp_sparse_vector_add_scaled(product,
                                              current_point_product,
                                              1.,
                                              -1. * direction_product_dot / bidir_product,
                                              zero_eps,
                                              block->inner_cache));

    SLEQP_CALL(sleqp_sparse_vector_add_scaled(block->inner_cache,
                                              current_damped_grad_diff,
                                              sizing_factor,
                                              direction_damped_grad_dot / block->damped_grad_point_diff_dots[j],
                                              zero_eps,
                                              product));
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
bfgs_initial_scale(BFGSBlock* block,
                   SleqpSparseVec* point_diff,
                   SleqpSparseVec* grad_diff,
                   double* initial_scale)
{
  double grad_point_dot;

  SLEQP_CALL(sleqp_sparse_vector_dot(grad_diff,
                                     point_diff,
                                     &grad_point_dot));

  const double point_diff_normsq = sleqp_sparse_vector_norm_sq(point_diff);

  assert(point_diff_normsq >= 0.);

  if(grad_point_dot == 0.)
  {
    (*initial_scale) = 1.;
    return SLEQP_OKAY;
  }

  (*initial_scale) = point_diff_normsq / grad_point_dot;

  (*initial_scale) = SLEQP_MAX((*initial_scale), initial_scale_min);

  if(block->damped)
  {
    // TODO: Find out if there is smoe better way
    // of applying the damping to the initial approximation
    (*initial_scale) = SLEQP_MIN((*initial_scale), damped_initial_scale_max);
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
bfgs_compute_sizing(BFGSBlock* block, const int val)
{
  const int begin = block->curr - block->len + 1;

  const int j = data_index(block, val);

  // First one...
  if(begin == val)
  {
    block->sizing_factors[j] = 1.;

    return SLEQP_OKAY;
  }

  block->sizing_factors[j] = 1.;

  assert(begin < val);

  const int i = data_index(block, val - 1);

  assert(i != j);

  if(block->sizing == SLEQP_BFGS_SIZING_CENTERED_OL)
  {
    const double factor = .5;

    const double numerator = (1. - factor) * (block->grad_point_diff_dots[i]) / (block->point_diff_inner_dots[i]) +
      (factor) * (block->grad_point_diff_dots[j]) / (block->point_diff_inner_dots[j]);

    const double denominator = (1. - factor) * (block->damped_grad_point_diff_dots[i]) / (block->point_diff_inner_dots[i]) +
      (factor) * (block->bidir_products[j]);

    block->sizing_factors[j] = numerator / denominator;

    block->sizing_factors[j] = SLEQP_MAX(block->sizing_factors[j], sizing_cutoff);

    block->sizing_factors[j] = SLEQP_MIN(block->sizing_factors[j], 1.);
  }

  assert(block->sizing_factors[j] >= 0.);

  return SLEQP_OKAY;
}

static
SLEQP_RETCODE
bfgs_compute_products(BFGSBlock* block,
                      const SleqpParams* params)
{
  const double eps = sleqp_params_get(params, SLEQP_PARAM_EPS);

  SLEQP_NUM_ASSERT_PARAM(eps);

  const double zero_eps = sleqp_params_get(params, SLEQP_PARAM_ZERO_EPS);

  assert(block->len > 0);

  {
    SLEQP_CALL(bfgs_initial_scale(block,
                                  block->point_diffs[block->curr],
                                  block->grad_diffs[block->curr],
                                  &block->initial_scale));

    sleqp_assert_is_pos(block->initial_scale, eps);
  }

  const int begin = block->curr - block->len + 1;

  for(int val = begin; val <= block->curr; ++val)
  {
    const int i = data_index(block, val);

    const SleqpSparseVec* current_point_diff = block->point_diffs[i];
    const SleqpSparseVec* current_grad_diff = block->grad_diffs[i];

    const SleqpSparseVec* direction = current_point_diff;
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

    SLEQP_CALL(sleqp_sparse_vector_dot(current_grad_diff,
                                       current_point_diff,
                                       &dot_product));

    SleqpSparseVec* current_point_product = block->point_products[i];
    SleqpSparseVec* current_damped_grad_diff = block->damped_grad_diffs[i];

    bool is_damped = false;

    if(block->damped && (dot_product < (damping_factor * bidir_product)))
    {
      double combination_factor = (1. - damping_factor) *
        bidir_product / (bidir_product - dot_product);

      sleqp_assert_is_pos(combination_factor, eps);
      sleqp_assert_is_lt(combination_factor, 1., eps);

      SLEQP_CALL(sleqp_sparse_vector_add_scaled(current_grad_diff,
                                                product,
                                                combination_factor,
                                                1. - combination_factor,
                                                zero_eps,
                                                current_damped_grad_diff));

      SLEQP_CALL(sleqp_sparse_vector_dot(current_damped_grad_diff,
                                         current_point_diff,
                                         &dot_product));

      is_damped = true;
    }
    else
    {
      SLEQP_CALL(sleqp_sparse_vector_copy(current_grad_diff, current_damped_grad_diff));
    }

    // save dot product
    {
      assert(dot_product > 0);

      block->damped_grad_point_diff_dots[i] = dot_product;
    }

    // set point product
    {
      SLEQP_CALL(sleqp_sparse_vector_copy(product, current_point_product));

      //SLEQP_CALL(sleqp_sparse_vector_scale(current_point_product, 1./bidir_product));

      block->bidir_products[i] = bidir_product;
    }

    if(block->sizing != SLEQP_BFGS_SIZING_NONE)
    {
      SLEQP_CALL(bfgs_compute_sizing(block, val));
    }

#if !defined(NDEBUG)
    if(!is_damped)
    {
      // check that secant equation is satisfied
      const int i = data_index(block, val);

      SLEQP_CALL(bfgs_hess_prod_range(block,
                                      params,
                                      block->point_diffs[i],
                                      block->prod_cache,
                                      val));

      sleqp_num_assert(sleqp_sparse_vector_eq(block->prod_cache,
                                              block->grad_diffs[i],
                                              eps));
    }
#endif
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
bfgs_block_push(BFGSBlock* block,
                const SleqpParams* params,
                const SleqpSparseVec* point_diff,
                const SleqpSparseVec* grad_diff)
{
  const int next = data_index(block, block->curr + 1);

  SLEQP_CALL(sleqp_sparse_vector_copy(point_diff,
                                      block->point_diffs[next]));

  SLEQP_CALL(sleqp_sparse_vector_copy(grad_diff,
                                      block->grad_diffs[next]));

  if(block->sizing != SLEQP_BFGS_SIZING_NONE)
  {
    block->point_diff_inner_dots[next] = sleqp_sparse_vector_norm_sq(point_diff);

    SLEQP_CALL(sleqp_sparse_vector_dot(point_diff,
                                       grad_diff,
                                       block->grad_point_diff_dots + next));
  }

  if(block->len < block->num)
  {
    ++block->len;
  }

  block->curr = next;

  SLEQP_CALL(bfgs_compute_products(block, params));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
bfgs_push(const SleqpIterate* previous_iterate,
          const SleqpIterate* current_iterate,
          const SleqpSparseVec* multipliers,
          void* data)
{
  BFGS* bfgs = (BFGS*) data;

  const double eps = sleqp_params_get(bfgs->params, SLEQP_PARAM_EPS);
  const double zero_eps = sleqp_params_get(bfgs->params, SLEQP_PARAM_ZERO_EPS);

  const int num_blocks = bfgs->num_blocks;

  // Compute gradient difference
  {
    SLEQP_CALL(sleqp_sparse_matrix_trans_vector_product(sleqp_iterate_get_cons_jac(previous_iterate),
                                                        multipliers,
                                                        zero_eps,
                                                        bfgs->prod_cache));

    SLEQP_CALL(sleqp_sparse_vector_add(bfgs->prod_cache,
                                       sleqp_iterate_get_func_grad(previous_iterate),
                                       zero_eps,
                                       bfgs->previous_grad));
  }

  {
    SLEQP_CALL(sleqp_sparse_matrix_trans_vector_product(sleqp_iterate_get_cons_jac(current_iterate),
                                                        multipliers,
                                                        zero_eps,
                                                        bfgs->prod_cache));

    SLEQP_CALL(sleqp_sparse_vector_add(bfgs->prod_cache,
                                       sleqp_iterate_get_func_grad(current_iterate),
                                       zero_eps,
                                       bfgs->current_grad));
  }

  {
    SLEQP_CALL(sleqp_sparse_vector_add_scaled(bfgs->previous_grad,
                                              bfgs->current_grad,
                                              -1.,
                                              1.,
                                              zero_eps,
                                              bfgs->grad_diff));
  }

  // Compute primal difference
  SLEQP_CALL(sleqp_sparse_vector_add_scaled(sleqp_iterate_get_primal(previous_iterate),
                                            sleqp_iterate_get_primal(current_iterate),
                                            -1.,
                                            1.,
                                            zero_eps,
                                            bfgs->point_diff));

  int k_point = 0;
  int k_grad = 0;

  int offset = 0;

  for(int i = 0; i < num_blocks; ++i)
  {
    BFGSBlock* block = bfgs->blocks + i;

    int next_offset = offset + block->dimension;

    SLEQP_CALL(sleqp_sparse_vector_clear(bfgs->block_grad_diff));
    SLEQP_CALL(sleqp_sparse_vector_clear(bfgs->block_point_diff));

    bfgs->block_grad_diff->dim = block->dimension;
    bfgs->block_point_diff->dim = block->dimension;

    while(k_point < bfgs->point_diff->nnz &&
          bfgs->point_diff->indices[k_point] < next_offset)
    {
      SLEQP_CALL(sleqp_sparse_vector_push(bfgs->block_point_diff,
                                          bfgs->point_diff->indices[k_point] - offset,
                                          bfgs->point_diff->data[k_point]));

      ++k_point;
    }

    while(k_grad < bfgs->grad_diff->nnz &&
          bfgs->grad_diff->indices[k_grad] < next_offset)
    {
      SLEQP_CALL(sleqp_sparse_vector_push(bfgs->block_grad_diff,
                                          bfgs->grad_diff->indices[k_grad] - offset,
                                          bfgs->grad_diff->data[k_grad]));

      ++k_grad;
    }

    const double point_normsq = sleqp_sparse_vector_norm_sq(bfgs->block_point_diff);

    assert(sleqp_sparse_vector_is_finite(bfgs->block_point_diff));
    assert(sleqp_sparse_vector_is_finite(bfgs->block_grad_diff));

    if(!sleqp_is_zero(point_normsq, eps))
    {
      SLEQP_CALL(bfgs_block_push(block,
                                 bfgs->params,
                                 bfgs->block_point_diff,
                                 bfgs->block_grad_diff));
    }


    offset = next_offset;
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
bfgs_reset(void* data)
{
  BFGS* bfgs = (BFGS*) data;

  const int num_blocks = bfgs->num_blocks;

  for (int i = 0; i < num_blocks; ++i)
  {
    BFGSBlock *block = bfgs->blocks + i;
    block->len = 0;
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
bfgs_hess_prod(const SleqpSparseVec* direction,
               SleqpSparseVec* product,
               void* data)
{
  BFGS* bfgs = (BFGS*) data;

  const int num_blocks = bfgs->num_blocks;

  int offset = 0;
  int k_direction = 0, k_product = 0;

  SLEQP_CALL(sleqp_sparse_vector_reserve(product, bfgs->num_variables));

  for(int i = 0; i < num_blocks; ++i)
  {
    BFGSBlock* block = bfgs->blocks + i;

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
      bfgs->block_direction->dim = block->dimension;
      bfgs->block_prod->dim = block->dimension;

      SLEQP_CALL(sleqp_sparse_vector_clear(bfgs->block_direction));
      SLEQP_CALL(sleqp_sparse_vector_clear(bfgs->block_prod));

      while(k_direction < direction->nnz &&
            direction->indices[k_direction] < next_offset)
      {
        SLEQP_CALL(sleqp_sparse_vector_push(bfgs->block_direction,
                                            direction->indices[k_direction] - offset,
                                            direction->data[k_direction]));

        ++k_direction;
      }

      SLEQP_CALL(bfgs_hess_prod_range(block,
                                      bfgs->params,
                                      bfgs->block_direction,
                                      bfgs->block_prod,
                                      block->curr));

      for(int k_block_prod = 0;
          k_block_prod < bfgs->block_prod->nnz;
          ++k_block_prod)
      {
        SLEQP_CALL(sleqp_sparse_vector_push(product,
                                            bfgs->block_prod->indices[k_block_prod] + offset,
                                            bfgs->block_prod->data[k_block_prod]));

        ++k_product;
      }
    }

    offset = next_offset;
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE bfgs_free(void* data)
{
  BFGS* bfgs = (BFGS*) data;

  if(!data)
  {
    return SLEQP_OKAY;
  }

  SLEQP_CALL(sleqp_sparse_vector_free(&bfgs->block_prod));
  SLEQP_CALL(sleqp_sparse_vector_free(&bfgs->block_direction));

  SLEQP_CALL(sleqp_sparse_vector_free(&bfgs->prod_cache));

  SLEQP_CALL(sleqp_sparse_vector_free(&bfgs->block_prod));
  SLEQP_CALL(sleqp_sparse_vector_free(&bfgs->block_prod));

  SLEQP_CALL(sleqp_sparse_vector_free(&bfgs->block_point_diff));
  SLEQP_CALL(sleqp_sparse_vector_free(&bfgs->block_grad_diff));

  SLEQP_CALL(sleqp_sparse_vector_free(&bfgs->current_grad));
  SLEQP_CALL(sleqp_sparse_vector_free(&bfgs->previous_grad));

  SLEQP_CALL(sleqp_sparse_vector_free(&bfgs->point_diff));
  SLEQP_CALL(sleqp_sparse_vector_free(&bfgs->grad_diff));

  int num_blocks = bfgs->num_blocks;

  for(int block = 0; block < num_blocks; ++block)
  {
    SLEQP_CALL(bfgs_block_free_at(bfgs->blocks + block));
  }

  sleqp_free(&bfgs->blocks);

  SLEQP_CALL(sleqp_options_release(&bfgs->options));

  SLEQP_CALL(sleqp_params_release(&bfgs->params));

  sleqp_free(&bfgs);

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_bfgs_create(SleqpQuasiNewton** star,
                                SleqpFunc* func,
                                SleqpParams* params,
                                SleqpOptions* options)
{
  SleqpQuasiNewtonCallbacks callbacks = {
    .push      = bfgs_push,
    .reset     = bfgs_reset,
    .hess_prod = bfgs_hess_prod,
    .free      = bfgs_free,
  };

  BFGS* bfgs;

  SLEQP_CALL(bfgs_create(&bfgs,
                         func,
                         params,
                         options));

  SLEQP_CALL(sleqp_quasi_newton_create(star,
                                       func,
                                       &callbacks,
                                       (void*) bfgs));

  SleqpQuasiNewton* quasi_newton = *star;

  SleqpFunc* bfgs_func = sleqp_quasi_newton_get_func(quasi_newton);

  SLEQP_CALL(sleqp_hessian_struct_copy(sleqp_func_get_hess_struct(func),
                                       sleqp_func_get_hess_struct(bfgs_func)));

  SLEQP_CALL(sleqp_func_set_psd_hessian(bfgs_func, true));

  return SLEQP_OKAY;
}
