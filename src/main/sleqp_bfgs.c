#include "sleqp_bfgs.h"

#include <assert.h>
#include <math.h>

#include "sleqp_assert.h"
#include "sleqp_cmp.h"
#include "sleqp_mem.h"

static const double damping_factor = 0.2;
static const double sizing_cutoff = 0.1;

/**
 * We follow the notation in "Numerical Optimization"
 **/

typedef struct BFGSBlock
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

struct SleqpBFGSData
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

  SleqpFunc* bfgs_func;
  SleqpFunc* func;

  SleqpTimer* update_timer;
};

static SLEQP_RETCODE
bfgs_func_set_value(SleqpFunc* func,
                    SleqpSparseVec* x,
                    SLEQP_VALUE_REASON reason,
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
bfgs_func_val(SleqpFunc* func,
              double* func_val,
              void* func_data)
{
  SleqpBFGSData* bfgs_data = (SleqpBFGSData*) func_data;

  SLEQP_CALL(sleqp_func_val(bfgs_data->func,
                            func_val));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
bfgs_func_grad(SleqpFunc* func,
               SleqpSparseVec* func_grad,
               void* func_data)
{
  SleqpBFGSData* bfgs_data = (SleqpBFGSData*) func_data;

  SLEQP_CALL(sleqp_func_grad(bfgs_data->func,
                             func_grad));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
bfgs_func_cons_val(SleqpFunc* func,
                   const SleqpSparseVec* cons_indices,
                   SleqpSparseVec* cons_val,
                   void* func_data)
{
  SleqpBFGSData* bfgs_data = (SleqpBFGSData*) func_data;

  SLEQP_CALL(sleqp_func_cons_val(bfgs_data->func,
                                 cons_indices,
                                 cons_val));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
bfgs_func_cons_jac(SleqpFunc* func,
                   const SleqpSparseVec* cons_indices,
                   SleqpSparseMatrix* cons_jac,
                   void* func_data)
{
  SleqpBFGSData* bfgs_data = (SleqpBFGSData*) func_data;

  SLEQP_CALL(sleqp_func_cons_jac(bfgs_data->func,
                                 cons_indices,
                                 cons_jac));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
bfgs_func_hess_prod(SleqpFunc* func,
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

static SLEQP_RETCODE bfgs_func_create(SleqpFunc** fstar,
                                      SleqpFunc* func,
                                      SleqpBFGSData* bfgs_data)
{
  const int num_variables = sleqp_func_get_num_variables(func);
  const int num_constraints = sleqp_func_get_num_constraints(func);

  SleqpFuncCallbacks callbacks = {
    .set_value = bfgs_func_set_value,
    .func_val = bfgs_func_val,
    .func_grad = bfgs_func_grad,
    .cons_val = bfgs_func_cons_val,
    .cons_jac = bfgs_func_cons_jac,
    .hess_prod = bfgs_func_hess_prod,
    .func_free = NULL
  };

  SLEQP_CALL(sleqp_func_create(fstar,
                               &callbacks,
                               num_variables,
                               num_constraints,
                               bfgs_data));

  SleqpFunc* bfgs_func = *fstar;

  SLEQP_CALL(sleqp_hessian_struct_copy(sleqp_func_get_hess_struct(func),
                                       sleqp_func_get_hess_struct(bfgs_func)));

  SLEQP_CALL(sleqp_func_set_psd_hessian(func, true));

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

SLEQP_RETCODE sleqp_bfgs_data_create(SleqpBFGSData** star,
                                     SleqpFunc* func,
                                     SleqpParams* params,
                                     SleqpOptions* options)
{
  SLEQP_CALL(sleqp_malloc(star));

  SleqpBFGSData* data = *star;

  *data = (SleqpBFGSData) {0};

  data->refcount = 1;

  SLEQP_CALL(sleqp_params_capture(params));
  data->params = params;

  SLEQP_CALL(sleqp_options_capture(options));
  data->options = options;

  const SLEQP_HESSIAN_EVAL hessian_eval = sleqp_options_get_hessian_eval(options);

  assert(hessian_eval == SLEQP_HESSIAN_EVAL_SIMPLE_BFGS ||
         hessian_eval == SLEQP_HESSIAN_EVAL_DAMPED_BFGS);

  const bool damped = (hessian_eval == SLEQP_HESSIAN_EVAL_DAMPED_BFGS);
  const int num = sleqp_options_get_quasi_newton_num_iterates(options);
  const SLEQP_BFGS_SIZING sizing = sleqp_options_get_bfgs_sizing(options);

  assert(num > 0);

  const int num_variables = sleqp_func_get_num_variables(func);

  SleqpHessianStruct* hessian_struct = sleqp_func_get_hess_struct(func);

  const int num_blocks = sleqp_hessian_struct_get_num_blocks(hessian_struct);

  data->num_variables = num_variables;
  data->num_blocks = num_blocks;
  data->func = func;

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

  SLEQP_CALL(bfgs_func_create(&data->bfgs_func,
                              data->func,
                              data));

  SLEQP_CALL(sleqp_timer_create(&(data->update_timer)));

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

static
SLEQP_RETCODE bfgs_initial_scale(BFGSBlock* block,
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

  if(block->damped)
  {
    // TODO: Find out if there is smoe better way
    // of applying the damping to the initial approximation
    (*initial_scale) = SLEQP_MAX((*initial_scale), 1.);
  }

  return SLEQP_OKAY;
}

static
SLEQP_RETCODE bfgs_compute_sizing(BFGSBlock* block, const int val)
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
SLEQP_RETCODE bfgs_compute_products(BFGSBlock* block,
                                    const SleqpParams* params)
{
  const double eps = sleqp_params_get(params, SLEQP_PARAM_EPS);
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

static
SLEQP_RETCODE bfgs_block_push(BFGSBlock* block,
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

SLEQP_RETCODE sleqp_bfgs_data_push(SleqpBFGSData* data,
                                   SleqpIterate* previous_iterate,
                                   SleqpIterate* current_iterate,
                                   SleqpSparseVec* multipliers)
{
  const double eps = sleqp_params_get(data->params, SLEQP_PARAM_EPS);
  const double zero_eps = sleqp_params_get(data->params, SLEQP_PARAM_ZERO_EPS);

  const int num_blocks = data->num_blocks;

  SLEQP_CALL(sleqp_timer_start(data->update_timer));

  // Compute gradient difference
  {
    SLEQP_CALL(sleqp_sparse_matrix_trans_vector_product(sleqp_iterate_get_cons_jac(previous_iterate),
                                                        multipliers,
                                                        zero_eps,
                                                        data->prod_cache));

    SLEQP_CALL(sleqp_sparse_vector_add(data->prod_cache,
                                       sleqp_iterate_get_func_grad(previous_iterate),
                                       zero_eps,
                                       data->previous_grad));
  }

  {
    SLEQP_CALL(sleqp_sparse_matrix_trans_vector_product(sleqp_iterate_get_cons_jac(current_iterate),
                                                        multipliers,
                                                        zero_eps,
                                                        data->prod_cache));

    SLEQP_CALL(sleqp_sparse_vector_add(data->prod_cache,
                                       sleqp_iterate_get_func_grad(current_iterate),
                                       zero_eps,
                                       data->current_grad));
  }

  {
    SLEQP_CALL(sleqp_sparse_vector_add_scaled(data->previous_grad,
                                              data->current_grad,
                                              -1.,
                                              1.,
                                              zero_eps,
                                              data->grad_diff));
  }

  // Compute primal difference
  SLEQP_CALL(sleqp_sparse_vector_add_scaled(sleqp_iterate_get_primal(previous_iterate),
                                            sleqp_iterate_get_primal(current_iterate),
                                            -1.,
                                            1.,
                                            zero_eps,
                                            data->point_diff));

  int k_point = 0;
  int k_grad = 0;

  int offset = 0;

  for(int i = 0; i < num_blocks; ++i)
  {
    BFGSBlock* block = data->blocks + i;

    int next_offset = offset + block->dimension;

    SLEQP_CALL(sleqp_sparse_vector_clear(data->block_grad_diff));
    SLEQP_CALL(sleqp_sparse_vector_clear(data->block_point_diff));

    data->block_grad_diff->dim = block->dimension;
    data->block_point_diff->dim = block->dimension;

    while(k_point < data->point_diff->nnz &&
          data->point_diff->indices[k_point] < next_offset)
    {
      SLEQP_CALL(sleqp_sparse_vector_push(data->block_point_diff,
                                          data->point_diff->indices[k_point] - offset,
                                          data->point_diff->data[k_point]));

      ++k_point;
    }

    while(k_grad < data->grad_diff->nnz &&
          data->grad_diff->indices[k_grad] < next_offset)
    {
      SLEQP_CALL(sleqp_sparse_vector_push(data->block_grad_diff,
                                          data->grad_diff->indices[k_grad] - offset,
                                          data->grad_diff->data[k_grad]));

      ++k_grad;
    }

    const double point_normsq = sleqp_sparse_vector_norm_sq(data->block_point_diff);

    assert(sleqp_sparse_vector_valid(data->block_point_diff));
    assert(sleqp_sparse_vector_valid(data->block_grad_diff));

    if(!sleqp_is_zero(point_normsq, eps))
    {
      SLEQP_CALL(bfgs_block_push(block,
                                 data->params,
                                 data->block_point_diff,
                                 data->block_grad_diff));
    }


    offset = next_offset;
  }

  SLEQP_CALL(sleqp_timer_stop(data->update_timer));

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_bfgs_data_hess_prod(const SleqpBFGSData* data,
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

SleqpTimer* sleqp_bfgs_update_timer(SleqpBFGSData* data)
{
  return data->update_timer;
}

static SLEQP_RETCODE bfgs_data_free(SleqpBFGSData** star)
{
  SleqpBFGSData* data = *star;

  if(!data)
  {
    return SLEQP_OKAY;
  }

  SLEQP_CALL(sleqp_timer_free(&(data->update_timer)));

  SLEQP_CALL(sleqp_func_release(&data->bfgs_func));

  SLEQP_CALL(sleqp_sparse_vector_free(&data->block_prod));
  SLEQP_CALL(sleqp_sparse_vector_free(&data->block_direction));

  SLEQP_CALL(sleqp_sparse_vector_free(&data->prod_cache));

  SLEQP_CALL(sleqp_sparse_vector_free(&data->block_prod));
  SLEQP_CALL(sleqp_sparse_vector_free(&data->block_prod));

  SLEQP_CALL(sleqp_sparse_vector_free(&data->block_point_diff));
  SLEQP_CALL(sleqp_sparse_vector_free(&data->block_grad_diff));

  SLEQP_CALL(sleqp_sparse_vector_free(&data->current_grad));
  SLEQP_CALL(sleqp_sparse_vector_free(&data->previous_grad));

  SLEQP_CALL(sleqp_sparse_vector_free(&data->point_diff));
  SLEQP_CALL(sleqp_sparse_vector_free(&data->grad_diff));

  int num_blocks = data->num_blocks;

  for(int block = 0; block < num_blocks; ++block)
  {
    SLEQP_CALL(bfgs_block_free_at(data->blocks + block));
  }

  sleqp_free(&data->blocks);

  SLEQP_CALL(sleqp_options_release(&data->options));

  SLEQP_CALL(sleqp_params_release(&data->params));

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
