#include "fact_ma57.h"

#include <string.h>

#include "hsl_ma57.h"
#include "hsl_matrix.h"

#include "cmp.h"
#include "defs.h"
#include "error.h"
#include "fail.h"
#include "log.h"
#include "mem.h"

// TODO: Find out how to derive estimation
//       number from the values provided
static const bool ma57_est_cond     = false;
static const bool ma57_solve_refine = false;

static const bool ma57_verbose = false;

static const double ma57_size_factor = 1.2;

static SLEQP_RETCODE
ma57_get_error_string(int value, const char** message)
{
  switch (value)
  {
  case MA57_ERROR_INSUFFICIENT_INTEGER_SPACE:
    *message = "Insufficient integer space";
    break;
  case MA57_ERROR_INSUFFICIENT_REAL_SPACE:
    *message = "Insufficient real space";
    break;
  case MA57_SUCCESS_ZERO_SOLUTION:
    *message = "Success - zero solution";
    break;
  case MA57_SUCCESS_INDEFINITE:
    *message = "Success - indefinite";
    break;
  case MA57_SUCCESS_RANK_DEFICIENT:
    *message = "Success - rank deficient";
    break;
  case MA57_SUCCESS_ENTRIES_IGNORED:
    *message = "Success - entries ignored";
    break;
  case MA57_SUCCESS_DUPLICATE_ENTRIES_IGNORED:
    *message = "Success - duplicate entries ignored";
    break;
  case MA57_SUCCESS_INVALID_ENTRIES_IGNORED:
    *message = "Success - invalid entries ignored";
    break;
  case MA57_SUCCESS:
    *message = "Success";
    break;
  case MA57_N_OUT_OF_RANGE:
    *message = "<N> out of range";
    break;
  case MA57_NE_OUT_OF_RANGE:
    *message = "<NE> out of range";
    break;
  case MA57_INSUFFICIENT_REAL_SPACE:
    *message = "Insufficient real space";
    break;
  case MA57_INSUFFICIENT_INTEGER_SPACE:
    *message = "Insufficient integer space";
    break;
  case MA57_TINY_PIVOT:
    *message = "Tiny pivot detected";
    break;
  case MA57_PIVOT_SIGN_CHANGE:
    *message = "Pivot sign change detected";
    break;
  case MA57_COPY_TO_SMALLER_ARRAY:
    *message = "Attempt to copy to a smaller array";
    break;
  case MA57_NO_REFINEMENT_CONVERGENCE:
    *message = "Iterative refinement does not converge";
    break;
  case MA57_PERMUTATION_ERROR:
    *message = "Error in permutation information";
    break;
  case MA57_PIVOTING_CONTROL_ERROR:
    *message = "Error in pivoting control";
    break;
  case MA57_LRHS_OUT_OF_RANGE:
    *message = "<LRHS> out of range";
    break;
  case MA57_JOB_OUT_OF_RANGE:
    *message = "<JOB> out of range";
    break;
  case MA57_REFINEMENT_CONTROL_ERROR:
    *message = "Error in refinement control";
    break;
  case MA57_CONDITION_ESTIMATE_ERROR:
    *message = "Error in condition estimate";
    break;
  case MA57_LKEEP_OUT_OF_RANGE:
    *message = "<LKEEP> out of range";
    break;
  case MA57_NRHS_OUT_OF_RANGE:
    *message = "<NRHS> out of range";
    break;
  case MA57_LWORK_OUT_OF_RANGE:
    *message = "<LWORK> out of range";
    break;
  case MA57_METIS_PACKAGE_MISSING:
    *message = "METIS package missing in MA57";
    break;
  default:
    *message = "(unknown error)";
    break;
  }

  return SLEQP_OKAY;
}

#define MA57_CHECK_ERROR(x)                                                    \
  do                                                                           \
  {                                                                            \
    int ma57_status = (x);                                                     \
                                                                               \
    if (ma57_status < 0)                                                       \
    {                                                                          \
      const char* ma57_error_string;                                           \
      SLEQP_CALL(ma57_get_error_string(ma57_status, &ma57_error_string));      \
                                                                               \
      sleqp_raise(SLEQP_INTERNAL_ERROR,                                        \
                  "Caught hsl_ma57 error <%d> (%s) in function %s",            \
                  ma57_status,                                                 \
                  ma57_error_string,                                           \
                  __func__);                                                   \
    }                                                                          \
  } while (0)

/** @brief MA57 integer control structure
 */
typedef struct MA57IControl
{
  int32_t error_message_stream;
  int32_t warning_message_stream;
  int32_t monitor_message_stream;
  int32_t stats_message_stream; /* not operative? */
  int32_t print_level;
  int32_t pivot_order;
  int32_t numerical_pivoting;
  int32_t flag_continue_factor;
  int32_t permitted_itref_steps;
  int32_t flag_estimate_cond;
  int32_t blas_block_size;
  int32_t merge_tree_nodes_threshold;
  int32_t level23_blas_threshold;
  int32_t dense_row_threshold;
  int32_t flag_scaling;
  int32_t flag_drop_small_pivots;
  int32_t unused[4];
} MA57IControl;

/** @brief MA57 real control structure
 */
typedef struct MA57RControl
{
  double pivoting_threshold;
  double zero_threshold;
  double refinement_contraction;
  double replacement_pivot;
  double static_pivoting_threshold;
} MA57RControl;

/** @brief MA57 integer information structure
 */
typedef struct MA57IInfo
{
  int32_t error; //< Success, error, or warning code
  int32_t error2;
  // ma57_symbolic
  int32_t n_out_of_range;
  int32_t n_off_diagonal_duplicates;
  int32_t n_forecast_reals_in_factor;
  int32_t n_forecast_ints_in_factor;
  int32_t forecast_frontal_matrix_size;
  int32_t n_nodes_in_tree;
  int32_t forecast_lfact_compress;
  int32_t forecast_lifact_compress;
  int32_t forecast_lfact_no_compress;
  int32_t forecast_lifact_no_compress;
  int32_t n_data_compresses;
  // ma57_factorize
  int32_t n_entries_in_factor;
  int32_t n_reals_in_factor;
  int32_t n_ints_in_factor;
  int32_t lfact_compress;
  int32_t lifact_compress;
  int32_t lfact_no_compress;
  int32_t lifact_no_compress;
  int32_t frontal_matrix_size;
  int32_t n_2x2_pivots;
  int32_t n_delayed_pivots;
  int32_t n_negative_eigenvalues;
  int32_t rank;
  int32_t n_pivot_sign_changes;
  int32_t pivot_step_commence;
  int32_t n_real_compresses;
  int32_t n_int_compresses;
  int32_t n_block_pivots;
  int32_t n_zeros_in_triangles;
  int32_t n_zeros_in_rectangle;
  int32_t n_zero_columns;
  int32_t n_static_pivots;
  // ma57_refine
  int32_t n_refinement_steps;
  // ma57_symbolic (again)
  int32_t ordering_used;
  // unused by ma57
  int32_t unused[4];
} MA57IInfo;

/** @brief MA57 real information structure
 */
typedef struct MA57RInfo
{
  // ma57_symbolic
  double n_forecast_assembly_fpadds;
  double n_forecast_elimination_flops;
  // ma57_factorize
  double n_assembly_fpadds;
  double n_elimination_flops;
  double n_gemm_extra_elimination_flops;
  // ma57_refine
  double error1;
  double error2;
  double matrix_inf_norm;
  double solution_inf_norm;
  double scaled_residuals_norm;
  double cond1;
  double cond2;
  double error_inf_norm_bound;
  // ma57_factorize (again)
  double max_pivot_change;
  double smallest_pivot;
  double min_scaling;
  double max_scaling;
  double max_modulus;
  // unused by ma57
  double unused[2];
} MA57RInfo;

/** @brief MA57 solver control and information structure
 */
typedef struct MA57ControlInfo
{
  union
  {
    double cntl_[5];
    MA57RControl cntl;
  };
  union
  {
    int32_t icntl_[20];
    MA57IControl icntl;
  };
  union
  {
    double rinfo_[20];
    MA57RInfo rinfo;
  };
  union
  {
    int32_t info_[40];
    MA57IInfo info;
  };
} MA57ControlInfo;

typedef struct MA57Factor
{
  double* factor;
  int32_t factor_size;

  int32_t* ifactor;
  int32_t ifactor_size;

} MA57Factor;

typedef struct MA57Workspace
{
  int32_t* keep;
  int32_t keep_size;

  int32_t* iwork;
  int32_t iwork_size;

  double* work;
  int32_t work_size;

  double* sol;
  double* rhs;
  double* res;
  int32_t rhs_sol_size;

} MA57Workspace;

typedef struct MA57Data
{
  MA57ControlInfo control_info;

  MA57Workspace workspace;

  MA57Factor factor;

  HSLMatrix matrix;

} MA57Data;

static SLEQP_RETCODE
ma57_symbolic(MA57Data* ma57_data)
{
  MA57ControlInfo* control_info = &(ma57_data->control_info);
  MA57Workspace* ma57_workspace = &(ma57_data->workspace);
  HSLMatrix* hsl_matrix         = &(ma57_data->matrix);

  const int32_t dim = hsl_matrix->dim;
  const int32_t nnz = hsl_matrix->nnz;

  const int32_t* rows = hsl_matrix->rows;
  const int32_t* cols = hsl_matrix->cols;

  const int32_t keep_size = ma57_workspace->keep_size;
  int32_t* keep           = ma57_workspace->keep;
  int32_t* iwork          = ma57_workspace->iwork;

  ma57ad_(&dim,
          &nnz,
          rows,
          cols,
          &keep_size,
          keep,
          iwork,
          control_info->icntl_,
          control_info->info_,
          control_info->rinfo_);

  MA57_CHECK_ERROR(control_info->info.error);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
ma57_increase_factor_size(MA57Data* ma57_data)
{
  MA57ControlInfo* control_info = &(ma57_data->control_info);
  MA57Factor* ma57_factor       = &(ma57_data->factor);
  MA57Workspace* ma57_workspace = &(ma57_data->workspace);
  HSLMatrix* hsl_matrix         = &(ma57_data->matrix);

  const int32_t dim = hsl_matrix->dim;
  int32_t* keep     = ma57_workspace->keep;

  const int32_t ic = 0; // copy real array

  double* factor            = ma57_factor->factor;
  const int32_t factor_size = ma57_factor->factor_size;

  double* new_factor;
  const int32_t new_factor_size = ma57_size_factor * factor_size;

  assert(new_factor_size > factor_size);

  SLEQP_CALL(sleqp_alloc_array(&new_factor, new_factor_size));

  int32_t empty_int = 0;

  ma57ed_(&dim,
          &ic,
          keep,
          factor,
          &factor_size,
          new_factor,
          &new_factor_size,
          &empty_int,
          &empty_int,
          &empty_int,
          &empty_int,
          control_info->info_);

  sleqp_free(&factor);

  ma57_factor->factor      = new_factor;
  ma57_factor->factor_size = new_factor_size;

  MA57_CHECK_ERROR(control_info->info.error);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
ma57_increase_ifactor_size(MA57Data* ma57_data)
{
  MA57ControlInfo* control_info = &(ma57_data->control_info);
  MA57Factor* ma57_factor       = &(ma57_data->factor);
  MA57Workspace* ma57_workspace = &(ma57_data->workspace);
  HSLMatrix* hsl_matrix         = &(ma57_data->matrix);

  const int32_t dim = hsl_matrix->dim;
  int32_t* keep     = ma57_workspace->keep;

  const int32_t ic = 1; // copy integer array

  int32_t* ifactor           = ma57_factor->ifactor;
  const int32_t ifactor_size = ma57_factor->ifactor_size;

  int32_t* new_ifactor;
  const int32_t new_ifactor_size = ma57_size_factor * ifactor_size;

  assert(new_ifactor_size > ifactor_size);

  SLEQP_CALL(sleqp_alloc_array(&new_ifactor, new_ifactor_size));

  int32_t empty_int = 0;
  double empty_dbl  = 0;

  ma57ed_(&dim,
          &ic,
          keep,
          &empty_dbl,
          &empty_int,
          &empty_dbl,
          &empty_int,
          ifactor,
          &ifactor_size,
          new_ifactor,
          &new_ifactor_size,
          control_info->info_);

  sleqp_free(&ifactor);

  ma57_factor->ifactor      = new_ifactor;
  ma57_factor->ifactor_size = new_ifactor_size;

  MA57_CHECK_ERROR(control_info->info.error);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
ma57_numeric(MA57Data* ma57_data)
{
  MA57ControlInfo* control_info = &(ma57_data->control_info);
  MA57Factor* ma57_factor       = &(ma57_data->factor);
  MA57Workspace* ma57_workspace = &(ma57_data->workspace);
  HSLMatrix* hsl_matrix         = &(ma57_data->matrix);

  const int32_t dim = hsl_matrix->dim;
  const int32_t nnz = hsl_matrix->nnz;

  const double* data = hsl_matrix->data;

  const int32_t keep_size = ma57_workspace->keep_size;
  int32_t* keep           = ma57_workspace->keep;
  int32_t* iwork          = ma57_workspace->iwork;

  bool insufficient_size = false;

  do
  {
    double* factor            = ma57_factor->factor;
    const int32_t factor_size = ma57_factor->factor_size;

    int32_t* ifactor           = ma57_factor->ifactor;
    const int32_t ifactor_size = ma57_factor->ifactor_size;

    insufficient_size = false;

    ma57bd_(&dim,
            &nnz,
            data,
            factor,
            &factor_size,
            ifactor,
            &ifactor_size,
            &keep_size,
            keep,
            iwork,
            control_info->icntl_,
            control_info->cntl_,
            control_info->info_,
            control_info->rinfo_);

    if (control_info->info.error == MA57_INSUFFICIENT_REAL_SPACE)
    {
      insufficient_size = true;

      SLEQP_CALL(ma57_increase_factor_size(ma57_data));
    }

    if (control_info->info.error == MA57_INSUFFICIENT_INTEGER_SPACE)
    {
      insufficient_size = true;

      SLEQP_CALL(ma57_increase_ifactor_size(ma57_data));
    }

  } while (insufficient_size);

  MA57_CHECK_ERROR(control_info->info.error);

  return SLEQP_OKAY;
}

static int
req_work_size(const int dim, MA57IControl* icntl)
{
  if (icntl->permitted_itref_steps == 1)
  {
    return dim;
  }
  else if (icntl->permitted_itref_steps > 1)
  {
    if (icntl->flag_estimate_cond > 0)
    {
      return 4 * dim;
    }

    return 3 * dim;
  }

  return 0;
}

static SLEQP_RETCODE
ma57_set_matrix(void* fact_data, SleqpSparseMatrix* matrix)
{
  MA57Data* ma57_data = (MA57Data*)fact_data;

  MA57ControlInfo* control_info = &(ma57_data->control_info);
  MA57Factor* ma57_factor       = &(ma57_data->factor);
  MA57Workspace* ma57_workspace = &(ma57_data->workspace);
  HSLMatrix* hsl_matrix         = &(ma57_data->matrix);

  const int num_rows = sleqp_sparse_matrix_num_cols(matrix);
  const int num_cols = sleqp_sparse_matrix_num_cols(matrix);
  assert(num_rows == num_cols);

  const int dim = num_cols;

  {
    const int required_work_size = req_work_size(dim, &(control_info->icntl));

    if (ma57_workspace->rhs_sol_size < dim)
    {
      SLEQP_CALL(sleqp_realloc(&(ma57_workspace->sol), dim));

      if (ma57_solve_refine)
      {
        SLEQP_CALL(sleqp_realloc(&(ma57_workspace->rhs), dim));
        SLEQP_CALL(sleqp_realloc(&(ma57_workspace->res), dim));
      }

      ma57_workspace->rhs_sol_size = dim;
    }

    if (ma57_workspace->work_size < required_work_size)
    {
      SLEQP_CALL(sleqp_realloc(&(ma57_workspace->work), required_work_size));
      ma57_workspace->work_size = required_work_size;
    }
  }

  // Convert matrix to HSL format
  {
    SLEQP_CALL(hsl_matrix_set(hsl_matrix, matrix));
  }

  // Prepare symbolic workspace
  {
    const int32_t nnz = hsl_matrix->nnz;

    const int32_t keep_size = 5 * dim + nnz + SLEQP_MAX(dim, nnz) + 42;

    if (ma57_workspace->keep_size < keep_size)
    {
      SLEQP_CALL(sleqp_realloc(&(ma57_workspace->keep), keep_size));

      memset(ma57_workspace->keep, 0, keep_size * sizeof(int32_t));

      ma57_workspace->keep_size = keep_size;
    }

    const int32_t iwork_size = 5 * dim;

    if (ma57_workspace->iwork_size < iwork_size)
    {
      SLEQP_CALL(sleqp_realloc(&(ma57_workspace->iwork), iwork_size));

      ma57_workspace->iwork_size = iwork_size;
    }
  }

  SLEQP_CALL(ma57_symbolic(ma57_data));

  {
    const int32_t factor_size
      = SLEQP_MAX(control_info->info.forecast_lfact_compress,
                  control_info->info.forecast_lfact_no_compress);

    if (ma57_factor->factor_size < factor_size)
    {
      SLEQP_CALL(sleqp_realloc(&(ma57_factor->factor), factor_size));
      ma57_factor->factor_size = factor_size;
    }

    const int32_t ifactor_size
      = SLEQP_MAX(control_info->info.forecast_lifact_compress,
                  control_info->info.forecast_lifact_no_compress);

    if (ma57_factor->ifactor_size < ifactor_size)
    {
      SLEQP_CALL(sleqp_realloc(&(ma57_factor->ifactor), ifactor_size));
      ma57_factor->ifactor_size = ifactor_size;
    }
  }

  SLEQP_CALL(ma57_numeric(ma57_data));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
ma57_solve(void* fact_data, const SleqpVec* rhs)
{
  MA57Data* ma57_data = (MA57Data*)fact_data;

  MA57ControlInfo* control_info = &(ma57_data->control_info);
  MA57Factor* ma57_factor       = &(ma57_data->factor);
  MA57Workspace* ma57_workspace = &(ma57_data->workspace);
  HSLMatrix* hsl_matrix         = &(ma57_data->matrix);

  const int32_t dim  = hsl_matrix->dim;
  const int32_t nnz  = hsl_matrix->nnz;
  const int32_t nrhs = 1;

  assert(dim == rhs->dim);

  double* factor            = ma57_factor->factor;
  const int32_t factor_size = ma57_factor->factor_size;

  int32_t* ifactor           = ma57_factor->ifactor;
  const int32_t ifactor_size = ma57_factor->ifactor_size;

  double* work            = ma57_workspace->work;
  const int32_t work_size = ma57_workspace->work_size;

  int32_t* iwork = ma57_workspace->iwork;

  if (ma57_solve_refine)
  {
    const int32_t job = 0; // Solve Ax=b

    double* rhs_dense = ma57_workspace->rhs;
    double* sol_dense = ma57_workspace->sol;
    double* res_dense = ma57_workspace->res;

    SLEQP_CALL(sleqp_vec_to_raw(rhs, rhs_dense));

    ma57dd_(&job,                  // const int32_t *job,
            &dim,                  // const int32_t *n,
            &nnz,                  // const int32_t *ne,
            hsl_matrix->data,      // const double *a,
            hsl_matrix->rows,      // const int32_t *irn,
            hsl_matrix->cols,      // const int32_t *jcn,
            factor,                // const double *fact,
            &factor_size,          // const int32_t *lfact,
            ifactor,               // const int32_t *ifact,
            &ifactor_size,         // const int32_t *lifact,
            rhs_dense,             // double *rhs,
            sol_dense,             // double *x,
            res_dense,             // double *resid,
            work,                  // double *work,
            iwork,                 // int32_t *iwork,
            control_info->icntl_,  // int32_t *icntl,
            control_info->cntl_,   // double *cntl,
            control_info->info_,   // int32_t *info,
            control_info->rinfo_); // double *rinfo
  }
  else
  {
    const int32_t job = 1; // Solve Ax=b

    double* rhs_sol = ma57_workspace->sol;

    SLEQP_CALL(sleqp_vec_to_raw(rhs, rhs_sol));

    ma57cd_(&job,
            &dim,
            factor,
            &factor_size,
            ifactor,
            &ifactor_size,
            &nrhs,
            rhs_sol,
            &dim,
            work,
            &work_size,
            iwork,
            control_info->icntl_,
            control_info->info_);
  }

  MA57_CHECK_ERROR(control_info->info.error);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
ma57_solution(void* fact_data,
              SleqpVec* sol,
              int begin,
              int end,
              double zero_eps)
{
  MA57Data* ma57_data = (MA57Data*)fact_data;

  MA57Workspace* ma57_workspace = &(ma57_data->workspace);

  SLEQP_CALL(sleqp_vec_set_from_raw(sol,
                                    ma57_workspace->sol + begin,
                                    end - begin,
                                    zero_eps));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
ma57_data_create(MA57Data** star)
{
  SLEQP_CALL(sleqp_malloc(star));

  MA57Data* ma57_data = *star;

  *ma57_data = (MA57Data){0};

  MA57ControlInfo* control_info = &(ma57_data->control_info);

  ma57id_(control_info->cntl_, control_info->icntl_);

  if (ma57_verbose)
  {
    // Print limited number of stats
    control_info->icntl.print_level = 3;
  }
  else
  {
    // Only print errors
    control_info->icntl.print_level = 1;
  }

  if (ma57_est_cond)
  {
    control_info->icntl.flag_estimate_cond = 1;
  }

  // AMD using MC47
  // There seem to be some issues with MeTiS integration
  control_info->icntl.pivot_order = 2;

  MA57_CHECK_ERROR(control_info->info.error);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
ma57_condition_estimate(void* fact_data, double* condition_estimate)
{
  // MA57Data* ma57_data = (MA57Data*) fact_data;

  (*condition_estimate) = SLEQP_NONE;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
ma57_free(void** star)
{
  MA57Data* ma57_data = (MA57Data*)(*star);

  MA57Factor* ma57_factor       = &(ma57_data->factor);
  MA57Workspace* ma57_workspace = &(ma57_data->workspace);
  HSLMatrix* hsl_matrix         = &(ma57_data->matrix);

  SLEQP_CALL(hsl_matrix_clear(hsl_matrix));

  sleqp_free(&(ma57_factor->factor));
  sleqp_free(&(ma57_factor->ifactor));

  sleqp_free(&(ma57_workspace->keep));
  sleqp_free(&(ma57_workspace->iwork));
  sleqp_free(&(ma57_workspace->work));

  sleqp_free(&(ma57_workspace->res));
  sleqp_free(&(ma57_workspace->rhs));
  sleqp_free(&(ma57_workspace->sol));

  sleqp_free(&ma57_data);

  *star = NULL;

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_fact_ma57_create(SleqpFact** star, SleqpParams* params)
{
  SleqpFactorizationCallbacks callbacks
    = {.set_matrix         = ma57_set_matrix,
       .solve              = ma57_solve,
       .solution           = ma57_solution,
       .condition_estimate = ma57_condition_estimate,
       .free               = ma57_free};

  MA57Data* ma57_data;

  SLEQP_CALL(ma57_data_create(&ma57_data));

  SLEQP_CALL(sleqp_fact_create(star,
                               SLEQP_FACT_MA57_NAME,
                               SLEQP_FACT_MA57_VERSION,
                               params,
                               &callbacks,
                               SLEQP_FACT_FLAGS_LOWER,
                               (void*)ma57_data));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_fact_create_default(SleqpFact** star, SleqpParams* params)
{
  SLEQP_CALL(sleqp_fact_ma57_create(star, params));

  return SLEQP_OKAY;
}
