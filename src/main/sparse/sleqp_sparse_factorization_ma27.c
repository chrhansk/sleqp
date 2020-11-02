#include "sleqp_sparse_factorization_ma27.h"

#include <assert.h>
#include <string.h>

#include "hsl_ma27.h"

#include "sleqp_cmp.h"
#include "sleqp_mem.h"

// 20 % increase as suggested in the manual
static const double ma27_alloc_factor = 1.2;

static const double ma27_fratio = 0.5;

static const bool ma27_verbose = false;


static SLEQP_RETCODE ma27_get_error_string(int value, const char** message)
{
  switch(value)
  {
  case MA27_NSTEPS_OUT_OF_RANGE:
    *message= "NSTEPS out of range";
    break;
  case MA27_PIVOT_SIGN_CHANGE:
    *message= "Pivot sign change detected";
    break;
  case MA27_SINGULAR_MATRIX:
    *message= "Singular matrix detected";
    break;
  case MA27_A_MEM_TOO_SMALL:
    *message= "Memory for A too small";
    break;
  case MA27_IW_MEM_TOO_SMALL:
    *message= "Memory for IW too small";
    break;
  case MA27_NZ_OUT_OF_RANGE:
    *message= "NZ out of range";
    break;
  case MA27_N_OUT_OF_RANGE:
    *message= "N out of range";
    break;
  case MA27_SUCCESS:
    *message= "Success";
    break;
  case MA27_WARN_IRN_ICN_OUT_OF_RANGE:
    *message= "Warning - IRN or ICN entry out of range";
    break;
  case MA27_WARN_INDEFINITE:
    *message= "Warning - Indefinite matrix";
    break;
  case MA27_WARN_RANK_DEFICIENT:
    *message= "Warning - Rank deficient matrix";
    break;
  default:
    *message= "<unknown hsl_ma27 error code>";
    break;
  }

  return SLEQP_OKAY;
}

#define MA27_CHECK_ERROR(x)                                   \
  do                                                          \
  {                                                           \
    int ma27_status = (x);                                    \
                                                              \
    if (ma27_status != MA27_SUCCESS)                          \
    {                                                         \
      const char* ma27_error_string;                          \
      SLEQP_CALL(ma27_get_error_string(ma27_status,           \
                                       &ma27_error_string));  \
                                                              \
      sleqp_log_error("Caught hsl_ma27 error <%d> (%s)",      \
                      ma27_status,                            \
                      ma27_error_string);                     \
                                                              \
      switch(ma27_status)                                     \
      {                                                       \
      case MA27_NSTEPS_OUT_OF_RANGE:                          \
      case MA27_NZ_OUT_OF_RANGE:                              \
      case MA27_N_OUT_OF_RANGE:                               \
      case MA27_WARN_IRN_ICN_OUT_OF_RANGE:                    \
        return SLEQP_ILLEGAL_ARGUMENT;                        \
      case MA27_A_MEM_TOO_SMALL:                              \
      case MA27_IW_MEM_TOO_SMALL:                             \
        return SLEQP_NOMEM;                                   \
      case MA27_PIVOT_SIGN_CHANGE:                            \
      case MA27_SINGULAR_MATRIX:                              \
      case MA27_WARN_INDEFINITE:                              \
      case MA27_WARN_RANK_DEFICIENT:                          \
        return SLEQP_MATH_ERROR;                              \
      default:                                                \
        return SLEQP_INTERNAL_ERROR;                          \
      }                                                       \
    }                                                         \
  } while(0)


typedef struct MA27IControl
{
  int32_t error_message_stream;
  int32_t diag_message_stream;
  int32_t print_level;
  int32_t not_of_interest[22];  /* as per MA27 documentation */
  int32_t unused[5];
} MA27IControl;

typedef struct MA27RControl
{
  double pivot_control; /* in [-0.5,+0.5] */
  double fratio;
  double pivtol;
  double unused[2];
} MA27RControl;

typedef struct MA27IInfo
{
  int32_t iflag;
  int32_t ierror;
  int32_t nrltot;
  int32_t nirtot;
  int32_t nrlnec;
  int32_t nirnec;
  int32_t nrladu;
  int32_t niradu;
  int32_t nrlbdu;
  int32_t nirbdu;
  int32_t ncmpa;
  int32_t ncmpbr;
  int32_t ncmpbi;
  int32_t ntwo;
  int32_t neig;
  int32_t unused[5];
} MA27IInfo;

typedef struct MA27ControlInfo
{
  union {
    double cntl_[5];
    MA27RControl cntl;
  };
  union {
    int32_t icntl_[30];
    MA27IControl icntl;
  };
  union {
    int32_t info_[20];
    MA27IInfo info;
  };
} MA27ControlInfo;


typedef struct MA27Matrix
{
  double* data;
  int32_t* rows;
  int32_t* cols;

  int32_t max_total_nnz;
  int32_t total_nnz;
  int32_t dim;

} MA27Matrix;

typedef struct MA27Factor
{
  double* factor;
  int32_t factor_size;

} MA27Factor;

typedef struct MA27Workspace
{
  int32_t* iw;
  int32_t iw_size;

  int32_t* iw1;
  int32_t iw1_size;

  int32_t* ikeep;
  int32_t ikeep_size;

  double* w;
  int32_t w_size;

} MA27Workspace;

typedef struct MA27State
{
  int32_t nsteps;

  int32_t maxfrt;

} MA27State;

typedef struct MA27Data
{
  MA27ControlInfo control_info;

  MA27Matrix matrix;

  MA27State state;

  MA27Workspace workspace;

  MA27Factor factor;

  double* rhs_sol;

  double ops;

} MA27Data;

static SLEQP_RETCODE ma27_matrix_reserve(MA27Matrix* ma27_matrix,
                                         int32_t total_nnz)
{
  if(ma27_matrix->max_total_nnz < total_nnz)
  {
    SLEQP_CALL(sleqp_realloc(&(ma27_matrix->cols), total_nnz));
    SLEQP_CALL(sleqp_realloc(&(ma27_matrix->rows), total_nnz));
    SLEQP_CALL(sleqp_realloc(&(ma27_matrix->data), total_nnz));

    ma27_matrix->max_total_nnz = total_nnz;
  }

  return SLEQP_OKAY;
}


static SLEQP_RETCODE ma27_symbolic(MA27Data* ma27_data)
{
  MA27ControlInfo* control_info = &(ma27_data->control_info);

  MA27Matrix* ma27_matrix = &(ma27_data->matrix);
  MA27Workspace* ma27_workspace = &(ma27_data->workspace);
  MA27State* ma27_state = &(ma27_data->state);

  const int32_t dim = ma27_matrix->dim;
  const int32_t total_nnz = ma27_matrix->total_nnz;

  const int32_t* rows = ma27_matrix->rows;
  const int32_t* cols = ma27_matrix->cols;

  int32_t iflag = 0;

  ma27ad_(&dim, &total_nnz, rows, cols, ma27_workspace->iw, &(ma27_workspace->iw_size),
          ma27_workspace->ikeep, ma27_workspace->iw1, &(ma27_state->nsteps), &iflag, control_info->icntl_, control_info->cntl_,
          control_info->info_, &ma27_data->ops);

  MA27_CHECK_ERROR(control_info->info.iflag);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE ma27_factor_reserve(MA27Factor* ma27_factor, int32_t factor_size)
{
  if(ma27_factor->factor_size < factor_size)
  {
    SLEQP_CALL(sleqp_realloc(&ma27_factor->factor, factor_size));
    ma27_factor->factor_size = factor_size;
  }

  return SLEQP_OKAY;
}

static SLEQP_RETCODE ma27_numeric(MA27Data* ma27_data)
{
  MA27ControlInfo* control_info = &(ma27_data->control_info);

  MA27Matrix* ma27_matrix = &(ma27_data->matrix);
  MA27Workspace* ma27_workspace = &(ma27_data->workspace);
  MA27State* ma27_state = &(ma27_data->state);
  MA27Factor* ma27_factor = &(ma27_data->factor);

  // reserve sufficient number of doubles
  {
    int32_t factor_size = ma27_alloc_factor*SLEQP_MAX(control_info->info.nrlnec, control_info->info.nrltot);

    factor_size = SLEQP_MAX(factor_size, ma27_matrix->total_nnz);

    ma27_factor_reserve(ma27_factor, factor_size);
  }

  // reserve sufficient number of integers
  {
    const int32_t iw_size = SLEQP_MAX(control_info->info.nirnec, control_info->info.nirtot);

    if(ma27_workspace->iw_size < iw_size)
    {
      SLEQP_CALL(sleqp_realloc(&(ma27_workspace->iw), iw_size));
      ma27_workspace->iw_size = iw_size;
    }
  }

  const int32_t dim = ma27_matrix->dim;

  bool insufficient_size;

  do
  {
    insufficient_size = false;

    const int32_t total_nnz = ma27_matrix->total_nnz;
    const int32_t factor_size = ma27_factor->factor_size;
    const int32_t* rows = ma27_matrix->rows;
    const int32_t* cols = ma27_matrix->cols;
    const double* matrix_data = ma27_matrix->data;
    double* factor_data = ma27_factor->factor;


    for(int32_t i = 0; i < total_nnz; ++i)
    {
      factor_data[i] = matrix_data[i];
    }

    ma27bd_(&dim, &total_nnz, rows, cols, factor_data,
            &factor_size, ma27_workspace->iw, &(ma27_workspace->iw_size), ma27_workspace->ikeep, &(ma27_state->nsteps),
            &(ma27_state->maxfrt), ma27_workspace->iw1, control_info->icntl_, control_info->cntl_, control_info->info_);

    if(control_info->info.iflag == MA27_A_MEM_TOO_SMALL)
    {
      insufficient_size = true;

      int32_t required_size = control_info->info.ierror;

      required_size = SLEQP_MAX(required_size, 1.2*ma27_factor->factor_size);

      ma27_factor_reserve(ma27_factor, required_size);
    }
  }
  while(insufficient_size);

  MA27_CHECK_ERROR(control_info->info.iflag);

  return SLEQP_OKAY;
}

// Switch to 1-based indexing, only keep upper diagonal
static SLEQP_RETCODE ma27_data_fill_matrix(MA27Matrix* ma27_matrix,
                                           SleqpSparseMatrix* matrix)
{
  const int dim = ma27_matrix->dim;

  double* ma27_data = ma27_matrix->data;
  int32_t* ma27_rows = ma27_matrix->rows;
  int32_t* ma27_cols = ma27_matrix->cols;

  const double* matrix_data = sleqp_sparse_matrix_get_data(matrix);
  const int* matrix_cols = sleqp_sparse_matrix_get_cols(matrix);
  const int* matrix_rows = sleqp_sparse_matrix_get_rows(matrix);
  const int matrix_nnz = sleqp_sparse_matrix_get_nnz(matrix);

  int32_t ma27_pos = 0;
  int32_t col = 0;

  for(int index = 0; index < matrix_nnz; ++index)
  {
    while(index >= matrix_cols[col + 1])
    {
      ++col;
    }

    const int32_t row = matrix_rows[index];
    const double entry = matrix_data[index];

    // Convert indices to be 1-based
    if(row <= col)
    {
      ma27_data[ma27_pos] = entry;
      ma27_cols[ma27_pos] = col + 1;
      ma27_rows[ma27_pos] = row + 1;

      ++ma27_pos;
    }
  }

  ma27_matrix->total_nnz = ma27_pos;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE ma27_data_set_matrix(void* factorization_data,
                                          SleqpSparseMatrix* matrix)
{
  MA27Data* ma27_data = (MA27Data*) factorization_data;

  MA27Matrix* ma27_matrix = &(ma27_data->matrix);
  MA27Workspace* ma27_workspace = &(ma27_data->workspace);

  const int32_t num_rows = sleqp_sparse_matrix_get_num_cols(matrix);
  const int32_t num_cols = sleqp_sparse_matrix_get_num_cols(matrix);
  assert(num_rows == num_cols);

  const int32_t dim = num_cols;

  const int32_t matrix_nnz = sleqp_sparse_matrix_get_nnz(matrix);
  const int32_t total_nnz = dim + matrix_nnz;

  // prepare rhs / sol
  {
    if(ma27_matrix->dim < dim)
    {
      SLEQP_CALL(sleqp_realloc(&ma27_data->rhs_sol, dim));
    }
  }

  // Convert matrix to MA27 format
  {
    SLEQP_CALL(ma27_matrix_reserve(ma27_matrix, total_nnz));

    ma27_matrix->dim = dim;

    SLEQP_CALL(ma27_data_fill_matrix(ma27_matrix, matrix));
  }

  //

  // Prepare symbolic workspace
  {
    const int32_t upper_nnz = total_nnz - dim;

    const int32_t iw_size = ma27_alloc_factor * SLEQP_MAX(2*upper_nnz + 3*dim + 1, upper_nnz + 3*dim + 1);

    if(ma27_workspace->iw_size < iw_size)
    {
      SLEQP_CALL(sleqp_realloc(&(ma27_workspace->iw), iw_size));
      ma27_workspace->iw_size = iw_size;
    }

    const int32_t ikeep_size = 3*dim;

    if(ma27_workspace->ikeep_size < ikeep_size)
    {
      SLEQP_CALL(sleqp_realloc(&(ma27_workspace->ikeep), ikeep_size));
      ma27_workspace->ikeep_size = ikeep_size;
    }

    const int32_t iw1_size = 2*dim;

    if(ma27_workspace->iw1_size < iw1_size)
    {
      SLEQP_CALL(sleqp_realloc(&(ma27_workspace->iw1), iw1_size));
      ma27_workspace->iw1_size = iw1_size;
    }
  }

  SLEQP_CALL(ma27_symbolic(ma27_data));

  SLEQP_CALL(ma27_numeric(ma27_data));

  return SLEQP_OKAY;
}


static MA27_ERROR ma27_solve(MA27Data* ma27_data)
{
  MA27Matrix* ma27_matrix = &(ma27_data->matrix);
  MA27Workspace* ma27_workspace = &(ma27_data->workspace);
  MA27State* ma27_state = &(ma27_data->state);
  MA27Factor* ma27_factor = &(ma27_data->factor);

  MA27ControlInfo* control_info = &(ma27_data->control_info);

  const int32_t dim = ma27_matrix->dim;
  const int32_t total_nnz = ma27_matrix->total_nnz;
  const int32_t max_total_nnz = ma27_matrix->max_total_nnz;
  const int32_t nsteps = ma27_state->nsteps;

  const double* factor_data = ma27_factor->factor;

  double* rhs_sol = ma27_data->rhs_sol;

  {
    const int32_t w_size = ma27_state->maxfrt;

    if(ma27_workspace->w_size < w_size)
    {
      SLEQP_CALL(sleqp_realloc(&(ma27_workspace->w), w_size));

      ma27_workspace->w_size = w_size;
    }
  }

  {
    const int32_t iw1_size = nsteps;

    if(ma27_workspace->iw1_size < iw1_size)
    {
      SLEQP_CALL(sleqp_realloc(&(ma27_workspace->iw1), iw1_size));

      ma27_workspace->iw1_size = iw1_size;
    }
  }



  ma27cd_(&dim, factor_data, &max_total_nnz, ma27_workspace->iw, &(ma27_workspace->iw_size), ma27_workspace->w,
          &ma27_state->maxfrt, rhs_sol, ma27_workspace->iw1, &nsteps, control_info->icntl_, control_info->info_);

  MA27_CHECK_ERROR(control_info->info.iflag);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE ma27_data_solve(void* factorization_data,
                                     SleqpSparseVec* rhs)
{
  MA27Data* ma27_data = (MA27Data*) factorization_data;
  MA27Matrix* ma27_matrix = &(ma27_data->matrix);

  assert(rhs->dim == ma27_matrix->dim);

  SLEQP_CALL(sleqp_sparse_vector_to_raw(rhs, ma27_data->rhs_sol));

  SLEQP_CALL(ma27_solve(ma27_data));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE ma27_data_get_sol(void* factorization_data,
                                       SleqpSparseVec* sol,
                                       int begin,
                                       int end,
                                       double zero_eps)
{
  MA27Data* ma27_data = (MA27Data*) factorization_data;

  SLEQP_CALL(sleqp_sparse_vector_from_raw(sol,
                                          ma27_data->rhs_sol + begin,
                                          end - begin,
                                          zero_eps));


  return SLEQP_OKAY;
}

static SLEQP_RETCODE ma27_data_get_condition_estimate(void* factorization_data,
                                                      double* condition_estimate)
{
  MA27Data* ma27_data = (MA27Data*) factorization_data;

  (*condition_estimate) = SLEQP_NONE;

  return SLEQP_OKAY;
}

// .solve = umfpack_sparse_factorization_solve,
// .get_sol = umfpack_sparse_factorization_get_sol,

static SLEQP_RETCODE ma27_data_free(void** star)
{
  MA27Data* ma27_data = (MA27Data*) (*star);

  MA27Matrix* ma27_matrix = &(ma27_data->matrix);
  MA27Workspace* ma27_workspace = &(ma27_data->workspace);
  MA27Factor* ma27_factor = &(ma27_data->factor);

  sleqp_free(&(ma27_matrix->cols));
  sleqp_free(&(ma27_matrix->rows));
  sleqp_free(&(ma27_matrix->data));

  sleqp_free(&(ma27_workspace->iw));
  sleqp_free(&(ma27_workspace->iw1));
  sleqp_free(&(ma27_workspace->ikeep));
  sleqp_free(&(ma27_workspace->w));

  sleqp_free(&ma27_data->rhs_sol);

  sleqp_free(&ma27_factor->factor);

  *star = NULL;

  sleqp_free(&ma27_data);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE ma27_data_create(MA27Data** star)
{
  SLEQP_CALL(sleqp_malloc(star));

  MA27Data* ma27_data = *star;

  *ma27_data = (MA27Data){0};

  MA27ControlInfo* control_info = &(ma27_data->control_info);

  ma27id_(control_info->icntl_, control_info->cntl_);

  control_info->cntl.fratio = ma27_fratio;

  if(ma27_verbose)
  {
    control_info->icntl.print_level = 2;
  }

  MA27_CHECK_ERROR(control_info->info.iflag);

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_sparse_factorization_ma27_create(SleqpSparseFactorization** star,
                                                     SleqpParams* params)
{
  SleqpSparseFactorizationCallbacks callbacks = {
    .set_matrix = ma27_data_set_matrix,
    .solve = ma27_data_solve,
    .get_sol = ma27_data_get_sol,
    .get_condition_estimate = ma27_data_get_condition_estimate,
    .free = ma27_data_free
  };

  MA27Data* ma27_data;

  SLEQP_CALL(ma27_data_create(&ma27_data));

  SLEQP_CALL(sleqp_sparse_factorization_create(star,
                                               params,
                                               &callbacks,
                                               (void*) ma27_data));

  return SLEQP_OKAY;
}
