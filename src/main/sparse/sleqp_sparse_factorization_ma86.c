#include "sleqp_sparse_factorization_ma86.h"

#include <assert.h>

#include <hsl_ma86d.h>
#include <hsl_mc68i.h>

#include "sleqp_mem.h"

static const bool ma86_verbose = false;

typedef enum {

  MA86_ERROR_ALLOCATION                = -1,
  MA86_ERROR_INVALID_ELIMINATION_ORDER = -2,
  MA86_ERROR_SINGULAR_MATRIX           = -3,
  MA86_ERROR_SIZE                      = -4,
  MA86_ERROR_INVALID_DATA              = -5,
  MA86_ERROR_JOB_OUT_OF_RANGE          = -6,
  MA86_ERROR_INVALID_SMALL_STATIC      = -7,
  MA86_ERROR_INVALID_MATRIX_TYPE       = -8,
  MA86_ERROR_INVALID_SCALE             = -9,

  MA86_SUCCESS                         = 0,

  MA86_WARNING_INSUFFICIENT_POOL_SIZE  = 1,
  MA86_WARNING_SINGULAR_MATRIX         = 2,

  MA86_WARNING = MA86_WARNING_INSUFFICIENT_POOL_SIZE | MA86_WARNING_SINGULAR_MATRIX
} MA866_STATUS;

typedef enum {
  MC68_ORDER_APX_MINDEG = 1,
  MC68_ORDER_MINDEG     = 2,
  MC68_ORDER_METIS      = 3,
  MC68_ORDER_MA47       = 4,
} MC68_ORDER;

static SLEQP_RETCODE ma86_get_error_string(int value, const char** message)
{

  switch(value)
  {
  case MA86_ERROR_ALLOCATION:
    *message = "Allocation error";
    break;
  case MA86_ERROR_INVALID_ELIMINATION_ORDER:
    *message = "Error in user-specified elimination order";
    break;
  case MA86_ERROR_SINGULAR_MATRIX:
    *message = "Matrix is singular";
    break;
  case MA86_ERROR_SIZE:
    *message = "Error in size of array";
    break;
  case MA86_ERROR_INVALID_DATA:
    *message = "Infinities found in matrix data";
    break;
  case MA86_ERROR_JOB_OUT_OF_RANGE:
    *message = "Job is out of range";
    break;
  case MA86_ERROR_INVALID_SMALL_STATIC:
    *message = "Invalid control.static_ size";
    break;
  case MA86_ERROR_INVALID_MATRIX_TYPE:
    *message = "Invalid matrix type";
    break;
  case MA86_ERROR_INVALID_SCALE:
    *message = "Invalid scale";
    break;
  case MA86_SUCCESS:
    *message = "Success";
    break;
  case MA86_WARNING_INSUFFICIENT_POOL_SIZE:
    *message = "Insufficient pool size";
    break;
  case MA86_WARNING_SINGULAR_MATRIX:
    *message = "Matrix is singular";
    break;
  default:
    *message = "(unknown)";
    break;
  }

  return SLEQP_OKAY;
}

#define MA86_IS_ERROR(value) (value < 0)

#define MA86_CHECK_ERROR(x)                                              \
  do                                                                     \
  {                                                                      \
    int ma86_status = (x);                                               \
                                                                         \
    if(ma86_status == MA86_SUCCESS)                                      \
    {                                                                    \
      break;                                                             \
    }                                                                    \
                                                                         \
    const char* ma86_error_string;                                       \
    SLEQP_CALL(ma86_get_error_string(ma86_status,                        \
                                     &ma86_error_string));               \
                                                                         \
    if(MA86_IS_ERROR(ma86_status))                                       \
    {                                                                    \
                                                                         \
      sleqp_log_error("Caught hsl_ma86 error <%d> (%s) in function %s",  \
                      ma86_status,                                       \
                      ma86_error_string,                                 \
                      __func__);                                         \
                                                                         \
      if(ma86_status == MA86_ERROR_ALLOCATION)                           \
      {                                                                  \
        return SLEQP_NOMEM;                                              \
      }                                                                  \
                                                                         \
      return SLEQP_INTERNAL_ERROR;                                       \
    }                                                                    \
    else                                                                 \
    {                                                                    \
      sleqp_log_warn("Caught hsl_ma86 warning <%d> (%s) in function %s", \
                     ma86_status,                                        \
                      ma86_error_string,                                 \
                     __func__);                                          \
                                                                         \
    }                                                                    \
  } while(0)

typedef struct MA86Data
{
  struct ma86_control control;
  struct ma86_info info;

  struct mc68_control control_c;
  struct mc68_info info_c;

  int dim;
  int max_dim;
  int* order;

  void* keep;
  SleqpSparseMatrix* matrix;
  double* rhs_sol;

} MA86Data;



static SLEQP_RETCODE ma86_data_create(MA86Data** star)
{
  SLEQP_CALL(sleqp_malloc(star));

  MA86Data* ma86_data = *star;

  *ma86_data = (MA86Data){0};

  ma86_default_control(&(ma86_data->control));

  // We expect all matrices here to be non-singular,
  // error out otherwise
  ma86_data->control.action = 0;

  if(ma86_verbose)
  {
    ma86_data->control.diagnostics_level = 2;
  }
  else
  {
    ma86_data->control.diagnostics_level = 0;
  }

  mc68_default_control(&(ma86_data->control_c));

  SLEQP_CALL(sleqp_sparse_matrix_create(&(ma86_data->matrix), 1, 1, 0));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE ma86_data_set_matrix(void* factorization_data,
                                          SleqpSparseMatrix* matrix)
{
  MA86Data* ma86_data = (MA86Data*) factorization_data;

  const int num_cols = sleqp_sparse_matrix_get_num_cols(matrix);
  const int num_rows = sleqp_sparse_matrix_get_num_rows(matrix);

  SLEQP_CALL(sleqp_sparse_matrix_resize(ma86_data->matrix, num_rows, num_cols));

  SLEQP_CALL(sleqp_sparse_lower_triangular(matrix, ma86_data->matrix));

  assert(num_cols == num_rows);

  const int dim = num_cols;

  if(ma86_data->max_dim < dim)
  {
    SLEQP_CALL(sleqp_realloc(&(ma86_data->order), dim));

    SLEQP_CALL(sleqp_realloc(&(ma86_data->rhs_sol), dim));

    ma86_data->max_dim = dim;
  }

  ma86_data->dim = dim;

  int* cols = sleqp_sparse_matrix_get_cols(ma86_data->matrix);
  int* rows = sleqp_sparse_matrix_get_rows(ma86_data->matrix);
  double* data = sleqp_sparse_matrix_get_data(ma86_data->matrix);

  mc68_order(MC68_ORDER_APX_MINDEG,
             dim,
             cols,
             rows,
             ma86_data->order,
             &(ma86_data->control_c),
             &(ma86_data->info_c));

  MA86_CHECK_ERROR(ma86_data->info_c.stat);

  // order

  ma86_analyse(dim,
               cols,
               rows,
               ma86_data->order,
               &(ma86_data->keep),
               &(ma86_data->control),
               &(ma86_data->info));

  MA86_CHECK_ERROR(ma86_data->info.stat);

  ma86_factor(dim,
              cols,
              rows,
              data,
              ma86_data->order,
              &(ma86_data->keep),
              &(ma86_data->control),
              &(ma86_data->info),
              NULL);

  MA86_CHECK_ERROR(ma86_data->info.stat);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE ma86_data_solve(void* factorization_data,
                                     SleqpSparseVec* rhs)
{
  MA86Data* ma86_data = (MA86Data*) factorization_data;

  const int job = 0;
  const int nrhs = 1;
  const int dim = ma86_data->dim;
  const int ldx = ma86_data->dim;

  assert(rhs->dim == dim);

  SLEQP_CALL(sleqp_sparse_vector_to_raw(rhs, ma86_data->rhs_sol));

  ma86_solve(job,
             nrhs,
             ldx,
             ma86_data->rhs_sol,
             ma86_data->order,
             &(ma86_data->keep),
             &(ma86_data->control),
             &(ma86_data->info),
             NULL);

  MA86_CHECK_ERROR(ma86_data->info.stat);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE ma86_data_get_sol(void* factorization_data,
                                       SleqpSparseVec* sol,
                                       int begin,
                                       int end,
                                       double zero_eps)
{
  MA86Data* ma86_data = (MA86Data*) factorization_data;

  SLEQP_CALL(sleqp_sparse_vector_from_raw(sol,
                                          ma86_data->rhs_sol + begin,
                                          end - begin,
                                          zero_eps));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE ma86_data_get_condition_estimate(void* factorization_data,
                                                      double* condition_estimate)
{
  MA86Data* ma86_data = (MA86Data*) factorization_data;

  (*condition_estimate) = SLEQP_NONE;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE ma86_data_free(void **star)
{
  MA86Data* ma86_data = (MA86Data*) (*star);

  ma86_finalise(&(ma86_data->keep), &(ma86_data->control));

  sleqp_free(&(ma86_data->rhs_sol));
  sleqp_free(&(ma86_data->order));

  SLEQP_CALL(sleqp_sparse_matrix_release(&(ma86_data->matrix)));

  sleqp_free(&ma86_data);

  *star = NULL;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_sparse_factorization_ma86_create(SleqpSparseFactorization** star,
                                                     SleqpParams* params)
{
  SleqpSparseFactorizationCallbacks callbacks = {
    .set_matrix = ma86_data_set_matrix,
    .solve = ma86_data_solve,
    .get_sol = ma86_data_get_sol,
    .get_condition_estimate = ma86_data_get_condition_estimate,
    .free = ma86_data_free
  };

  MA86Data* ma86_data;

  SLEQP_CALL(ma86_data_create(&ma86_data));

  SLEQP_CALL(sleqp_sparse_factorization_create(star,
                                               params,
                                               &callbacks,
                                               (void*) ma86_data));

  return SLEQP_OKAY;
}