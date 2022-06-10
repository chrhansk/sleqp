#include "factorization_ma97.h"

#include <assert.h>

#include <hsl_ma97d.h>
#include <hsl_mc68i.h>

#include "defs.h"
#include "log.h"
#include "mem.h"

static const bool ma97_verbose = true;

static const bool ma97_check_matrix = false;

typedef enum
{

  MA97_ERROR_WRONG_CALL_SEQUENCE           = -1,
  MA97_ERROR_NEGATIVE_N                    = -2,
  MA97_ERROR_INVALID_PTR                   = -3,
  MA97_ERROR_COLUMN_OUT_OF_RANGE           = -4,
  MA97_ERROR_MATRIX_TYPE_OUT_OF_RANGE      = -5,
  MA97_ERROR_MATRIX_NOT_HERMITIAN          = -6,
  MA97_ERROR_MATRIX_SINGULAR               = -7,
  MA97_ERROR_MATRIX_NOT_POS_DEF            = -8,
  MA97_ERROR_MATRIX_INVALID_DATA           = -9,
  MA97_ERROR_MATRIX_MISSING_PTRS           = -10,
  MA97_ERROR_MATRIX_INVALID_ORDERING       = -11,
  MA97_ERROR_MATRIX_INVALID_ARRAY_SIZE     = -12,
  MA97_ERROR_JOB_OUT_OF_RANGE              = -13,
  MA97_ERROR_INVALID_MATRIX_TYPE_POSDEF    = -14,
  MA97_ERROR_INVALID_MATRIX_TYPE_HERMITIAN = -15,
  MA97_ERROR_ALLOCATION                    = -16,
  MA97_ERROR_METIS_NOT_AVAILABLE           = -17,
  MA97_ERROR_MC68                          = -18,
  MA97_ERROR_MC77                          = -19,
  MA97_ERROR_MISSING_VALUES                = -20,

  MA97_SUCCESS = 0,

  MA97_WARN_OUT_OF_RANGE_INDICES                         = 1,
  MA97_WARN_DUPLICATE_INDICES                            = 2,
  MA97_WARN_INVALID_INDICES                              = 3,
  MA97_WARN_MISSING_DIAGONAL_ENTRIES                     = 4,
  MA97_WARN_MISSING_DIAGONAL_ENTRIES_AND_INVALID_INDICES = 5,
  MA97_WARN_STRUCTURALLY_SINGULAR                        = 6,
  MA97_WARN_SINGULAR                                     = 7,
  MA97_WARN_MISSING_SCALING                              = 8,
} MA97_STATUS;

static SLEQP_RETCODE
ma97_get_error_string(int value, const char** message)
{
  switch (value)
  {
  case MA97_ERROR_WRONG_CALL_SEQUENCE:
    *message = "Wrong calling sequence";
    break;
  case MA97_ERROR_NEGATIVE_N:
    *message = "Negative dimension";
    break;
  case MA97_ERROR_INVALID_PTR:
    *message = "Invalid column entries";
    break;
  case MA97_ERROR_COLUMN_OUT_OF_RANGE:
    *message = "All column entries out of range";
    break;
  case MA97_ERROR_MATRIX_TYPE_OUT_OF_RANGE:
    *message = "Matrix type out of range";
    break;
  case MA97_ERROR_MATRIX_NOT_HERMITIAN:
    *message = "Matrix is not Hermitian";
    break;
  case MA97_ERROR_MATRIX_SINGULAR:
    *message = "Matrix is singular";
    break;
  case MA97_ERROR_MATRIX_NOT_POS_DEF:
    *message = "Matrix is not positive definite";
    break;
  case MA97_ERROR_MATRIX_INVALID_DATA:
    *message = "Infinities found in matrix data";
    break;
  case MA97_ERROR_MATRIX_MISSING_PTRS:
    *message = "Missing row / column pointers";
    break;
  case MA97_ERROR_MATRIX_INVALID_ORDERING:
    *message = "Error in user-specified elimination order";
    break;
  case MA97_ERROR_MATRIX_INVALID_ARRAY_SIZE:
    *message = "Invalid array size";
    break;
  case MA97_ERROR_JOB_OUT_OF_RANGE:
    *message = "Job out of range";
    break;
  case MA97_ERROR_INVALID_MATRIX_TYPE_POSDEF:
    *message = "Invalid matrix type";
    break;
  case MA97_ERROR_INVALID_MATRIX_TYPE_HERMITIAN:
    *message = "Invalid matrix type";
    break;
  case MA97_ERROR_ALLOCATION:
    *message = "Allocation error";
    break;
  case MA97_ERROR_METIS_NOT_AVAILABLE:
    *message = "METIS not available";
    break;
  case MA97_ERROR_MC68:
    *message = "Error in MC68";
    break;
  case MA97_ERROR_MC77:
    *message = "Error in MC77";
    break;
  case MA97_ERROR_MISSING_VALUES:
    *message = "Missing values for ordering";
    break;

  case MA97_SUCCESS:
    *message = "Success";
    break;

  case MA97_WARN_OUT_OF_RANGE_INDICES:
    *message = "Matrix indices out of range";
    break;
  case MA97_WARN_DUPLICATE_INDICES:
    *message = "Dupliate matrix indices";
    break;
  case MA97_WARN_INVALID_INDICES:
    *message = "Invalid matrix indices";
    break;
  case MA97_WARN_MISSING_DIAGONAL_ENTRIES:
    *message = "Missing diagonal entries";
    break;
  case MA97_WARN_MISSING_DIAGONAL_ENTRIES_AND_INVALID_INDICES:
    *message = "Missing diagonal entries and invalid indices";
    break;
  case MA97_WARN_STRUCTURALLY_SINGULAR:
    *message = "Matrix structurally singular";
    break;
  case MA97_WARN_SINGULAR:
    *message = "Matrix singular";
    break;
  case MA97_WARN_MISSING_SCALING:
    *message = "Scaling for mathcing-based ordering not found";
    break;
  default:
    *message = "(unknown)";
    break;
  }

  return SLEQP_OKAY;
}

#define MA97_IS_ERROR(value) (value < 0)

#define MA97_CHECK_ERROR(x)                                                    \
  do                                                                           \
  {                                                                            \
    int ma97_status = (x);                                                     \
                                                                               \
    if (ma97_status == MA97_SUCCESS)                                           \
    {                                                                          \
      break;                                                                   \
    }                                                                          \
                                                                               \
    const char* ma97_error_string;                                             \
    SLEQP_CALL(ma97_get_error_string(ma97_status, &ma97_error_string));        \
                                                                               \
    if (MA97_IS_ERROR(ma97_status))                                            \
    {                                                                          \
                                                                               \
      sleqp_raise(SLEQP_INTERNAL_ERROR,                                        \
                  "Caught hsl_ma97 error <%d> (%s) in function %s",            \
                  ma97_status,                                                 \
                  ma97_error_string,                                           \
                  __func__);                                                   \
    }                                                                          \
    else                                                                       \
    {                                                                          \
      sleqp_log_warn("Caught hsl_ma97 warning <%d> (%s) in function %s",       \
                     ma97_status,                                              \
                     ma97_error_string,                                        \
                     __func__);                                                \
    }                                                                          \
  } while (0)

typedef enum
{
  MC68_ORDER_USER_SUPPLIED = 0,
  MC68_ORDER_APX_MINDEG    = 1,
  MC68_ORDER_MINDEG        = 2,
  MC68_ORDER_METIS         = 3,
  MC68_ORDER_MA47          = 4,
} MC68_ORDER;

typedef enum
{
  MA97_REAL_POS_DEF = 3,
  MA97_REAL_INDEF   = 4,
  MA97_HERM_POS_DEF = -3,
  MA97_HERM_INDEF   = -4,

} MA_97_MATRIX_TYPE;

typedef struct MA97Data
{
  struct ma97_control control;
  struct ma97_info info;

  struct mc68_control control_c;
  struct mc68_info info_c;

  int dim;
  int max_dim;
  int* order;

  void* akeep;
  void* fkeep;
  double* rhs_sol;

} MA97Data;

static SLEQP_RETCODE
ma97_data_create(MA97Data** star)
{
  SLEQP_CALL(sleqp_malloc(star));

  MA97Data* ma97_data = *star;

  *ma97_data = (MA97Data){0};

  ma97_default_control(&(ma97_data->control));

  // We expect all matrices here to be non-singular,
  // error out otherwise
  ma97_data->control.action = 0;

  ma97_data->control.ordering = MC68_ORDER_USER_SUPPLIED;

  if (ma97_verbose)
  {
    ma97_data->control.print_level = 2;
  }
  else
  {
    ma97_data->control.print_level = 0;
  }

  mc68_default_control(&(ma97_data->control_c));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
ma97_data_set_matrix(void* factorization_data, SleqpSparseMatrix* matrix)
{
  MA97Data* ma97_data = (MA97Data*)factorization_data;

  const int num_cols = sleqp_sparse_matrix_num_cols(matrix);
  const int num_rows = sleqp_sparse_matrix_num_rows(matrix);

  assert(num_cols == num_rows);

  const int dim = num_cols;

  if (ma97_data->max_dim < dim)
  {
    SLEQP_CALL(sleqp_realloc(&(ma97_data->order), dim));

    SLEQP_CALL(sleqp_realloc(&(ma97_data->rhs_sol), dim));

    ma97_data->max_dim = dim;
  }

  ma97_data->dim = dim;

  int* cols    = sleqp_sparse_matrix_cols(matrix);
  int* rows    = sleqp_sparse_matrix_rows(matrix);
  double* data = sleqp_sparse_matrix_data(matrix);

  mc68_order(MC68_ORDER_APX_MINDEG,
             dim,
             cols,
             rows,
             ma97_data->order,
             &(ma97_data->control_c),
             &(ma97_data->info_c));

  MA97_CHECK_ERROR(ma97_data->info_c.stat);

  // order

  ma97_analyse(ma97_check_matrix,
               dim,
               cols,
               rows,
               data,
               &(ma97_data->akeep),
               &(ma97_data->control),
               &(ma97_data->info),
               ma97_data->order);

  MA97_CHECK_ERROR(ma97_data->info.stat);

  ma97_factor(MA97_REAL_INDEF,
              cols,
              rows,
              data,
              &(ma97_data->akeep),
              &(ma97_data->fkeep),
              &(ma97_data->control),
              &(ma97_data->info),
              NULL); // No scaling

  /*

  */

  MA97_CHECK_ERROR(ma97_data->info.stat);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
ma97_data_solve(void* factorization_data, const SleqpVec* rhs)
{
  MA97Data* ma97_data = (MA97Data*)factorization_data;

  const int nrhs = 1;
  const int dim  = ma97_data->dim;

  assert(rhs->dim == dim);

  SLEQP_CALL(sleqp_vec_to_raw(rhs, ma97_data->rhs_sol));

  ma97_solve(0,
             nrhs,
             ma97_data->rhs_sol,
             dim,
             &(ma97_data->akeep),
             &(ma97_data->fkeep),
             &(ma97_data->control),
             &(ma97_data->info));

  MA97_CHECK_ERROR(ma97_data->info.stat);

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
ma97_data_solution(void* factorization_data,
                   SleqpVec* sol,
                   int begin,
                   int end,
                   double zero_eps)
{
  MA97Data* ma97_data = (MA97Data*)factorization_data;

  SLEQP_CALL(sleqp_vec_set_from_raw(sol,
                                    ma97_data->rhs_sol + begin,
                                    end - begin,
                                    zero_eps));

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
ma97_data_condition_estimate(void* factorization_data,
                             double* condition_estimate)
{
  // MA97Data* ma97_data = (MA97Data*) factorization_data;

  (*condition_estimate) = SLEQP_NONE;

  return SLEQP_OKAY;
}

static SLEQP_RETCODE
ma97_data_free(void** star)
{
  MA97Data* ma97_data = (MA97Data*)(*star);

  // ma97_finalise(&(ma97_data->keep), &(ma97_data->control));

  sleqp_free(&(ma97_data->rhs_sol));
  sleqp_free(&(ma97_data->order));

  sleqp_free(&ma97_data);

  *star = NULL;

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_factorization_ma97_create(SleqpFactorization** star, SleqpParams* params)
{
  SleqpFactorizationCallbacks callbacks
    = {.set_matrix         = ma97_data_set_matrix,
       .solve              = ma97_data_solve,
       .solution           = ma97_data_solution,
       .condition_estimate = ma97_data_condition_estimate,
       .free               = ma97_data_free};

  MA97Data* ma97_data;

  SLEQP_CALL(ma97_data_create(&ma97_data));

  SLEQP_CALL(sleqp_factorization_create(star,
                                        SLEQP_FACT_MA97_NAME,
                                        SLEQP_FACT_MA97_VERSION,
                                        params,
                                        &callbacks,
                                        SLEQP_FACTORIZATION_LOWER,
                                        (void*)ma97_data));

  return SLEQP_OKAY;
}

SLEQP_RETCODE
sleqp_factorization_create_default(SleqpFactorization** star,
                                   SleqpParams* params)
{
  SLEQP_CALL(sleqp_factorization_ma97_create(star, params));

  return SLEQP_OKAY;
}
