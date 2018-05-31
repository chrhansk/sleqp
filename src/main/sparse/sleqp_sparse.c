#include "sleqp_sparse.h"

#include <assert.h>

#include "sleqp_mem.h"

#define MAX(a,b) (((a)>(b))?(a):(b))

SLEQP_RETCODE sleqp_sparse_vector_create(SleqpSparseVec** vstar,
                                         size_t dim,
                                         size_t nnz_max)
{
  assert(nnz_max <= dim);

  sleqp_malloc(vstar);

  SleqpSparseVec *vec = *vstar;

  vec->nnz = 0;
  vec->dim = dim;
  vec->nnz_max = nnz_max;

  sleqp_calloc(&vec->data, nnz_max);
  sleqp_calloc(&vec->indices, nnz_max);

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_sparse_vector_free(SleqpSparseVec** vstar)
{
  SleqpSparseVec *vec = *vstar;

  sleqp_free(vec->indices);
  sleqp_free(vec->data);

  sleqp_free(vec);

  *vstar = NULL;

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_sparse_matrix_create(SleqpSparseMatrix** mstar,
                                         size_t num_rows,
                                         size_t num_cols,
                                         size_t nnz_max)
{
  sleqp_malloc(mstar);

  SleqpSparseMatrix* matrix = *mstar;

  matrix->nnz = 0;
  matrix->nnz_max = nnz_max;

  matrix->num_cols = num_cols;
  matrix->num_rows = num_rows;

  sleqp_calloc(&matrix->data, nnz_max);
  sleqp_calloc(&matrix->cols, (MAX(nnz_max, num_cols) + 1));
  sleqp_calloc(&matrix->rows, nnz_max);

  return SLEQP_OKAY;
}

SLEQP_RETCODE sleqp_sparse_matrix_free(SleqpSparseMatrix** mstar)
{
  SleqpSparseMatrix* matrix = *mstar;

  sleqp_free(matrix->rows);
  sleqp_free(matrix->cols);
  sleqp_free(matrix->data);

  sleqp_free(matrix);

  *mstar = NULL;

  return SLEQP_OKAY;
}
