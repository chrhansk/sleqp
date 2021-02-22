#ifndef HSL_MATRIX_H
#define HSL_MATRIX_H

#include <stdint.h>

#include "sleqp_types.h"
#include "sparse/sleqp_sparse_matrix.h"

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct HSLMatrix
  {
    double* data;
    int32_t* rows;
    int32_t* cols;

    int32_t max_nnz;
    int32_t nnz;
    int32_t dim;

  } HSLMatrix;

  SLEQP_NODISCARD
  SLEQP_RETCODE hsl_matrix_set(HSLMatrix* hsl_matrix,
                               SleqpSparseMatrix* matrix);

  SLEQP_NODISCARD
  SLEQP_RETCODE hsl_matrix_clear(HSLMatrix* hsl_matrix);

#ifdef __cplusplus
}
#endif

#endif /* HSL_MATRIX_H */
