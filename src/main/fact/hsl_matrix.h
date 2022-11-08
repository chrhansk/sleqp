#ifndef HSL_MATRIX_H
#define HSL_MATRIX_H

#include <stdint.h>

#include "sparse/mat.h"
#include "types.h"

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
SLEQP_RETCODE
hsl_matrix_set(HSLMatrix* hsl_matrix, SleqpMat* matrix);

SLEQP_NODISCARD
SLEQP_RETCODE
hsl_matrix_clear(HSLMatrix* hsl_matrix);

#endif /* HSL_MATRIX_H */
