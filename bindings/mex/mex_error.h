#ifndef SLEQP_MEX_ERROR_H
#define SLEQP_MEX_ERROR_H

#include <mex.h>

#include "sleqp.h"

#define MEX_EXPECT(array, what, message, ...)                                  \
  do                                                                           \
  {                                                                            \
    if (!(array && (what)(array)))                                             \
    {                                                                          \
      sleqp_raise(SLEQP_FAILED_ASSERTION, message, ##__VA_ARGS__);             \
    }                                                                          \
  } while (false)

#define MEX_EXPECT_DOUBLE(array)                                               \
  MEX_EXPECT(array, mxIsDouble, "Expected type double")

#define MEX_EXPECT_SPARSE(array)                                               \
  MEX_EXPECT(array, mxIsSparse, "Expected sparse")

#define MEX_EXPECT_STRUCT(array)                                               \
  MEX_EXPECT(array, mxIsStruct, "Expected struct")

#define MEX_EXPECT_SCALAR(array)                                               \
  MEX_EXPECT(array, mxIsScalar, "Expected scalar")

#define MEX_EXPECT_LOGICAL_SCALAR(array)                                       \
  MEX_EXPECT(array, mxIsLogicalScalar, "Expected logical scalar")

#define MEX_EXPECT_FUNCTION_HANDLE(array)                                      \
  MEX_EXPECT(array, mxIsFunctionHandle, "Expected function handle")

#define MEX_EXPECT_SHAPE(array, rows, cols)                                    \
  do                                                                           \
  {                                                                            \
    if (mxGetM(array) != (rows))                                               \
    {                                                                          \
      sleqp_raise(SLEQP_FAILED_ASSERTION,                                      \
                  "Invalid number of rows: expected %d, actual %zu",           \
                  (rows),                                                      \
                  mxGetM(array));                                              \
    }                                                                          \
    if (mxGetN(array) != (cols))                                               \
    {                                                                          \
      sleqp_raise(SLEQP_FAILED_ASSERTION,                                      \
                  "Invalid number of cols: expected %d, actual %zu",           \
                  (cols),                                                      \
                  mxGetN(array));                                              \
    }                                                                          \
  } while (false)

#define MEX_EXPECT_NUM_ELEMENTS(array, elements)                               \
  do                                                                           \
  {                                                                            \
    if (mxGetNumberOfElements(array) != (elements))                            \
    {                                                                          \
      sleqp_raise(SLEQP_FAILED_ASSERTION,                                      \
                  "Invalid number of elements: expected %d, actual %zu",       \
                  (elements),                                                  \
                  mxGetNumberOfElements(array));                               \
    }                                                                          \
  } while (false)

#endif /* SLEQP_MEX_ERROR_H */
