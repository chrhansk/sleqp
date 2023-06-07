#ifndef SLEQP_MEX_FUNC_COMMON_H
#define SLEQP_MEX_FUNC_COMMON_H

#include <assert.h>
#include <mex.h>
#include <threads.h>

#include "sleqp.h"

#define MEX_MSG_BUF_SIZE 512

#define MEX_CALL_SIMPLE(x)                                                     \
  do                                                                           \
  {                                                                            \
    mxArray* exception = (x);                                                  \
                                                                               \
    if (exception)                                                             \
    {                                                                          \
      sleqp_raise(SLEQP_FUNC_EVAL_ERROR, "Unknown exception in Matlab call");  \
    }                                                                          \
                                                                               \
  } while (0)

#define MEX_CALL(x)                                                            \
  do                                                                           \
  {                                                                            \
    mxArray* exception = (x);                                                  \
                                                                               \
    if (exception)                                                             \
    {                                                                          \
      char msg_buf[MEX_MSG_BUF_SIZE];                                          \
      mxArray* lhs;                                                            \
      MEX_CALL_SIMPLE(                                                         \
        mexCallMATLABWithTrap(1, &lhs, 1, &exception, MATLAB_FUNC_DISP));      \
      assert(mxIsChar(lhs));                                                   \
      mxGetString(lhs, msg_buf, MEX_MSG_BUF_SIZE);                             \
                                                                               \
      sleqp_raise(SLEQP_FUNC_EVAL_ERROR,                                       \
                  "Exception '%s' in Matlab call",                             \
                  msg_buf);                                                    \
    }                                                                          \
  } while (0)

#define MEX_ARRAY_LEN(array) sizeof(array) / sizeof(array[0])

#define MEX_EVAL(lhs, rhs)                                                     \
  do                                                                           \
  {                                                                            \
    const int nlhs = MEX_ARRAY_LEN((lhs));                                     \
    const int nrhs = MEX_ARRAY_LEN((rhs));                                     \
                                                                               \
    MEX_CALL(mexCallMATLABWithTrap(nlhs,                                       \
                                   (mxArray**)(lhs),                           \
                                   nrhs,                                       \
                                   (mxArray**)(rhs),                           \
                                   MATLAB_FUNC_FEVAL));                        \
  } while (false)

SLEQP_RETCODE
mex_callback_has_field(const mxArray* mex_callbacks,
                       const char* name,
                       bool* has_field);

SLEQP_RETCODE
mex_callback_from_struct(const mxArray* mex_callbacks,
                         const char* name,
                         mxArray** star);

SLEQP_RETCODE
mex_array_to_real(mxArray* array, double* value);

SLEQP_RETCODE
mex_array_to_vec(const mxArray* array, SleqpSettings* settings, SleqpVec* vec);

SLEQP_RETCODE
mex_eval_into_real(int nrhs, mxArray** rhs, double* value);

SLEQP_RETCODE
mex_eval_into_bool(int nrhs, mxArray** rhs, bool* value);

SLEQP_RETCODE
mex_eval_into_vec(int nrhs,
                  mxArray** rhs,
                  SleqpSettings* settings,
                  SleqpVec* vec);

SLEQP_RETCODE
mex_eval_into_sparse_matrix(int nrhs,
                            mxArray** rhs,
                            SleqpSettings* settings,
                            SleqpMat* matrix);

#define MEX_EVAL_INTO_REAL(array, value)                                       \
  SLEQP_CALL(mex_eval_into_real(MEX_ARRAY_LEN(array), array, value))

#define MEX_EVAL_INTO_BOOL(array, value)                                       \
  SLEQP_CALL(mex_eval_into_bool(MEX_ARRAY_LEN(array), array, value))

#define MEX_EVAL_INTO_VEC(array, settings, vec)                                \
  SLEQP_CALL(mex_eval_into_vec(MEX_ARRAY_LEN(array), array, settings, vec))

#define MEX_EVAL_INTO_SPARSE_MATRIX(array, settings, matrix)                   \
  SLEQP_CALL(mex_eval_into_sparse_matrix(MEX_ARRAY_LEN(array),                 \
                                         array,                                \
                                         settings,                             \
                                         matrix))

#endif /* SLEQP_MEX_FUNC_COMMON_H */
