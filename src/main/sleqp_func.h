#ifndef SLEQP_FUNC_H
#define SLEQP_FUNC_H

#include "sleqp_types.h"
#include "sparse/sleqp_sparse.h"

#ifdef __cplusplus
extern "C" {
#endif

  // TODO: add SLEQPProblem* at some point

  //typedef struct SleqpFunc SleqpFunc;

  /**
   * Sets the current iterate value.
   **/
  typedef SLEQP_RETCODE (*SLEQP_FUNC_SET)(SleqpSparseVec* x,
                                          size_t num_variables,
                                          void* func_data);

  // evaluate functions (obj + cons)
  /**
   * Evaluates the function values / gradients
   **/
  typedef SLEQP_RETCODE (*SLEQP_FUNC_EVAL)(size_t num_variables,
                                           int* indices,
                                           SleqpSparseVec* fvals,
                                           SleqpSparseMatrix* grad,
                                           void* func_data);

  // lambda_0 * d^{T} * H(f) * d - d^{T} * (lambda^{T} * H(g)) * d
  typedef SLEQP_RETCODE (*SLEQP_HESS_EVAL_BILINEAR)(size_t num_variables,
                                                    double* fval,
                                                    SleqpSparseVec* direction,
                                                    SleqpSparseVec* multipliers,
                                                    void* func_data);

  // objective:



  /*
  typedef SLEQP_RETCODE (*SLEQP_JAC_EVAL)(double* x,
                                          size_t num_variables,
                                          SleqpSparseMatrix* jacobi,
                                           void* func_data);

  SLEQP_RETCODE sleqp_func_create(SleqpFunc** fstar,
                                  SLEQP_FUNC_EVAL func_eval,
                                  SLEQP_GRAD_EVAL grad_eval,
                                  void* func_data);

  SLEQP_RETCODE sleqp_func_eval(SleqpFunc* func,
                                double* x,
                                size_t num_variables,
                                double* v);

  SLEQP_RETCODE sleqp_grad_eval(SleqpFunc* func,
                                double* x,
                                size_t num_variables,
                                double* g);

  SLEQP_RETCODE sleqp_func_free(SleqpFunc** fstar);
  */

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_FUNC_H */
