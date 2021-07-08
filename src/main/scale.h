#ifndef SLEQP_SCALE_H
#define SLEQP_SCALE_H

#include "pub_scale.h"

#ifdef __cplusplus
extern "C" {
#endif

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_scaling_set_func(SleqpScaling* scaling,
                                       SleqpFunc* func);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_scale_point(SleqpScaling* scaling,
                                  SleqpSparseVec* point);

  double sleqp_scale_func_val(SleqpScaling* scaling,
                              double func_val);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_scale_func_grad(SleqpScaling* scaling,
                                      SleqpSparseVec* func_grad);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_scale_cons_val(SleqpScaling* scaling,
                                     SleqpSparseVec* cons_val);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_scale_cons_general(SleqpScaling* scaling,
                                         SleqpSparseVec* general_val);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_scale_cons_linear(SleqpScaling* scaling,
                                        SleqpSparseVec* linear_val);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_scale_cons_jac(SleqpScaling* scaling,
                                     SleqpSparseMatrix* cons_jac);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_scale_linear_coeffs(SleqpScaling* scaling,
                                          SleqpSparseMatrix* linear_coeffs);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_scale_cons_duals(SleqpScaling* scaling,
                                       SleqpSparseVec* cons_duals);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_scale_var_duals(SleqpScaling* scaling,
                                      SleqpSparseVec* var_duals);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_scale_hessian_product(SleqpScaling* scaling,
                                            SleqpSparseVec* product);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_scale_iterate(SleqpScaling* scaling,
                                    SleqpIterate* iterate);


  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_unscale_point(SleqpScaling* scaling,
                                    SleqpSparseVec* scaled_point);

  double sleqp_unscale_func_val(SleqpScaling* scaling,
                                double unscaled_func_val);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_unscale_func_grad(SleqpScaling* scaling,
                                        SleqpSparseVec* scaled_func_grad);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_unscale_cons_val(SleqpScaling* scaling,
                                       SleqpSparseVec* scaled_cons_val);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_unscale_cons_jac(SleqpScaling* scaling,
                                       SleqpSparseMatrix* scaled_cons_jac);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_unscale_cons_duals(SleqpScaling* scaling,
                                         SleqpSparseVec* scaled_cons_duals);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_unscale_var_duals(SleqpScaling* scaling,
                                        SleqpSparseVec* scaled_var_duals);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_unscale_hessian_direction(SleqpScaling* scaling,
                                                SleqpSparseVec* direction,
                                                SleqpSparseVec* cons_duals);

  SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_unscale_iterate(SleqpScaling* scaling,
                                      SleqpIterate* scaled_iterate);


#ifdef __cplusplus
}
#endif

#endif /* SLEQP_SCALE_H */
