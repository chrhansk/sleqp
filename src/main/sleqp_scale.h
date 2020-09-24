#ifndef SLEQP_SCALE_H
#define SLEQP_SCALE_H

/**
 * @file sleqp_scale.h
 * @brief Definition of the problem scaling.
 **/

#include "sleqp_problem.h"
#include "sleqp_iterate.h"

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * @defgroup scaling NLP scaling
   *
   * The scaled version of the basic NLP is given by
   *
   * \f[
   * \begin{aligned}
   * \min \: & f'(x')                         \\
   * \text{s.t. } \: & l' \leq c'(x') \leq u' \\
   * & l'_x \leq x' \leq u'_x,
   * \end{aligned}
   * \f]
   *
   * defined by scaling weights \f$ \lambda \in \mathbb{R} \f$,
   * \f$ \alpha \in \mathbb{N}^{m} \f$,
   * \f$ \beta \in \mathbb{N}^{n} \f$.
   *
   * The weights \f$ \alpha, \beta \f$ yield scaling
   * factors \f$ a, b \f$:
   *
   * \f[
   * \begin{aligned}
   * a_i := 2^{\alpha_i}      \\
   * b_i := 2^{\beta_i}       \\
   * \end{aligned}
   * \f]
   *
   * The constraint scales \f$ a \f$ and variable scales
   * \f$ b \f$ yield matrices
   * \f$ A = \operatorname{diag}(a) \f$,
   * \f$ B = \operatorname{diag}(b) \f$.
   *
   * \f[
   * \begin{aligned}
   * f'(\cdot) &:= \lambda f(B \cdot)      \\
   * c'(\cdot) &:= A c(B \cdot)            \\
   * l' &:= A l                            \\
   * u' &:= A u                            \\
   * l_x' &:= B l_x                        \\
   * x'   &:= B x                          \\
   * u_x' &:= B u_x                        \\
   * \end{aligned}
   * \f]
   *
   * Note that the scaling is exact (apart from over- / underflows)
   * in the sense that the unscaling is inverse to the scaling
   * even on floating points.
   *
   * @see Functions
   *
   *
   * @{
   **/

  typedef struct SleqpScalingData SleqpScalingData;

  SLEQP_RETCODE sleqp_scaling_create(SleqpScalingData** scaling,
                                     int num_variables,
                                     int num_constraints);

  int sleqp_scaling_get_num_variables(SleqpScalingData* scaling);
  int sleqp_scaling_get_num_constraints(SleqpScalingData* scaling);

  int sleqp_scaling_get_func_weight(SleqpScalingData* scaling);

  SLEQP_RETCODE sleqp_scaling_set_func_weight(SleqpScalingData* scaling,
                                              int weight);

  SLEQP_RETCODE sleqp_scaling_set_func_weight_from_nominal(SleqpScalingData* scaling,
                                                           double nominal_value);

  SLEQP_RETCODE sleqp_scaling_set_var_weight(SleqpScalingData* scaling,
                                             int index,
                                             int weight);

  /**
   * Sets variable scaling weights in order for the scaling of the primal values
   * to map all of the given nominal values to [.5, 1.)
   **/
  SLEQP_RETCODE sleqp_scaling_set_var_weights_from_nominal(SleqpScalingData* scaling,
                                                           double* nominal_values);

  SLEQP_RETCODE sleqp_scaling_set_var_weight_from_nominal(SleqpScalingData* scaling,
                                                          int index,
                                                          double nominal_value);

  SLEQP_RETCODE sleqp_scaling_set_cons_weight(SleqpScalingData* scaling,
                                              int index,
                                              int weight);

  SLEQP_RETCODE sleqp_scaling_set_cons_weights_from_nominal(SleqpScalingData* scaling,
                                                            double* nominal_values);

  SLEQP_RETCODE sleqp_scaling_set_cons_weight_from_nominal(SleqpScalingData* scaling,
                                                           int index,
                                                           double nominal_value);

  int* sleqp_scaling_get_var_weights(SleqpScalingData* scaling);

  int* sleqp_scaling_get_cons_weights(SleqpScalingData* scaling);

  SLEQP_RETCODE sleqp_scaling_set_func(SleqpScalingData* scaling,
                                       SleqpFunc* func);

  /** @name Scaling
   *  Functions to perform scaling.
   */
  /**@{*/

  SLEQP_RETCODE sleqp_scale_point(SleqpScalingData* scaling,
                                  SleqpSparseVec* point);

  double sleqp_scale_func_val(SleqpScalingData* scaling,
                              double func_val);

  SLEQP_RETCODE sleqp_scale_func_grad(SleqpScalingData* scaling,
                                      SleqpSparseVec* func_grad);

  SLEQP_RETCODE sleqp_scale_cons_val(SleqpScalingData* scaling,
                                     SleqpSparseVec* cons_val);

  SLEQP_RETCODE sleqp_scale_cons_jac(SleqpScalingData* scaling,
                                     SleqpSparseMatrix* cons_jac);

  SLEQP_RETCODE sleqp_scale_cons_duals(SleqpScalingData* scaling,
                                       SleqpSparseVec* cons_duals);

  SLEQP_RETCODE sleqp_scale_var_duals(SleqpScalingData* scaling,
                                      SleqpSparseVec* var_duals);

  SLEQP_RETCODE sleqp_scale_hessian_product(SleqpScalingData* scaling,
                                            SleqpSparseVec* product);

  SLEQP_RETCODE sleqp_scale_iterate(SleqpScalingData* scaling,
                                    SleqpIterate* iterate);

  /**@}*/

  /** @name Unscaling
   *  Functions to perform unscaling.
   */
  /**@{*/

  SLEQP_RETCODE sleqp_unscale_point(SleqpScalingData* scaling,
                                    SleqpSparseVec* scaled_point);

  double sleqp_unscale_func_val(SleqpScalingData* scaling,
                                double unscaled_func_val);

  SLEQP_RETCODE sleqp_unscale_func_grad(SleqpScalingData* scaling,
                                        SleqpSparseVec* scaled_func_grad);

  SLEQP_RETCODE sleqp_unscale_cons_val(SleqpScalingData* scaling,
                                       SleqpSparseVec* scaled_cons_val);

  SLEQP_RETCODE sleqp_unscale_cons_jac(SleqpScalingData* scaling,
                                       SleqpSparseMatrix* scaled_cons_jac);

  SLEQP_RETCODE sleqp_unscale_cons_duals(SleqpScalingData* scaling,
                                         SleqpSparseVec* scaled_cons_duals);

  SLEQP_RETCODE sleqp_unscale_var_duals(SleqpScalingData* scaling,
                                        SleqpSparseVec* scaled_var_duals);

  SLEQP_RETCODE sleqp_unscale_hessian_direction(SleqpScalingData* scaling,
                                                SleqpSparseVec* direction,
                                                SleqpSparseVec* cons_duals);

  SLEQP_RETCODE sleqp_unscale_iterate(SleqpScalingData* scaling,
                                      SleqpIterate* scaled_iterate);

  /**@}*/

  SLEQP_RETCODE sleqp_func_scaling_from_gradient(SleqpScalingData* scaling,
                                                 SleqpSparseVec* gradient,
                                                 double eps);

  SLEQP_RETCODE sleqp_scaling_from_cons_jac(SleqpScalingData* scaling,
                                            SleqpSparseMatrix* cons_jac,
                                            double eps);

  SLEQP_RETCODE sleqp_scaling_capture(SleqpScalingData* scaling);

  SLEQP_RETCODE sleqp_scaling_release(SleqpScalingData** scaling);

  /**
   * @}
   **/

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_SCALE_H */
