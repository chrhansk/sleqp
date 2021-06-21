#ifndef SLEQP_SCALE_H
#define SLEQP_SCALE_H

/**
 * @file sleqp_scale.h
 * @brief Definition of the problem scaling.
 **/

#include "sleqp_export.h"
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

  typedef struct SleqpScaling SleqpScaling;

  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_scaling_create(SleqpScaling** scaling,
                                     int num_variables,
                                     int num_constraints);

  SLEQP_RETCODE sleqp_scaling_reset(SleqpScaling* scaling);

  SLEQP_EXPORT int sleqp_scaling_get_num_variables(SleqpScaling* scaling);
  SLEQP_EXPORT int sleqp_scaling_get_num_constraints(SleqpScaling* scaling);

  SLEQP_EXPORT int sleqp_scaling_get_func_weight(SleqpScaling* scaling);

  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_scaling_set_func_weight(SleqpScaling* scaling,
                                              int weight);

  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_scaling_set_func_weight_from_nominal(SleqpScaling* scaling,
                                                           double nominal_value);

  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_scaling_set_var_weight(SleqpScaling* scaling,
                                             int index,
                                             int weight);

  /**
   * Sets variable scaling weights in order for the scaling of the primal values
   * to map all of the given nominal values to [.5, 1.)
   **/
  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_scaling_set_var_weights_from_nominal(SleqpScaling* scaling,
                                                           double* nominal_values);

  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_scaling_set_var_weight_from_nominal(SleqpScaling* scaling,
                                                          int index,
                                                          double nominal_value);

  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_scaling_set_cons_weight(SleqpScaling* scaling,
                                              int index,
                                              int weight);

  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_scaling_set_cons_weights_from_nominal(SleqpScaling* scaling,
                                                            double* nominal_values);

  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_scaling_set_cons_weight_from_nominal(SleqpScaling* scaling,
                                                           int index,
                                                           double nominal_value);

  SLEQP_EXPORT int* sleqp_scaling_get_var_weights(SleqpScaling* scaling);

  SLEQP_EXPORT int* sleqp_scaling_get_cons_weights(SleqpScaling* scaling);

  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_scaling_set_func(SleqpScaling* scaling,
                                       SleqpFunc* func);

  /** @name Scaling
   *  Functions to perform scaling.
   */
  /**@{*/

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

  /**@}*/

  /** @name Unscaling
   *  Functions to perform unscaling.
   */
  /**@{*/

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

  /**@}*/

  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_func_scaling_from_gradient(SleqpScaling* scaling,
                                                 SleqpSparseVec* gradient,
                                                 double eps);

  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_scaling_from_cons_jac(SleqpScaling* scaling,
                                            SleqpSparseMatrix* cons_jac,
                                            double eps);

  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_scaling_capture(SleqpScaling* scaling);

  SLEQP_EXPORT SLEQP_NODISCARD
  SLEQP_RETCODE sleqp_scaling_release(SleqpScaling** star);

  /**
   * @}
   **/

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_SCALE_H */
