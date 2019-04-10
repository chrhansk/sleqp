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
   * & l'_x \leq x' \leq u_x,
   * \end{aligned}
   * \f]
   *
   * defined by scaling factors \f$ \lambda \in \mathbb{R} \f$,
   * \f$ a \in \mathbb{R}^{m} \f$,
   * \f$ b \in \mathbb{R}^{n} \f$.
   *
   * The constraint scales \f$ a \f$ and variable scales
   * \f$ b \f$ yield matrices
   * \f$ A = \operatorname{diag}(a) \f$,
   * \f$ B = \operatorname{diag}(b) \f$.
   *
   * \f[
   * \begin{aligned}
   * f'(\cdot) &:= \lambda f(B^{-1} \cdot) \\
   * c'(\cdot) &:= A f(B^{-1} \cdot)       \\
   * l' &:= A l                            \\
   * u' &:= A u                            \\
   * l_x' &:= B l_x                        \\
   * u_x' &:= B u_x                        \\
   * \end{aligned}
   * \f]
   *
   *
   * @see Functions
   *
   *
   * @{
   **/

  typedef struct SleqpScalingData SleqpScalingData;

  SLEQP_RETCODE sleqp_scaling_create(SleqpScalingData** scaling,
                                     SleqpProblem* problem,
                                     SleqpParams* params);

  SleqpProblem* sleqp_scaling_get_scaled_problem(SleqpScalingData* scaling);

  SLEQP_RETCODE sleqp_scaling_set_func_scale(SleqpScalingData* scaling,
                                             int scale);

  SLEQP_RETCODE sleqp_scaling_set_var_scale(SleqpScalingData* scaling,
                                            int index,
                                            int scale);

  SLEQP_RETCODE sleqp_scaling_set_cons_scale(SleqpScalingData* scaling,
                                             int index,
                                             int scale);

  SLEQP_RETCODE sleqp_scaling_flush(SleqpScalingData* scaling);

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

  SLEQP_RETCODE sleqp_scale_iterate(SleqpScalingData* scaling,
                                    SleqpIterate* iterate);

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

  SLEQP_RETCODE sleqp_unscale_iterate(SleqpScalingData* scaling,
                                      SleqpIterate* scaled_iterate);

  SLEQP_RETCODE sleqp_scaling_free(SleqpScalingData** scaling);

  /**
   * @}
   **/

#ifdef __cplusplus
}
#endif

#endif /* SLEQP_SCALE_H */
