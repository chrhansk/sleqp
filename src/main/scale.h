#ifndef SLEQP_SCALE_H
#define SLEQP_SCALE_H

#include "pub_scale.h"

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_scaling_set_func(SleqpScaling* scaling, SleqpFunc* func);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_scale_point(SleqpScaling* scaling, SleqpVec* point);

double
sleqp_scale_obj_val(SleqpScaling* scaling, double obj_val);

double
sleqp_scale_lsq_obj_val(SleqpScaling* scaling, double obj_val);

SLEQP_RETCODE
sleqp_scale_lsq_residuals(SleqpScaling* scaling, SleqpVec* lsq_residuals);

SLEQP_RETCODE
sleqp_scale_lsq_forward_direction(SleqpScaling* scaling,
                                  SleqpVec* forward_direction);

SLEQP_RETCODE
sleqp_scale_lsq_adjoint_direction(SleqpScaling* scaling,
                                  SleqpVec* adjoint_direction);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_scale_obj_grad(SleqpScaling* scaling, SleqpVec* obj_grad);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_scale_cons_val(SleqpScaling* scaling, SleqpVec* cons_val);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_scale_cons_general(SleqpScaling* scaling, SleqpVec* general_val);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_scale_cons_linear(SleqpScaling* scaling, SleqpVec* linear_val);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_scale_cons_jac(SleqpScaling* scaling, SleqpSparseMatrix* cons_jac);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_scale_linear_coeffs(SleqpScaling* scaling,
                          SleqpSparseMatrix* linear_coeffs);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_scale_cons_duals(SleqpScaling* scaling, SleqpVec* cons_duals);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_scale_var_duals(SleqpScaling* scaling, SleqpVec* var_duals);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_scale_hessian_product(SleqpScaling* scaling, SleqpVec* product);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_scale_iterate(SleqpScaling* scaling, SleqpIterate* iterate, bool lsq);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_unscale_point(SleqpScaling* scaling, SleqpVec* scaled_point);

double
sleqp_unscale_obj_val(SleqpScaling* scaling, double unscaled_obj_val);

double
sleqp_unscale_lsq_obj_val(SleqpScaling* scaling, double obj_val);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_unscale_obj_grad(SleqpScaling* scaling, SleqpVec* scaled_obj_grad);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_unscale_cons_val(SleqpScaling* scaling, SleqpVec* scaled_cons_val);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_unscale_cons_jac(SleqpScaling* scaling,
                       SleqpSparseMatrix* scaled_cons_jac);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_unscale_cons_duals(SleqpScaling* scaling, SleqpVec* scaled_cons_duals);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_unscale_var_duals(SleqpScaling* scaling, SleqpVec* scaled_var_duals);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_unscale_hessian_direction(SleqpScaling* scaling,
                                SleqpVec* direction,
                                SleqpVec* cons_duals);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_unscale_iterate(SleqpScaling* scaling,
                      SleqpIterate* scaled_iterate,
                      bool lsq);

#endif /* SLEQP_SCALE_H */
