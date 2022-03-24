#ifndef SLEQP_PROBLEM_H
#define SLEQP_PROBLEM_H

#include "pub_problem.h"

#include "func.h"
#include "params.h"
#include "types.h"

bool
sleqp_problem_has_nonlinear_cons(SleqpProblem* problem);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_problem_set_value(SleqpProblem* problem,
                        SleqpVec* x,
                        SLEQP_VALUE_REASON reason,
                        bool* reject,
                        int* obj_grad_nnz,
                        int* cons_val_nnz,
                        int* cons_jac_nnz);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_problem_eval(SleqpProblem* problem,
                   double* obj_val,
                   SleqpVec* obj_grad,
                   SleqpVec* cons_val,
                   SleqpSparseMatrix* cons_jac);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_problem_obj_val(SleqpProblem* problem, double* obj_val);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_problem_obj_grad(SleqpProblem* problem, SleqpVec* obj_grad);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_problem_cons_val(SleqpProblem* problem, SleqpVec* cons_val);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_problem_cons_jac(SleqpProblem* problem, SleqpSparseMatrix* cons_jac);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_problem_hess_prod(SleqpProblem* problem,
                        const double* obj_dual,
                        const SleqpVec* direction,
                        const SleqpVec* cons_duals,
                        SleqpVec* product);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_problem_hess_bilinear(SleqpProblem* problem,
                            const double* obj_dual,
                            const SleqpVec* direction,
                            const SleqpVec* cons_duals,
                            double* bilinear_prod);

bool
sleqp_problem_is_unconstrained(SleqpProblem* problem);

#endif /* SLEQP_PROBLEM_H */
