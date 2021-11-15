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
                        SleqpSparseVec* x,
                        SLEQP_VALUE_REASON reason,
                        bool* reject,
                        int* func_grad_nnz,
                        int* cons_val_nnz,
                        int* cons_jac_nnz);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_problem_eval(SleqpProblem* problem,
                   const SleqpSparseVec* cons_indices,
                   double* func_val,
                   SleqpSparseVec* func_grad,
                   SleqpSparseVec* cons_val,
                   SleqpSparseMatrix* cons_jac);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_problem_val(SleqpProblem* problem, double* func_val);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_problem_grad(SleqpProblem* problem, SleqpSparseVec* func_grad);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_problem_cons_val(SleqpProblem* problem,
                       const SleqpSparseVec* cons_indices,
                       SleqpSparseVec* cons_val);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_problem_cons_jac(SleqpProblem* problem,
                       const SleqpSparseVec* cons_indices,
                       SleqpSparseMatrix* cons_jac);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_problem_hess_prod(SleqpProblem* problem,
                        const double* func_dual,
                        const SleqpSparseVec* direction,
                        const SleqpSparseVec* cons_duals,
                        SleqpSparseVec* product);

SLEQP_NODISCARD
SLEQP_RETCODE
sleqp_problem_hess_bilinear(SleqpProblem* problem,
                            const double* func_dual,
                            const SleqpSparseVec* direction,
                            const SleqpSparseVec* cons_duals,
                            double* bilinear_prod);

bool
sleqp_problem_is_unconstrained(SleqpProblem* problem);

#endif /* SLEQP_PROBLEM_H */
