#!/usr/bin/env python

from math import sqrt
import unittest

import numpy as np

import sleqp

num_constraints = 0
num_variables = 2

class LSQFunc(sleqp.Func):
    """
    Simple LSQ function to solve the Rosenbrock problem.
    Uses an inexact Hessian based on first-order information.

    To adapt to any least squares problem, it is sufficient
    to replace the `eval_lsq_func` and `eval_lsq_jac` functions.
    """
    def __init__(self,
                 num_vars,
                 num_cons,
                 *args,
                 **keywords):
        self.a = 1
        self.b = 100

    def set_value(self, x, reason):
        self.x = x

        self.lsq_val = self.eval_lsq_func()

        self.lsq_jac = self.eval_lsq_jac()

    def func_val(self):
        return .5 * np.inner(self.lsq_val, self.lsq_val)

    def eval_lsq_func(self):
        x0 = self.x[0]
        x1 = self.x[1]

        a = self.a
        b = self.b

        return np.array([(a - x0), sqrt(b)*(x1 - (x0*x0))])

    def eval_lsq_jac(self):
        x0 = self.x[0]
        x1 = self.x[1]

        a = self.a
        b = self.b

        return np.array([[-1., 0.],
                         [-2*sqrt(b)*x0, sqrt(b)]])

    def func_grad_nnz(self):
        return self.num_variables

    def func_grad(self):
        return np.dot(self.lsq_val, self.lsq_jac)

    def cons_vals(self):
        return np.zeros((num_constraints,))

    def cons_jac(self):
        return np.zeros((num_constraints, num_variables))

    def hess_prod(self, func_dual, direction, cons_dual):

        if not func_dual: return np.zeros((self.num_vars,))

        # use the first-order approximation to the Hessian
        return func_dual * np.dot(np.transpose(self.lsq_jac),
                                  np.dot(self.lsq_jac, direction))

class LSQImplicitFunc(sleqp.Func):
    """
    Simple LSQ function to solve the Rosenbrock problem.
    Uses an inexact Hessian based on first-order information.
    This function only computes Jacobian product, never
    the acutal Jacobian (well, or at least it could in principle).

    To adapt to any least squares problem, it is sufficient
    to replace the `eval_lsq_func`, `eval_lsq_jac_forward`,
    and `eval_lsq_jac_adjoint` functions.

    """
    def __init__(self,
                 num_vars,
                 num_cons,
                 *args,
                 **keywords):
        self.a = 1
        self.b = 100

    def set_value(self, x, reason):
        self.x = x

        self.lsq_val = self.eval_lsq_func()

    def func_val(self):
        return .5 * np.inner(self.lsq_val, self.lsq_val)

    def eval_lsq_func(self):
        x0 = self.x[0]
        x1 = self.x[1]

        a = self.a
        b = self.b

        return np.array([(a - x0), sqrt(b)*(x1 - (x0*x0))])

    def _eval_jac(self):
        x0 = self.x[0]
        x1 = self.x[1]

        a = self.a
        b = self.b

        jac = np.array([[-1., 0.],
                        [-2*sqrt(b)*x0, sqrt(b)]])

        return jac

    def eval_lsq_jac_forward(self, direction):
        jac = self._eval_jac()
        return np.dot(jac, direction)

    def eval_lsq_jac_adjoint(self, direction):
        jac = self._eval_jac()
        return np.dot(np.transpose(jac), direction)

    def cons_vals(self):
        return np.zeros((num_constraints,))

    def cons_jac(self):
        return np.zeros((num_constraints, num_variables))

    def func_grad_nnz(self):
        return self.num_variables

    def func_grad(self):
        return self.eval_lsq_jac_adjoint(self.lsq_val)

    def hess_prod(self, func_dual, direction, cons_dual):

        if not func_dual: return np.zeros((self.num_vars,))

        # use the first-order approximation to the Hessian
        return func_dual * self.eval_lsq_jac_adjoint(self.eval_lsq_jac_forward(direction))


class LSQTest(unittest.TestCase):

    def setUp(self):
        self.num_vars = 2
        self.num_cons = 0

        self.var_lb = np.full((self.num_vars,), -5.)
        self.var_ub = np.full((self.num_vars,), 5.)

        self.cons_lb = np.full((self.num_cons,), 0.)
        self.cons_ub = np.full((self.num_cons,), 0.)

        self.params = sleqp.Params()

        self.options = sleqp.Options(deriv_check=sleqp.DerivCheck.First)

        self.initial_sol = np.full((self.num_vars,), 0.)

        self.target_sol = np.array([1., 1.])


    def test_simple(self):

        func = LSQFunc(self.num_vars, self.num_cons)

        problem = sleqp.Problem(func,
                                self.params,
                                self.var_lb,
                                self.var_ub,
                                self.cons_lb,
                                self.cons_ub)

        solver = sleqp.Solver(problem,
                              self.params,
                              self.options,
                              self.initial_sol)

        solver.solve(100, 3600)

        self.assertEqual(solver.status, sleqp.Status.Optimal)

        self.assertTrue(np.allclose(self.target_sol, solver.solution.primal))


    def test_implicit(self):

        func = LSQImplicitFunc(self.num_vars, self.num_cons)

        problem = sleqp.Problem(func,
                                self.params,
                                self.var_lb,
                                self.var_ub,
                                self.cons_lb,
                                self.cons_ub)

        solver = sleqp.Solver(problem,
                              self.params,
                              self.options,
                              self.initial_sol)

        solver.solve(100, 3600)

        self.assertEqual(solver.status, sleqp.Status.Optimal)

        self.assertTrue(np.allclose(self.target_sol, solver.solution.primal))

if __name__ == "__main__":
    unittest.main()
