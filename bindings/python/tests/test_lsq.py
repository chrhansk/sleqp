#!/usr/bin/env python

import math
import unittest

import numpy as np

import sleqp

num_variables = 2
num_residuals = 2
num_constraints = 0

class Func(sleqp.LSQFunc):

  def __init__(self,
               num_variables,
               num_constraints,
               num_residuals,
               levenberg_marquardt,
               params,
               *args,
               **keywords):
    self.a = 1
    self.b = 1
    self.values = np.zeros((num_variables,))


  def set_value(self, values, reason):
    self.values[:] = values

  def func_val(self):
    # No additional function value beside lsq residuals
    return 0.


  def lsq_residuals(self):
    x0 = self.values[0]
    x1 = self.values[1]

    return np.array([self.a - x0,
                     math.sqrt(self.b)*(x1 - (x0*x0))])


  def lsq_jac_forward(self, forward_direction):
    x0 = self.values[0]
    x1 = self.values[1]

    d0 = forward_direction[0]
    d1 = forward_direction[1]

    return np.array([-1. * d0,
                     math.sqrt(self.b)*(-2.*x0*d0 + d1)])


  def lsq_jac_adjoint(self, adjoint_direction):
    x0 = self.values[0]
    x1 = self.values[1]

    d0 = adjoint_direction[0]
    d1 = adjoint_direction[1]

    return np.array([-1.*d0 - 2*math.sqrt(self.b)*x0*d1,
                     math.sqrt(self.b)*(d1)])


  def func_grad_nnz(self):
    return 0


  def cons_val_nnz(self):
    return 0


  def cons_jac_nnz(self):
    return 0


  def func_grad(self):
    return np.zeros((num_variables,))


  def cons_vals(self):
    return np.zeros((num_constraints,))


  def cons_jac(self):
    return np.zeros((num_constraints, num_variables))


  def hess_prod(self, func_dual, direction, cons_dual):
    return np.zeros((num_variables,))


class LSQTest(unittest.TestCase):
  def setUp(self):
    inf = sleqp.inf()

    var_lb = np.array([-inf, -inf])
    var_ub = np.array([inf, inf])

    cons_lb = np.array([])
    cons_ub = np.array([])

    self.initial_sol = np.array([0., 0.])

    self.params = sleqp.Params()

    self.options = sleqp.Options()

    func = Func(num_variables, num_constraints, num_residuals, 0., self.params)

    problem = sleqp.Problem(func,
                            self.params,
                            var_lb,
                            var_ub,
                            cons_lb,
                            cons_ub)

    self.solver = sleqp.Solver(problem,
                               self.params,
                               self.options,
                               self.initial_sol)

    self.expected_sol = np.array([1., 1.])


  def test_solve(self):
    self.solver.solve(100, 3600)

    self.assertEqual(self.solver.status, sleqp.Status.Optimal)

    self.assertTrue(np.allclose(self.expected_sol,
                                self.solver.solution.primal))


  def test_solve_nogil(self):
    sleqp.set_release_gil(True)

    try:
      self.solver.solve(100, 3600)

      self.assertEqual(self.solver.status, sleqp.Status.Optimal)

      self.assertTrue(np.allclose(self.expected_sol,
                                  self.solver.solution.primal))

    finally:
      sleqp.set_release_gil(False)


if __name__ == '__main__':
  unittest.main()
