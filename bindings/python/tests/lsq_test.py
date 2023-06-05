#!/usr/bin/env python

import math
import unittest

import numpy as np

import sleqp

num_variables = 2
num_residuals = 2
num_constraints = 0

class Func:

  def __init__(self):
    self.a = 1
    self.b = 1
    self.values = np.zeros((num_variables,))


  def set_value(self, values, reason):
    self.values[:] = values

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

  def cons_vals(self):
    return np.zeros((num_constraints,))


  def cons_jac(self):
    return np.zeros((num_constraints, num_variables))


class LSQTest(unittest.TestCase):
  def setUp(self):
    inf = sleqp.inf()

    self.var_lb = np.array([-inf, -inf])
    self.var_ub = np.array([inf, inf])

    self.cons_lb = np.array([])
    self.cons_ub = np.array([])

    self.initial_sol = np.array([0., 0.])

    self.settings = sleqp.Settings()

    self.expected_sol = np.array([1., 1.])


  def test_solve(self):
    func = Func()

    problem = sleqp.LSQProblem(func,
                               self.settings,
                               self.var_lb,
                               self.var_ub,
                               self.cons_lb,
                               self.cons_ub,
                               num_residuals,
                               regularization=1e-4)

    solver = sleqp.Solver(problem,
                          self.settings,
                          self.initial_sol)

    solver.solve(100, 3600)

    self.assertEqual(solver.status, sleqp.Status.Optimal)

    self.assertTrue(np.allclose(self.expected_sol,
                                solver.solution.primal))

  def test_solve_ml(self):
    func = Func()

    problem = sleqp.LSQProblem(func,
                               self.settings,
                               self.var_lb,
                               self.var_ub,
                               self.cons_lb,
                               self.cons_ub,
                               num_residuals)

    solver = sleqp.Solver(problem,
                          self.settings,
                          self.initial_sol)

    solver.solve(100, 3600)

    self.assertEqual(solver.status, sleqp.Status.Optimal)

    self.assertTrue(np.allclose(self.expected_sol,
                                solver.solution.primal))

  def test_solve_nogil(self):
    func = Func()

    problem = sleqp.LSQProblem(func,
                               self.settings,
                               self.var_lb,
                               self.var_ub,
                               self.cons_lb,
                               self.cons_ub,
                               num_residuals)

    solver = sleqp.Solver(problem,
                          self.settings,
                          self.initial_sol)

    sleqp.set_release_gil(True)

    try:
      solver.solve(100, 3600)

      self.assertEqual(solver.status, sleqp.Status.Optimal)

      self.assertTrue(np.allclose(self.expected_sol,
                                  solver.solution.primal))

    finally:
      sleqp.set_release_gil(False)


if __name__ == '__main__':
  unittest.main()
