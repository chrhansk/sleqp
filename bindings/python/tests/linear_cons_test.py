#!/usr/bin/env python

import numpy as np
import unittest

import sleqp

num_variables = 2
num_linear = 1
num_general = 0
num_constraints = num_linear + num_general

class Func:
    def set_value(self, v, reason):
      self.v = v

    def obj_val(self):
      n = np.linalg.norm(self.v)
      return .5 * n*n

    def obj_grad(self):
      return self.v

    def hess_prod(self, direction, cons_dual):
      return direction


class LinearConsTest(unittest.TestCase):

  def setUp(self):
    self.func = Func()

    var_lb = np.array([-np.inf, -np.inf])
    var_ub = np.array([np.inf, np.inf])

    cons_lb = np.zeros((num_general,))
    cons_ub = np.zeros((num_general,))

    self.initial_sol = np.array([1., 5.])

    self.linear_lb = np.array([1.])
    self.linear_ub = np.array([1.])

    self.linear_coeffs = np.zeros((num_linear, num_variables))

    self.linear_coeffs[0, 0] = 1.
    self.linear_coeffs[0, 1] = 1.

    self.problem = sleqp.Problem(self.func,
                                 var_lb,
                                 var_ub,
                                 cons_lb,
                                 cons_ub,
                                 linear_coeffs=self.linear_coeffs,
                                 linear_lb=self.linear_lb,
                                 linear_ub=self.linear_ub)

    self.optimal_sol = np.array([.5, .5])


  def test_solve(self):
    solver = sleqp.Solver(self.problem,
                          self.initial_sol)

    solver.solve(max_num_iterations=100)

    self.assertEqual(solver.status, sleqp.Status.Optimal)

    solution = solver.solution

    self.assertTrue(np.allclose(self.optimal_sol, solution.primal))

    self.assertEqual(solution.working_set.cons_state(0),
                     sleqp.ActiveState.ActiveBoth)


  def test_scaled_solve(self):
    scaling = sleqp.Scaling(num_variables, num_constraints)

    scaling.set_var_weight(0, -1)
    scaling.set_var_weight(1, -2)

    scaling.set_cons_weight(0, 2)

    solver = sleqp.Solver(self.problem,
                          self.initial_sol)

    solver.solve(max_num_iterations=100)

    self.assertEqual(solver.status, sleqp.Status.Optimal)

    solution = solver.solution

    self.assertTrue(np.allclose(self.optimal_sol, solution.primal))

    self.assertEqual(solution.working_set.cons_state(0),
                     sleqp.ActiveState.ActiveBoth)


if __name__ == '__main__':
  import logging
  logging.basicConfig(level=logging.INFO)
  unittest.main()
