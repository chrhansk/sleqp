#!/usr/bin/env python

import numpy as np
import unittest

import sleqp

from tests.rosenbrock_fixture import RosenbrockFunc

num_variables = 2
num_constraints = 0

class UnconstrainedTest(unittest.TestCase):

  def setUp(self):
    inf = sleqp.inf()

    var_lb = np.array([-inf, -inf])
    var_ub = np.array([inf, inf])

    cons_lb = np.array([])
    cons_ub = np.array([])

    self.initial_sol = np.array([0., 0.])

    self.params = sleqp.Params()

    self.options = sleqp.Options()

    func = RosenbrockFunc()

    problem = sleqp.Problem(func,
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
