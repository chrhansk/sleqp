#!/usr/bin/env python

import numpy as np
import unittest

import sleqp

from tests.rosenbrock_fixture import *


class UnconstrainedTest(unittest.TestCase):

  def setUp(self):
    self.settings = sleqp.Settings()

    func = RosenbrockFunc()

    problem = sleqp.Problem(func,
                            self.settings,
                            var_lb,
                            var_ub,
                            cons_lb,
                            cons_ub)

    self.solver = sleqp.Solver(problem,
                               self.settings,
                               initial_sol)

  def test_solve(self):
    self.solver.solve(100, 3600)

    self.assertEqual(self.solver.status, sleqp.Status.Optimal)

    self.assertTrue(np.allclose(expected_sol,
                                self.solver.solution.primal))

  def test_solve_nogil(self):
    sleqp.set_release_gil(True)

    try:
      self.solver.solve(100, 3600)

      self.assertEqual(self.solver.status, sleqp.Status.Optimal)

      self.assertTrue(np.allclose(expected_sol,
                                  self.solver.solution.primal))

    finally:
      sleqp.set_release_gil(False)


if __name__ == '__main__':
    unittest.main()
