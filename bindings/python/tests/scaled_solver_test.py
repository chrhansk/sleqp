#!/usr/bin/env python

import numpy as np
import unittest

import sleqp

from .constrained_fixture import *

class ScaledSolverTest(unittest.TestCase):

  def setUp(self):
    self.params = sleqp.Params()

    self.func = ConstrainedFunc()

    self.problem = sleqp.Problem(self.func,
                                 var_lb,
                                 var_ub,
                                 cons_lb,
                                 cons_ub)
    self.options = sleqp.Options()

  def test_scaled_solve(self):
    scaling = sleqp.Scaling(num_variables, num_constraints)

    scaling.func_weight = -3

    scaling.variable_weights = np.array([2, -2, 10, 1],
                                        dtype=int)

    scaling.constraint_weights = np.array([0, -1],
                                          dtype=int)

    solver = sleqp.Solver(self.problem,
                          self.params,
                          self.options,
                          initial_sol,
                          scaling)

    solver.solve(max_num_iterations=100)

    self.assertEqual(solver.status, sleqp.Status.Optimal)

    solution = solver.solution

    self.assertTrue(np.allclose(expected_sol, solution.primal))
