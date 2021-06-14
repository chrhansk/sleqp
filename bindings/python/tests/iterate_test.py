#!/usr/bin/env python

import numpy as np
import unittest

import sleqp

from .constrained_fixture import *

class IterateTest(unittest.TestCase):

  def test_iterate(self):
    params = sleqp.Params()

    options = sleqp.Options()

    func = ConstrainedFunc()

    problem = sleqp.Problem(func,
                            params,
                            var_lb,
                            var_ub,
                            cons_lb,
                            cons_ub)

    solver = sleqp.Solver(problem,
                          params,
                          options,
                          initial_sol)

    solver.solve(max_num_iterations=1000)

    self.assertEqual(solver.status, sleqp.Status.Optimal)

    solution = solver.solution

    func.set_value(solution.primal, sleqp.ValueReason.NoReason)

    expected_func_val = func.func_val()

    self.assertTrue(np.allclose(np.array([expected_func_val]),
                                np.array([solution.func_val])))

    expected_cons_vals = func.cons_vals()

    self.assertTrue(np.allclose(expected_cons_vals,
                                solution.cons_val))

    expected_cons_jac = func.cons_jac()

    actual_cons_jac = solution.cons_jac.toarray()

    self.assertTrue(np.allclose(expected_cons_jac,
                                actual_cons_jac))
