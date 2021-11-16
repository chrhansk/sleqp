#!/usr/bin/env python

import numpy as np
import unittest

import sleqp

from .constrained_fixture import *

class SolverTest(unittest.TestCase):

  def setUp(self):
    self.params = sleqp.Params()

    self.func = ConstrainedFunc()

    self.problem = sleqp.Problem(self.func,
                                 self.params,
                                 var_lb,
                                 var_ub,
                                 cons_lb,
                                 cons_ub)

  def get_solver(self, options=None, params=None):

    if options is None:
      options = sleqp.Options()

    if params is None:
      params = self.params

    return sleqp.Solver(self.problem,
                        params,
                        options,
                        initial_sol)

  def test_solve(self):
    solver = self.get_solver()

    solver.solve(max_num_iterations=100)

    self.assertEqual(solver.status, sleqp.Status.Optimal)

    solution = solver.solution

    self.assertTrue(np.allclose(expected_sol, solution.primal))

  def test_solve_linesearch(self):
    options = sleqp.Options(linesearch=sleqp.LineSearch.Exact)

    solver = self.get_solver(options=options)

    solver.solve(max_num_iterations=1000)

    self.assertEqual(solver.status, sleqp.Status.Optimal)

    solution = solver.solution

    self.assertTrue(np.allclose(expected_sol, solution.primal))

  def test_solve_lp_duals(self):
    options = sleqp.Options(dual_estimation_type=sleqp.DualEstimationType.LP)

    solver = self.get_solver(options=options)

    solver.solve(max_num_iterations=1000)

    self.assertEqual(solver.status, sleqp.Status.Optimal)

    solution = solver.solution

    self.assertTrue(np.allclose(expected_sol, solution.primal))

  def test_solve_no_newton(self):
    tol = 1e-4

    options = sleqp.Options(perform_newton_step=False)
    params = sleqp.Params(stationarity_tolerance=tol)

    solver = self.get_solver(options=options,
                             params=params)

    solver.solve(max_num_iterations=1000)

    self.assertEqual(solver.status, sleqp.Status.Optimal)

    solution = solver.solution

    self.assertTrue(np.allclose(expected_sol,
                                solution.primal,
                                rtol=tol))

  def test_solve_linear(self):
    options = sleqp.Options(use_quadratic_model=False)

    solver = self.get_solver(options=options)

    solver.solve(max_num_iterations=1000)

    self.assertEqual(solver.status, sleqp.Status.Optimal)

    solution = solver.solution

    self.assertTrue(np.allclose(expected_sol, solution.primal))

  def test_solve_linear_no_newton(self):
    options = sleqp.Options(perform_newton_step=False,
                            use_quadratic_model=False)

    solver = self.get_solver(options=options)

    solver.solve(max_num_iterations=1000)

    self.assertEqual(solver.status, sleqp.Status.Optimal)

    solution = solver.solution

    self.assertTrue(np.allclose(expected_sol, solution.primal))

  def test_iterate(self):
    solver = self.get_solver()

    solver.solve(max_num_iterations=1000)

    self.assertEqual(solver.status, sleqp.Status.Optimal)

    solution = solver.solution

    self.func.set_value(solution.primal, sleqp.ValueReason.NoReason)

    expected_obj_val = self.func.obj_val()

    self.assertTrue(np.allclose(np.array([expected_obj_val]),
                                np.array([solution.obj_val])))

    expected_cons_vals = self.func.cons_vals()

    self.assertTrue(np.allclose(expected_cons_vals,
                                solution.cons_val))

    expected_cons_jac = self.func.cons_jac()

    actual_cons_jac = solution.cons_jac.toarray()

    self.assertTrue(np.allclose(expected_cons_jac,
                                actual_cons_jac))


if __name__ == '__main__':
  import logging
  logging.basicConfig(level=logging.INFO)
  unittest.main()
