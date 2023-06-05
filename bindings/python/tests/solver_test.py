#!/usr/bin/env python

import numpy as np
import unittest

import sleqp

from .constrained_fixture import *

class SolverTest(unittest.TestCase):

  def setUp(self):
    self.settings = sleqp.Settings()

    self.func = ConstrainedFunc()

    self.problem = sleqp.Problem(self.func,
                                 self.settings,
                                 var_lb,
                                 var_ub,
                                 cons_lb,
                                 cons_ub)

  def get_solver(self, settings=None):

    if settings is None:
      settings = self.settings

    return sleqp.Solver(self.problem,
                        settings,
                        initial_sol)

  def test_solve(self):
    solver = self.get_solver()

    solver.solve(max_num_iterations=100)

    self.assertEqual(solver.status, sleqp.Status.Optimal)

    solution = solver.solution

    self.assertTrue(np.allclose(expected_sol, solution.primal))

  def test_solve_linesearch(self):
    settings = sleqp.Settings(linesearch=sleqp.LineSearch.Exact)

    solver = self.get_solver(settings=settings)

    solver.solve(max_num_iterations=1000)

    self.assertEqual(solver.status, sleqp.Status.Optimal)

    solution = solver.solution

    self.assertTrue(np.allclose(expected_sol, solution.primal))

  def test_solve_lp_duals(self):
    settings = sleqp.Settings(dual_estimation_type=sleqp.DualEstimationType.LP)

    solver = self.get_solver(settings=settings)

    solver.solve(max_num_iterations=1000)

    self.assertEqual(solver.status, sleqp.Status.Optimal)

    solution = solver.solution

    self.assertTrue(np.allclose(expected_sol, solution.primal))

  def test_solve_no_newton(self):
    tol = 1e-4

    settings = sleqp.Settings(stat_tol=tol, perform_newton_step=False)

    solver = self.get_solver(settings=settings)

    solver.solve(max_num_iterations=1000)

    self.assertEqual(solver.status, sleqp.Status.Optimal)

    solution = solver.solution

    self.assertTrue(np.allclose(expected_sol,
                                solution.primal,
                                rtol=tol))

  def test_solve_linear(self):
    settings = sleqp.Settings(use_quadratic_model=False)

    solver = self.get_solver(settings=settings)

    solver.solve(max_num_iterations=1000)

    self.assertEqual(solver.status, sleqp.Status.Optimal)

    solution = solver.solution

    self.assertTrue(np.allclose(expected_sol, solution.primal))

  def test_solve_linear_no_newton(self):
    settings = sleqp.Settings(perform_newton_step=False,
                              use_quadratic_model=False)

    solver = self.get_solver(settings=settings)

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
