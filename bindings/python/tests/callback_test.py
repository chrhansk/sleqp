#!/usr/bin/env python

import numpy as np
import unittest

import sleqp

from tests.rosenbrock_fixture import RosenbrockFunc

num_variables = 2
num_constraints = 0

class CallbackTest(unittest.TestCase):

  class AcceptedIterateHandler:
    def __init__(self):
      self.called = False

    def __call__(self, solver, iterate):
      self.called = True

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


  def test_callback_func(self):

    callback_called = [False]

    def accepted_iterate(solver, iterate):
      callback_called[0] = True

    self.solver.add_callback(sleqp.SolverEvent.AcceptedIterate,
                             accepted_iterate)

    self.solver.solve(100, 3600)

    self.assertEqual(self.solver.status, sleqp.Status.Optimal)

    self.assertTrue(callback_called[0])

  def test_callback_error(self):

    def accepted_iterate(solver, iterate):
      raise Exception("Error")

    self.solver.add_callback(sleqp.SolverEvent.AcceptedIterate,
                             accepted_iterate)

    self.assertRaises(sleqp.SLEQPError,
                      self.solver.solve, 100, 3600)

  def test_callback_object(self):

    accepted_iterate = self.AcceptedIterateHandler()

    self.solver.add_callback(sleqp.SolverEvent.AcceptedIterate,
                             accepted_iterate)

    self.solver.solve(100, 3600)

    self.assertEqual(self.solver.status, sleqp.Status.Optimal)

    self.assertTrue(accepted_iterate.called)

  def test_callback_object_nogil(self):

    sleqp.set_release_gil(True)

    try:
      accepted_iterate = self.AcceptedIterateHandler()

      self.solver.add_callback(sleqp.SolverEvent.AcceptedIterate,
                               accepted_iterate)

      self.solver.solve(100, 3600)

      self.assertEqual(self.solver.status, sleqp.Status.Optimal)

      self.assertTrue(accepted_iterate.called)
    finally:
      sleqp.set_release_gil(False)

  def test_remove_callback_handler(self):

    accepted_iterate = self.AcceptedIterateHandler()

    handle = self.solver.add_callback(sleqp.SolverEvent.AcceptedIterate,
                                      accepted_iterate)

    handle.unregister()

    self.solver.solve(100, 3600)

    self.assertEqual(self.solver.status, sleqp.Status.Optimal)

    self.assertFalse(accepted_iterate.called)

if __name__ == '__main__':
    unittest.main()
