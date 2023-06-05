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

    func = RosenbrockFunc()

    problem = sleqp.Problem(func,
                            var_lb,
                            var_ub,
                            cons_lb,
                            cons_ub)

    self.solver = sleqp.Solver(problem,
                               self.initial_sol)


  def test_callback_func(self):

    callback_called = [False]

    def accepted_iterate(solver, iterate):
      callback_called[0] = True

    self.solver.add_callback(sleqp.SolverEvent.AcceptedIterate,
                             accepted_iterate)

    self.solver.solve(max_num_iterations=100)

    self.assertEqual(self.solver.status, sleqp.Status.Optimal)

    self.assertTrue(callback_called[0])

  def test_callback_error(self):

    exception = Exception("Error")

    def accepted_iterate(solver, iterate):
      raise exception

    self.solver.add_callback(sleqp.SolverEvent.AcceptedIterate,
                             accepted_iterate)

    fail = False

    try:
      self.solver.solve(max_num_iterations=100)
    except Exception as exc:
      cause = exc.__cause__
      orig = cause.__cause__

      self.assertTrue(isinstance(cause, sleqp.CallbackError))
      self.assertTrue(orig is exception)
      fail = True

    self.assertTrue(fail)

  def test_callback_abort(self):

    def performed_iteration(solver):
      solver.abort()

    self.solver.add_callback(sleqp.SolverEvent.PerformedIteration,
                             performed_iteration)

    self.solver.solve(max_num_iterations=100)

    self.assertEqual(self.solver.iterations, 1)

  def test_callback_object(self):

    accepted_iterate = self.AcceptedIterateHandler()

    self.solver.add_callback(sleqp.SolverEvent.AcceptedIterate,
                             accepted_iterate)

    self.solver.solve(max_num_iterations=100)

    self.assertEqual(self.solver.status, sleqp.Status.Optimal)

    self.assertTrue(accepted_iterate.called)

  def test_callback_object_nogil(self):

    sleqp.set_release_gil(True)

    try:
      accepted_iterate = self.AcceptedIterateHandler()

      self.solver.add_callback(sleqp.SolverEvent.AcceptedIterate,
                               accepted_iterate)

      self.solver.solve(max_num_iterations=100)

      self.assertEqual(self.solver.status, sleqp.Status.Optimal)

      self.assertTrue(accepted_iterate.called)
    finally:
      sleqp.set_release_gil(False)

  def test_remove_callback_handler(self):

    accepted_iterate = self.AcceptedIterateHandler()

    handle = self.solver.add_callback(sleqp.SolverEvent.AcceptedIterate,
                                      accepted_iterate)

    handle.unregister()

    self.solver.solve(max_num_iterations=100)

    self.assertEqual(self.solver.status, sleqp.Status.Optimal)

    self.assertFalse(accepted_iterate.called)

  def test_finished_callback(self):

    did_finish = 0

    def finished(solver, iterate):
      nonlocal did_finish
      did_finish += 1

    self.solver.add_callback(sleqp.SolverEvent.Finished,
                             finished)

    self.solver.solve(max_num_iterations=100)

    self.assertEqual(did_finish, 1)


if __name__ == '__main__':
    unittest.main()
