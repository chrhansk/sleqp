#!/usr/bin/env python

import numpy as np
import unittest

import sleqp

from .rosenbrock_fixture import *


class RejectFunc:
  def __init__(self):
    self.func = RosenbrockFunc()
    self.count = 0
    self.rejected = False

  def set_value(self, v, reason):
    if reason == sleqp.ValueReason.TryingIterate:
      self.count += 1
      if self.count % 4 == 1:
        # reject
        self.rejected = True
        return True

    self.rejected = False

    return self.func.set_value(v, reason)

  def obj_val(self):
    assert not self.rejected
    return self.func.obj_val()

  def obj_grad_nnz(self):
    assert not self.rejected
    return self.func.obj_grad_nnz()

  def obj_grad(self):
    assert not self.rejected
    return self.func.obj_grad()

  def cons_vals(self):
    assert not self.rejected
    return self.func.cons_vals()

  def cons_jac(self):
    assert not self.rejected
    return self.func.cons_jac()

  def hess_prod(self, *args):
    assert not self.rejected
    return self.func.hess_prod(*args)


class RejectTest(unittest.TestCase):

  def setUp(self):
    self.settings = sleqp.Settings()

    self.func = RejectFunc()

    self.problem = sleqp.Problem(self.func,
                                 self.settings,
                                 var_lb,
                                 var_ub,
                                 cons_lb,
                                 cons_ub)

  def test_solve(self):
    solver = sleqp.Solver(self.problem,
                          self.settings,
                          initial_sol)

    solver.solve(max_num_iterations=100)

    self.assertEqual(solver.status, sleqp.Status.Optimal)

    solution = solver.solution

    self.assertTrue(np.allclose(expected_sol, solution.primal))
