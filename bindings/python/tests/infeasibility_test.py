#!/usr/bin/env python

import math
import numpy as np
import unittest

import sleqp

num_variables = 1
num_constraints = 1

class Func:
  def set_value(self, v, reason):
    [self.x] = v

  def func_val(self):
    return 0.

  def func_grad(self):
    return np.array([0.])

  def cons_vals(self):
    x = self.x
    return np.array([x**2])

  def cons_jac(self):
    x = self.x
    return np.array([[2*x]])

  def hess_prod(self, func_dual, direction, cons_duals):
    return 2*cons_duals.item()*direction


class SecondOrderTest(unittest.TestCase):

  def setUp(self):
    var_lb = np.array([-np.inf])
    var_ub = np.array([np.inf])

    cons_lb = np.array([-1.])
    cons_ub = np.array([-1.])

    self.params = sleqp.Params()

    self.options = sleqp.Options()

    func = Func()

    problem = sleqp.Problem(func,
                            var_lb,
                            var_ub,
                            cons_lb,
                            cons_ub)

    self.initial_sol = np.array([10.])

    self.solver = sleqp.Solver(problem,
                               self.params,
                               self.options,
                               self.initial_sol)

    self.expected_sol = np.array([0.])

  def test_solve(self):
    self.solver.solve(100, 3600)

    # Want to detect the infeasibility
    # at least before the iterations run out
    self.assertTrue(self.solver.iterations < 100)

    self.assertTrue(np.allclose(self.solver.solution.primal,
                                self.expected_sol))
