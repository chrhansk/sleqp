#!/usr/bin/env python

import math
import numpy as np
import unittest

import sleqp

num_variables = 2
num_constraints = 1

class Func:
  def set_value(self, v, reason):
    self.v = v

  def func_val(self):
    [x, y] = self.v

    return 2*(x**2 + y**2 - 1) - x

  def func_grad(self):
    [x, y] = self.v

    return np.array([4*x - 1, 4*y])

  def cons_vals(self):
    [x, y] = self.v

    return np.array([x**2 + y**2 - 1])

  def cons_jac(self):
    [x, y] = self.v
    return np.array([[2*x, 2*y]])

  def hess_prod(self, func_dual, direction, cons_duals):
    return 4*func_dual*direction + 2*cons_duals.item()*direction


class SecondOrderTest(unittest.TestCase):

  def setUp(self):
    var_lb = np.array([-np.inf, -np.inf])
    var_ub = np.array([np.inf, np.inf])

    cons_lb = np.array([0.])
    cons_ub = np.array([0.])

    self.params = sleqp.Params()

    self.options = sleqp.Options()

    func = Func()

    problem = sleqp.Problem(func,
                            var_lb,
                            var_ub,
                            cons_lb,
                            cons_ub)

    theta = math.pi/20.

    self.initial_sol = np.array([np.cos(theta),
                                 np.sin(theta)])

    self.solver = sleqp.Solver(problem,
                               self.params,
                               self.options,
                               self.initial_sol)

    self.expected_sol = np.array([1., 0.])

  def test_soc_step(self):
    self.solver.solve(1)

    self.assertEqual(self.solver.states[sleqp.SolverState.LastStepType],
                     sleqp.StepType.AcceptedSOC)

  def test_solve(self):
    self.solver.solve()

    self.assertTrue(np.allclose(self.solver.solution.primal,
                                self.expected_sol))
