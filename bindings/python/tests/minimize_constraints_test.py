#!/usr/bin/env python

import numpy as np
import unittest

from scipy.optimize import (LinearConstraint,
                            NonlinearConstraint)

import sleqp


class ConstrainedMinimizeTest(unittest.TestCase):

  def setUp(self):
    self.obj = lambda x: (x[0] - 1)**2 + (x[1] - 2.5)**2

    self.grad = lambda x: [2.*(x[0] - 1), 2.*(x[1] - 2.5)]

    self.hessp = lambda x, d: 2*d

    cons_matrix = np.array([[ 1., -2.],
                            [-1., -2.],
                            [-1.,  2.]])

    cons_lb = np.array([-2.,
                        -6.,
                         2.])

    cons_ub = np.array([np.inf,
                        np.inf,
                        np.inf])

    self.linear_cons = LinearConstraint(cons_matrix,
                                        cons_lb,
                                        cons_ub)

    nonlinear_fun = lambda x: np.dot(cons_matrix, x)
    nonlinear_jac = lambda x: cons_matrix

    self.nonlinear_cons = NonlinearConstraint(nonlinear_fun,
                                              cons_lb,
                                              cons_ub,
                                              jac=nonlinear_jac,
                                              hess='2-point')

    self.initial_sol = np.array([2, 0])
    self.expected_sol = np.array([1.4, 1.7])
    self.bounds = ((0, None), (0, None))

  def test_linear_constraint(self):

    res = sleqp.minimize(self.obj,
                         self.initial_sol,
                         bounds=self.bounds,
                         constraints=self.linear_cons)

    self.assertTrue(res.success)
    self.assertTrue(np.allclose(res.x, self.expected_sol))

  def test_nonlinear_constraint(self):

    res = sleqp.minimize(self.obj,
                         self.initial_sol,
                         jac=self.grad,
                         hessp=self.hessp,
                         bounds=self.bounds,
                         constraints=self.nonlinear_cons)

    self.assertTrue(res.success)
    self.assertTrue(np.allclose(res.x, self.expected_sol))
