#!/usr/bin/env python

import numpy as np
import unittest

import sleqp


class ConstrainedMinimizeTest(unittest.TestCase):

  def setUp(self):
    self.obj = lambda x: (x[0] - 1)**2 + (x[1] - 2.5)**2

    self.cons_funcs = [
      lambda x:  x[0] - 2 * x[1] + 2,
      lambda x: -x[0] - 2 * x[1] + 6,
      lambda x: -x[0] + 2 * x[1] + 2
    ]

    self.cons_jacs = [
      lambda x: [1., -2.],
      lambda x: [-1., -2.],
      lambda x: [-1., 2.]
    ]

    self.initial_sol = (2, 0)
    self.expected_sol = np.array([1.4, 1.7])

  def test_constrained(self):
    cons = ({'type': 'ineq', 'fun': self.cons_funcs[0]},
            {'type': 'ineq', 'fun': self.cons_funcs[1]},
            {'type': 'ineq', 'fun': self.cons_funcs[2]})

    bnds = ((0, None), (0, None))

    res = sleqp.minimize(self.obj, self.initial_sol, bounds=bnds, constraints=cons)

    self.assertTrue(res.success)
    self.assertTrue(np.allclose(res.x, self.expected_sol))

  def test_constrained_bounds(self):
    from scipy.optimize import Bounds

    cons = ({'type': 'ineq', 'fun': self.cons_funcs[0]},
            {'type': 'ineq', 'fun': self.cons_funcs[1]},
            {'type': 'ineq', 'fun': self.cons_funcs[2]})

    bnds = Bounds([0., 0.], [np.inf, np.inf])

    res = sleqp.minimize(self.obj, self.initial_sol, bounds=bnds, constraints=cons)

    self.assertTrue(res.success)
    self.assertTrue(np.allclose(res.x, self.expected_sol))

  def test_constrained_jac(self):
    cons = ({'type': 'ineq', 'fun': self.cons_funcs[0], 'jac': self.cons_jacs[0]},
            {'type': 'ineq', 'fun': self.cons_funcs[1], 'jac': self.cons_jacs[1]},
            {'type': 'ineq', 'fun': self.cons_funcs[2], 'jac': self.cons_jacs[2]})

    bnds = ((0, None), (0, None))

    res = sleqp.minimize(self.obj, self.initial_sol, bounds=bnds, constraints=cons)

    self.assertTrue(res.success)
    self.assertTrue(np.allclose(res.x, self.expected_sol))

  def test_constrained_jac_findiff_twopoint(self):
    cons = ({'type': 'ineq', 'fun': self.cons_funcs[0], 'jac': '2-point'},
            {'type': 'ineq', 'fun': self.cons_funcs[1], 'jac': '2-point'},
            {'type': 'ineq', 'fun': self.cons_funcs[2], 'jac': '2-point'})

    bnds = ((0, None), (0, None))

    res = sleqp.minimize(self.obj, self.initial_sol, bounds=bnds, constraints=cons)

    self.assertTrue(res.success)
    self.assertTrue(np.allclose(res.x, self.expected_sol))

  def test_constrained_jac_findiff_threepoint(self):
    cons = ({'type': 'ineq', 'fun': self.cons_funcs[0], 'jac': '3-point'},
            {'type': 'ineq', 'fun': self.cons_funcs[1], 'jac': '3-point'},
            {'type': 'ineq', 'fun': self.cons_funcs[2], 'jac': '3-point'})

    bnds = ((0, None), (0, None))

    res = sleqp.minimize(self.obj, self.initial_sol, bounds=bnds, constraints=cons)

    self.assertTrue(res.success)
    self.assertTrue(np.allclose(res.x, self.expected_sol))

  def test_constrained_jac_findiff_cs(self):
    cons = ({'type': 'ineq', 'fun': self.cons_funcs[0], 'jac': 'cs'},
            {'type': 'ineq', 'fun': self.cons_funcs[1], 'jac': 'cs'},
            {'type': 'ineq', 'fun': self.cons_funcs[2], 'jac': 'cs'})

    bnds = ((0, None), (0, None))

    res = sleqp.minimize(self.obj, self.initial_sol, bounds=bnds, constraints=cons)

    self.assertTrue(res.success)
    self.assertTrue(np.allclose(res.x, self.expected_sol))
