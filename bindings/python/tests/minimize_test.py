#!/usr/bin/env python

import numpy as np
import unittest

import sleqp


class MinimizeTest(unittest.TestCase):

  def test_unconstrained(self):
    from scipy.optimize import rosen, rosen_der, rosen_hess_prod

    x0 = [1.3, 0.7, 0.8, 1.9, 1.2]

    res = sleqp.minimize(rosen, x0, grad=rosen_der, hessp=rosen_hess_prod)

    self.assertTrue(res.success)
    self.assertTrue(np.allclose(res.x,
                                np.array([ 1.,  1.,  1.,  1.,  1.])))

  def test_constrained(self):

    fun = lambda x: (x[0] - 1)**2 + (x[1] - 2.5)**2

    cons = ({'type': 'ineq', 'fun': lambda x:  x[0] - 2 * x[1] + 2},
            {'type': 'ineq', 'fun': lambda x: -x[0] - 2 * x[1] + 6},
            {'type': 'ineq', 'fun': lambda x: -x[0] + 2 * x[1] + 2})

    bnds = ((0, None), (0, None))

    res = sleqp.minimize(fun, (2, 0), bounds=bnds, constraints=cons)

    self.assertTrue(res.success)
    self.assertTrue(np.allclose(res.x, np.array([1.4 ,1.7])))
