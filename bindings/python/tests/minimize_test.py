#!/usr/bin/env python

import numpy as np
import unittest

from scipy.optimize import rosen, rosen_der, rosen_hess_prod

import sleqp

class MinimizeTest(unittest.TestCase):

  def setUp(self):
    self.x0 = [1.3, 0.7, 0.8, 1.9, 1.2]
    self.xopt = np.array([ 1.,  1.,  1.,  1.,  1.])

  def test_unconstrained(self):
    res = sleqp.minimize(rosen, self.x0, grad=rosen_der, hessp=rosen_hess_prod)

    self.assertTrue(res.success)
    self.assertTrue(np.allclose(res.x, self.xopt))

  def test_callback(self):
    def callback(x):
      assert x.shape == np.array(self.x0).shape

    res = sleqp.minimize(rosen, self.x0, grad=rosen_der, hessp=rosen_hess_prod, callback=callback)

    self.assertTrue(res.success)
    self.assertTrue(np.allclose(res.x, self.xopt))

  def test_callback_abort(self):
    def callback(x):
      return True

    res = sleqp.minimize(rosen, self.x0, grad=rosen_der, hessp=rosen_hess_prod, callback=callback)

    self.assertEqual(res['nit'], 1)
