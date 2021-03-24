#!/usr/bin/env python

import numpy as np
import unittest

import sleqp


class MinimizeTest(unittest.TestCase):

  def test_unconstrained(self):
    from scipy.optimize import rosen, rosen_der, rosen_hess_prod

    x0 = [1.3, 0.7, 0.8, 1.9, 1.2]

    res = sleqp.minimize(rosen, x0, grad=rosen_der, hessp=rosen_hess_prod)

    self.assertTrue(np.allclose(res.x,
                                np.array([ 1.,  1.,  1.,  1.,  1.])))
