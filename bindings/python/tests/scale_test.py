#!/usr/bin/env python

import numpy as np
import unittest

import sleqp

from .zero_func import *

class ScaleTest(unittest.TestCase):

  def setUp(self):
    inf = sleqp.inf()

    var_lb = np.array([-inf]*num_variables)
    var_ub = np.array([inf]*num_variables)

    cons_lb = np.array([-inf]*num_constraints)
    cons_ub = np.array([inf]*num_constraints)

    self.settings = sleqp.Settings()

    self.func = ZeroFunc()

    self.problem = sleqp.Problem(self.func,
                                 self.settings,
                                 var_lb,
                                 var_ub,
                                 cons_lb,
                                 cons_ub)

    self.scaling = sleqp.Scaling(num_variables, num_constraints)

  def test_nominal_obj_weight(self):

    values = [(50., 6), (-50., 6), (0.75, 0), (0.3, -1)]

    for (nominal, expected_weight) in values:
      self.scaling.set_obj_weight_from_nominal(nominal)

      self.assertEqual(self.scaling.obj_weight, expected_weight)

  def test_nominal_var_weights(self):

    self.scaling.set_var_weights_from_nominal(np.array([1.9,
                                                        0.5,
                                                        -3.9,
                                                        60.]))

    expected_weights = np.array([1, 0, 2, 6],dtype=int)
    actual_weights = self.scaling.var_weights

    self.assertTrue((expected_weights == actual_weights).all())

  def test_nominal_cons_weights(self):

    self.scaling.set_cons_weights_from_nominal(np.array([1.9,
                                                         -60.]))

    expected_weights = np.array([1, 6],dtype=int)
    actual_weights = self.scaling.cons_weights

    self.assertTrue((expected_weights == actual_weights).all())
