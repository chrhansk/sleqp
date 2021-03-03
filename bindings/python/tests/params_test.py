#!/usr/bin/env python

import unittest

import sleqp

class ParamsTest(unittest.TestCase):

  def test_invalid_construction(self):

    def test():
      params = sleqp.Params(invalid_value=1)

    self.assertRaises(AttributeError, test)

  def test_roundtrip(self):
    params = sleqp.Params(eps=1e-4)

    self.assertEqual(params.eps, 1e-4)

    params.eps = 1e-6

    self.assertEqual(params.eps, 1e-6)
