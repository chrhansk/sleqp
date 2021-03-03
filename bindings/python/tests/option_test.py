#!/usr/bin/env python

import unittest

import sleqp

class OptionTest(unittest.TestCase):

  def test_invalid_construction(self):

    def test():
      options = sleqp.Options(invalid_value=1)

      import pdb; pdb.set_trace()

    self.assertRaises(Exception, test)

  def test_roundtrip(self):
    options = sleqp.Options(dual_estimation_type=sleqp.DualEstimationType.LP)

    self.assertEqual(options.dual_estimation_type,
                     sleqp.DualEstimationType.LP)

    options.dual_estimation_type = sleqp.DualEstimationType.LSQ

    self.assertEqual(options.dual_estimation_type,
                     sleqp.DualEstimationType.LSQ)
