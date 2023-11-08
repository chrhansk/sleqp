#!/usr/bin/env python

import unittest

import sleqp


class SettingsTest(unittest.TestCase):

    def test_invalid_construction(self):

        def test():
            settings = sleqp.Settings(invalid_value=1)

        self.assertRaises(AttributeError, test)

    def test_roundtrip(self):
        settings = sleqp.Settings(eps=1e-4)

        self.assertEqual(settings.eps, 1e-4)

        settings.eps = 1e-6

        self.assertEqual(settings.eps, 1e-6)

    def test_set_enum(self):
        settings = sleqp.Settings(deriv_check=sleqp.DerivCheck.FirstCons)

        self.assertEqual(settings.deriv_check, sleqp.DerivCheck.FirstCons)

    def test_set_enum_str(self):
        settings = sleqp.Settings(deriv_check="FirstCons")

        self.assertEqual(settings.deriv_check, sleqp.DerivCheck.FirstCons)
