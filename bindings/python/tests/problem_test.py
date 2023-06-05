#!/usr/bin/env python

import numpy as np
import unittest

import sleqp

num_variables = 4
num_constraints = 2

class Func:

    def set_value(self, values, reason):
        assert((values == initial_sol).all())

    def obj_val(self):
        return 0.

    def cons_vals(self):
        return np.zeros((num_constraints,))

    def cons_jac(self):
        return np.zeros((num_constraints, num_variables))

    def hess_prod(self, direction, cons_dual):
        return np.zeros((num_variables,))


class ProblemTest(unittest.TestCase):

    def setUp(self):

        self.var_lb = np.array([-1.]*num_variables)
        self.var_ub = np.array([1.]*num_variables)

        self.cons_lb = np.array([-2.]*num_constraints)
        self.cons_ub = np.array([2.]*num_constraints)

        self.func = Func()

        self.problem = sleqp.Problem(self.func,
                                     self.var_lb,
                                     self.var_ub,
                                     self.cons_lb,
                                     self.cons_ub)

    def test_var_lb(self):
        self.assertTrue((self.problem.var_lb == self.var_lb).all())

    def test_var_ub(self):
        self.assertTrue((self.problem.var_ub == self.var_ub).all())

    def test_cons_lb(self):
        self.assertTrue((self.problem.cons_lb == self.cons_lb).all())

    def test_cons_ub(self):
        self.assertTrue((self.problem.cons_ub == self.cons_ub).all())

    def test_invalid_cons_bounds(self):

        cons_lb = self.cons_lb[:]
        cons_ub = self.cons_ub[:]

        cons_lb[0] = 1.
        cons_ub[0] = 0.

        with self.assertRaises(ValueError):
          sleqp.Problem(self.func,
                        self.var_lb,
                        self.var_ub,
                        cons_lb,
                        cons_ub)

    def test_invalid_var_bounds(self):

        var_lb = self.var_lb[:]
        var_ub = self.var_ub[:]

        var_lb[0] = 1.
        var_ub[0] = 0.

        with self.assertRaises(ValueError):
          sleqp.Problem(self.func,
                        var_lb,
                        var_ub,
                        self.cons_lb,
                        self.cons_ub)
