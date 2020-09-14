#!/usr/bin/env python

import numpy as np
import unittest

import sleqp

num_variables = 4
num_constraints = 2

class Func(sleqp.Func):

    def set_value(self, values, reason):
        assert((values == initial_sol).all())

    def func_val(self):
        return 0.

    def cons_vals(self):
        return np.zeros((num_constraints,))

    def cons_jac(self):
        return np.zeros((num_constraints, num_variables))

    def hess_prod(self, func_dual, direction, cons_dual):
        return np.zeros((num_variables,))


class ProblemTest(unittest.TestCase):

    def setUp(self):

        self.var_lb = np.array([-1.]*num_variables)
        self.var_ub = np.array([1.]*num_variables)

        self.cons_lb = np.array([-2.]*num_constraints)
        self.cons_ub = np.array([2.]*num_constraints)

        self.params = sleqp.Params()

        self.func = Func(num_variables, num_constraints)

        self.problem = sleqp.Problem(self.func,
                                     self.params,
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
