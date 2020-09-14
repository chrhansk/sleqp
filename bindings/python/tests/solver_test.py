#!/usr/bin/env python

import numpy as np
import unittest

import sleqp

num_variables = 4
num_constraints = 2

initial_sol = np.array([1., 0., 2., 0.])

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

class SolverTest(unittest.TestCase):

    def setUp(self):
        inf = sleqp.inf()

        var_lb = np.array([-inf]*num_variables)
        var_ub = np.array([inf]*num_variables)

        cons_lb = np.array([-inf]*num_constraints)
        cons_ub = np.array([inf]*num_constraints)

        self.params = sleqp.Params()

        self.func = Func(num_variables, num_constraints)

        self.problem = sleqp.Problem(self.func,
                                     self.params,
                                     var_lb,
                                     var_ub,
                                     cons_lb,
                                     cons_ub)

        self.options = sleqp.Options()

    # Solution round trip array -> sparse vec -> array
    def test_set_solution(self):

        solver = sleqp.Solver(self.problem,
                              self.params,
                              self.options,
                              initial_sol)

        solver.solve(max_num_iterations=0)

        sol = solver.solution.primal

        self.assertTrue((sol == initial_sol).all())
