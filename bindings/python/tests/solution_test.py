#!/usr/bin/env python

import numpy as np
import unittest

import sleqp

from .zero_func import *

class SolutionTest(unittest.TestCase):

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

    # Solution round trip array -> sparse vec -> array
    def test_set_solution(self):

        solver = sleqp.Solver(self.problem,
                              self.settings,
                              initial_sol)

        solver.solve(max_num_iterations=0)

        sol = solver.solution.primal

        self.assertTrue((sol == initial_sol).all())
