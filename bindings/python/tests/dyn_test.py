#!/usr/bin/env python

import numpy as np
import unittest

import sleqp

from tests.rosenbrock_fixture import *


def error_func(*args):
    raise Exception()


class DynRosenbrockFunc:
    def __init__(self):
        self.func = RosenbrockFunc()
        self.rng = np.random.default_rng(seed=42)

    def set_value(self, v, reason):
        self.func.set_value(v, reason)

    def obj_grad(self):
        return self.func.obj_grad()

    def hess_prod(self, direction, cons_duals):
        return self.func.hess_prod(direction,
                                   cons_duals)

    def eval(self):
        assert self.error_bound is not None
        assert self.obj_weight is not None
        val = self.func.obj_val()

        noise = self.rng.uniform(-1., 1., 1).item()

        factor = self.error_bound / self.obj_weight

        return (val + factor), (abs(noise) * self.error_bound)

    def set_obj_weight(self, obj_weight):
        assert obj_weight > 0.
        self.obj_weight = obj_weight

    def set_error_bound(self, error_bound):
        assert error_bound > 0.
        self.error_bound = error_bound


class DynTest(unittest.TestCase):

    def setUp(self):
        self.func = DynRosenbrockFunc()

        self.problem = sleqp.DynProblem(self.func,
                                        var_lb,
                                        var_ub,
                                        cons_lb,
                                        cons_ub)

        self.solver = sleqp.Solver(self.problem,
                                   initial_sol)

    def test_solve(self):
        self.solver.solve(100, 3600.)

        self.assertEqual(self.solver.status, sleqp.Status.Optimal)

        self.assertTrue(np.allclose(expected_sol,
                                    self.solver.solution.primal))

    def test_eval_error(self):
        self.func.eval = error_func
        with self.assertRaises(Exception):
            self.solver.solve(100, 3600.)

    def test_weight_error(self):
        self.func.set_obj_weight = error_func
        with self.assertRaises(Exception):
            solver = sleqp.Solver(self.problem,
                                  initial_sol)
            solver.solve(100, 3600.)

    def test_bound_error(self):
        self.func.set_error_bound = error_func
        with self.assertRaises(Exception):
            self.solver.solve(100, 3600.)

    def test_solve_nogil(self):
        sleqp.set_release_gil(True)

        try:
            self.solver.solve(100, 3600.)

            self.assertEqual(self.solver.status, sleqp.Status.Optimal)

            self.assertTrue(np.allclose(expected_sol,
                                        self.solver.solution.primal))

        finally:
            sleqp.set_release_gil(False)


if __name__ == '__main__':
    unittest.main()
