#!/usr/bin/env python

import numpy as np
import unittest

import sleqp

num_variables = 2
num_constraints = 0

class Func(sleqp.Func):

  def __init__(self, num_variables, num_constraints):
    self.a = 1
    self.b = 1

  def set_value(self, v, reason):
    self.v = v

  def func_val(self):
    [x, y] = self.v
    (a, b) = (self.a, self.b)

    xsq = x**2

    return (a - x)**2 + b*(y - xsq)**2

  def func_grad_nnz(self):
    return 2

  def func_grad(self):
    [x, y] = self.v
    (a, b) = (self.a, self.b)

    xsq = x**2

    g = np.array([(4*b*x*(xsq - y)) + 2*x - 2*a,
                  -2*b*(xsq - y)])

    return g

  def hess_prod(self, func_dual, direction, _):
    [x, y] = self.v
    (a, b) = (self.a, self.b)
    [dx, dy] = direction

    xsq = x**2

    product = np.array([((8.*b*xsq + 4.*b*(xsq - y) + 2.)*dx - (4.*b*x)*dy)*func_dual,
                        ((-4.*b*x)*dx + (2.*b)*dy)*func_dual])

    return product


class UnconstrainedTest(unittest.TestCase):

  def setUp(self):
    inf = sleqp.inf()

    var_lb = np.array([-inf, -inf])
    var_ub = np.array([inf, inf])

    cons_lb = np.array([])
    cons_ub = np.array([])

    x = np.array([0., 0.])

    self.params = sleqp.Params()

    self.options = sleqp.Options()

    #a = 1.
    #b = 100.

    self.func = Func(num_variables, num_constraints)

    self.problem = sleqp.Problem(self.func,
                                 self.params,
                                 var_lb,
                                 var_ub,
                                 cons_lb,
                                 cons_ub)

    self.solver = sleqp.Solver(self.problem,
                               self.params,
                               self.options,
                               x)

  def test_solve(self):
    self.solver.solve(100, 3600)

    self.assertEqual(self.solver.status, sleqp.Status.Optimal)

    s = np.array([1., 1.])

    self.assertTrue(np.allclose(s, self.solver.solution.primal))

if __name__ == '__main__':
    unittest.main()
