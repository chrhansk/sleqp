#!/usr/bin/env python

import numpy as np
import unittest

import sleqp


class PolishingFunc:

  def set_value(self, v, reason):
    self.v = np.copy(v)

  def obj_val(self):
    [x, y] = self.v

    return 0.5 * (x - 1.)**2 + (y - 2.5)**2

  def obj_grad(self):
    [x, y] = self.v

    return np.array([x - 1, y - 2.5])

  def hess_prod(self, direction, _):
    return direction

# The first linear constraint is active,
# the second is not. Correspondingly, its
# dual must be zero. However, for a sufficiently
# large trust region, the second constraint is
# active at the LP step and therefore in the
# working set. The polishing should
# remove it.
class PolishingTest(unittest.TestCase):

  def setUp(self):

    inf = np.inf

    self.var_lb = np.array([-inf, -inf])
    self.var_ub = np.array([inf, inf])

    self.cons_lb = np.array([])
    self.cons_ub = np.array([])

    self.func = PolishingFunc()

    self.params = sleqp.Params()

    self.options = sleqp.Options()

    self.linear_coeffs = np.array([[ 1., -2.],
                                   [-1., -2.]])

    self.linear_lb = np.array([-2.,
                               -6.])

    self.linear_ub = np.array([inf,
                               inf])

    self.problem = sleqp.Problem(self.func,
                                 self.params,
                                 self.var_lb,
                                 self.var_ub,
                                 self.cons_lb,
                                 self.cons_ub,
                                 linear_coeffs=self.linear_coeffs,
                                 linear_lb=self.linear_lb,
                                 linear_ub=self.linear_ub)

    self.initial_sol = np.array([2., 0.])

    self.expected_primal = np.array([1.4, 1.7])

    self.expected_cons_dual = np.array([-0.4, 0.])

  def check_sol(self, solution):
    primal = solution.primal
    cons_dual = solution.cons_dual

    self.assertTrue(np.allclose(self.expected_primal,
                                primal))

    self.assertTrue(np.allclose(self.expected_cons_dual,
                                cons_dual))

    num_cons = self.problem.cons_lb.size
    working_set = solution.working_set

    for i in range(num_cons):
      state = working_set.cons_state(i)

      if cons_dual[i] ==  0.:
        self.assertEqual(state, sleqp.ActiveState.Inactive)

  def test_polishing_zero_dual(self):

    self.options.polishing_type = sleqp.PolishingType.ZeroDual

    solver = sleqp.Solver(self.problem,
                          self.params,
                          self.options,
                          self.initial_sol)

    solver.solve(100, 3600)

    self.assertEqual(solver.status, sleqp.Status.Optimal)

    self.check_sol(solver.solution)

  def test_polishing_inactive(self):

    self.options.polishing_type = sleqp.PolishingType.Inactive

    solver = sleqp.Solver(self.problem,
                          self.params,
                          self.options,
                          self.initial_sol)

    solver.solve(100, 3600)

    self.assertEqual(solver.status, sleqp.Status.Optimal)

    self.check_sol(solver.solution)


if __name__ == '__main__':
  unittest.main()
