#!/usr/bin/env python

import numpy as np
import unittest

import sleqp

num_variables = 4
num_constraints = 2

def sq(x):
  return x*x

class Func(sleqp.Func):

  def set_value(self, v):
    print("set_value")
    self.v= v
    print(v)

  def func_val(self):
    x = self.v

    return x[0]*x[3]*(x[0] + x[1]+ x[2]) + x[2]

  def func_grad(self):
    x = self.v

    return np.array([(x[0] + x[1] + x[2])*x[3] + x[0]*x[3],
                     x[0]*x[3],
                     x[0]*x[3] + 1,
                     (x[0] + x[1] + x[2])*x[0]])

  def func_grad_nnz(self):
    return num_variables

  def cons_val_nnz(self):
    return num_constraints

  def cons_jac_nnz(self):
    return num_variables * num_constraints

  def cons_vals(self):
    x = self.v

    return np.array([x[0]*x[1]*x[2]*x[3],
                     sq(x[0]) + sq(x[1]) + sq(x[2]) +  sq(x[3])])

  def cons_jac(self):
    x = self.v

    return np.array([[x[1]*x[2]*x[3], x[0]*x[2]*x[3], x[0]*x[1]*x[3], x[0]*x[1]*x[2]],
                     [2*x[0], 2*x[1], 2*x[2], 2*x[3]]])

  def hess_prod(self, func_dual, direction, cons_dual):
    x = self.v

    p = np.array([0]*num_variables, dtype=np.float64)

    p0 = 0.

    p0 += (2*x[3]*func_dual + 2*cons_dual[1])*direction[0]
    p0 += (x[3]*func_dual + x[2]*x[3]*cons_dual[0])*direction[1]
    p0 += (x[3]*func_dual + x[1]*x[3]*cons_dual[0])*direction[2]
    p0 += ((2*x[0] + x[1] + x[2])*func_dual + x[1]*x[2]*cons_dual[0])*direction[3]

    p[0] = p0

    p1 = 0.

    p1 += (x[3]*func_dual + x[2]*x[3]*cons_dual[0])*direction[0]
    p1 += (2*cons_dual[1])*direction[1]
    p1 += (x[0]*x[3]*cons_dual[0])*direction[2]
    p1 += (x[0]*func_dual + (x[0]*x[2])*cons_dual[0])*direction[3]

    p[1] = p1

    p2 = 0.

    p2 += (x[3]*func_dual + x[1]*x[3]*cons_dual[0])*direction[0]
    p2 += (x[0]*x[3]*cons_dual[0])*direction[1]
    p2 += (2*cons_dual[1])*direction[2]
    p2 += (x[0]*func_dual + x[0]*x[1]*cons_dual[0])*direction[3]

    p[2] = p2

    p3 = 0.

    p3 += ((2*x[0] + x[1] + x[2])*func_dual + x[1]*x[2]*cons_dual[0])*direction[0]
    p3 += (x[0]*func_dual + x[0]*x[2]*cons_dual[0])*direction[1]
    p3 += (x[0]*func_dual + x[0]*x[1]*cons_dual[0])*direction[2]
    p3 += (2*cons_dual[1])*direction[3]

    p[3] = p3

    return p


class ConstrainedTest(unittest.TestCase):

  def setUp(self):
    inf = sleqp.inf()

    var_lb = np.array([1.]*num_variables)
    var_ub = np.array([5.]*num_variables)

    cons_lb = np.array([25., 40.])
    cons_ub = np.array([inf, 40.])

    x = np.array([1., 5., 5., 1.])


    self.params = sleqp.Params()

    self.options = sleqp.Options()

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
    self.solver.solve(10, 3600)

    self.assertEqual(self.solver.status, sleqp.Status.Optimal)

    s = np.array([1., 4.742999, 3.821151, 1.379408])

    self.assertTrue(np.allclose(s, self.solver.primal))


if __name__ == '__main__':
    unittest.main()
