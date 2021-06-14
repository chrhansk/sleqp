#!/usr/bin/env python

import numpy as np
import unittest

import sleqp

num_variables = 2
num_constraints = 0

inner = Exception("Error in set_value")

class ErrorFunc:
  def set_value(self, v, reason):
    raise inner

class TypeErrorFunc:
  def set_value(self, v):
    pass

  def func_val(self):
    return 0

  def func_grad(self):
    return "wrong"

  def hess_prod(self, func_dual, direction, cons_dual):
    return "wrong"

class MatrixErrorFunc:
  def set_value(self, v, reason):
    pass

  def set_matrix_value(self, m):
    self.m = m

  def func_val(self):
    return 0

  def func_grad_nnz(self):
    return 1

  def cons_jac_nnz(self):
    return 1

  def func_grad(self):
    return np.array([0])

  def hess_prod(self, func_dual, direction, cons_dual):
    return np.array([0])

  def cons_jac(self):
    return self.m



class FuncErrorTest(unittest.TestCase):

  def setUp(self):
    inf = sleqp.inf()

    self.var_lb = np.array([-inf, -inf])
    self.var_ub = np.array([inf, inf])

    self.cons_lb = np.array([])
    self.cons_ub = np.array([])

    self.x = np.array([0., 0.])

    self.params = sleqp.Params()
    self.options = sleqp.Options()

  def test_error_func(self):
    func = ErrorFunc()

    self.problem = sleqp.Problem(func,
                                 self.params,
                                 self.var_lb,
                                 self.var_ub,
                                 self.cons_lb,
                                 self.cons_ub)

    self.solver = sleqp.Solver(self.problem,
                               self.params,
                               self.options,
                               self.x)

    with self.assertRaises(sleqp.SLEQPError):
      self.solver.solve()

  def test_error_chain(self):
    func = ErrorFunc()

    self.problem = sleqp.Problem(func,
                                 self.params,
                                 self.var_lb,
                                 self.var_ub,
                                 self.cons_lb,
                                 self.cons_ub)

    self.solver = sleqp.Solver(self.problem,
                               self.params,
                               self.options,
                               self.x)

    try:
      self.solver.solve()
    except sleqp.SLEQPError as err:
      self.assertEqual(err.__cause__, inner)

  def test_type_error_func(self):
    func = TypeErrorFunc()

    problem = sleqp.Problem(func,
                            self.params,
                            self.var_lb,
                            self.var_ub,
                            self.cons_lb,
                            self.cons_ub)

    solver = sleqp.Solver(problem,
                          self.params,
                          self.options,
                          self.x)

    with self.assertRaises(sleqp.SLEQPError):
      solver.solve()

  def test_matrix_error_func(self):
    func = MatrixErrorFunc()

    problem = sleqp.Problem(func,
                            self.params,
                            self.var_lb,
                            self.var_ub,
                            self.cons_lb,
                            self.cons_ub)

    solver = sleqp.Solver(problem,
                          self.params,
                          self.options,
                          self.x)

    with self.assertRaises(sleqp.SLEQPError):
      solver.solve()

if __name__ == '__main__':
    unittest.main()
