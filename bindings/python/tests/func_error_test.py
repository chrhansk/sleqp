#!/usr/bin/env python

import numpy as np
import unittest

import sleqp

num_variables = 2
num_constraints = 0

inner = Exception("Inner exception")


class ErrorFunc:
  def set_value(self, v, reason):
    raise inner


class TypeErrorFunc:
  def set_value(self, v):
    pass

  def obj_val(self):
    return 0

  def obj_grad(self):
    return "wrong"

  def hess_prod(self, obj_dual, direction, cons_dual):
    return "wrong"


class MatrixErrorFunc:
  def set_value(self, v, reason):
    pass

  def set_matrix_value(self, m):
    self.m = m

  def obj_val(self):
    return 0

  def obj_grad_nnz(self):
    return 1

  def cons_jac_nnz(self):
    return 1

  def obj_grad(self):
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

    with self.assertRaises(Exception):
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

    failed = False

    try:
      self.solver.solve()
    except Exception as err:
      eval_error = err.__cause__
      orig_error = eval_error.__cause__
      self.assertTrue(isinstance(eval_error, sleqp.EvaluationError))
      self.assertEqual(orig_error, inner)

      failed = True

    self.assertTrue(failed)

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

    with self.assertRaises(Exception):
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

    with self.assertRaises(Exception):
      solver.solve()

if __name__ == '__main__':
    unittest.main()
