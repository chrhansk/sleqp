#!/usr/bin/env python

import numpy as np
import unittest

import scipy.sparse

import sleqp

num_variables = 2
num_constraints = 1

class MatrixErrorFunc(sleqp.Func):
  def set_value(self, v):
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
    return np.array([0]*num_variables)

  def hess_prod(self, func_dual, direction, cons_dual):
    return np.array([0]*num_variables)

  def cons_jac(self):
    return self.m


class MatrixErrorTest(unittest.TestCase):
  def setUp(self):
    inf = sleqp.inf()

    self.var_lb = np.array([-inf]*num_variables)
    self.var_ub = np.array([inf]*num_variables)

    self.cons_lb = np.array([0]*num_constraints)
    self.cons_ub = np.array([])*num_constraints

    self.x = np.array([0.]*num_variables)

    self.params = sleqp.Params()

  def test_string_error(self):
    func = MatrixErrorFunc(num_variables,
                           num_constraints)

    func.set_matrix_value("asd")

    problem = sleqp.Problem(func,
                            self.params,
                            self.var_lb,
                            self.var_ub,
                            self.cons_lb,
                            self.cons_ub)

    solver = sleqp.Solver(problem,
                          self.params,
                          self.x)

    with self.assertRaises(sleqp.SLEQPError):
      solver.solve(1, 3600)

  def test_wrong_shape(self):
    func = MatrixErrorFunc(num_variables,
                           num_constraints)

    func.set_matrix_value(np.array((2, 2, 2)))

    problem = sleqp.Problem(func,
                            self.params,
                            self.var_lb,
                            self.var_ub,
                            self.cons_lb,
                            self.cons_ub)

    solver = sleqp.Solver(problem,
                          self.params,
                          self.x)

    with self.assertRaises(sleqp.SLEQPError):
        solver.solve(1, 3600)

  def test_sparse_type(self):
    func = MatrixErrorFunc(num_variables,
                           num_constraints)

    m = scipy.sparse.lil_matrix((num_constraints, num_variables))

    func.set_matrix_value(m)

    problem = sleqp.Problem(func,
                            self.params,
                            self.var_lb,
                            self.var_ub,
                            self.cons_lb,
                            self.cons_ub)

    solver = sleqp.Solver(problem,
                          self.params,
                          self.x)

    solver.solve(1, 3600)
