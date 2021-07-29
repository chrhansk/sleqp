#!/usr/bin/env python

import numpy as np
import unittest

from scipy.optimize import rosen, rosen_hess, rosen_der, rosen_hess_prod

from scipy.sparse import dok_matrix
from scipy.sparse.linalg import LinearOperator

import sleqp


def convert_to_sparse(array):
  (m, n) = array.shape
  mat = dok_matrix((m, n))

  for i in range(m):
    for j in range(n):
      mat[i, j] = array[i, j]

  return mat

def sparse_hessian(hess):
  def eval(x, *args):
    return convert_to_sparse(hess(x, *args))

  return eval


class MinimizeTest(unittest.TestCase):

  def setUp(self):
    self.x0 = [1.3, 0.7, 0.8, 1.9, 1.2]
    self.xopt = np.array([ 1.,  1.,  1.,  1.,  1.])

  def test_unconstrained(self):
    res = sleqp.minimize(rosen, self.x0, jac=rosen_der, hessp=rosen_hess_prod)

    self.assertTrue(res.success)
    self.assertTrue(np.allclose(res.x, self.xopt))

  def test_unconstrained_findiff(self):
    res = sleqp.minimize(rosen, self.x0, jac='3-point', hessp=rosen_hess_prod)

    self.assertTrue(res.success)
    self.assertTrue(np.allclose(res.x, self.xopt))


  def test_dense_hessian(self):
    res = sleqp.minimize(rosen, self.x0, jac=rosen_der, hess=rosen_hess)

    self.assertTrue(res.success)
    self.assertTrue(np.allclose(res.x, self.xopt))

  def test_sparse_hessian(self):
    res = sleqp.minimize(rosen, self.x0, jac=rosen_der, hess=sparse_hessian(rosen_hess))

    self.assertTrue(res.success)
    self.assertTrue(np.allclose(res.x, self.xopt))

  def test_findiff_hessian(self):
    res = sleqp.minimize(rosen, self.x0, jac=rosen_der, hess='3-point')

    self.assertTrue(res.success)
    self.assertTrue(np.allclose(res.x, self.xopt))

  def test_hessian_linear_operator(self):

    n = len(self.x0)

    def hess_op(x, *args):
      prod = lambda d: rosen_hess_prod(x, d, *args)
      return LinearOperator((n, n), matvec=prod)

    res = sleqp.minimize(rosen, self.x0, jac=rosen_der, hess=hess_op)

    self.assertTrue(res.success)
    self.assertTrue(np.allclose(res.x, self.xopt))

  def test_callback(self):
    def callback(x):
      assert x.shape == np.array(self.x0).shape

    res = sleqp.minimize(rosen, self.x0, jac=rosen_der, hessp=rosen_hess_prod, callback=callback)

    self.assertTrue(res.success)
    self.assertTrue(np.allclose(res.x, self.xopt))

  def test_callback_abort(self):
    def callback(x):
      return True

    res = sleqp.minimize(rosen, self.x0, jac=rosen_der, hessp=rosen_hess_prod, callback=callback)

    self.assertEqual(res['nit'], 1)
