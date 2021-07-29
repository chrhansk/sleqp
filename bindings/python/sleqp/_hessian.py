import numpy as np
import scipy.sparse as sps

from scipy.sparse.linalg import LinearOperator

from sleqp._derivative import create_derivative


class Hessian:
  def set_value(self, x):
    pass


  def product(self, args, d):
    pass


class HessianProduct:
  def __init__(self, hessp):
    self.hessp = hessp

  def set_value(self, x):
    self.x = x

  def product(self, args, d):
    return self.hessp(self.x, d, *args)


class SparseHessian:
  def __init__(self, hess):
    self.hess = hess
    self.hval = None

  def set_value(self, x):
    self.x = x
    self.hval = None

  def _eval(self, args):
    if self.hval is None:
      self.hval = self.hess(self.x, *args)

  def product(self, args, d):
    self._eval(args)
    return self.hval.dot(d)


class DenseHessian:
  def __init__(self, hess):
    self.hess = hess
    self.hval = None

  def set_value(self, x):
    self.x = x
    self.hval = None

  def _eval(self, args):
    if self.hval is None:
      self.hval = np.atleast_2d(self.hess(self.x, *args))

  def product(self, args, d):
    self._eval(args)
    return np.dot(self.hval, d)


class HessianOperator:
  def __init__(self, hess):
    self.hess = hess

  def set_value(self, x):
    self.x = x
    self.hess_op = None

  def _eval(self, args):
    if self.hess_op is None:
      self.hess_op = self.hess(self.x, *args)

  def product(self, args, d):
    self._eval(args)
    return self.hess_op(d)


def create_hessian(x0, args, grad, hess, hessp):

  if hess and hessp:
    raise ValueError("Specified both Hessian and Hessian product")

  if hessp:
    return HessianProduct(hessp)

  elif callable(hess):

    hessian = hess(x0, *args)

    if sps.issparse(hessian):
      return SparseHessian(hess)

    elif isinstance(hessian, LinearOperator):
      return HessianOperator(hess)
    else:
      return DenseHessian(hess)

  raise ValueError("Invalid Hessian")
