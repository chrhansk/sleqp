import numpy as np
import scipy.sparse as sps

from scipy.sparse.linalg import LinearOperator

from sleqp._derivative import FD_METHODS, create_derivative


class Hessian:
  def set_value(self, x):
    pass

  def product(self, args, grad, d):
    pass


class HessianProduct(Hessian):
  def __init__(self, hessp):
    self.hessp = hessp

  def set_value(self, x):
    self.x = x

  def product(self, args, direction, grad=None):
    return self.hessp(self.x, direction, *args)


# Callable Hessian
class SparseHessian(Hessian):
  def __init__(self, x, hess, hess_val):
    self.hess = hess
    self.x = x
    self.hval = hess_val

  def set_value(self, x):
    self.x = x
    self.hval = None

  def _eval(self, args):
    if self.hval is None:
      self.hval = self.hess(self.x, *args)

  def product(self, args, direction, grad=None):
    self._eval(args)
    return self.hval.dot(direction)


# Callable Hessian
class DenseHessian(Hessian):
  def __init__(self, x, hess, hess_val):
    self.hess = hess
    self.x = x
    self.hval = self._prepare_hess_val(hess_val)

  def set_value(self, x):
    self.x = x
    self.hval = None

  def _prepare_hess_val(self, hess_val):
    return np.atleast_2d(hess_val)

  def _eval(self, args):
    if self.hval is None:
      self.hval = self._prepare_hess_val(self.hess(self.x, *args))

  def product(self, args, direction, grad=None):
    self._eval(args)
    return np.dot(self.hval, direction)


# Callable Hessian
class HessianOperator(Hessian):
  def __init__(self, x, hess, hess_val):
    self.hess = hess
    self.x = x
    self.hess_op = hess_val

  def set_value(self, x):
    self.x = x
    self.hess_op = None

  def _eval(self, args):
    if self.hess_op is None:
      self.hess_op = self.hess(self.x, *args)

  def product(self, args, direction, grad=None):
    self._eval(args)
    return self.hess_op(direction)


class FindiffHessian(Hessian):
  def __init__(self, deriv):
    self.deriv = deriv
    self.hess = None

  def set_value(self, x):
    self.x = x
    self.hess = None

  def _eval(self, grad, args):
    if self.hess is None:
      self.hess = self.deriv(self.x, args, grad)

  def product(self, args, direction, grad=None):
    self._eval(grad, args)
    return np.dot(self.hess, direction)


class GenericHessian(Hessian):
  def __init__(self, hess):
    self.hess = hess
    self.impl = None

  def product(self, args, direction, grad=None):
    hess_val = self.hess(self.x, *args)

    self.impl = _create_callable_hessian(self.x,
                                         self.hess,
                                         hess_val)

    def set_value(x):
      self.impl.set_value(x)

    def product(args, direction, grad=None):
      return self.impl.product(args, direction, grad)

    self.set_value = set_value
    self.product = product

    return self.impl.product(args, direction, grad)

  def set_value(self, x):
    self.x = x


def _create_callable_hessian(x, hess, hess_val):
  if sps.issparse(hess_val):
    return SparseHessian(x, hess, hess_val)

  elif isinstance(hess_val, LinearOperator):
    return HessianOperator(x, hess, hess_val)
  else:
    return DenseHessian(x, hess, hess_val)

def create_hessian(grad, hess, hessp):

  if hess and hessp:
    raise ValueError("Specified both Hessian and Hessian product")

  if hessp:
    return HessianProduct(hessp)

  if hess in FD_METHODS:
    deriv = create_derivative(grad, hess)
    return FindiffHessian(deriv)

  elif callable(hess):
    return GenericHessian(hess)

  raise ValueError("Invalid Hessian")
