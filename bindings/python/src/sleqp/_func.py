import numpy as np

from sleqp._derivative import create_derivative
from sleqp._hessian import create_hessian

class PointFunc:
  def __init__(self, fun, deriv, hessian, dimension=1):
    self.x = None
    self.fun = fun
    self.deriv = deriv
    self.hessian = hessian
    self.dimension = dimension

  def set_value(self, x):
    self.x = x
    self.fval = None
    self.gval = None

    if self.hessian:
      self.hessian.set_value(x)

  def val(self, args=()):
    if self.fval is None:
      self.fval = self.fun(self.x, *args)
    return self.fval

  def grad(self, args=()):
    if self.gval is None:
      self.gval = self.deriv(self.x, args, self.fval)
      self.gval = np.atleast_1d(self.gval)
    return self.gval

  def hess_prod(self, direction, args=()):
    return self.hessian.product(args, direction)


class CombinedFunc:
  def __init__(self, fun, hessian, dimension=1):
    self.x = None
    self.fun = fun
    self.hessian = hessian
    self.dimension = dimension

  def set_value(self, x):
    self.x = x
    self.fval = None
    self.gval = None

    if self.hessian:
      self.hessian.set_value(x)

  def _eval(self, x, args=()):
    (self.fval, self.gval) = self.fun(x, *args)
    self.gval = np.atleast_1d(self.gval)

  def val(self, args=()):
    if self.fval is None:
      self._eval(self.x, *args)

    return self.fval

  def grad(self, args=()):
    if self.gval is None:
      self._eval(self,x, *args)

    return self.gval

  def hess_prod(self, direction, args=()):
    prod = self.hessian.product(args, direction)
    return np.atleast_1d(prod)


def actual_gradient(fun, grad):
  if grad is True:
    def _grad(x, *args):
      return fun(x, *args)[1]

    return _grad

  return grad

def create_func(fun, grad, hess=None, hessp=None, dimension=1):

  actual_grad = actual_gradient(fun, grad)
  hessian = None

  if hess is not None or hessp is not None:
    hessian = create_hessian(actual_grad, hess, hessp)

  if grad is True:
    return CombinedFunc(fun, hessian, dimension)

  deriv = create_derivative(fun, grad)

  return PointFunc(fun, deriv, hessian, dimension)
