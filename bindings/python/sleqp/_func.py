import numpy as np

from sleqp._derivative import create_derivative


class PointFunc:
  def __init__(self, fun, deriv):
    self.x = None
    self.fun = fun
    self.deriv = deriv

  def set_value(self, x):
    self.x = x
    self.fval = None
    self.gval = None

  def val(self, args=()):
    if self.fval is None:
      self.fval = self.fun(self.x, *args)
    return self.fval

  def grad(self, args=()):
    if self.gval is None:
      self.gval = self.deriv(self.x, args, self.fval)
    return self.gval

class CombinedFunc:
  def __init__(self, fun):
    self.x = None
    self.fun = fun

  def set_value(self, x):
    self.x = x
    self.fval = None
    self.gval = None

  def _eval(self, x, args=()):
    (self.fval, self.gval) = self.fun(x, *args)

  def val(self, args=()):
    if self.fval is None:
      self._eval(self.x, *args)

    return self.fval

  def grad(self, args=()):
    if self.gval is None:
      self._eval(self,x, *args)

    return self.gval


def create_func(fun, jac):
  if jac is True:
    return CombinedFunc(fun)

  deriv = create_derivative(fun, jac)

  return PointFunc(fun, deriv)
