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


def create_simple_func(fun, jac):
  if jac is True:
    return CombinedFunc(fun)

  deriv = create_derivative(fun, jac)

  return PointFunc(fun, deriv)


class ConstraintFunc:

  def __init__(self, num_variables, constraints):
    self.num_constraints = len(constraints)
    self.num_variables = num_variables
    self.funcs = self._create_simple_funcs(constraints)

    self.shape = (self.num_constraints,
                  self.num_variables)


  def _create_simple_funcs(self, constraints):
    funcs = []

    for constraint in constraints:
      fun = constraint['fun']
      jac = constraint.get('jac')

      funcs.append(create_simple_func(fun, jac))

    return funcs

  def set_value(self, x):
    self.x = np.asarray(x)
    for func in self.funcs:
      func.set_value(self.x)

  def val(self, args=()):
    values = np.empty((self.num_constraints,))

    for i, func in enumerate(self.funcs):
      values[i] = func.val(args)

    return values

  def jac(self, args=()):
    jac = np.empty(self.shape)

    for i, func in enumerate(self.funcs):
      jac[i, :] = func.grad(args)

    return jac


def create_constraint_func(num_variables,
                           constraints):
  return ConstraintFunc(num_variables, constraints)
