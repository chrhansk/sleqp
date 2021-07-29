import numpy as np

from sleqp._func import create_func


class ConstraintFunc:

  def __init__(self, num_variables, constraints):
    self.num_constraints = len(constraints)
    self.num_variables = num_variables
    self.funcs = self._create_funcs(constraints)

    self.shape = (self.num_constraints,
                  self.num_variables)


  def _create_funcs(self, constraints):
    funcs = []

    for constraint in constraints:
      fun = constraint['fun']
      jac = constraint.get('jac')

      funcs.append(create_func(fun, jac))

    return funcs

  def set_value(self, x):
    self.x = x

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
