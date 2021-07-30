import numpy as np

from sleqp._func import create_func


class ConstraintFunc:

  def __init__(self, num_variables, constraints):
    self.num_variables = num_variables
    self.funcs = self._create_funcs(constraints)

    self.num_constraints = sum(func.dimension for func in self.funcs)

    self.shape = (self.num_constraints,
                  self.num_variables)

  def _create_func(self, constraint):

      try:
        fun = constraint.fun
        jac = constraint.jac
        hess = constraint.hess
        dim = len(constraint.lb)
        return create_func(fun, jac, hess=hess, dimension=dim)

      except AttributeError:
        pass

      fun = constraint['fun']
      jac = constraint.get('jac')
      hess = constraint.get('hess')
      hessp = constraint.get('hessp')

      return create_func(fun, jac, hess, hessp)

  def _create_funcs(self, constraints):
    funcs = []

    for constraint in constraints:
      funcs.append(self._create_func(constraint))

    return funcs

  def set_value(self, x):
    self.x = x

    for func in self.funcs:
      func.set_value(self.x)

  def val(self, args=()):
    values = []

    for i, func in enumerate(self.funcs):
      values.append(func.val(args))

    if values:
      values = np.hstack(values)
    else:
      values = np.array([])

    assert values.shape == (self.num_constraints,)

    return values

  def jac(self, args=()):
    jac = []

    for i, func in enumerate(self.funcs):
      func_grad = func.grad(args)
      assert func_grad.shape == (func.dimension, self.num_variables)
      jac.append(func_grad)

    if jac:
      jac = np.hstack(jac)
    else:
      return np.zeros(self.shape)

    return jac

  def hess_prod(self, direction, duals, args=()):
    assert direction.shape == (self.num_variables,)
    assert duals.shape == (self.num_constraints,)

    product = np.zeros_like(direction)

    for i, func in enumerate(self.funcs):
      if duals[i] == 0.:
        continue

      product += duals[i] * func.hess_prod(direction, args)

    return product


def create_constraint_func(num_variables,
                           constraints):
  return ConstraintFunc(num_variables, constraints)
