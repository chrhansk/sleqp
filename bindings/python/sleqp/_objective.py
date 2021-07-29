import numpy as np

from sleqp._func import create_func


class Objective:

  def __init__(self, num_variables, fun, jac):
    self.func = create_func(fun, jac)
    self.num_variables = num_variables

  def set_value(self, x):
    self.func.set_value(x)

  def val(self, args=()):
    return self.func.val(args)

  def grad(self, args=()):
    grad = self.func.grad(args)
    return grad.reshape((self.num_variables,))


def create_objective(num_variables, fun, jac):
  return Objective(num_variables, fun, jac)
