#!/usr/bin/env python

import numpy as np

import sleqp

num_variables = 2
num_constraints = 0

class RosenbrockFunc(sleqp.Func):

  def __init__(self, num_variables, num_constraints):
    self.a = 1
    self.b = 1

  def set_value(self, v, reason):
    self.v = v

  def func_val(self):
    [x, y] = self.v
    (a, b) = (self.a, self.b)

    xsq = x**2

    return (a - x)**2 + b*(y - xsq)**2

  def func_grad_nnz(self):
    return 2

  def func_grad(self):
    [x, y] = self.v
    (a, b) = (self.a, self.b)

    xsq = x**2

    g = np.array([(4*b*x*(xsq - y)) + 2*x - 2*a,
                  -2*b*(xsq - y)])

    return g

  def cons_vals(self):
    return np.zeros((num_constraints,))

  def cons_jac(self):
    return np.zeros((num_constraints, num_variables))

  def hess_prod(self, func_dual, direction, _):
    [x, y] = self.v
    (a, b) = (self.a, self.b)
    [dx, dy] = direction

    xsq = x**2

    product = np.array([((8.*b*xsq + 4.*b*(xsq - y) + 2.)*dx - (4.*b*x)*dy)*func_dual,
                        ((-4.*b*x)*dx + (2.*b)*dy)*func_dual])

    return product
