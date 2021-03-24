import numpy as np
import scipy.optimize

import sleqp

_base_pert = np.sqrt(np.finfo(float).eps)

class _MinFunc:

  def __init__(self, fun, grad, hessp, args, num_variables):
    self.x = np.zeros((num_variables,))
    self.fun = fun
    self.grad = grad
    self.hessp = hessp
    self.args = args
    self.num_variables = num_variables
    self._func_val = None

  def set_value(self, v, reason):
    self.x[:] = v
    self._func_val = None

  def func_val(self):
    if self._func_val is None:
      self._func_val = self.fun(self.x, *self.args)

    return self._func_val

  def func_grad(self):
    if self.grad is not None:
      return self.grad(self.x, *self.args)

    grad = np.empty_like(self.x)

    x = self.x
    xd = np.copy(self.x)

    orig_func_val = self.func_val()

    for i in range(self.num_variables):

      pert = _base_pert = max(abs(x[i]), 1.)
      xd[i] += pert

      pert_func_val = self.fun(xd, *self.args)

      grad[i] = (pert_func_val - orig_func_val) / pert

      xd[i] = self.x[i]

    return grad

  def hess_prod(self, func_dual, direction, _):
    prod = self.hessp(self.x,
                      direction,
                      *self.args)

    return func_dual * prod


class Result:
  def __init__(self, x):
    self.x = x

def _create_variable_bounds(num_variables, bounds):
  inf = sleqp.inf()

  var_lb = np.full((num_variables,), -inf)
  var_ub = np.full((num_variables,), inf)

  if bounds is None:
    return (var_lb, var_ub)

  if isinstance(bounds, scipy.optimize.Bounds):
    var_lb = np.array(bounds.lb)
    var_ub = np.array(bounds.ub)
    return (var_lb, var_ub)

  for i, (lb, ub) in enumerate(bounds):
    if lb is not None:
      var_lb[i] = lb
    if ub is not None:
      var_ub[i] = ub

  return (var_lb, var_ub)

def minimize(fun, x0, args=(), grad=None, hessp=None, bounds=None):

  num_variables = len(x0)
  num_constraints = 0

  inf = sleqp.inf()

  initial_sol = np.array(x0)

  (var_lb, var_ub) = _create_variable_bounds(num_variables,
                                             bounds)

  cons_lb = np.full((num_constraints,), -inf)
  cons_ub = np.full((num_constraints,), inf)

  min_func = _MinFunc(fun, grad, hessp, args, num_variables)

  problem = sleqp.Problem(min_func,
                          var_lb,
                          var_ub,
                          cons_lb,
                          cons_ub)

  params = sleqp.Params()

  options = sleqp.Options()

  if hessp is None:
    options.hessian_eval = sleqp.HessianEval.SR1

  solver = sleqp.Solver(problem,
                        params,
                        options,
                        initial_sol)

  solver.solve(100, 3600)

  return Result(solver.solution.primal)
