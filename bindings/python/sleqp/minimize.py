import numpy as np

import sleqp

from sleqp._bounds import (create_constraint_bounds,
                           create_variable_bounds)

from sleqp._cons_func import create_constraint_func
from sleqp._objective import create_objective
from sleqp._hessian import create_hessian

_base_pert = np.sqrt(np.finfo(float).eps)

class _MinFunc:

  def __init__(self, objective, cons_func, hessian, args, num_variables):
    self.x = np.zeros((num_variables,))
    self.objective = objective
    self.cons_func = cons_func
    self.hessian = hessian
    self.args = args
    self.num_variables = num_variables

    self.obj_val = None
    self.obj_grad = None

  def set_value(self, v, reason):
    if (self.x == v).all():
      return

    self.x[:] = v

    self.objective.set_value(self.x)

    if self.hessian is not None:
      self.hessian.set_value(self.x)

    if self.cons_func:
      self.cons_func.set_value(self.x)

    self.obj_val = None
    self.obj_grad = None

  def cons_vals(self):
    return self.cons_func.val(self.args)

  def cons_jac(self):
    return self.cons_func.jac(self.args)

  def func_val(self):
    self._eval_obj_val()
    return self.obj_val

  def _eval_obj_val(self):
    if self.obj_val is None:
      self.obj_val = self.objective.val(self.args)

  def _eval_obj_grad(self):
    self._eval_obj_val()

    if self.obj_grad is None:
      self.obj_grad = self.objective.grad(self.args)

  def func_grad(self):
    self._eval_obj_grad()
    return self.obj_grad

  def hess_prod(self, func_dual, direction, _):
    self._eval_obj_grad()
    prod = self.hessian.product(self.args,
                                direction,
                                self.obj_grad)

    return func_dual * prod


class OptimizeResult:
  def __init__(self, solver):
    self.x = solver.solution.primal
    self.status = solver.status
    self.success = (self.status == sleqp.Status.Optimal)
    self.nit = solver.iterations
    self.fun = solver.solution.func_val

  def __getitem__(self, key):
    return getattr(self, key)


def _add_solver_callback(solver, callback):

  def accepted_iterate(solver, iterate):
    abort = callback(iterate.primal)

    if abort is True:
      solver.abort()

  solver.add_callback(sleqp.SolverEvent.AcceptedIterate,
                      accepted_iterate)


def minimize(fun, x0, args=(), jac=None, hess=None, hessp=None, bounds=None, constraints=None, callback=None):
  """
  A drop-in replacement for :func:`scipy.optimize.minimize`, minimizing a scalar function
  of one or more variables subjecting to constraints.

  :param fun:
      The objective function to be minimized.
      ``fun(x, *args) -> float``
      where ``x`` is an 1-D array with shape (n,) and ``args``
      is a tuple of the fixed parameters needed to completely
      specify the function.
  :type fun: callable
  :param x0:
        Initial guess. Array of real elements of size (n,),
        where 'n' is the number of independent variables.
  :type x0: :class:`numpy.ndarray`, shape (n,)
  :param args:
        Extra arguments passed to the objective function and its
        derivatives (`fun`, `jac` and `hess` functions).
  :type args: tuple, optional
  :param bounds:
        Bounds on variables. There are two ways to specify the bounds:
        1. Instance of :class:`scipy.optimize.Bounds` class.
        2. Sequence of ``(min, max)`` pairs for each element in `x`. ``None`` is used to specify no bound.
  :type bounds: sequence or :class:`scipy.optimize.Bounds`, optional
  :param callback:
        Called after each iteration. If callback returns True
        the algorithm execution is terminated. The signature is: ``callback(xk)``
        where ``xk`` is the current guess.
  :type callback: callable, optional
  :return: An improved guess
  :rtype: :class:`numpy.ndarray`, shape (n,)
  """
  if not isinstance(args, tuple):
    args = (args,)

  num_variables = len(x0)
  num_constraints = 0

  cons_lb = np.zeros((0,))
  cons_ub = np.zeros((0,))
  cons_func = None

  if constraints is not None:
    num_constraints = len(constraints)

    (cons_lb, cons_ub) = create_constraint_bounds(constraints)
    cons_func = create_constraint_func(num_variables, constraints)

  objective = create_objective(num_variables, fun, jac)

  options = sleqp.Options()

  hessian = None

  if hessp is None and hess is None:
    options.hessian_eval = sleqp.HessianEval.DampedBFGS
  else:
    _jac = jac
    if jac is True:
      def _jac(x, *args):
        return fun(x, *args)[1]

    hessian = create_hessian(_jac, hess, hessp)

  initial_sol = np.array(x0)

  (var_lb, var_ub) = create_variable_bounds(num_variables, bounds)

  min_func = _MinFunc(objective, cons_func, hessian, args, num_variables)

  params = sleqp.Params()

  problem = sleqp.Problem(min_func,
                          params,
                          var_lb,
                          var_ub,
                          cons_lb,
                          cons_ub)

  solver = sleqp.Solver(problem,
                        params,
                        options,
                        initial_sol)

  if callback is not None:
    _add_solver_callback(solver, callback)

  solver.solve(100, 3600)

  return OptimizeResult(solver)
