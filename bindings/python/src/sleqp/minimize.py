import logging

import numpy as np
from scipy.optimize import OptimizeResult

import sleqp

from sleqp._bounds import create_variable_bounds
from sleqp._cons import create_constraints
from sleqp._func import create_func
from sleqp._hessian import create_hessian


class _MinFunc:

    def __init__(self, objective, constraints, args, num_variables):
        self.x = np.zeros((num_variables,))
        self.objective = objective
        self.constraints = constraints
        self.args = args
        self.num_variables = num_variables

        self._obj_val = None
        self._obj_grad = None

    def set_value(self, v, reason):
        self.x[:] = v

        self.objective.set_value(self.x)

        if self.constraints:
            self.constraints.set_value(self.x)

        self._obj_val = None
        self._obj_grad = None

    def cons_vals(self):
        return self.constraints.val(self.args)

    def cons_jac(self):
        return self.constraints.jac(self.args)

    def obj_val(self):
        self._eval_obj_val()
        return self._obj_val

    def _eval_obj_val(self):
        if self._obj_val is None:
            self._obj_val = self.objective.val(self.args)

    def _eval_obj_grad(self):
        self._eval_obj_val()

        if self._obj_grad is None:
            self._obj_grad = self.objective.grad(self.args)

            if self._obj_grad.ndim == 2:
                self._obj_grad = self._obj_grad[0, :]

    def obj_grad(self):
        self._eval_obj_grad()
        return self._obj_grad

    def hess_prod(self, direction, cons_duals):
        self._eval_obj_grad()

        prod = self.objective.hess_prod(direction,
                                        self.args)

        if self.constraints:
            cons_prod = self.constraints.hess_prod(direction,
                                                   cons_duals,
                                                   self.args)
            prod += cons_prod

        return prod


def _create_result(solver):
    solution = solver.solution

    result = OptimizeResult()

    result["mult_g"] = solution.cons_dual
    result["mult_x"] = solution.vars_dual
    result["x"] = solution.primal
    result["success"] = (solver.status == sleqp.Status.Optimal)
    result["status"] = solver.status.value
    result["message"] = solver.status.desc

    result["fun"] = solution.obj_val
    result["jac"] = solution.obj_grad
    result["nit"] = solver.iterations

    result["maxcv"] = solver.states[sleqp.SolverState.ScaledFeasRes]

    return result


def _add_solver_callback(solver, callback):

    def accepted_iterate(solver, iterate):
        abort = callback(iterate.primal)

        if abort is True:
            solver.abort()

    solver.add_callback(sleqp.SolverEvent.AcceptedIterate,
                        accepted_iterate)


_handler = None
_attached = False


def _attach_handler():
    global _handler, _attached

    if _handler is None:
        _handler = logging.StreamHandler()
        _handler.setLevel(logging.INFO)

    if not (_attached):
        sleqp.logger.addHandler(_handler)
        _attached = True


def _detach_handler():
    global _handler, _attached

    if _handler is None:
        return

    if _attached:
        sleqp.logger.removeHandler(_handler)
        _attached = False


class LogHelper:
    def __init__(self, verbose):
        self.verbose = verbose
        self.old_level = None

    def __enter__(self):
        if not (self.verbose):
            return

        # Set level to at least INFO
        if sleqp.logger.level < logging.INFO:
            self.old_level = sleqp.logger.level
            sleqp.logger.setLevel(logging.INFO)

        _attach_handler()

    def __exit__(self, exc_type, exc_value, exc_tb):
        if self.verbose:
            _detach_handler()

        # Restore old level
        if self.old_level is not None:
            sleqp.logger.setLevel(self.old_level)


def minimize(fun, x0, args=(), jac=None, hess=None, hessp=None, bounds=None, constraints=None, callback=None, **options):
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

    :return: The optimization result
    :rtype: :class:`scipy.optimize.Optimizeresult`
    """
    if not isinstance(args, tuple):
        args = (args,)

    num_variables = len(x0)

    cons_lb = np.zeros((0,))
    cons_ub = np.zeros((0,))
    cons_func = None

    constraints = create_constraints(num_variables, constraints)

    objective = create_func(fun, jac, hess, hessp)

    verbose = options.pop('verbose', False)

    settings = sleqp.Settings(**options)

    if hessp is None and hess is None:
        settings.hessian_eval = sleqp.HessianEval.DampedBFGS

    initial_sol = np.array(x0)

    (var_lb, var_ub) = create_variable_bounds(num_variables, bounds)

    min_func = _MinFunc(objective,
                        constraints.general.func,
                        args,
                        num_variables)

    problem = sleqp.Problem(min_func,
                            var_lb,
                            var_ub,
                            constraints.general.lb,
                            constraints.general.ub,
                            settings,
                            linear_coeffs=constraints.linear.coeffs,
                            linear_lb=constraints.linear.lb,
                            linear_ub=constraints.linear.ub)

    solver = sleqp.Solver(problem,
                          initial_sol)

    if callback is not None:
        _add_solver_callback(solver, callback)

    with LogHelper(verbose):
        solver.solve()

    return _create_result(solver)
