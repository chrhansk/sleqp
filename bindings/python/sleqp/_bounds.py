import numpy as np

import sleqp


def create_variable_bounds(num_variables, bounds):
  inf = sleqp.inf()

  var_lb = np.full((num_variables,), -inf)
  var_ub = np.full((num_variables,), inf)

  if bounds is None:
    return (var_lb, var_ub)

  try:
    var_lb = np.array(bounds.lb)
    var_ub = np.array(bounds.ub)
    return (var_lb, var_ub)
  except AttributeError:
    pass

  for i, (lb, ub) in enumerate(bounds):
    if lb is not None:
      var_lb[i] = lb
    if ub is not None:
      var_ub[i] = ub

  return (var_lb, var_ub)


def create_constraint_bounds(constraints):
  num_constraints = len(constraints)

  inf = sleqp.inf()
  cons_lb = np.full((num_constraints,), -inf)
  cons_ub = np.full((num_constraints,), inf)

  for i, constraint in enumerate(constraints):
    if constraint['type'] == 'ineq':
      cons_lb[i] = 0.
      cons_ub[i] = inf
    else:
      assert constraint['type'] == 'eq'
      cons_lb[i] = 0.
      cons_ub[i] = 0.

  return (cons_lb, cons_ub)
