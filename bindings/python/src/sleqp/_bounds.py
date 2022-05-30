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

  inf = sleqp.inf()
  cons_lb = []
  cons_ub = []

  for i, constraint in enumerate(constraints):

    try:
      cons_lb.append(constraint.lb)
      cons_ub.append(constraint.ub)
      continue
    except AttributeError:
      pass

    if constraint['type'] == 'ineq':
      cons_lb.append(np.array([0.]))
      cons_ub.append(np.array([inf]))
    else:
      assert constraint['type'] == 'eq'
      cons_lb.append(np.array([0.]))
      cons_ub.append(np.array([0.]))

  if cons_lb:
    cons_lb = np.hstack(cons_lb)
  else:
    cons_lb = np.array([])

  if cons_ub:
    cons_ub = np.hstack(cons_ub)
  else:
    cons_ub = np.array([])

  assert cons_lb.shape == cons_ub.shape

  assert (cons_lb <= cons_ub).all()

  return (cons_lb, cons_ub)
