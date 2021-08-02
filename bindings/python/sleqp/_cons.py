from dataclasses import dataclass

import numpy as np

from scipy.optimize import (LinearConstraint,
                            NonlinearConstraint)

from sleqp._bounds import create_constraint_bounds
from sleqp._cons_func import create_constraint_func

from sleqp._linear import (create_linear_constraints, LinearConstraints)


@dataclass
class GeneralConstraints:
  lb: object = np.zeros((0,))
  func: object = None
  ub: object = np.zeros((0,))


@dataclass
class Constraints:
  general: object = GeneralConstraints()
  linear: object = LinearConstraints()


def split_linear_constraints(constraints):
  general = []
  linear = []

  for constraint in constraints:
    if isinstance(constraint, LinearConstraint):
      linear.append(constraint)
    else:
      general.append(constraint)

  return (general, linear)


def create_general_constraints(num_variables, constraints):
  (cons_lb, cons_ub) = create_constraint_bounds(constraints)

  cons_func = create_constraint_func(num_variables, constraints)

  return GeneralConstraints(cons_lb, cons_func, cons_ub)


def create_constraints(num_variables, constraints):

  if constraints is None:
    return Constraints()

  if not (isinstance(constraints, list) or
          isinstance(constraints, tuple)):
    constraints = [constraints]

  (general, linear) = split_linear_constraints(constraints)

  general = create_general_constraints(num_variables, general)
  linear = create_linear_constraints(linear)

  return Constraints(general=general, linear=linear)
