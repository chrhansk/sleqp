from dataclasses import dataclass

import numpy as np
import scipy.sparse


@dataclass
class LinearConstraints:
  lb: object = None
  coeffs: object = None
  ub: object = None


def convert_to_sparse(array):
  (m, n) = array.shape

  mat = scipy.sparse.dok_matrix((m, n))

  for i in range(m):
    for j in range(n):
      value = array[i, j]

      if value:
        mat[i, j] = value

  return mat


def create_linear_constraints(constraints):

  if not(constraints):
    return LinearConstraints()

  linear_lb = np.hstack([constraint.lb for constraint in constraints])

  linear_ub = np.hstack([constraint.ub for constraint in constraints])

  matrices = []

  for constraint in constraints:
    matrix = constraint.A
    if not(scipy.sparse.issparse(matrix)):
      matrix = convert_to_sparse(matrix)

    matrices.append([matrix])

  matrix = scipy.sparse.bmat(matrices)

  return LinearConstraints(lb=linear_lb, coeffs=matrix, ub=linear_ub)
