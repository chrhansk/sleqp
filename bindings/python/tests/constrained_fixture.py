import numpy as np

import sleqp

num_variables = 4
num_constraints = 2

inf = sleqp.inf()

def sq(x):
  return x*x

class ConstrainedFunc:

  def set_value(self, v, reason):
    self.v = v

  def obj_val(self):
    x = self.v

    return x[0]*x[3]*(x[0] + x[1]+ x[2]) + x[2]

  def obj_grad(self):
    x = self.v

    return np.array([(x[0] + x[1] + x[2])*x[3] + x[0]*x[3],
                     x[0]*x[3],
                     x[0]*x[3] + 1,
                     (x[0] + x[1] + x[2])*x[0]])

  def obj_grad_nnz(self):
    return num_variables

  def cons_val_nnz(self):
    return num_constraints

  def cons_jac_nnz(self):
    return num_variables * num_constraints

  def cons_vals(self):
    x = self.v

    return np.array([x[0]*x[1]*x[2]*x[3],
                     sq(x[0]) + sq(x[1]) + sq(x[2]) +  sq(x[3])])

  def cons_jac(self):
    x = self.v

    return np.array([[x[1]*x[2]*x[3], x[0]*x[2]*x[3], x[0]*x[1]*x[3], x[0]*x[1]*x[2]],
                     [2*x[0], 2*x[1], 2*x[2], 2*x[3]]])

  def hess_prod(self, direction, cons_dual):
    x = self.v

    p = np.array([0]*num_variables, dtype=np.float64)

    p0 = 0.

    p0 += (2*x[3] + 2*cons_dual[1])*direction[0]
    p0 += (x[3] + x[2]*x[3]*cons_dual[0])*direction[1]
    p0 += (x[3] + x[1]*x[3]*cons_dual[0])*direction[2]
    p0 += ((2*x[0] + x[1] + x[2]) + x[1]*x[2]*cons_dual[0])*direction[3]

    p[0] = p0

    p1 = 0.

    p1 += (x[3] + x[2]*x[3]*cons_dual[0])*direction[0]
    p1 += (2*cons_dual[1])*direction[1]
    p1 += (x[0]*x[3]*cons_dual[0])*direction[2]
    p1 += (x[0] + (x[0]*x[2])*cons_dual[0])*direction[3]

    p[1] = p1

    p2 = 0.

    p2 += (x[3] + x[1]*x[3]*cons_dual[0])*direction[0]
    p2 += (x[0]*x[3]*cons_dual[0])*direction[1]
    p2 += (2*cons_dual[1])*direction[2]
    p2 += (x[0] + x[0]*x[1]*cons_dual[0])*direction[3]

    p[2] = p2

    p3 = 0.

    p3 += ((2*x[0] + x[1] + x[2]) + x[1]*x[2]*cons_dual[0])*direction[0]
    p3 += (x[0] + x[0]*x[2]*cons_dual[0])*direction[1]
    p3 += (x[0] + x[0]*x[1]*cons_dual[0])*direction[2]
    p3 += (2*cons_dual[1])*direction[3]

    p[3] = p3

    return p



var_lb = np.array([1.]*num_variables)
var_ub = np.array([5.]*num_variables)

cons_lb = np.array([25., 40.])
cons_ub = np.array([inf, 40.])

initial_sol = np.array([1., 5., 5., 1.])
expected_sol = np.array([1., 4.742999, 3.821151, 1.379408])
