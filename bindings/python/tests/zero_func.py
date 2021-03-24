import numpy as np

import sleqp

num_variables = 4
num_constraints = 2

initial_sol = np.array([1., 0., 2., 0.])

class ZeroFunc:

    def set_value(self, values, reason):
        assert((values == initial_sol).all())

    def func_val(self):
        return 0.

    def func_grad(self):
        return np.zeros((num_variables,))

    def cons_vals(self):
        return np.zeros((num_constraints,))

    def cons_jac(self):
        return np.zeros((num_constraints, num_variables))

    def hess_prod(self, func_dual, direction, cons_dual):
        return np.zeros((num_variables,))
