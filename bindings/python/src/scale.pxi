#cython: language_level=3

cdef class Scaling:
  cdef dict __dict__
  cdef csleqp.SleqpScalingData* scaling
  cdef csleqp.SleqpSparseVec* gradient
  cdef csleqp.SleqpSparseMatrix* cons_jac

  cdef int num_variables
  cdef int num_constraints

  def __cinit__(self,
                Problem problem,
                Params params):
    csleqp_call(csleqp.sleqp_scaling_create(&self.scaling,
                                            problem.problem,
                                            params.params))

    csleqp_call(csleqp.sleqp_sparse_vector_create(&self.gradient,
                                                  problem.num_variables,
                                                  0))

    csleqp_call(csleqp.sleqp_sparse_matrix_create(&self.cons_jac,
                                                  problem.num_constraints,
                                                  problem.num_variables,
                                                  0))

    self.num_variables = problem.num_variables
    self.num_constraints = problem.num_constraints


  def __dealloc__(self):
    csleqp_call(csleqp.sleqp_sparse_matrix_release(&self.cons_jac))
    csleqp_call(csleqp.sleqp_sparse_vector_free(&self.gradient))
    csleqp_call(csleqp.sleqp_scaling_free(&self.scaling))

  def set_func_weight(self, int weight):
    csleqp_call(csleqp.sleqp_scaling_set_func_weight(self.scaling,
                                                     weight))

  def set_func_weight_from_nominal(self, nominal_value):
    csleqp_call(csleqp.sleqp_scaling_set_func_weight_from_nominal(self.scaling,
                                                                  nominal_value))

  def set_variable_weight(self, int index, int weight):
    csleqp_call(csleqp.sleqp_scaling_set_var_weight(self.scaling,
                                                    index,
                                                    weight))

  def set_variable_weights_from_nominal(self, nominal_array):
    assert nominal_array.shape == (self.num_variables,)
    cdef double[:] nominal_values = nominal_array
    csleqp_call(csleqp.sleqp_scaling_set_var_weights_from_nominal(self.scaling,
                                                                  &nominal_values[0]))

  def set_constraint_weight(self, int index, int weight):
    csleqp_call(csleqp.sleqp_scaling_set_cons_weight(self.scaling,
                                                     index,
                                                     weight))

  def set_constraint_weights_from_nominal(self, nominal_array):
    assert nominal_array.shape == (self.num_constraints,)
    cdef double[:] nominal_values = nominal_array
    csleqp_call(csleqp.sleqp_scaling_set_cons_weights_from_nominal(self.scaling,
                                                                   &nominal_values[0]))

  def set_from_gradient(self,
                        np.ndarray grad_array):
    csleqp_call(array_to_sleqp_sparse_vec(grad_array,
                                          self.gradient))

    csleqp_call(csleqp.sleqp_func_scaling_from_gradient(self.scaling,
                                                        self.gradient))

  def set_from_cons_jac(self,
                        object mat):
    csleqp_call(matrix_to_sleqp_sparse_matrix(mat,
                                              self.cons_jac))

    csleqp_call(csleqp.sleqp_scaling_from_cons_jac(self.scaling,
                                                   self.cons_jac))
