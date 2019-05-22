#cython: language_level=3

cdef class Scaling:
  cdef csleqp.SleqpScalingData* scaling
  cdef csleqp.SleqpSparseVec* gradient
  cdef csleqp.SleqpSparseMatrix* cons_jac

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


  def __dealloc__(self):
    csleqp_call(csleqp.sleqp_sparse_matrix_free(&self.cons_jac))
    csleqp_call(csleqp.sleqp_sparse_vector_free(&self.gradient))
    csleqp_call(csleqp.sleqp_scaling_free(&self.scaling))

  def set_func_weight(self, int weight):
    csleqp_call(csleqp.sleqp_scaling_set_func_weight(self.scaling,
                                                     weight))

  def set_variable_weight(self, int index, int weight):
    csleqp_call(csleqp.sleqp_scaling_set_var_weight(self.scaling,
                                                    index,
                                                    weight))

  def set_constraint_weight(self, int index, int weight):
    csleqp_call(csleqp.sleqp_scaling_set_cons_weight(self.scaling,
                                                     index,
                                                     weight))

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
