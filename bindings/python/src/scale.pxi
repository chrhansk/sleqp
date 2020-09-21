#cython: language_level=3

class Array(np.ndarray):
  """
  Thin wrapper around a numpy array, storing a
  reference to prevent the destruction of
  an underlying object
  """
  def __new__(cls, input_array, token=None):
    obj = np.asarray(input_array).view(cls)
    obj._token = token
    return obj

  def __array_finalize__(self, obj):
    if obj is None:
      return
    self._token = getattr(obj, '_token', None)

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
    csleqp_call(csleqp.sleqp_scaling_release(&self.scaling))

  @property
  def func_weight(self):
    return csleqp.sleqp_scaling_get_func_weight(self.scaling)

  @func_weight.setter
  def func_weight(self, value):
    csleqp_call(csleqp.sleqp_scaling_set_func_weight(self.scaling,
                                                     value))

  @property
  def variable_weights(self):
    length = self.num_variables
    cdef int[:] values = <int[:length]> csleqp.sleqp_scaling_get_var_weights(self.scaling)

    array = np.asarray(values)
    array.flags.writeable = False

    return Array(array, self)

  @property
  def constraint_weights(self):
    length = self.num_constraints
    cdef int[:] values = <int[:length]> csleqp.sleqp_scaling_get_cons_weights(self.scaling)

    array = np.asarray(values)
    array.flags.writeable = False

    return Array(array, self)

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
