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

  def __cinit__(self,
                int num_variables,
                int num_constraints,
                *args,
                **keywords):

    csleqp_call(csleqp.sleqp_scaling_create(&self.scaling,
                                            num_variables,
                                            num_constraints))

    csleqp_call(csleqp.sleqp_sparse_vector_create_empty(&self.gradient,
                                                        num_variables))

    csleqp_call(csleqp.sleqp_sparse_matrix_create(&self.cons_jac,
                                                  num_constraints,
                                                  num_variables,
                                                  0))

  @property
  def num_variables(self):
    return csleqp.sleqp_scaling_get_num_variables(self.scaling)

  @property
  def num_constraints(self):
    return csleqp.sleqp_scaling_get_num_constraints(self.scaling)

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

  @variable_weights.setter
  def variable_weights(self, values):
    assert values.shape == (self.num_variables,)
    assert values.dtype == np.int64

    for i in range(self.num_variables):
      self.set_variable_weight(i, values[i])

  @property
  def constraint_weights(self):
    length = self.num_constraints
    cdef int[:] values = <int[:length]> csleqp.sleqp_scaling_get_cons_weights(self.scaling)

    array = np.asarray(values)
    array.flags.writeable = False

    return Array(array, self)

  @constraint_weights.setter
  def constraint_weights(self, values):
    assert values.shape == (self.num_constraints,)
    assert values.dtype == np.int64

    for i in range(self.num_constraints):
      self.set_constraint_weight(i, values[i])

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

  def set_variable_weight_from_nominal(self, int index, float nominal_value):
    csleqp_call(csleqp.sleqp_scaling_set_var_weight_from_nominal(self.scaling,
                                                                 index,
                                                                 nominal_value))

  def set_constraint_weight(self, int index, int weight):
    csleqp_call(csleqp.sleqp_scaling_set_cons_weight(self.scaling,
                                                     index,
                                                     weight))

  def set_constraint_weights_from_nominal(self, nominal_array):
    assert nominal_array.shape == (self.num_constraints,)
    cdef double[:] nominal_values = nominal_array
    csleqp_call(csleqp.sleqp_scaling_set_cons_weights_from_nominal(self.scaling,
                                                                   &nominal_values[0]))

  def set_constraint_weight_from_nominal(self, int index, float nominal_value):
    csleqp_call(csleqp.sleqp_scaling_set_cons_weight_from_nominal(self.scaling,
                                                                  index,
                                                                  nominal_value))

  def set_from_gradient(self,
                        np.ndarray grad_array,
                        float eps):
    csleqp_call(array_to_sleqp_sparse_vec(grad_array,
                                          self.gradient))

    csleqp_call(csleqp.sleqp_func_scaling_from_gradient(self.scaling,
                                                        self.gradient,
                                                        eps))

  def set_from_cons_jac(self,
                        object mat,
                        float eps):
    csleqp_call(matrix_to_sleqp_sparse_matrix(mat,
                                              self.cons_jac))

    csleqp_call(csleqp.sleqp_scaling_from_cons_jac(self.scaling,
                                                   self.cons_jac,
                                                   eps))
