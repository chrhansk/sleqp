#cython: language_level=3

class Array(np.ndarray):
  """
  Thin wrapper around a :class:`numpy.ndarray`, storing a
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
  cdef csleqp.SleqpScaling* scaling
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
    return csleqp.sleqp_scaling_num_vars(self.scaling)

  @property
  def num_constraints(self):
    return csleqp.sleqp_scaling_num_cons(self.scaling)

  def __dealloc__(self):
    csleqp_call(csleqp.sleqp_sparse_matrix_release(&self.cons_jac))
    csleqp_call(csleqp.sleqp_sparse_vector_free(&self.gradient))
    csleqp_call(csleqp.sleqp_scaling_release(&self.scaling))

  def __str__(self):
    val = "Scaling(num_variables={0}, num_constraints={1})\n".format(self.num_variables,
                                                                     self.num_constraints)

    val += "Function weight: {0}\n".format(self.func_weight)

    val += "Variable weights: {0}\n".format(self.variable_weights)

    val += "Constraint weights: {0}\n".format(self.variable_weights)

    return val

  @staticmethod
  def weights_to_nominal_values(weights):
    return np.ldexp(np.ones(weights.shape), weights)

  @staticmethod
  def nominal_values_to_weights(nominal_values):
    return np.frexp(nominal_values)[1]

  @staticmethod
  def nominal_value_to_weight(nominal_value):
    return np.frexp(np.array([nominal_value]))[1].item()

  @staticmethod
  def weight_to_nominal_value(weight):
    return np.ldexp(np.ones(1), np.array([weight])).item()

  @property
  def func_weight(self):
    return csleqp.sleqp_scaling_func_weight(self.scaling)

  @func_weight.setter
  def func_weight(self, value):
    csleqp_call(csleqp.sleqp_scaling_set_func_weight(self.scaling,
                                                     value))

  @property
  def variable_weights(self):
    length = self.num_variables
    cdef int[:] values = <int[:length]> csleqp.sleqp_scaling_var_weights(self.scaling)

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
    cdef int[:] values = <int[:length]> csleqp.sleqp_scaling_cons_weights(self.scaling)

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
