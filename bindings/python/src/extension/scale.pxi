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
  """
  Scaling of a problem, consisting of weights (powers of 2) for
  variables, constraints, and objective. Scaling weights can also be
  given as nominal values, i.e., values which are transformed to ~1 by
  the scaling.

  """
  cdef dict __dict__
  cdef csleqp.SleqpScaling* scaling
  cdef csleqp.SleqpVec* gradient
  cdef csleqp.SleqpMat* cons_jac

  def __cinit__(self,
                int num_vars,
                int num_cons,
                *args,
                **keywords):

    csleqp_call(csleqp.sleqp_scaling_create(&self.scaling,
                                            num_vars,
                                            num_cons))

    csleqp_call(csleqp.sleqp_vec_create_empty(&self.gradient,
                                                        num_vars))

    csleqp_call(csleqp.sleqp_mat_create(&self.cons_jac,
                                                  num_cons,
                                                  num_vars,
                                                  0))

  @property
  def num_vars(self):
    """
    Number of variables in the scaling
    """
    return csleqp.sleqp_scaling_num_vars(self.scaling)

  @property
  def num_cons(self):
    """
    Number of constraints in the scaling
    """
    return csleqp.sleqp_scaling_num_cons(self.scaling)

  def __dealloc__(self):
    csleqp_call(csleqp.sleqp_mat_release(&self.cons_jac))
    csleqp_call(csleqp.sleqp_vec_free(&self.gradient))
    csleqp_call(csleqp.sleqp_scaling_release(&self.scaling))

  def __str__(self):
    val = "Scaling(num_variables={0}, num_constraints={1})\n".format(self.num_vars,
                                                                     self.num_cons)

    val += "Objective weight: {0}\n".format(self.obj_weight)

    val += "Variable weights: {0}\n".format(self.variable_weights)

    val += "Constraint weights: {0}\n".format(self.variable_weights)

    return val

  @staticmethod
  def weights_to_nominal_values(weights):
    """
    Convert weights to nominal values
    """
    return np.ldexp(np.ones(weights.shape), weights)

  @staticmethod
  def nominal_values_to_weights(nominal_values):
    """
    Convert nominal values to weights
    """
    return np.frexp(nominal_values)[1]

  @staticmethod
  def nominal_value_to_weight(nominal_value):
    """
    Convert single nominal value to weight
    """
    return np.frexp(np.array([nominal_value]))[1].item()

  @staticmethod
  def weight_to_nominal_value(weight):
    """
    Convert single weight to nominal value
    """
    return np.ldexp(np.ones(1), np.array([weight])).item()

  @property
  def obj_weight(self):
    """
    Weight of objective
    """
    return csleqp.sleqp_scaling_obj_weight(self.scaling)

  @obj_weight.setter
  def obj_weight(self, value):
    csleqp_call(csleqp.sleqp_scaling_set_obj_weight(self.scaling,
                                                    value))

  @property
  def var_weights(self):
    """
    Weights of variables
    """
    length = self.num_vars
    cdef int[:] values = <int[:length]> csleqp.sleqp_scaling_var_weights(self.scaling)

    array = np.asarray(values)
    array.flags.writeable = False

    return Array(array, self)

  @var_weights.setter
  def var_weights(self, values):
    assert values.shape == (self.num_vars,)
    assert values.dtype == np.int64

    for i in range(self.num_vars):
      self.set_variable_weight(i, values[i])

  @property
  def cons_weights(self):
    """
    Weights of constraints
    """
    length = self.num_cons
    cdef int[:] values = <int[:length]> csleqp.sleqp_scaling_cons_weights(self.scaling)

    array = np.asarray(values)
    array.flags.writeable = False

    return Array(array, self)

  @cons_weights.setter
  def cons_weights(self, values):
    assert values.shape == (self.num_cons,)
    assert values.dtype == np.int64

    for i in range(self.num_cons):
      self.set_constraint_weight(i, values[i])

  def set_obj_weight_from_nominal(self, nominal_value):
    csleqp_call(csleqp.sleqp_scaling_set_obj_weight_from_nominal(self.scaling,
                                                                 nominal_value))

  def set_var_weight(self, int index, int weight):
    csleqp_call(csleqp.sleqp_scaling_set_var_weight(self.scaling,
                                                    index,
                                                    weight))

  def set_var_weights_from_nominal(self, nominal_array):
    assert nominal_array.shape == (self.num_vars,)
    cdef double[:] nominal_values = nominal_array
    csleqp_call(csleqp.sleqp_scaling_set_var_weights_from_nominal(self.scaling,
                                                                  &nominal_values[0]))

  def set_var_weight_from_nominal(self, int index, float nominal_value):
    csleqp_call(csleqp.sleqp_scaling_set_var_weight_from_nominal(self.scaling,
                                                                 index,
                                                                 nominal_value))

  def set_cons_weight(self, int index, int weight):
    csleqp_call(csleqp.sleqp_scaling_set_cons_weight(self.scaling,
                                                     index,
                                                     weight))

  def set_cons_weights_from_nominal(self, nominal_array):
    assert nominal_array.shape == (self.num_cons,)
    cdef double[:] nominal_values = nominal_array
    csleqp_call(csleqp.sleqp_scaling_set_cons_weights_from_nominal(self.scaling,
                                                                   &nominal_values[0]))

  def set_cons_weight_from_nominal(self, int index, float nominal_value):
    csleqp_call(csleqp.sleqp_scaling_set_cons_weight_from_nominal(self.scaling,
                                                                  index,
                                                                  nominal_value))

  def set_from_obj_grad(self,
                        np.ndarray grad_array,
                        float eps):
    csleqp_call(array_to_sleqp_sparse_vec(grad_array,
                                          self.gradient))

    csleqp_call(csleqp.sleqp_obj_scaling_from_grad(self.scaling,
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
