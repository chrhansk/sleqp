#cython: language_level=3

cdef csleqp.SLEQP_RETCODE sleqp_lsq_eval(int num_variables,
                                         csleqp.SleqpSparseVec* residual,
                                         void* func_data):
  try:
    func = (<object> func_data)

    array = func.lsq_residuals()
    csleqp_call(array_to_sleqp_sparse_vec(array, residual))

  except BaseException as exception:
    func.call_exception = exception
    return csleqp.SLEQP_INTERNAL_ERROR

  return csleqp.SLEQP_OKAY

cdef csleqp.SLEQP_RETCODE sleqp_lsq_jac_forward(int num_variables,
                                                const csleqp.SleqpSparseVec* forward_direction,
                                                csleqp.SleqpSparseVec* product,
                                                void* func_data):
  try:
      func = (<object> func_data)

      forward_direction_array = sleqp_sparse_vec_to_array(forward_direction)

      product_array = func.lsq_jac_forward(forward_direction_array)

      csleqp_call(array_to_sleqp_sparse_vec(product_array, product))

  except BaseException as exception:
    func.call_exception = exception
    return csleqp.SLEQP_INTERNAL_ERROR

  return csleqp.SLEQP_OKAY


cdef csleqp.SLEQP_RETCODE sleqp_lsq_jac_adjoint(int num_variables,
                                                const csleqp.SleqpSparseVec* adjoint_direction,
                                                csleqp.SleqpSparseVec* product,
                                                void* func_data):
  try:
      func = (<object> func_data)

      adjoint_direction_array = sleqp_sparse_vec_to_array(adjoint_direction)

      product_array = func.lsq_jac_adjoint(adjoint_direction_array)

      csleqp_call(array_to_sleqp_sparse_vec(product_array, product))
  except BaseException as exception:
    func.call_exception = exception
    return csleqp.SLEQP_INTERNAL_ERROR

  return csleqp.SLEQP_OKAY


cdef class LSQFunc:
  cdef csleqp.SleqpFunc* func
  cdef int num_variables
  cdef int num_constraints
  cdef int num_residuals

  cdef public object call_exception

  def __cinit__(self,
                int num_variables,
                int num_constraints,
                int num_residuals,
                double levenberg_marquardt,
                Params params,
                *args,
                **keywords):

    cdef csleqp.SleqpLSQCallbacks callbacks

    callbacks.set_value            = &sleqp_func_set
    callbacks.lsq_eval             = &sleqp_lsq_eval
    callbacks.lsq_jac_forward      = &sleqp_lsq_jac_forward
    callbacks.lsq_jac_adjoint      = &sleqp_lsq_jac_adjoint
    callbacks.eval_additional      = &sleqp_func_eval
    callbacks.hess_prod_additional = &sleqp_func_hess_product
    callbacks.func_free            = &sleqp_func_free

    csleqp_call(csleqp.sleqp_lsq_func_create(&self.func,
                                             &callbacks,
                                             num_variables,
                                             num_residuals,
                                             levenberg_marquardt,
                                             params.params,
                                             <void*> self))

    self.num_variables = num_variables
    self.num_constraints = num_constraints

    self.call_exception = None

    assert(self.func)


  cpdef void set_value(self, x):
    pass

  cpdef double func_val(self):
    pass

  cpdef object lsq_residuals(self):
    pass

  cpdef object lsq_jac_forward(self, forward_direction: np.array):
    pass

  cpdef object lsq_jac_adjoint(self, adjoint_direction: np.array):
    pass

  cpdef int func_grad_nnz(self):
    return 0

  cpdef int cons_val_nnz(self):
    return 0

  cpdef int cons_jac_nnz(self):
    return 0

  cpdef object func_grad(self):
    pass

  cpdef object cons_vals(self):
    pass

  cpdef object cons_jac(self):
    pass

  cpdef object hess_prod(self,
                         func_dual: float,
                         direction: np.array,
                         cons_dual: np.array):
    pass

  @property
  def num_variables(self) -> int:
      return self.num_variables

  @property
  def num_constraints(self) -> int:
      return self.num_constraints

  def __dealloc__(self):
    csleqp_call(csleqp.sleqp_func_release(&self.func))
