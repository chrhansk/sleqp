#cython: language_level=3

cdef csleqp.SLEQP_RETCODE sleqp_lsq_eval(csleqp.SleqpFunc* func,
                                         csleqp.SleqpSparseVec* residual,
                                         void* func_data):
  try:
    func_obj = (<object> func_data)

    array = func_obj.lsq_residuals()
    csleqp_call(array_to_sleqp_sparse_vec(array, residual))

  except BaseException as exception:
    func_obj.call_exception = exception
    return csleqp.SLEQP_INTERNAL_ERROR

  return csleqp.SLEQP_OKAY


cdef csleqp.SLEQP_RETCODE sleqp_lsq_eval_nogil(csleqp.SleqpFunc* func,
                                               csleqp.SleqpSparseVec* residual,
                                               void* func_data) nogil:
  with gil:
    return sleqp_lsq_eval(func,
                          residual,
                          func_data)


cdef csleqp.SLEQP_RETCODE sleqp_lsq_jac_forward(csleqp.SleqpFunc* func,
                                                const csleqp.SleqpSparseVec* forward_direction,
                                                csleqp.SleqpSparseVec* product,
                                                void* func_data):
  try:
      func_obj = (<object> func_data)

      forward_direction_array = sleqp_sparse_vec_to_array(forward_direction)

      product_array = func_obj.lsq_jac_forward(forward_direction_array)

      assert product_array is not None, "lsq_jac_forward(...) returned 'None'"

      csleqp_call(array_to_sleqp_sparse_vec(product_array, product))

  except BaseException as exception:
    func_obj.call_exception = exception
    return csleqp.SLEQP_INTERNAL_ERROR

  return csleqp.SLEQP_OKAY


cdef csleqp.SLEQP_RETCODE sleqp_lsq_jac_forward_nogil(csleqp.SleqpFunc* func,
                                                      const csleqp.SleqpSparseVec* forward_direction,
                                                      csleqp.SleqpSparseVec* product,
                                                      void* func_data) nogil:
  with gil:
    return sleqp_lsq_jac_forward(func,
                                 forward_direction,
                                 product,
                                 func_data)


cdef csleqp.SLEQP_RETCODE sleqp_lsq_jac_adjoint(csleqp.SleqpFunc* func,
                                                const csleqp.SleqpSparseVec* adjoint_direction,
                                                csleqp.SleqpSparseVec* product,
                                                void* func_data):
  try:
      func_obj = (<object> func_data)

      adjoint_direction_array = sleqp_sparse_vec_to_array(adjoint_direction)

      product_array = func_obj.lsq_jac_adjoint(adjoint_direction_array)

      assert product_array is not None, "lsq_jac_adjoint(...) returned 'None'"

      csleqp_call(array_to_sleqp_sparse_vec(product_array, product))
  except BaseException as exception:
    func_obj.call_exception = exception
    return csleqp.SLEQP_INTERNAL_ERROR

  return csleqp.SLEQP_OKAY


cdef object lsq_funcs = weakref.WeakSet()


cdef csleqp.SLEQP_RETCODE sleqp_lsq_jac_adjoint_nogil(csleqp.SleqpFunc* func,
                                                      const csleqp.SleqpSparseVec* adjoint_direction,
                                                      csleqp.SleqpSparseVec* product,
                                                      void* func_data) nogil:
  with gil:
    return sleqp_lsq_jac_adjoint(func,
                                 adjoint_direction,
                                 product,
                                 func_data)


cdef set_lsq_func_callbacks(csleqp.SleqpLSQCallbacks* callbacks):
  if release_gil:
    callbacks.set_value            = &sleqp_func_set_nogil
    callbacks.lsq_eval             = &sleqp_lsq_eval_nogil
    callbacks.lsq_jac_forward      = &sleqp_lsq_jac_forward_nogil
    callbacks.lsq_jac_adjoint      = &sleqp_lsq_jac_adjoint_nogil
    callbacks.additional_func_val  = &sleqp_func_val_nogil
    callbacks.additional_func_grad = &sleqp_func_grad_nogil
    callbacks.additional_cons_val  = &sleqp_func_cons_val_nogil
    callbacks.additional_cons_jac  = &sleqp_func_cons_jac_nogil
    callbacks.additional_hess_prod = &sleqp_func_hess_product_nogil
    callbacks.additional_hess_prod = &sleqp_func_hess_product_nogil
  else:
    callbacks.set_value            = &sleqp_func_set
    callbacks.lsq_eval             = &sleqp_lsq_eval
    callbacks.lsq_jac_forward      = &sleqp_lsq_jac_forward
    callbacks.lsq_jac_adjoint      = &sleqp_lsq_jac_adjoint
    callbacks.additional_func_val  = &sleqp_func_val
    callbacks.additional_func_grad = &sleqp_func_grad
    callbacks.additional_cons_val  = &sleqp_func_cons_val
    callbacks.additional_cons_jac  = &sleqp_func_cons_jac
    callbacks.additional_hess_prod = &sleqp_func_hess_product

  callbacks.func_free = &sleqp_func_free


cdef update_lsq_func_callbacks():
  cdef LSQFunc func
  cdef csleqp.SleqpLSQCallbacks callbacks

  set_lsq_func_callbacks(&callbacks)

  for obj in lsq_funcs:
    func = <LSQFunc> obj

    csleqp_call(csleqp.sleqp_lsq_func_set_callbacks(func.func,
                                                    &callbacks))


cdef class LSQFunc:

  cdef dict __dict__

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

    set_lsq_func_callbacks(&callbacks)

    csleqp_call(csleqp.sleqp_lsq_func_create(&self.func,
                                             &callbacks,
                                             num_variables,
                                             num_constraints,
                                             num_residuals,
                                             levenberg_marquardt,
                                             params.params,
                                             <void*> self))

    self.num_variables = num_variables
    self.num_constraints = num_constraints

    self.call_exception = None

    lsq_funcs.add(self)

    assert(self.func)


  cpdef void set_value(self, value: np.array, reason: ValueReason):
    pass

  cpdef double func_val(self):
    pass

  cpdef object lsq_residuals(self):
    return None

  cpdef object lsq_jac_forward(self, forward_direction: np.array):
    return None

  cpdef object lsq_jac_adjoint(self, adjoint_direction: np.array):
    return None

  cpdef int func_grad_nnz(self):
    return 0

  cpdef int cons_val_nnz(self):
    return 0

  cpdef int cons_jac_nnz(self):
    return 0

  cpdef object func_grad(self):
    return None

  cpdef object cons_vals(self):
    return None

  cpdef object cons_jac(self):
    return None

  cpdef object hess_prod(self,
                         func_dual: float,
                         direction: np.array,
                         cons_dual: np.array):
    return None

  @property
  def num_variables(self) -> int:
      return self.num_variables

  @property
  def num_constraints(self) -> int:
      return self.num_constraints

  def __dealloc__(self):
    csleqp_call(csleqp.sleqp_func_release(&self.func))
