#cython: language_level=3


cdef csleqp.SLEQP_RETCODE sleqp_func_set(csleqp.SleqpFunc* func,
                                         csleqp.SleqpSparseVec* x,
                                         csleqp.SLEQP_VALUE_REASON reason,
                                         int* func_grad_nnz,
                                         int* cons_val_nnz,
                                         int* cons_jac_nnz,
                                         void* func_data):

  try:
    x_array = sleqp_sparse_vec_to_array(x)

    func_obj = (<object> func_data)

    func_obj.set_value(x_array, ValueReason(reason))

    func_grad_nnz[0] = func_obj.func_grad_nnz()
    cons_val_nnz[0] = func_obj.cons_val_nnz()
    cons_jac_nnz[0] =  func_obj.cons_jac_nnz()

  except BaseException as exception:
    func_obj.call_exception = exception
    return csleqp.SLEQP_INTERNAL_ERROR

  return csleqp.SLEQP_OKAY


cdef csleqp.SLEQP_RETCODE sleqp_func_set_nogil(csleqp.SleqpFunc* func,
                                               csleqp.SleqpSparseVec* x,
                                               csleqp.SLEQP_VALUE_REASON reason,
                                               int* func_grad_nnz,
                                               int* cons_val_nnz,
                                               int* cons_jac_nnz,
                                               void* func_data) nogil:
  with gil:
    return sleqp_func_set(func,
                          x,
                          reason,
                          func_grad_nnz,
                          cons_val_nnz,
                          cons_jac_nnz,
                          func_data)


cdef csleqp.SLEQP_RETCODE sleqp_func_eval(csleqp.SleqpFunc* func,
                                          const csleqp.SleqpSparseVec* cons_indices,
                                          double* func_val,
                                          csleqp.SleqpSparseVec* func_grad,
                                          csleqp.SleqpSparseVec* cons_vals,
                                          csleqp.SleqpSparseMatrix* cons_jac,
                                          void* func_data):
  cdef int num_variables
  cdef int num_constraints
  try:

    func_obj = (<object> func_data)

    num_variables = csleqp.sleqp_func_get_num_variables(func)
    num_constraints = csleqp.sleqp_func_get_num_constraints(func)

    if func_val:
      result = func_obj.func_val()
      assert result is not None, "func_val() returned 'None'"
      func_val[0] = result

    if func_grad:
      grad_array = func_obj.func_grad()

      assert grad_array is not None, "func_grad() returned 'None'"
      assert grad_array.ndim == 1, "Gradient must be a vector"
      assert grad_array.size == num_variables, "Gradient has wrong size"

      csleqp_call(array_to_sleqp_sparse_vec(grad_array, func_grad))

    if cons_vals:
      cons_array = func_obj.cons_vals()

      assert cons_array is not None, "cons_vals() returned 'None'"
      assert cons_array.ndim == 1, "Constraint values must be a vector"
      assert cons_array.size == num_constraints, "Constraint values have wrong size"

      csleqp_call(array_to_sleqp_sparse_vec(cons_array, cons_vals))

    if cons_jac:
      cons_jac_mat = func_obj.cons_jac()

      expected_shape = (num_constraints, num_variables)

      assert cons_jac_mat is not None, "cons_jac() returned 'None'"
      assert cons_jac_mat.ndim == 2, "Constraint Jacobian must be a matrix"
      assert cons_jac_mat.shape == expected_shape, "Constraint Jacobian has wrong shape"

      csleqp_call(matrix_to_sleqp_sparse_matrix(cons_jac_mat, cons_jac))


  except BaseException as exception:
    func_obj.call_exception = exception
    return csleqp.SLEQP_INTERNAL_ERROR

  return csleqp.SLEQP_OKAY


cdef csleqp.SLEQP_RETCODE sleqp_func_eval_nogil(csleqp.SleqpFunc* func,
                                                const csleqp.SleqpSparseVec* cons_indices,
                                                double* func_val,
                                                csleqp.SleqpSparseVec* func_grad,
                                                csleqp.SleqpSparseVec* cons_vals,
                                                csleqp.SleqpSparseMatrix* cons_jac,
                                                void* func_data) nogil:
  with gil:
    return sleqp_func_eval(func,
                           cons_indices,
                           func_val,
                           func_grad,
                           cons_vals,
                           cons_jac,
                           func_data)


cdef csleqp.SLEQP_RETCODE sleqp_func_hess_product(csleqp.SleqpFunc* func,
                                                  const double* func_dual,
                                                  const csleqp.SleqpSparseVec* direction,
                                                  const csleqp.SleqpSparseVec* cons_dual,
                                                  csleqp.SleqpSparseVec* product,
                                                  void* func_data):

  cdef int num_variables

  try:

    assert cons_dual
    assert direction
    assert product

    func_obj = (<object> func_data)

    num_variables = csleqp.sleqp_func_get_num_variables(func)

    f_dual = func_dual[0] if func_dual else 0.
    direction_array = sleqp_sparse_vec_to_array(direction)
    cons_dual_array = sleqp_sparse_vec_to_array(cons_dual)

    product_array = func_obj.hess_prod(f_dual, direction_array, cons_dual_array)

    assert product_array is not None, "hess_prod(...) returned 'None'"
    assert product_array.ndim == 1, "Hessian product must be a vector"
    assert product_array.size == num_variables, "Hessian product has wrong size"

    csleqp_call(array_to_sleqp_sparse_vec(product_array, product))

  except BaseException as exception:
    func_obj.call_exception = exception
    return csleqp.SLEQP_INTERNAL_ERROR

  return csleqp.SLEQP_OKAY


cdef csleqp.SLEQP_RETCODE sleqp_func_hess_product_nogil(csleqp.SleqpFunc* func,
                                                        const double* func_dual,
                                                        const csleqp.SleqpSparseVec* direction,
                                                        const csleqp.SleqpSparseVec* cons_dual,
                                                        csleqp.SleqpSparseVec* product,
                                                        void* func_data) nogil:
  with gil:
    return sleqp_func_hess_product(func,
                                   func_dual,
                                   direction,
                                   cons_dual,
                                   product,
                                   func_data)


cdef csleqp.SLEQP_RETCODE sleqp_func_free(void* func_data):
  return csleqp.SLEQP_OKAY


cdef object funcs = weakref.WeakSet()


cdef set_func_callbacks(csleqp.SleqpFuncCallbacks* callbacks):
  if release_gil:
    callbacks[0].set_value = &sleqp_func_set_nogil
    callbacks[0].func_eval = &sleqp_func_eval_nogil
    callbacks[0].hess_prod = &sleqp_func_hess_product_nogil
  else:
    callbacks[0].set_value = &sleqp_func_set
    callbacks[0].func_eval = &sleqp_func_eval
    callbacks[0].hess_prod = &sleqp_func_hess_product

  callbacks.func_free = &sleqp_func_free


cdef update_func_callbacks():
  cdef Func func
  cdef csleqp.SleqpFuncCallbacks callbacks

  set_func_callbacks(&callbacks)

  for obj in funcs:
    func = <Func> obj

    csleqp_call(csleqp.sleqp_func_set_callbacks(func.func,
                                                &callbacks))


cdef class Func:

  cdef dict __dict__

  cdef csleqp.SleqpFunc* func

  cdef public object call_exception

  def __cinit__(self,
                int num_variables,
                int num_constraints,
                *args,
                **keywords):

    cdef csleqp.SleqpFuncCallbacks callbacks

    set_func_callbacks(&callbacks)

    csleqp_call(csleqp.sleqp_func_create(&self.func,
                                         &callbacks,
                                         num_variables,
                                         num_constraints,
                                         <void*> self))

    self.call_exception = None

    funcs.add(self)

    assert(self.func)

  cpdef void set_value(self, value, reason):
    pass

  cpdef double func_val(self):
    pass

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
    pass

  @property
  def hess_struct(self) -> HessianStruct:
    return HessianStruct(self)

  @property
  def num_variables(self) -> int:
    return csleqp.sleqp_func_get_num_variables(self.func)

  @property
  def num_constraints(self) -> int:
    return csleqp.sleqp_func_get_num_constraints(self.func)

  def __dealloc__(self):
    csleqp_call(csleqp.sleqp_func_release(&self.func))
