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

    def try_call(f):
      try:
        return f()
      except AttributeError:
        return 0

    func_grad_nnz[0] = try_call(lambda : func_obj.func_grad_nnz())
    cons_val_nnz[0] = try_call(lambda : func_obj.cons_val_nnz())
    cons_jac_nnz[0] = try_call(lambda : func_obj.cons_jac_nnz())

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


cdef csleqp.SLEQP_RETCODE sleqp_func_val(csleqp.SleqpFunc* func,
                                         double* func_val,
                                         void* func_data):
  cdef int num_variables
  cdef int num_constraints
  try:

    func_obj = (<object> func_data)

    num_variables = csleqp.sleqp_func_get_num_variables(func)
    num_constraints = csleqp.sleqp_func_get_num_constraints(func)

    result = func_obj.func_val()
    assert result is not None, "func_val() returned 'None'"
    func_val[0] = result

  except BaseException as exception:
    func_obj.call_exception = exception
    return csleqp.SLEQP_INTERNAL_ERROR

  return csleqp.SLEQP_OKAY


cdef csleqp.SLEQP_RETCODE sleqp_func_val_nogil(csleqp.SleqpFunc* func,
                                                double* func_val,
                                                void* func_data) nogil:
  with gil:
    return sleqp_func_val(func,
                          func_val,
                          func_data)


cdef csleqp.SLEQP_RETCODE sleqp_func_grad(csleqp.SleqpFunc* func,
                                          csleqp.SleqpSparseVec* func_grad,
                                          void* func_data):
  cdef int num_variables
  cdef int num_constraints
  try:

    func_obj = (<object> func_data)

    num_variables = csleqp.sleqp_func_get_num_variables(func)
    num_constraints = csleqp.sleqp_func_get_num_constraints(func)

    grad_array = func_obj.func_grad()

    assert grad_array is not None, "func_grad() returned 'None'"
    assert grad_array.ndim == 1, "Gradient must be a vector"
    assert grad_array.size == num_variables, "Gradient has wrong size"

    csleqp_call(array_to_sleqp_sparse_vec(grad_array, func_grad))

  except BaseException as exception:
    func_obj.call_exception = exception
    return csleqp.SLEQP_INTERNAL_ERROR

  return csleqp.SLEQP_OKAY


cdef csleqp.SLEQP_RETCODE sleqp_func_grad_nogil(csleqp.SleqpFunc* func,
                                                csleqp.SleqpSparseVec* func_grad,
                                                void* func_data) nogil:
  with gil:
    return sleqp_func_grad(func,
                           func_grad,
                           func_data)


cdef csleqp.SLEQP_RETCODE sleqp_func_cons_val(csleqp.SleqpFunc* func,
                                              const csleqp.SleqpSparseVec* cons_indices,
                                              csleqp.SleqpSparseVec* cons_vals,
                                              void* func_data):
  cdef int num_variables
  cdef int num_constraints
  try:

    func_obj = (<object> func_data)

    num_variables = csleqp.sleqp_func_get_num_variables(func)
    num_constraints = csleqp.sleqp_func_get_num_constraints(func)

    cons_array = func_obj.cons_vals()

    assert cons_array is not None, "cons_vals() returned 'None'"
    assert cons_array.ndim == 1, "Constraint values must be a vector"
    assert cons_array.size == num_constraints, "Constraint values have wrong size"

    csleqp_call(array_to_sleqp_sparse_vec(cons_array, cons_vals))

  except BaseException as exception:
    func_obj.call_exception = exception
    return csleqp.SLEQP_INTERNAL_ERROR

  return csleqp.SLEQP_OKAY


cdef csleqp.SLEQP_RETCODE sleqp_func_cons_val_nogil(csleqp.SleqpFunc* func,
                                                    const csleqp.SleqpSparseVec* cons_indices,
                                                    csleqp.SleqpSparseVec* cons_vals,
                                                    void* func_data) nogil:
  with gil:
    return sleqp_func_cons_val(func,
                               cons_indices,
                               cons_vals,
                               func_data)


cdef csleqp.SLEQP_RETCODE sleqp_func_cons_jac(csleqp.SleqpFunc* func,
                                              const csleqp.SleqpSparseVec* cons_indices,
                                              csleqp.SleqpSparseMatrix* cons_jac,
                                              void* func_data):
  cdef int num_variables
  cdef int num_constraints
  try:

    func_obj = (<object> func_data)

    num_variables = csleqp.sleqp_func_get_num_variables(func)
    num_constraints = csleqp.sleqp_func_get_num_constraints(func)

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


cdef csleqp.SLEQP_RETCODE sleqp_func_cons_jac_nogil(csleqp.SleqpFunc* func,
                                                    const csleqp.SleqpSparseVec* cons_indices,
                                                    csleqp.SleqpSparseMatrix* cons_jac,
                                                    void* func_data) nogil:
  with gil:
    return sleqp_func_cons_jac(func,
                               cons_indices,
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
    callbacks[0].func_val  = &sleqp_func_val_nogil
    callbacks[0].func_grad = &sleqp_func_grad_nogil
    callbacks[0].cons_val  = &sleqp_func_cons_val_nogil
    callbacks[0].cons_jac  = &sleqp_func_cons_jac_nogil
    callbacks[0].hess_prod = &sleqp_func_hess_product_nogil
  else:
    callbacks[0].set_value = &sleqp_func_set
    callbacks[0].func_val  = &sleqp_func_val
    callbacks[0].func_grad = &sleqp_func_grad
    callbacks[0].cons_val  = &sleqp_func_cons_val
    callbacks[0].cons_jac  = &sleqp_func_cons_jac
    callbacks[0].hess_prod = &sleqp_func_hess_product

  callbacks.func_free = &sleqp_func_free


cdef update_func_callbacks():
  cdef _Func func
  cdef csleqp.SleqpFuncCallbacks callbacks

  set_func_callbacks(&callbacks)

  for obj in funcs:
    func = <_Func> obj

    if not func.cfunc:
      continue

    csleqp_call(csleqp.sleqp_func_set_callbacks(func.cfunc,
                                                &callbacks))


cdef class _Func:
  cdef csleqp.SleqpFunc* cfunc
  cdef object __weakref__

  def __cinit__(self):
    self.cfunc = NULL

  cdef void set_func(self, csleqp.SleqpFunc* cfunc):
    self.release()
    self.cfunc = cfunc
    self.capture()

  cdef void capture(self):
    if self.cfunc != NULL:
      csleqp_call(csleqp.sleqp_func_capture(self.cfunc))

  cdef void release(self):
    if self.cfunc != NULL:
      csleqp_call(csleqp.sleqp_func_release(&self.cfunc))

  def __dealloc__(self):
    self.release()


cdef csleqp.SLEQP_RETCODE create_func(csleqp.SleqpFunc** cfunc,
                                      object func,
                                      int num_variables,
                                      int num_constraints):
  cdef csleqp.SleqpFuncCallbacks callbacks
  cdef _Func ofunc = _Func()
  cdef csleqp.SLEQP_RETCODE retcode

  assert func is not None

  set_func_callbacks(&callbacks)

  retcode = csleqp.sleqp_func_create(cfunc,
                                     &callbacks,
                                     num_variables,
                                     num_constraints,
                                     <void*> func)

  if cfunc[0] != NULL:
    ofunc.set_func(cfunc[0])
    funcs.add(ofunc)

  return retcode
