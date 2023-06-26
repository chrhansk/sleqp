#cython: language_level=3

cdef store_func_exc(func_obj, exception):
  func_obj.call_exception = exception

cdef csleqp.SLEQP_RETCODE sleqp_func_set(csleqp.SleqpFunc* func,
                                         csleqp.SleqpVec* x,
                                         csleqp.SLEQP_VALUE_REASON reason,
                                         csleqp.bool* reject,
                                         void* func_data):

  try:
    x_array = _readonly(sleqp_sparse_vec_to_array(x))

    func_obj = (<object> func_data)

    do_reject = func_obj.set_value(x_array, ValueReason(reason))

    if do_reject is True:
      reject[0] = True
      return csleqp.SLEQP_OKAY
    else:
      reject[0] = False

  except BaseException as exception:
    store_func_exc(func_obj, exception)
    return csleqp.SLEQP_ERROR

  return csleqp.SLEQP_OKAY


cdef csleqp.SLEQP_RETCODE sleqp_func_set_nogil(csleqp.SleqpFunc* func,
                                               csleqp.SleqpVec* x,
                                               csleqp.SLEQP_VALUE_REASON reason,
                                               csleqp.bool* reject,
                                               void* func_data) nogil:
  with gil:
    return sleqp_func_set(func,
                          x,
                          reason,
                          reject,
                          func_data)


cdef csleqp.SLEQP_RETCODE sleqp_func_nonzeros(csleqp.SleqpFunc* func,
                                              int* obj_grad_nnz,
                                              int* cons_val_nnz,
                                              int* cons_jac_nnz,
                                              int* hess_prod_nnz,
                                              void* func_data):

  try:
    func_obj = (<object> func_data)

    def try_call(f):
      try:
        return f()
      except AttributeError:
        return csleqp.SLEQP_NONE

    obj_grad_nnz[0] = try_call(lambda : func_obj.obj_grad_nnz())
    cons_val_nnz[0] = try_call(lambda : func_obj.cons_val_nnz())
    cons_jac_nnz[0] = try_call(lambda : func_obj.cons_jac_nnz())
    hess_prod_nnz[0] = try_call(lambda : func_obj.hess_prod_nnz())

  except BaseException as exception:
    store_func_exc(func_obj, exception)
    return csleqp.SLEQP_ERROR

  return csleqp.SLEQP_OKAY


cdef csleqp.SLEQP_RETCODE sleqp_func_nonzeros_nogil(csleqp.SleqpFunc* func,
                                                    int* obj_grad_nnz,
                                                    int* cons_val_nnz,
                                                    int* cons_jac_nnz,
                                                    int* hess_prod_nnz,
                                                    void* func_data) nogil:
  with gil:
    return sleqp_func_nonzeros(func,
                               obj_grad_nnz,
                               cons_val_nnz,
                               cons_jac_nnz,
                               hess_prod_nnz,
                               func_data)


cdef csleqp.SLEQP_RETCODE sleqp_func_obj_val(csleqp.SleqpFunc* func,
                                             double* obj_val,
                                             void* func_data):
  try:

    func_obj = (<object> func_data)

    result = func_obj.obj_val()
    assert result is not None, "obj_val() returned 'None'"
    obj_val[0] = result

  except BaseException as exception:
    store_func_exc(func_obj, exception)
    return csleqp.SLEQP_ERROR

  return csleqp.SLEQP_OKAY


cdef csleqp.SLEQP_RETCODE sleqp_func_obj_val_nogil(csleqp.SleqpFunc* func,
                                                   double* obj_val,
                                                   void* func_data) nogil:
  with gil:
    return sleqp_func_obj_val(func,
                              obj_val,
                              func_data)


cdef csleqp.SLEQP_RETCODE sleqp_func_obj_grad(csleqp.SleqpFunc* func,
                                              csleqp.SleqpVec* obj_grad,
                                              void* func_data):
  cdef int num_vars
  try:

    func_obj = (<object> func_data)

    num_vars = csleqp.sleqp_func_num_vars(func)

    grad_array = func_obj.obj_grad()

    assert grad_array is not None, "obj_grad() returned 'None'"
    assert grad_array.ndim == 1, "Gradient must be a vector"
    assert grad_array.size == num_vars, "Gradient has wrong size"

    csleqp_call(array_to_sleqp_sparse_vec(grad_array, obj_grad))

  except BaseException as exception:
    store_func_exc(func_obj, exception)
    return csleqp.SLEQP_ERROR

  return csleqp.SLEQP_OKAY


cdef csleqp.SLEQP_RETCODE sleqp_func_obj_grad_nogil(csleqp.SleqpFunc* func,
                                                    csleqp.SleqpVec* obj_grad,
                                                    void* func_data) nogil:
  with gil:
    return sleqp_func_obj_grad(func,
                               obj_grad,
                               func_data)


cdef csleqp.SLEQP_RETCODE sleqp_func_cons_val(csleqp.SleqpFunc* func,
                                              csleqp.SleqpVec* cons_vals,
                                              void* func_data):
  cdef int num_cons
  try:

    func_obj = (<object> func_data)

    num_cons = csleqp.sleqp_func_num_cons(func)

    cons_array = func_obj.cons_vals()

    assert cons_array is not None, "cons_vals() returned 'None'"
    assert cons_array.ndim == 1, "Constraint values must be a vector"
    assert cons_array.size == num_cons, "Constraint values have wrong size"

    csleqp_call(array_to_sleqp_sparse_vec(cons_array, cons_vals))

  except BaseException as exception:
    store_func_exc(func_obj, exception)
    return csleqp.SLEQP_ERROR

  return csleqp.SLEQP_OKAY


cdef csleqp.SLEQP_RETCODE sleqp_func_cons_val_nogil(csleqp.SleqpFunc* func,
                                                    csleqp.SleqpVec* cons_vals,
                                                    void* func_data) nogil:
  with gil:
    return sleqp_func_cons_val(func,
                               cons_vals,
                               func_data)


cdef csleqp.SLEQP_RETCODE sleqp_func_cons_jac(csleqp.SleqpFunc* func,
                                              csleqp.SleqpMat* cons_jac,
                                              void* func_data):
  cdef int num_vars
  cdef int num_cons
  try:

    func_obj = (<object> func_data)

    num_vars = csleqp.sleqp_func_num_vars(func)
    num_cons = csleqp.sleqp_func_num_cons(func)

    cons_jac_mat = func_obj.cons_jac()

    expected_shape = (num_cons, num_vars)

    assert cons_jac_mat is not None, "cons_jac() returned 'None'"
    assert cons_jac_mat.ndim == 2, "Constraint Jacobian must be a matrix"
    assert cons_jac_mat.shape == expected_shape, "Constraint Jacobian has wrong shape"

    csleqp_call(matrix_to_sleqp_sparse_matrix(cons_jac_mat, cons_jac))


  except BaseException as exception:
    store_func_exc(func_obj, exception)
    return csleqp.SLEQP_ERROR

  return csleqp.SLEQP_OKAY


cdef csleqp.SLEQP_RETCODE sleqp_func_cons_jac_nogil(csleqp.SleqpFunc* func,
                                                    csleqp.SleqpMat* cons_jac,
                                                    void* func_data) nogil:
  with gil:
    return sleqp_func_cons_jac(func,
                               cons_jac,
                               func_data)


cdef csleqp.SLEQP_RETCODE sleqp_func_hess_prod(csleqp.SleqpFunc* func,
                                               const csleqp.SleqpVec* direction,
                                               const csleqp.SleqpVec* cons_dual,
                                               csleqp.SleqpVec* product,
                                               void* func_data):

  cdef int num_vars

  try:

    assert cons_dual
    assert direction
    assert product

    func_obj = (<object> func_data)

    num_vars = csleqp.sleqp_func_num_vars(func)

    direction_array = _readonly(sleqp_sparse_vec_to_array(direction))
    cons_dual_array = _readonly(sleqp_sparse_vec_to_array(cons_dual))

    product_array = func_obj.hess_prod(direction_array, cons_dual_array)

    assert product_array is not None, "hess_prod(...) returned 'None'"
    assert product_array.ndim == 1, "Hessian product must be a vector"
    assert product_array.size == num_vars, "Hessian product has wrong size"

    csleqp_call(array_to_sleqp_sparse_vec(product_array, product))

  except BaseException as exception:
    store_func_exc(func_obj, exception)
    return csleqp.SLEQP_ERROR

  return csleqp.SLEQP_OKAY


cdef csleqp.SLEQP_RETCODE sleqp_func_hess_prod_nogil(csleqp.SleqpFunc* func,
                                                     const csleqp.SleqpVec* direction,
                                                     const csleqp.SleqpVec* cons_dual,
                                                     csleqp.SleqpVec* product,
                                                     void* func_data) nogil:
  with gil:
    return sleqp_func_hess_prod(func,
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
    callbacks[0].nonzeros  = &sleqp_func_nonzeros_nogil
    callbacks[0].obj_val   = &sleqp_func_obj_val_nogil
    callbacks[0].obj_grad  = &sleqp_func_obj_grad_nogil
    callbacks[0].cons_val  = &sleqp_func_cons_val_nogil
    callbacks[0].cons_jac  = &sleqp_func_cons_jac_nogil
    callbacks[0].hess_prod = &sleqp_func_hess_prod_nogil
  else:
    callbacks[0].set_value = &sleqp_func_set
    callbacks[0].nonzeros  = &sleqp_func_nonzeros
    callbacks[0].obj_val   = &sleqp_func_obj_val
    callbacks[0].obj_grad  = &sleqp_func_obj_grad
    callbacks[0].cons_val  = &sleqp_func_cons_val
    callbacks[0].cons_jac  = &sleqp_func_cons_jac
    callbacks[0].hess_prod = &sleqp_func_hess_prod

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
                                      int num_vars,
                                      int num_cons):
  cdef csleqp.SleqpFuncCallbacks callbacks
  cdef _Func ofunc = _Func()
  cdef csleqp.SLEQP_RETCODE retcode

  assert func is not None

  set_func_callbacks(&callbacks)

  retcode = csleqp.sleqp_func_create(cfunc,
                                     &callbacks,
                                     num_vars,
                                     num_cons,
                                     <void*> func)

  if cfunc[0] != NULL:
    ofunc.set_func(cfunc[0])
    funcs.add(ofunc)

  return retcode
