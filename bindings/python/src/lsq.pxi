#cython: language_level=3

cdef csleqp.SLEQP_RETCODE sleqp_lsq_residuals(csleqp.SleqpFunc* func,
                                              csleqp.SleqpSparseVec* residual,
                                              void* func_data):
  try:
    func_obj = (<object> func_data)

    array = func_obj.lsq_residuals()
    csleqp_call(array_to_sleqp_sparse_vec(array, residual))

  except BaseException as exception:
    store_func_exc(func_obj, exception)
    return csleqp.SLEQP_ERROR

  return csleqp.SLEQP_OKAY


cdef csleqp.SLEQP_RETCODE sleqp_lsq_residuals_nogil(csleqp.SleqpFunc* func,
                                                    csleqp.SleqpSparseVec* residual,
                                                    void* func_data) nogil:
  with gil:
    return sleqp_lsq_residuals(func,
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
    store_func_exc(func_obj, exception)
    return csleqp.SLEQP_ERROR

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
    store_func_exc(func_obj, exception)
    return csleqp.SLEQP_ERROR

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
    callbacks.lsq_residuals        = &sleqp_lsq_residuals_nogil
    callbacks.lsq_jac_forward      = &sleqp_lsq_jac_forward_nogil
    callbacks.lsq_jac_adjoint      = &sleqp_lsq_jac_adjoint_nogil
    callbacks.cons_val             = &sleqp_func_cons_val_nogil
    callbacks.cons_jac             = &sleqp_func_cons_jac_nogil
  else:
    callbacks.set_value            = &sleqp_func_set
    callbacks.lsq_residuals        = &sleqp_lsq_residuals
    callbacks.lsq_jac_forward      = &sleqp_lsq_jac_forward
    callbacks.lsq_jac_adjoint      = &sleqp_lsq_jac_adjoint
    callbacks.cons_val             = &sleqp_func_cons_val
    callbacks.cons_jac             = &sleqp_func_cons_jac

  callbacks.func_free = &sleqp_func_free


cdef update_lsq_func_callbacks():
  cdef _Func func
  cdef csleqp.SleqpLSQCallbacks callbacks

  set_lsq_func_callbacks(&callbacks)

  for obj in lsq_funcs:
    func = <_Func> obj

    csleqp_call(csleqp.sleqp_lsq_func_set_callbacks(func.cfunc,
                                                    &callbacks))

cdef csleqp.SLEQP_RETCODE create_lsq_func(csleqp.SleqpFunc** cfunc,
                                          object func,
                                          int num_variables,
                                          int num_constraints,
                                          int num_residuals,
                                          double levenberg_marquardt,
                                          csleqp.SleqpParams* params):
  cdef csleqp.SleqpLSQCallbacks callbacks
  cdef csleqp.SLEQP_RETCODE retcode

  assert func is not None

  set_lsq_func_callbacks(&callbacks)

  retcode = csleqp.sleqp_lsq_func_create(cfunc,
                                         &callbacks,
                                         num_variables,
                                         num_constraints,
                                         num_residuals,
                                         levenberg_marquardt,
                                         params,
                                         <void*> func)
  
  return retcode
    
