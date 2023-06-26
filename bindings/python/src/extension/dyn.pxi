#cython: language_level=3

cdef csleqp.SLEQP_RETCODE sleqp_dyn_eval(csleqp.SleqpFunc* func,
                                         double* obj_val,
                                         csleqp.SleqpVec* cons_val,
                                         double* error,
                                         void* func_data):
  try:
    func_obj = (<object> func_data)

    result = func_obj.eval()

    assert result is not None, "eval() returned 'None'"

    num_cons = csleqp.sleqp_func_num_cons(func)

    obj_val[0] = result[0]
    error[0] = result[1]

    if num_cons > 0:
      csleqp_call(array_to_sleqp_sparse_vec(result[2], cons_val))

  except BaseException as exception:
    store_func_exc(func_obj, exception)
    return csleqp.SLEQP_ERROR

  return csleqp.SLEQP_OKAY

cdef csleqp.SLEQP_RETCODE sleqp_dyn_eval_nogil(csleqp.SleqpFunc* func,
                                               double* obj_val,
                                               csleqp.SleqpVec* cons_val,
                                               double* error,
                                               void* func_data) nogil:
  with gil:
    return sleqp_dyn_eval(func,
                          obj_val,
                          cons_val,
                          error,
                          func_data)

cdef csleqp.SLEQP_RETCODE sleqp_dyn_set_error_bound(csleqp.SleqpFunc* func,
                                                    double error_bound,
                                                    void* func_data):
  try:
    func_obj = (<object> func_data)

    func_obj.set_error_bound(error_bound)

  except BaseException as exception:
    store_func_exc(func_obj, exception)
    return csleqp.SLEQP_ERROR

cdef csleqp.SLEQP_RETCODE sleqp_dyn_set_error_bound_nogil(csleqp.SleqpFunc* func,
                                                          double error_bound,
                                                          void* func_data) nogil:
  with gil:
    return sleqp_dyn_set_error_bound(func,
                                     error_bound,
                                     func_data)

cdef csleqp.SLEQP_RETCODE sleqp_dyn_set_obj_weight(csleqp.SleqpFunc* func,
                                                   double obj_weight,
                                                   void* func_data):
  try:
    func_obj = (<object> func_data)

    func_obj.set_obj_weight(obj_weight)

  except BaseException as exception:
    store_func_exc(func_obj, exception)
    return csleqp.SLEQP_ERROR

cdef csleqp.SLEQP_RETCODE sleqp_dyn_set_obj_weight_nogil(csleqp.SleqpFunc* func,
                                                         double obj_weight,
                                                         void* func_data) nogil:
  with gil:
    return sleqp_dyn_set_obj_weight(func,
                                    obj_weight,
                                    func_data)

cdef csleqp.SLEQP_RETCODE sleqp_dyn_set_cons_weights(csleqp.SleqpFunc* func,
                                                     const double* cons_weights,
                                                     void* func_data):
  cdef double[:] cons_view
  try:
    func_obj = (<object> func_data)

    num_cons = csleqp.sleqp_func_num_cons(func)
    cons_view = _readonly(<double[:num_cons]> cons_weights)

    func_obj.set_cons_weights(np.copy(cons_view))

  except BaseException as exception:
    store_func_exc(func_obj, exception)
    return csleqp.SLEQP_ERROR

cdef csleqp.SLEQP_RETCODE sleqp_dyn_set_cons_weights_nogil(csleqp.SleqpFunc* func,
                                                           const double* cons_weights,
                                                           void* func_data) nogil:
  with gil:
    return sleqp_dyn_set_cons_weights(func,
                                      cons_weights,
                                      func_data)


cdef set_dyn_func_callbacks(csleqp.SleqpDynFuncCallbacks* callbacks):
  if release_gil:
    callbacks[0].set_value        = &sleqp_func_set_nogil
    callbacks[0].nonzeros         = &sleqp_func_nonzeros_nogil
    callbacks[0].set_error_bound  = &sleqp_dyn_set_error_bound_nogil
    callbacks[0].set_obj_weight   = &sleqp_dyn_set_obj_weight_nogil
    callbacks[0].set_cons_weights = &sleqp_dyn_set_cons_weights_nogil
    callbacks[0].eval             = &sleqp_dyn_eval_nogil
    callbacks[0].obj_grad         = &sleqp_func_obj_grad_nogil
    callbacks[0].cons_jac         = &sleqp_func_cons_jac_nogil
    callbacks[0].hess_prod        = &sleqp_func_hess_prod_nogil
  else:
    callbacks[0].set_value        = &sleqp_func_set
    callbacks[0].nonzeros         = &sleqp_func_nonzeros
    callbacks[0].set_error_bound  = &sleqp_dyn_set_error_bound
    callbacks[0].set_obj_weight   = &sleqp_dyn_set_obj_weight
    callbacks[0].set_cons_weights = &sleqp_dyn_set_cons_weights
    callbacks[0].eval             = &sleqp_dyn_eval
    callbacks[0].obj_grad         = &sleqp_func_obj_grad
    callbacks[0].cons_jac         = &sleqp_func_cons_jac
    callbacks[0].hess_prod        = &sleqp_func_hess_prod

  callbacks[0].func_free        = &sleqp_func_free


cdef object dyn_funcs = weakref.WeakSet()

cdef update_dyn_func_callbacks():
  cdef _Func func
  cdef csleqp.SleqpDynFuncCallbacks callbacks

  set_dyn_func_callbacks(&callbacks)

  for obj in dyn_funcs:
    func = <_Func> obj
    csleqp_call(csleqp.sleqp_dyn_func_set_callbacks(func.cfunc,
                                                    &callbacks))

cdef csleqp.SLEQP_RETCODE create_dyn_func(csleqp.SleqpFunc** cfunc,
                                          object func,
                                          int num_variables,
                                          int num_constraints):
  cdef csleqp.SleqpDynFuncCallbacks callbacks

  assert func is not None

  set_dyn_func_callbacks(&callbacks)

  return csleqp.sleqp_dyn_func_create(cfunc,
                                      &callbacks,
                                      num_variables,
                                      num_constraints,
                                      <void*> func)
