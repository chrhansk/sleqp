#cython: language_level=3

cdef csleqp.SLEQP_RETCODE sleqp_func_set(csleqp.SleqpSparseVec* x,
                                         csleqp.SLEQP_VALUE_REASON reason,
                                         int num_variables,
                                         int* func_grad_nnz,
                                         int* cons_val_nnz,
                                         int* cons_jac_nnz,
                                         void* func_data):

  try:
    x_array = sleqp_sparse_vec_to_array(x)

    func = (<object> func_data)

    func.set_value(x_array, ValueReason(reason))

    func_grad_nnz[0] = func.func_grad_nnz()
    cons_val_nnz[0] = func.cons_val_nnz()
    cons_jac_nnz[0] =  func.cons_jac_nnz()

  except BaseException as exception:
    func.call_exception = exception
    return csleqp.SLEQP_INTERNAL_ERROR

  return csleqp.SLEQP_OKAY

cdef csleqp.SLEQP_RETCODE sleqp_func_eval(int num_variables,
                                          const csleqp.SleqpSparseVec* cons_indices,
                                          double* func_val,
                                          csleqp.SleqpSparseVec* func_grad,
                                          csleqp.SleqpSparseVec* cons_vals,
                                          csleqp.SleqpSparseMatrix* cons_jac,
                                          void* func_data):
  try:

    func = (<object> func_data)

    if func_val:
      func_val[0] = func.func_val()

    if func_grad:
      grad_array = func.func_grad()
      csleqp_call(array_to_sleqp_sparse_vec(grad_array, func_grad))

    if cons_vals:
      cons_array = func.cons_vals()
      csleqp_call(array_to_sleqp_sparse_vec(cons_array, cons_vals))

    if cons_jac:
      cons_jac_mat = func.cons_jac()
      csleqp_call(matrix_to_sleqp_sparse_matrix(cons_jac_mat, cons_jac))


  except BaseException as exception:
    func.call_exception = exception
    return csleqp.SLEQP_INTERNAL_ERROR

  return csleqp.SLEQP_OKAY

cdef csleqp.SLEQP_RETCODE sleqp_func_hess_product(int num_variables,
                                                  const double* func_dual,
                                                  const csleqp.SleqpSparseVec* direction,
                                                  const csleqp.SleqpSparseVec* cons_dual,
                                                  csleqp.SleqpSparseVec* product,
                                                  void* func_data):

  try:

    assert cons_dual
    assert direction
    assert product

    func = (<object> func_data)

    f_dual = func_dual[0] if func_dual else 0.
    direction_array = sleqp_sparse_vec_to_array(direction)
    cons_dual_array = sleqp_sparse_vec_to_array(cons_dual)

    product_array = func.hess_prod(f_dual, direction_array, cons_dual_array)

    csleqp_call(array_to_sleqp_sparse_vec(product_array, product))

  except BaseException as exception:
    func.call_exception = exception
    return csleqp.SLEQP_INTERNAL_ERROR

  return csleqp.SLEQP_OKAY

cdef csleqp.SLEQP_RETCODE sleqp_func_free(void* func_data):
  return csleqp.SLEQP_OKAY


cdef class Func:
  cdef csleqp.SleqpFunc* func
  cdef int num_variables
  cdef int num_constraints

  cdef public object call_exception

  def __cinit__(self,
                int num_variables,
                int num_constraints,
                *args,
                **keywords):

    cdef csleqp.SleqpFuncCallbacks callbacks

    callbacks.set_value = &sleqp_func_set
    callbacks.func_eval = &sleqp_func_eval
    callbacks.hess_prod = &sleqp_func_hess_product
    callbacks.func_free = &sleqp_func_free

    csleqp_call(csleqp.sleqp_func_create(&self.func,
                                         &callbacks,
                                         num_variables,
                                         <void*> self))

    self.num_variables = num_variables
    self.num_constraints = num_constraints

    self.call_exception = None

    assert(self.func)

  cpdef set_value(self, x, reason) -> None:
    pass

  cpdef func_val(self) -> float:
    pass

  cpdef func_grad_nnz(self) -> int:
    return 0

  cpdef cons_val_nnz(self) -> int:
    return 0

  cpdef cons_jac_nnz(self) -> int:
    return 0

  cpdef func_grad(self) -> np.array:
    pass

  cpdef cons_vals(self) -> np.array:
    pass

  cpdef cons_jac(self) -> typing.Union[np.array, scipy.sparse.csc_matrix]:
    pass

  cpdef hess_prod(self,
                  func_dual: float,
                  direction: np.array,
                  cons_dual: np.array) -> np.array:
    pass

  @property
  def hess_struct(self) -> HessianStruct:
      return HessianStruct(self)

  @property
  def num_variables(self) -> int:
      return self.num_variables

  @property
  def num_constraints(self) -> int:
      return self.num_constraints

  def __dealloc__(self):
    csleqp_call(csleqp.sleqp_func_release(&self.func))
