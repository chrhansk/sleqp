cdef csleqp.SLEQP_RETCODE sleqp_func_set(csleqp.SleqpSparseVec* x,
                                         int num_variables,
                                         int* func_grad_nnz,
                                         int* cons_val_nnz,
                                         int* cons_jac_nnz,
                                         void* func_data):

  try:
    x_array = sleqp_sparse_vec_to_array(x)
    (<object> func_data).set_value(x_array)

    func_grad_nnz[0] = (<object> func_data).func_grad_nnz()
    cons_val_nnz[0] = (<object> func_data).cons_val_nnz()
    cons_jac_nnz[0] =  (<object> func_data).cons_jac_nnz()
  except Exception:
    traceback.print_exc()
    return csleqp.SLEQP_INTERNAL_ERROR

  return csleqp.SLEQP_OKAY

cdef csleqp.SLEQP_RETCODE sleqp_func_eval(int num_variables,
                                          csleqp.SleqpSparseVec* cons_indices,
                                          double* func_val,
                                          csleqp.SleqpSparseVec* func_grad,
                                          csleqp.SleqpSparseVec* cons_vals,
                                          csleqp.SleqpSparseMatrix* cons_jac,
                                          void* func_data):
  try:

    if func_val:
      func_val[0] = (<object> func_data).func_val()

    if func_grad:
      grad_array = (<object> func_data).func_grad()
      array_to_sleqp_sparse_vec(grad_array, func_grad)

    if cons_vals:
      cons_array = (<object> func_data).cons_vals()
      array_to_sleqp_sparse_vec(cons_array, cons_vals)

    if cons_jac:
      cons_jac_mat = (<object> func_data).cons_jac()
      matrix_to_sleqp_sparse_matrix(cons_jac_mat, cons_jac)


  except Exception:
    traceback.print_exc()
    return csleqp.SLEQP_INTERNAL_ERROR

  return csleqp.SLEQP_OKAY

cdef csleqp.SLEQP_RETCODE sleqp_func_hess_product(int num_variables,
                                                  double* func_dual,
                                                  csleqp.SleqpSparseVec* direction,
                                                  csleqp.SleqpSparseVec* cons_dual,
                                                  csleqp.SleqpSparseVec* product,
                                                  void* func_data):

  try:

    assert cons_dual
    assert direction
    assert product

    f_dual = func_dual[0] if func_dual else 0.
    direction_array = sleqp_sparse_vec_to_array(direction)
    cons_dual_array = sleqp_sparse_vec_to_array(cons_dual)

    product_array = (<object> func_data).hess_prod(f_dual, direction_array, cons_dual_array)

    array_to_sleqp_sparse_vec(product_array, product)

  except Exception:
    traceback.print_exc()
    return csleqp.SLEQP_INTERNAL_ERROR

  return csleqp.SLEQP_OKAY


cdef class Func:
  cdef csleqp.SleqpFunc* func
  cdef int num_variables
  cdef int num_constraints

  def __cinit__(self, int num_variables, int num_constraints):
    csleqp_call(csleqp.sleqp_func_create(&self.func,
                                         &sleqp_func_set,
                                         &sleqp_func_eval,
                                         &sleqp_func_hess_product,
                                         num_variables,
                                         <void*> self))

    self.num_variables = num_variables
    self.num_constraints = num_constraints

    assert(self.func)

  cpdef set_value(self, x):
    pass

  cpdef func_val(self):
    pass

  cpdef func_grad_nnz(self):
    return 0

  cpdef cons_val_nnz(self):
    return 0

  cpdef cons_jac_nnz(self):
    return 0

  cpdef func_grad(self):
    pass

  cpdef cons_vals(self):
    pass

  cpdef cons_jac(self):
    pass

  cpdef hess_prod(self, func_dual, direction, cons_dual):
    pass

  @property
  def num_variables(self):
      return self.num_variables

  @property
  def num_constraints(self):
      return self.num_constraints

  def __dealloc__(self):
    csleqp_call(csleqp.sleqp_func_free(&self.func))
