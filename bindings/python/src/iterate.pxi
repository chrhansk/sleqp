#cython: language_level=3

cdef class Iterate:
  cdef csleqp.SleqpIterate* iterate

  cdef dict __dict__

  cdef _release(self):
    if self.iterate:
      csleqp_call(csleqp.sleqp_iterate_release(&self.iterate))

  cdef _set_iterate(self, csleqp.SleqpIterate* iterate):
    self._release()

    csleqp_call(csleqp.sleqp_iterate_capture(iterate))

    self.iterate = iterate

  def __dealloc__(self):
    self._release()

  @property
  def primal(self) -> np.array:
    assert self.iterate

    return sleqp_sparse_vec_to_array(csleqp.sleqp_iterate_get_primal(self.iterate))

  @property
  def vars_dual(self) -> np.array:
    assert self.iterate

    return sleqp_sparse_vec_to_array(csleqp.sleqp_iterate_get_vars_dual(self.iterate))

  @property
  def cons_dual(self) -> np.array:
    assert self.iterate

    return sleqp_sparse_vec_to_array(csleqp.sleqp_iterate_get_cons_dual(self.iterate))

  @property
  def working_set(self) -> WorkingSet:
    cdef csleqp.SleqpWorkingSet* _working_set
    cdef WorkingSet working_set = WorkingSet()

    assert self.iterate

    _working_set = csleqp.sleqp_iterate_get_working_set(self.iterate)

    working_set._set_working_set(_working_set)

    return working_set

  @property
  def cons_jac(self) -> scipy.sparse.csc_matrix:
    assert self.iterate

    return sleqp_sparse_matrix_to_scipy(csleqp.sleqp_iterate_get_cons_jac(self.iterate))

  @property
  def func_val(self) -> float:
    assert self.iterate

    return csleqp.sleqp_iterate_get_func_val(self.iterate)

  @property
  def cons_val(self) -> np.array:
    assert self.iterate

    return sleqp_sparse_vec_to_array(csleqp.sleqp_iterate_get_cons_val(self.iterate))
