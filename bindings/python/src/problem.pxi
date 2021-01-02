#cython: language_level=3

cdef class Problem:

  cdef dict __dict__

  cdef csleqp.SleqpProblem* problem

  cdef csleqp.SleqpSparseVec* var_lb
  cdef csleqp.SleqpSparseVec* var_ub

  cdef csleqp.SleqpSparseVec* cons_lb
  cdef csleqp.SleqpSparseVec* cons_ub

  cdef object _func
  cdef Params params

  def __cinit__(self,
                object func,
                Params params,
                np.ndarray var_lb,
                np.ndarray var_ub,
                np.ndarray cons_lb,
                np.ndarray cons_ub):

    cdef csleqp.SleqpFunc* cfunc = NULL

    num_constraints = cons_lb.shape[0]
    num_variables = var_lb.shape[0]

    if isinstance(func, Func):
        cfunc = (<Func> func).func
    elif isinstance(func, LSQFunc):
        cfunc = (<LSQFunc> func).func

    assert cfunc != NULL, "Invalid function type"

    csleqp_call(csleqp.sleqp_sparse_vector_create_empty(&self.var_lb,
                                                        num_variables))

    csleqp_call(csleqp.sleqp_sparse_vector_create_empty(&self.var_ub,
                                                        num_variables))

    csleqp_call(csleqp.sleqp_sparse_vector_create_empty(&self.cons_lb,
                                                        num_constraints))

    csleqp_call(csleqp.sleqp_sparse_vector_create_empty(&self.cons_ub,
                                                        num_constraints))

    array_to_sleqp_sparse_vec(var_lb, self.var_lb)
    array_to_sleqp_sparse_vec(var_ub, self.var_ub)
    array_to_sleqp_sparse_vec(cons_lb, self.cons_lb)
    array_to_sleqp_sparse_vec(cons_ub, self.cons_ub)

    self._func = func
    self.params = params

    csleqp_call(csleqp.sleqp_problem_create(&self.problem,
                                            cfunc,
                                            params.params,
                                            self.var_lb,
                                            self.var_ub,
                                            self.cons_lb,
                                            self.cons_ub))

  @property
  def func(self):
      return self._func

  @property
  def num_variables(self) -> int:
    return self.problem.num_variables

  @property
  def num_constraints(self) -> int:
    return self.problem.num_constraints

  @property
  def var_lb(self) -> np.array:
    return sleqp_sparse_vec_to_array(self.problem.var_lb)

  @property
  def var_ub(self) -> np.array:
    return sleqp_sparse_vec_to_array(self.problem.var_ub)

  @property
  def cons_lb(self) -> np.array:
    return sleqp_sparse_vec_to_array(self.problem.cons_lb)

  @property
  def cons_ub(self) -> np.array:
    return sleqp_sparse_vec_to_array(self.problem.cons_ub)

  def __dealloc__(self):
    csleqp_call(csleqp.sleqp_problem_free(&self.problem))

    csleqp_call(csleqp.sleqp_sparse_vector_free(&self.cons_ub))
    csleqp_call(csleqp.sleqp_sparse_vector_free(&self.cons_lb))

    csleqp_call(csleqp.sleqp_sparse_vector_free(&self.var_ub))
    csleqp_call(csleqp.sleqp_sparse_vector_free(&self.var_lb))
