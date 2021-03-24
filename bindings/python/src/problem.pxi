#cython: language_level=3


cdef csleqp.SLEQP_RETCODE create_problem(csleqp.SleqpProblem** problem,
                                         csleqp.SleqpFunc* cfunc,
                                         np.ndarray var_lb,
                                         np.ndarray var_ub,
                                         np.ndarray cons_lb,
                                         np.ndarray cons_ub):

  cdef csleqp.SleqpSparseVec* var_lb_vec
  cdef csleqp.SleqpSparseVec* var_ub_vec

  cdef csleqp.SleqpSparseVec* cons_lb_vec
  cdef csleqp.SleqpSparseVec* cons_ub_vec

  cdef int num_variables = var_lb.shape[0]
  cdef int num_constraints = cons_lb.shape[0]

  assert cfunc != NULL

  try:

    csleqp_call(csleqp.sleqp_sparse_vector_create_empty(&var_lb_vec,
                                                        num_variables))

    csleqp_call(csleqp.sleqp_sparse_vector_create_empty(&var_ub_vec,
                                                        num_variables))

    csleqp_call(csleqp.sleqp_sparse_vector_create_empty(&cons_lb_vec,
                                                        num_constraints))

    csleqp_call(csleqp.sleqp_sparse_vector_create_empty(&cons_ub_vec,
                                                        num_constraints))

    array_to_sleqp_sparse_vec(var_lb, var_lb_vec)
    array_to_sleqp_sparse_vec(var_ub, var_ub_vec)
    array_to_sleqp_sparse_vec(cons_lb, cons_lb_vec)
    array_to_sleqp_sparse_vec(cons_ub, cons_ub_vec)


    csleqp_call(csleqp.sleqp_problem_create(problem,
                                            cfunc,
                                            var_lb_vec,
                                            var_ub_vec,
                                            cons_lb_vec,
                                            cons_ub_vec))

    return csleqp.SLEQP_OKAY

  except SLEQPError as error:
    return error.code

  except BaseException as exception:
    return csleqp.SLEQP_INTERNAL_ERROR

  finally:
    csleqp_call(csleqp.sleqp_sparse_vector_free(&cons_ub_vec))
    csleqp_call(csleqp.sleqp_sparse_vector_free(&cons_lb_vec))

    csleqp_call(csleqp.sleqp_sparse_vector_free(&var_ub_vec))
    csleqp_call(csleqp.sleqp_sparse_vector_free(&var_lb_vec))


cdef class _Problem:
  cdef csleqp.SleqpProblem* cproblem
  cdef _Func funcref

  def __cinit__(self):
    self.cproblem = NULL

  @staticmethod
  cdef _Problem create(csleqp.SleqpFunc* cfunc,
                       np.ndarray var_lb,
                       np.ndarray var_ub,
                       np.ndarray cons_lb,
                       np.ndarray cons_ub):

    cdef _Problem _problem = _Problem()

    cdef int num_variables = var_lb.shape[0]
    cdef int num_constraints = cons_lb.shape[0]

    csleqp_call(create_problem(&_problem.cproblem,
                               cfunc,
                               var_lb,
                               var_ub,
                               cons_lb,
                               cons_ub))

    _problem.funcref = _Func()
    _problem.funcref.set_func(cfunc)

    return _problem

  @property
  def num_variables(self) -> int:
    return self.cproblem.num_variables

  @property
  def num_constraints(self) -> int:
    return self.cproblem.num_constraints

  @property
  def var_lb(self) -> np.array:
    return sleqp_sparse_vec_to_array(self.cproblem.var_lb)

  @property
  def var_ub(self) -> np.array:
    return sleqp_sparse_vec_to_array(self.cproblem.var_ub)

  @property
  def cons_lb(self) -> np.array:
    return sleqp_sparse_vec_to_array(self.cproblem.cons_lb)

  @property
  def cons_ub(self) -> np.array:
    return sleqp_sparse_vec_to_array(self.cproblem.cons_ub)

  @property
  def hess_struct(self) -> HessianStruct:
    return HessianStruct(self.funcref)

  def __dealloc__(self):
    csleqp_call(csleqp.sleqp_problem_free(&self.cproblem))


cdef class Problem:
  cdef dict __dict__

  cdef _Problem problem
  cdef _Func funcref

  cdef object _func

  def __cinit__(self,
                object func,
                np.ndarray var_lb,
                np.ndarray var_ub,
                np.ndarray cons_lb,
                np.ndarray cons_ub):

    cdef int num_variables = var_lb.shape[0]
    cdef int num_constraints = cons_lb.shape[0]
    cdef csleqp.SleqpFunc* cfunc

    self.funcref = _Func()

    csleqp_call(create_func(&cfunc,
                            func,
                            num_variables,
                            num_constraints))

    assert cfunc != NULL

    self._func = func

    try:
      self.problem = _Problem.create(cfunc,
                                     var_lb,
                                     var_ub,
                                     cons_lb,
                                     cons_ub)

      self.funcref.set_func(cfunc)
      funcs.add(self.funcref)

    finally:
      csleqp_call(csleqp.sleqp_func_release(&cfunc))

    try:
      func.set_hessian_struct(self.hess_struct)
    except AttributeError:
      pass

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
    return self.problem.var_lb

  @property
  def var_ub(self) -> np.array:
    return self.problem.var_ub

  @property
  def cons_lb(self) -> np.array:
    return self.problem.cons_lb

  @property
  def cons_ub(self) -> np.array:
    return self.problem.cons_ub

  @property
  def hess_struct(self) -> HessianStruct:
    return self.problem.hess_struct

  def _get_problem(self):
    return self.problem

cdef class LSQProblem:
  cdef dict __dict__

  cdef _Problem problem
  cdef _Func funcref

  cdef object _func

  def __cinit__(self,
                object func,
                np.ndarray var_lb,
                np.ndarray var_ub,
                np.ndarray cons_lb,
                np.ndarray cons_ub,
                num_residuals,
                Params params,
                **properties):

    cdef int num_variables = var_lb.shape[0]
    cdef int num_constraints = cons_lb.shape[0]
    cdef csleqp.SleqpFunc* cfunc

    self.funcref = _Func()

    csleqp_call(create_lsq_func(&cfunc,
                                func,
                                num_variables,
                                num_constraints,
                                num_residuals,
                                properties.get('regularization', 0.),
                                params.params))

    assert cfunc != NULL

    self._func = func

    try:
      self.problem = _Problem.create(cfunc,
                                     var_lb,
                                     var_ub,
                                     cons_lb,
                                     cons_ub)

      self.funcref.set_func(cfunc)
      lsq_funcs.add(self.funcref)

    finally:
      csleqp_call(csleqp.sleqp_func_release(&cfunc))

    try:
      func.set_hessian_struct(self.hess_struct)
    except AttributeError:
      pass

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
    return self.problem.var_lb

  @property
  def var_ub(self) -> np.array:
    return self.problem.var_ub

  @property
  def cons_lb(self) -> np.array:
    return self.problem.cons_lb

  @property
  def cons_ub(self) -> np.array:
    return self.problem.cons_ub

  @property
  def hess_struct(self) -> HessianStruct:
    return self.problem.hess_struct

  def _get_problem(self):
    return self.problem
