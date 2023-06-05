#cython: language_level=3


cdef csleqp.SLEQP_RETCODE create_problem(csleqp.SleqpProblem** problem,
                                         csleqp.SleqpFunc* cfunc,
                                         np.ndarray var_lb,
                                         np.ndarray var_ub,
                                         np.ndarray general_lb,
                                         np.ndarray general_ub,
                                         csleqp.SleqpSettings* csettings,
                                         object linear_coeffs = None,
                                         np.ndarray linear_lb = None,
                                         np.ndarray linear_ub = None):

  assert var_lb is not None
  assert var_ub is not None

  assert general_lb is not None
  assert general_ub is not None

  cdef csleqp.SleqpVec* var_lb_vec
  cdef csleqp.SleqpVec* var_ub_vec

  cdef csleqp.SleqpVec* general_lb_vec
  cdef csleqp.SleqpVec* general_ub_vec

  cdef csleqp.SleqpVec* linear_lb_vec
  cdef csleqp.SleqpVec* linear_ub_vec

  cdef csleqp.SleqpMat* linear_coeffs_mat

  cdef int num_vars = var_lb.shape[0]
  cdef int num_cons = general_lb.shape[0]

  assert cfunc != NULL

  cdef int num_linear_constraints = 0

  if linear_coeffs is not None:
    num_linear_constraints = linear_coeffs.shape[0]

  try:

    csleqp_call(csleqp.sleqp_vec_create_empty(&var_lb_vec,
                                              num_vars))

    csleqp_call(csleqp.sleqp_vec_create_empty(&var_ub_vec,
                                              num_vars))

    csleqp_call(csleqp.sleqp_vec_create_empty(&general_lb_vec,
                                              num_cons))

    csleqp_call(csleqp.sleqp_vec_create_empty(&general_ub_vec,
                                              num_cons))

    csleqp_call(csleqp.sleqp_mat_create(&linear_coeffs_mat,
                                        num_linear_constraints,
                                        num_vars,
                                        0))

    csleqp_call(csleqp.sleqp_vec_create_empty(&linear_lb_vec,
                                              num_linear_constraints))

    csleqp_call(csleqp.sleqp_vec_create_empty(&linear_ub_vec,
                                              num_linear_constraints))

    array_to_sleqp_sparse_vec(var_lb, var_lb_vec)
    array_to_sleqp_sparse_vec(var_ub, var_ub_vec)
    array_to_sleqp_sparse_vec(general_lb, general_lb_vec)
    array_to_sleqp_sparse_vec(general_ub, general_ub_vec)

    if linear_lb is not None:
      array_to_sleqp_sparse_vec(linear_lb, linear_lb_vec)

    if linear_ub is not None:
      array_to_sleqp_sparse_vec(linear_ub, linear_ub_vec)

    if linear_coeffs is not None:
      csleqp_call(matrix_to_sleqp_sparse_matrix(linear_coeffs,
                                                linear_coeffs_mat))

    csleqp_call(csleqp.sleqp_problem_create(problem,
                                            cfunc,
                                            var_lb_vec,
                                            var_ub_vec,
                                            general_lb_vec,
                                            general_ub_vec,
                                            linear_coeffs_mat,
                                            linear_lb_vec,
                                            linear_ub_vec,
                                            csettings))

    return csleqp.SLEQP_OKAY

  except BaseException as exception:
    return csleqp.SLEQP_ERROR

  finally:
    csleqp_call(csleqp.sleqp_vec_free(&linear_ub_vec))
    csleqp_call(csleqp.sleqp_vec_free(&linear_lb_vec))

    csleqp_call(csleqp.sleqp_mat_release(&linear_coeffs_mat))

    csleqp_call(csleqp.sleqp_vec_free(&general_ub_vec))
    csleqp_call(csleqp.sleqp_vec_free(&general_lb_vec))

    csleqp_call(csleqp.sleqp_vec_free(&var_ub_vec))
    csleqp_call(csleqp.sleqp_vec_free(&var_lb_vec))


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
                       np.ndarray cons_ub,
                       csleqp.SleqpSettings* csettings,
                       object linear_coeffs = None,
                       np.ndarray linear_lb = None,
                       np.ndarray linear_ub = None):

    cdef _Problem _problem = _Problem()

    cdef int num_vars = var_lb.shape[0]
    cdef int num_cons = cons_lb.shape[0]

    csleqp_call(create_problem(&_problem.cproblem,
                               cfunc,
                               var_lb,
                               var_ub,
                               cons_lb,
                               cons_ub,
                               csettings,
                               linear_coeffs,
                               linear_lb,
                               linear_ub))

    _problem.funcref = _Func()
    _problem.funcref.set_func(cfunc)

    return _problem

  @property
  def num_vars(self) -> int:
    return csleqp.sleqp_problem_num_vars(self.cproblem)

  @property
  def num_cons(self) -> int:
    return csleqp.sleqp_problem_num_cons(self.cproblem)

  @property
  def var_lb(self) -> np.array:
    return sleqp_sparse_vec_to_array(csleqp.sleqp_problem_vars_lb(self.cproblem))

  @property
  def var_ub(self) -> np.array:
    return sleqp_sparse_vec_to_array(csleqp.sleqp_problem_vars_ub(self.cproblem))

  @property
  def cons_lb(self) -> np.array:
    return sleqp_sparse_vec_to_array(csleqp.sleqp_problem_cons_lb(self.cproblem))

  @property
  def cons_ub(self) -> np.array:
    return sleqp_sparse_vec_to_array(csleqp.sleqp_problem_cons_ub(self.cproblem))

  @property
  def hess_struct(self) -> HessianStruct:
    return HessianStruct(self.funcref)

  def __dealloc__(self):
    csleqp_call(csleqp.sleqp_problem_release(&self.cproblem))


cdef class Problem:
  cdef dict __dict__

  cdef _Problem problem
  cdef _Func funcref

  cdef Settings _settings

  cdef object _func

  def __cinit__(self,
                object func,
                np.ndarray var_lb,
                np.ndarray var_ub,
                np.ndarray cons_lb,
                np.ndarray cons_ub,
                Settings settings=None,
                **properties):

    cdef int num_vars = var_lb.shape[0]
    cdef int num_cons = cons_lb.shape[0]
    cdef csleqp.SleqpFunc* cfunc
    cdef csleqp.SleqpSettings* csettings

    if settings is None:
      self._settings = Settings()
    else:
      self._settings = settings

    csettings = <csleqp.SleqpSettings*> self._settings.settings

    self.funcref = _Func()

    csleqp_call(create_func(&cfunc,
                            func,
                            num_vars,
                            num_cons))

    assert cfunc != NULL

    self._func = func

    try:
      self.problem = _Problem.create(cfunc,
                                     var_lb,
                                     var_ub,
                                     cons_lb,
                                     cons_ub,
                                     csettings,
                                     properties.get('linear_coeffs', None),
                                     properties.get('linear_lb', None),
                                     properties.get('linear_ub', None))

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
  def num_vars(self) -> int:
    return self.problem.num_vars

  @property
  def num_cons(self) -> int:
    return self.problem.num_cons

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

  cdef Settings _settings

  cdef object _func

  def __cinit__(self,
                object func,
                np.ndarray var_lb,
                np.ndarray var_ub,
                np.ndarray cons_lb,
                np.ndarray cons_ub,
                num_residuals,
                Settings settings=None,
                **properties):

    cdef int num_vars = var_lb.shape[0]
    cdef int num_cons = cons_lb.shape[0]
    cdef csleqp.SleqpFunc* cfunc
    cdef csleqp.SleqpSettings* csettings

    if settings is None:
      self._settings = Settings()
    else:
      self._settings = settings

    csettings = <csleqp.SleqpSettings*> self._settings.settings

    self.funcref = _Func()

    csleqp_call(create_lsq_func(&cfunc,
                                func,
                                num_vars,
                                num_cons,
                                num_residuals,
                                properties.get('regularization', 0.),
                                csettings))

    assert cfunc != NULL

    self._func = func

    try:
      self.problem = _Problem.create(cfunc,
                                     var_lb,
                                     var_ub,
                                     cons_lb,
                                     cons_ub,
                                     csettings,
                                     properties.get('linear_coeffs', None),
                                     properties.get('linear_lb', None),
                                     properties.get('linear_ub', None))

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
  def num_vars(self) -> int:
    return self.problem.num_vars

  @property
  def num_cons(self) -> int:
    return self.problem.num_cons

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


cdef class DynProblem:
  cdef dict __dict__

  cdef _Problem problem
  cdef _Func funcref

  cdef Settings _settings

  cdef object _func

  def __cinit__(self,
                object func,
                np.ndarray var_lb,
                np.ndarray var_ub,
                np.ndarray cons_lb,
                np.ndarray cons_ub,
                Settings settings=None,
                **properties):
    cdef int num_vars = var_lb.shape[0]
    cdef int num_cons = cons_lb.shape[0]
    cdef csleqp.SleqpFunc* cfunc
    cdef csleqp.SleqpSettings* csettings

    if settings is None:
      self._settings = Settings()
    else:
      self._settings = settings

    csettings = <csleqp.SleqpSettings*> self._settings.settings

    self.funcref = _Func()

    csleqp_call(create_dyn_func(&cfunc,
                                func,
                                num_vars,
                                num_cons))

    assert cfunc != NULL
    self._func = func

    try:
      self.problem = _Problem.create(cfunc,
                                     var_lb,
                                     var_ub,
                                     cons_lb,
                                     cons_ub,
                                     csettings,
                                     properties.get('linear_coeffs', None),
                                     properties.get('linear_lb', None),
                                     properties.get('linear_ub', None))

      self.funcref.set_func(cfunc)
      dyn_funcs.add(self.funcref)

      pass
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
  def num_vars(self) -> int:
    return self.problem.num_vars

  @property
  def num_cons(self) -> int:
    return self.problem.num_cons

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
