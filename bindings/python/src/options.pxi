#cython: language_level=3

cdef class Options:
  cdef csleqp.SleqpOptions* options
  cdef dict __dict__

  def __cinit__(self):
    csleqp_call(csleqp.sleqp_options_create(&self.options))

  def __init__(self, **values):
    self.props = ['perform_newton_step',
                  'perform_soc',
                  'use_quadratic_model',
                  'deriv_check',
                  'hessian_eval',
                  'dual_estimation_type',
                  'quasi_newton_num_iterates',
                  'max_newton_iterations',
                  'bfgs_sizing',
                  'tr_solver']

    for key, value in values.items():
      self._set_prop(key, value)

  cdef _set_prop(self, name, value):
    if not name in self.props:
      raise AttributeError("Invalid property {0}".format(name))

    setattr(self, name, value)

  def __dealloc__(self):
    csleqp_call(csleqp.sleqp_options_release(&self.options))

  @property
  def perform_newton_step(self) -> bool:
    return csleqp.sleqp_options_get_perform_newton_step(self.options)

  @property
  def perform_soc(self) -> bool:
    return csleqp.sleqp_options_get_perform_soc(self.options)

  @property
  def use_quadratic_model(self) -> bool:
    return csleqp.sleqp_options_get_use_quadratic_model(self.options)

  @property
  def deriv_check(self) -> DerivCheck:
    return DerivCheck(csleqp.sleqp_options_get_deriv_check(self.options))

  @property
  def hessian_eval(self) -> HessianEval:
    return HessianEval(csleqp.sleqp_options_get_hessian_eval(self.options))

  @property
  def dual_estimation_type(self) -> DualEstimationType:
    return DualEstimationType(csleqp.sleqp_options_get_dual_estimation_type(self.options))

  @property
  def quasi_newton_num_iterates(self) -> int:
    return csleqp.sleqp_options_get_quasi_newton_num_iterates(self.options)

  @property
  def max_newton_iterations(self):
      cdef int iter = csleqp.sleqp_options_get_max_newton_iterations(self.options)

      return iter if iter != csleqp.SLEQP_NONE else None

  @property
  def bfgs_sizing(self) -> Sizing:
      cdef int sizing = csleqp.sleqp_options_get_bfgs_sizing(self.options)
      return Sizing(sizing)

  @property
  def tr_solver(self) -> TRSolver:
      cdef int tr_solver = csleqp.sleqp_options_get_tr_solver(self.options)
      return TRSolver(tr_solver)

  @perform_newton_step.setter
  def perform_newton_step(self, value: bool) -> None:
    csleqp_call(csleqp.sleqp_options_set_perform_newton_step(self.options, value))

  @perform_soc.setter
  def perform_soc(self, value: bool) -> None:
    csleqp_call(csleqp.sleqp_options_set_perform_soc(self.options, value))

  @use_quadratic_model.setter
  def use_quadratic_model(self, value: bool) -> None:
    csleqp_call(csleqp.sleqp_options_set_use_quadratic_model(self.options, value))

  @deriv_check.setter
  def deriv_check(self, value) -> None:
    csleqp_call(csleqp.sleqp_options_set_deriv_check(self.options, value.value))

  @hessian_eval.setter
  def hessian_eval(self, value) -> None:
    csleqp_call(csleqp.sleqp_options_set_hessian_eval(self.options, value.value))

  @dual_estimation_type.setter
  def dual_estimation_type(self, value) -> None:
    csleqp_call(csleqp.sleqp_options_set_dual_estimation_type(self.options, value.value))

  @quasi_newton_num_iterates.setter
  def quasi_newton_num_iterates(self, value: int) -> None:
    csleqp_call(csleqp.sleqp_options_set_quasi_newton_num_iterates(self.options, value))

  @max_newton_iterations.setter
  def max_newton_iterations(self, value):
    if value is None:
      value = csleqp.SLEQP_NONE

    csleqp_call(csleqp.sleqp_options_set_max_newton_iterations(self.options, value))

  @bfgs_sizing.setter
  def bfgs_sizing(self, value):
    csleqp_call(csleqp.sleqp_options_set_bfgs_sizing(self.options, value.value))

  @tr_solver.setter
  def tr_solver(self, value):
    csleqp_call(csleqp.sleqp_options_set_tr_solver(self.options, value.value))

  def values(self) -> set:
    return {key: getattr(self, key) for key in self.props}

  def __str__(self) -> str:
    return 'Options: {0}'.format(self.values())

  def __repr__(self) -> str:
    return 'Options({0})'.format(repr(self.values()))
