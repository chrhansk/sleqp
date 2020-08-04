#cython: language_level=3

cdef class Options:
  cdef csleqp.SleqpOptions* options
  cdef dict __dict__

  def __cinit__(self):
    csleqp_call(csleqp.sleqp_options_create(&self.options))

  def __init__(self, **values):
    self.props = ['perform_newton_step',
                  'perform_soc',
                  'deriv_check',
                  'hessian_eval',
                  'dual_estimation_type',
                  'quasi_newton_num_iterates']

    for key, value in values.items():
      setattr(self, key, value)

  def __dealloc__(self):
    csleqp_call(csleqp.sleqp_options_free(&self.options))

  @property
  def perform_newton_step(self):
    return csleqp.sleqp_options_get_perform_newton_step(self.options)

  @property
  def perform_soc(self):
    return csleqp.sleqp_options_get_perform_soc(self.options)

  @property
  def deriv_check(self):
    return DerivCheck(csleqp.sleqp_options_get_deriv_check(self.options))

  @property
  def hessian_eval(self):
    return HessianEval(csleqp.sleqp_options_get_hessian_eval(self.options))

  @property
  def dual_estimation_type(self):
    return DualEstimationType(csleqp.sleqp_options_get_dual_estimation_type(self.options))

  @property
  def quasi_newton_num_iterates(self):
    return csleqp.sleqp_options_get_quasi_newton_num_iterates(self.options)

  @perform_newton_step.setter
  def perform_newton_step(self, value):
    csleqp.sleqp_options_set_perform_newton_step(self.options, value)

  @perform_soc.setter
  def perform_soc(self, value):
    csleqp.sleqp_options_set_perform_soc(self.options, value)

  @deriv_check.setter
  def deriv_check(self, value):
    csleqp_call(csleqp.sleqp_options_set_deriv_check(self.options, value.value))

  @hessian_eval.setter
  def hessian_eval(self, value):
    csleqp_call(csleqp.sleqp_options_set_hessian_eval(self.options, value.value))

  @dual_estimation_type.setter
  def dual_estimation_type(self, value):
    csleqp_call(csleqp.sleqp_options_set_dual_estimation_type(self.options, value.value))

  @quasi_newton_num_iterates.setter
  def quasi_newton_num_iterates(self, int value):
    csleqp_call(csleqp.sleqp_options_set_quasi_newton_num_iterates(self.options, value))

  def values(self):
    return {key: getattr(self, key) for key in self.props}

  def __str__(self):
    return 'Options: {0}'.format(self.values())

  def __repr__(self):
    return 'Options({0})'.format(repr(self.values()))
