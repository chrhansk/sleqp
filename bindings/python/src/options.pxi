#cython: language_level=3

cdef class Options:
  cdef csleqp.SleqpOptions* options
  cdef dict __dict__

  def __cinit__(self):
    csleqp_call(csleqp.sleqp_options_create(&self.options))

  def __init__(self, **values):
    self.props = ['deriv_check',
                  'hessian_eval']

    for key, value in values.items():
      setattr(self, key, value)

  def __dealloc__(self):
    csleqp_call(csleqp.sleqp_options_free(&self.options))

  @property
  def deriv_check(self):
    return DerivCheck(csleqp.sleqp_options_get_deriv_check(self.options))

  @property
  def hessian_eval(self):
    return HessianEval(csleqp.sleqp_options_get_hessian_eval(self.options))

  @deriv_check.setter
  def deriv_check(self, value):
    csleqp_call(csleqp.sleqp_options_set_deriv_check(self.options, value.value))

  @hessian_eval.setter
  def hessian_eval(self, value):
    csleqp_call(csleqp.sleqp_options_set_hessian_eval(self.options, value.value))

  def values(self):
    return {key: getattr(self, key) for key in self.props}

  def __str__(self):
    return 'Options: {0}'.format(self.values())

  def __repr__(self):
    return 'Options({0})'.format(repr(self.values()))
