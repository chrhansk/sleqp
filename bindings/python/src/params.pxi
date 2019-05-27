#cython: language_level=3

cdef class Params:
  cdef csleqp.SleqpParams* params
  cdef dict __dict__

  def __cinit__(self):
    csleqp_call(csleqp.sleqp_params_create(&self.params))

  def __init__(self, **values):
    self.props = ['zero_eps',
                  'eps',
                  'deriv_perturbation',
                  'deriv_tolerance',
                  'cauchy_tau',
                  'cauchy_eta',
                  'linesearch_tau',
                  'linesearch_eta',
                  'linesearch_cutoff',
                  'optimality_tolerance',
                  'accepted_reduction']

    for key, value in values.items():
      setattr(self, key, value)

  def __dealloc__(self):
    csleqp_call(csleqp.sleqp_params_free(&self.params))

  @property
  def zero_eps(self):
    return csleqp.sleqp_params_get_zero_eps(self.params)

  @property
  def eps(self):
    return csleqp.sleqp_params_get_eps(self.params)

  @property
  def deriv_perturbation(self):
    return csleqp.sleqp_params_get_deriv_perturbation(self.params)

  @property
  def deriv_tolerance(self):
    return csleqp.sleqp_params_get_deriv_tolerance(self.params)

  @property
  def cauchy_tau(self):
    return csleqp.sleqp_params_get_cauchy_tau(self.params)

  @property
  def cauchy_eta(self):
    return csleqp.sleqp_params_get_cauchy_eta(self.params)

  @property
  def linesearch_tau(self):
    return csleqp.sleqp_params_get_linesearch_tau(self.params)

  @property
  def linesearch_eta(self):
    return csleqp.sleqp_params_get_linesearch_eta(self.params)

  @property
  def linesearch_cutoff(self):
    return csleqp.sleqp_params_get_linesearch_cutoff(self.params)

  @property
  def optimality_tolerance(self):
    return csleqp.sleqp_params_get_optimality_tolerance(self.params)

  @property
  def accepted_reduction(self):
    return csleqp.sleqp_params_get_accepted_reduction(self.params)

  @zero_eps.setter
  def zero_eps(self, value):
    csleqp_call(csleqp.sleqp_params_set_zero_eps(self.params, value))

  @eps.setter
  def eps(self, value):
    csleqp_call(csleqp.sleqp_params_set_eps(self.params, value))

  @deriv_perturbation.setter
  def deriv_perturbation(self, value):
    csleqp_call(csleqp.sleqp_params_set_deriv_perturbation(self.params, value))

  @deriv_tolerance.setter
  def deriv_tolerance(self, value):
    csleqp_call(csleqp.sleqp_params_set_deriv_tolerance(self.params, value))

  @cauchy_tau.setter
  def cauchy_tau(self, value):
    csleqp_call(csleqp.sleqp_params_set_cauchy_tau(self.params, value))

  @cauchy_eta.setter
  def cauchy_eta(self, value):
    csleqp_call(csleqp.sleqp_params_set_cauchy_eta(self.params, value))

  @linesearch_tau.setter
  def linesearch_tau(self, value):
    csleqp_call(csleqp.sleqp_params_set_linesearch_tau(self.params, value))

  @linesearch_eta.setter
  def linesearch_eta(self, value):
    csleqp_call(csleqp.sleqp_params_set_linesearch_eta(self.params, value))

  @linesearch_cutoff.setter
  def linesearch_cutoff(self, value):
    csleqp_call(csleqp.sleqp_params_set_linesearch_cutoff(self.params, value))

  @optimality_tolerance.setter
  def optimality_tolerance(self, value):
    csleqp_call(csleqp.sleqp_params_set_optimality_tolerance(self.params, value))

  @accepted_reduction.setter
  def accepted_reduction(self, value):
    csleqp_call(csleqp.sleqp_params_set_accepted_reduction(self.params, value))

  def values(self):
    return {key: getattr(self, key) for key in self.props}

  def __str__(self):
    return 'Params: {0}'.format(self.values())

  def __repr__(self):
    return 'Params({0})'.format(repr(self.values()))
