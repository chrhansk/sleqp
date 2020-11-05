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
                  'accepted_reduction',
                  'deadpoint_bound',
                  'newton_relative_tolerance']

    for key, value in values.items():
      setattr(self, key, value)

  def __dealloc__(self):
    csleqp_call(csleqp.sleqp_params_release(&self.params))

  @property
  def zero_eps(self) -> float:
    return csleqp.sleqp_params_get_zero_eps(self.params)

  @property
  def eps(self) -> float:
    return csleqp.sleqp_params_get_eps(self.params)

  @property
  def deriv_perturbation(self) -> float:
    return csleqp.sleqp_params_get_deriv_perturbation(self.params)

  @property
  def deriv_tolerance(self) -> float:
    return csleqp.sleqp_params_get_deriv_tolerance(self.params)

  @property
  def cauchy_tau(self) -> float:
    return csleqp.sleqp_params_get_cauchy_tau(self.params)

  @property
  def cauchy_eta(self) -> float:
    return csleqp.sleqp_params_get_cauchy_eta(self.params)

  @property
  def linesearch_tau(self) -> float:
    return csleqp.sleqp_params_get_linesearch_tau(self.params)

  @property
  def linesearch_eta(self) -> float:
    return csleqp.sleqp_params_get_linesearch_eta(self.params)

  @property
  def linesearch_cutoff(self) -> float:
    return csleqp.sleqp_params_get_linesearch_cutoff(self.params)

  @property
  def optimality_tolerance(self) -> float:
    return csleqp.sleqp_params_get_optimality_tolerance(self.params)

  @property
  def accepted_reduction(self) -> float:
    return csleqp.sleqp_params_get_accepted_reduction(self.params)

  @property
  def deadpoint_bound(self) -> float:
    return csleqp.sleqp_params_get_deadpoint_bound(self.params)

  @property
  def newton_relative_tolerance(self) -> float:
    return csleqp.sleqp_params_get_newton_relative_tolerance(self.params)

  @zero_eps.setter
  def zero_eps(self, value: float) -> None:
    csleqp_call(csleqp.sleqp_params_set_zero_eps(self.params, value))

  @eps.setter
  def eps(self, value: float) -> None:
    csleqp_call(csleqp.sleqp_params_set_eps(self.params, value))

  @deriv_perturbation.setter
  def deriv_perturbation(self, value: float) -> None:
    csleqp_call(csleqp.sleqp_params_set_deriv_perturbation(self.params, value))

  @deriv_tolerance.setter
  def deriv_tolerance(self, value: float) -> None:
    csleqp_call(csleqp.sleqp_params_set_deriv_tolerance(self.params, value))

  @cauchy_tau.setter
  def cauchy_tau(self, value: float) -> None:
    csleqp_call(csleqp.sleqp_params_set_cauchy_tau(self.params, value))

  @cauchy_eta.setter
  def cauchy_eta(self, value: float) -> None:
    csleqp_call(csleqp.sleqp_params_set_cauchy_eta(self.params, value))

  @linesearch_tau.setter
  def linesearch_tau(self, value: float) -> None:
    csleqp_call(csleqp.sleqp_params_set_linesearch_tau(self.params, value))

  @linesearch_eta.setter
  def linesearch_eta(self, value: float) -> None:
    csleqp_call(csleqp.sleqp_params_set_linesearch_eta(self.params, value))

  @linesearch_cutoff.setter
  def linesearch_cutoff(self, value: float) -> None:
    csleqp_call(csleqp.sleqp_params_set_linesearch_cutoff(self.params, value))

  @optimality_tolerance.setter
  def optimality_tolerance(self, value: float) -> None:
    csleqp_call(csleqp.sleqp_params_set_optimality_tolerance(self.params, value))

  @accepted_reduction.setter
  def accepted_reduction(self, value: float) -> None:
    csleqp_call(csleqp.sleqp_params_set_accepted_reduction(self.params, value))

  @deadpoint_bound.setter
  def deadpoint_bound(self, value: float) -> None:
    csleqp_call(csleqp.sleqp_params_set_deadpoint_bound(self.params, value))

  @newton_relative_tolerance.setter
  def newton_relative_tolerance(self, value: float) -> None:
    csleqp_call(csleqp.sleqp_params_set_newton_relative_tolerance(self.params, value))

  def values(self) -> set:
    return {key: getattr(self, key) for key in self.props}

  def __str__(self) -> str:
    return 'Params: {0}'.format(self.values())

  def __repr__(self) -> str:
    return 'Params({0})'.format(repr(self.values()))
