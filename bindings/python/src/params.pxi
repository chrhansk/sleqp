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
                  'feasibility_tolerance',
                  'slackness_tolerance',
                  'stationarity_tolerance',
                  'accepted_reduction',
                  'deadpoint_bound']

    for key, value in values.items():
      self._set_prop(key, value)

  cdef _set_prop(self, name, value):
    if not name in self.props:
      raise AttributeError("Invalid property {0}".format(name))

    setattr(self, name, value)

  def __dealloc__(self):
    csleqp_call(csleqp.sleqp_params_release(&self.params))

  @property
  def zero_eps(self) -> float:
    return csleqp.sleqp_params_get(self.params,
                                   csleqp.SLEQP_PARAM_ZERO_EPS)

  @property
  def eps(self) -> float:
    return csleqp.sleqp_params_get(self.params,
                                   csleqp.SLEQP_PARAM_EPS)

  @property
  def deriv_perturbation(self) -> float:
    return csleqp.sleqp_params_get(self.params,
                                   csleqp.SLEQP_PARAM_DERIV_PERTURBATION)

  @property
  def deriv_tolerance(self) -> float:
    return csleqp.sleqp_params_get(self.params,
                                   csleqp.SLEQP_PARAM_DERIV_TOL)

  @property
  def cauchy_tau(self) -> float:
    return csleqp.sleqp_params_get(self.params,
                                   csleqp.SLEQP_PARAM_CAUCHY_TAU)

  @property
  def cauchy_eta(self) -> float:
    return csleqp.sleqp_params_get(self.params,
                                   csleqp.SLEQP_PARAM_CAUCHY_ETA)

  @property
  def linesearch_tau(self) -> float:
    return csleqp.sleqp_params_get(self.params,
                                   csleqp.SLEQP_PARAM_LINESEARCH_TAU)

  @property
  def linesearch_eta(self) -> float:
    return csleqp.sleqp_params_get(self.params,
                                   csleqp.SLEQP_PARAM_LINESEARCH_ETA)

  @property
  def linesearch_cutoff(self) -> float:
    return csleqp.sleqp_params_get(self.params,
                                   csleqp.SLEQP_PARAM_LINESEARCH_CUTOFF)

  @property
  def feasibility_tolerance(self) -> float:
    return csleqp.sleqp_params_get(self.params,
                                   csleqp.SLEQP_PARAM_FEASIBILITY_TOL)

  @property
  def slackness_tolerance(self) -> float:
    return csleqp.sleqp_params_get(self.params,
                                   csleqp.SLEQP_PARAM_SLACKNESS_TOL)

  @property
  def stationarity_tolerance(self) -> float:
    return csleqp.sleqp_params_get(self.params,
                                   csleqp.SLEQP_PARAM_STATIONARITY_TOL)

  @property
  def accepted_reduction(self) -> float:
    return csleqp.sleqp_params_get(self.params,
                                   csleqp.SLEQP_PARAM_ACCEPTED_REDUCTION)

  @property
  def deadpoint_bound(self) -> float:
    return csleqp.sleqp_params_get(self.params,
                                   csleqp.SLEQP_PARAM_DEADPOINT_BOUND)

  @zero_eps.setter
  def zero_eps(self, value: float) -> None:
    csleqp_call(csleqp.sleqp_params_set(self.params,
                                        csleqp.SLEQP_PARAM_ZERO_EPS,
                                        value))

  @eps.setter
  def eps(self, value: float) -> None:
    csleqp_call(csleqp.sleqp_params_set(self.params,
                                        csleqp.SLEQP_PARAM_EPS,
                                        value))

  @deriv_perturbation.setter
  def deriv_perturbation(self, value: float) -> None:
    csleqp_call(csleqp.sleqp_params_set(self.params,
                                        csleqp.SLEQP_PARAM_DERIV_PERTURBATION,
                                        value))

  @deriv_tolerance.setter
  def deriv_tolerance(self, value: float) -> None:
    csleqp_call(csleqp.sleqp_params_set(self.params,
                                        csleqp.SLEQP_PARAM_DERIV_TOL,
                                        value))

  @cauchy_tau.setter
  def cauchy_tau(self, value: float) -> None:
    csleqp_call(csleqp.sleqp_params_set(self.params,
                                        csleqp.SLEQP_PARAM_CAUCHY_TAU,
                                        value))

  @cauchy_eta.setter
  def cauchy_eta(self, value: float) -> None:
    csleqp_call(csleqp.sleqp_params_set(self.params,
                                        csleqp.SLEQP_PARAM_CAUCHY_ETA,
                                        value))

  @linesearch_tau.setter
  def linesearch_tau(self, value: float) -> None:
    csleqp_call(csleqp.sleqp_params_set(self.params,
                                        csleqp.SLEQP_PARAM_LINESEARCH_TAU,
                                        value))

  @linesearch_eta.setter
  def linesearch_eta(self, value: float) -> None:
    csleqp_call(csleqp.sleqp_params_set(self.params,
                                        csleqp.SLEQP_PARAM_LINESEARCH_ETA,
                                        value))

  @linesearch_cutoff.setter
  def linesearch_cutoff(self, value: float) -> None:
    csleqp_call(csleqp.sleqp_params_set(self.params,
                                        csleqp.SLEQP_PARAM_LINESEARCH_CUTOFF,
                                        value))

  @feasibility_tolerance.setter
  def feasibility_tolerance(self, value: float) -> None:
    csleqp_call(csleqp.sleqp_params_set(self.params,
                                        csleqp.SLEQP_PARAM_FEASIBILITY_TOL,
                                        value))

  @slackness_tolerance.setter
  def slackness_tolerance(self, value: float) -> None:
    csleqp_call(csleqp.sleqp_params_set(self.params,
                                        csleqp.SLEQP_PARAM_SLACKNESS_TOL,
                                        value))

  @stationarity_tolerance.setter
  def stationarity_tolerance(self, value: float) -> None:
    csleqp_call(csleqp.sleqp_params_set(self.params,
                                        csleqp.SLEQP_PARAM_STATIONARITY_TOL,
                                        value))

  @accepted_reduction.setter
  def accepted_reduction(self, value: float) -> None:
    csleqp_call(csleqp.sleqp_params_set(self.params,
                                        csleqp.SLEQP_PARAM_ACCEPTED_REDUCTION,
                                        value))

  @deadpoint_bound.setter
  def deadpoint_bound(self, value: float) -> None:
    csleqp_call(csleqp.sleqp_params_set(self.params,
                                        csleqp.SLEQP_PARAM_DEADPOINT_BOUND,
                                        value))

  def values(self) -> set:
    return {key: getattr(self, key) for key in self.props}

  def __str__(self) -> str:
    return 'Params: {0}'.format(self.values())

  def __repr__(self) -> str:
    return 'Params({0})'.format(repr(self.values()))
