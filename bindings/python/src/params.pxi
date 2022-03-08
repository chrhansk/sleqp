#cython: language_level=3

# global constant properties
cdef dict prop_map = {
  'zero_eps':           csleqp.SLEQP_PARAM_ZERO_EPS,
  'eps':                csleqp.SLEQP_PARAM_EPS,
  'deriv_perturbation': csleqp.SLEQP_PARAM_DERIV_PERTURBATION,
  'deriv_tol':          csleqp.SLEQP_PARAM_DERIV_TOL,
  'cauchy_tau':         csleqp.SLEQP_PARAM_CAUCHY_TAU,
  'cauchy_eta':         csleqp.SLEQP_PARAM_CAUCHY_ETA,
  'linesearch_tau':     csleqp.SLEQP_PARAM_LINESEARCH_TAU,
  'linesearch_eta':     csleqp.SLEQP_PARAM_LINESEARCH_ETA,
  'linesearch_cutoff':  csleqp.SLEQP_PARAM_LINESEARCH_CUTOFF,
  'feas_tol':           csleqp.SLEQP_PARAM_FEAS_TOL,
  'slack_tol':          csleqp.SLEQP_PARAM_SLACK_TOL,
  'stat_tol':           csleqp.SLEQP_PARAM_STAT_TOL,
  'accepted_reduction': csleqp.SLEQP_PARAM_ACCEPTED_REDUCTION,
  'deadpoint_bound':    csleqp.SLEQP_PARAM_DEADPOINT_BOUND
}

cdef class Params:
  cdef csleqp.SleqpParams* params

  def __cinit__(self):
    csleqp_call(csleqp.sleqp_params_create(&self.params))

  def __init__(self, **values):
    for key, value in values.items():
      setattr(self, key, value)

  def __dealloc__(self):
    csleqp_call(csleqp.sleqp_params_release(&self.params))

  def props(self) -> dict:
    return {key: getattr(self, key) for key in prop_map.keys()}

  def __str__(self) -> str:
    return 'Params: {0}'.format(self.values())

  def __repr__(self) -> str:
    return 'Params({0})'.format(repr(self.values()))

  cdef _prop_val(self, name):
    prop_val = prop_map.get(name)
    if prop_val is None:
      raise AttributeError('Params has no attribute \'{0}\''.format(name))
    return prop_val

  def __getattr__(self, name):
    return csleqp.sleqp_params_value(self.params,
                                     self._prop_val(name))

  def __setattr__(self, name, value):
    csleqp_call(csleqp.sleqp_params_set_value(self.params,
                                              self._prop_val(name),
                                              value))
