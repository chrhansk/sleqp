#cython: language_level=3

cdef class Params:
  cdef csleqp.SleqpParams* params

  def __cinit__(self):
    csleqp_call(csleqp.sleqp_params_create(&self.params))

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
    return csleqp.sleqp_params_get_deriv_pertubation(self.params)

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
