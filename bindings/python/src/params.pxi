cdef class Params:
  cdef csleqp.SleqpParams* params

  def __cinit__(self):
    csleqp_call(csleqp.sleqp_params_create(&self.params))

  def __dealloc__(self):
    csleqp_call(csleqp.sleqp_params_free(&self.params))
