#cython: language_level=3

# SLEQP_FAILED_ASSERTION,
# SLEQP_NOMEM,
# SLEQP_INTERNAL_ERROR,
# SLEQP_FUNC_EVAL_ERROR,
# SLEQP_MATH_ERROR,
# SLEQP_INVALID_DERIV,
# SLEQP_ILLEGAL_ARGUMENT

cdef _raise_exception():
  cdef csleqp.SLEQP_ERROR_TYPE error_type = csleqp.sleqp_error_type()
  cdef const char* error_msg = csleqp.sleqp_error_msg()
  raise Exception(error_msg.decode('UTF-8'))

cdef csleqp_call(csleqp.SLEQP_RETCODE retcode):
  if retcode != csleqp.SLEQP_OKAY:
    _raise_exception()
