class SLEQPError(Exception):
  def __init__(self, code):
    assert code != csleqp.SLEQP_OKAY
    self.code = code

  def __str__(self):
    messages = {
      csleqp.SLEQP_NOMEM: "Out of memory",
      csleqp.SLEQP_ILLEGAL_ARGUMENT: "Illegal argument",
      csleqp.SLEQP_INVALID_DERIV: "Invalid derivative",
      csleqp.SLEQP_INTERNAL_ERROR: "Internal error"
    }

    assert self.code in messages

    return messages[self.code]

cpdef csleqp_call(csleqp.SLEQP_RETCODE retcode):
  if retcode != csleqp.SLEQP_OKAY:
    raise SLEQPError(retcode)
