# cython: language_level=3

from collections import defaultdict


class CallbackError(Exception):
    def __init__(self, *args, **kwds):
        super().__init__(*args, **kwds)


class EvaluationError(Exception):
    def __init__(self, *args, **kwds):
        super().__init__(*args, **kwds)


cdef object _exception_map = {
    csleqp.SLEQP_FAILED_ASSERTION: AssertionError,
    csleqp.SLEQP_NOMEM: MemoryError,
    csleqp.SLEQP_INTERNAL_ERROR: Exception,
    csleqp.SLEQP_FUNC_EVAL_ERROR: EvaluationError,
    csleqp.SLEQP_CALLBACK_ERROR: CallbackError,
    csleqp.SLEQP_MATH_ERROR: ArithmeticError,
    csleqp.SLEQP_INVALID_DERIV: Exception,
    csleqp.SLEQP_ILLEGAL_ARGUMENT: ValueError
}

cdef _get_exception():
    cdef csleqp.SLEQP_ERROR_TYPE error_type = csleqp.sleqp_error_type()
    cdef const char * error_msg = csleqp.sleqp_error_msg()
    msg = error_msg.decode('UTF-8')
    if error_type in _exception_map:
        exception = _exception_map[error_type]
        return exception(msg)
    return Exception(error_msg.decode('UTF-8'))

cdef _raise_exception():
    raise _get_exception()

cdef csleqp_call(csleqp.SLEQP_RETCODE retcode):
    if retcode != csleqp.SLEQP_OKAY:
        _raise_exception()
