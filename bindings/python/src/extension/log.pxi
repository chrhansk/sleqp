import logging

cimport libc.time

cdef void sleqp_python_handler(csleqp.SLEQP_LOG_LEVEL level,
                               libc.time.time_t time,
                               const char* message) noexcept:
  levels = {
    csleqp.SLEQP_LOG_DEBUG: logging.DEBUG,
    csleqp.SLEQP_LOG_ERROR: logging.ERROR,
    csleqp.SLEQP_LOG_WARN: logging.WARNING,
    csleqp.SLEQP_LOG_INFO: logging.INFO,
  }

  assert level in levels

  msg = message.decode("UTF-8", "replace")

  sleqp_logger.log(levels[level], msg)

cdef void sleqp_python_handler_nogil(csleqp.SLEQP_LOG_LEVEL level,
                                     libc.time.time_t time,
                                     const char* message) noexcept nogil:
  with gil:
    sleqp_python_handler(level, time, message)


cdef update_log_handler():
  if release_gil:
    csleqp.sleqp_log_set_handler(sleqp_python_handler_nogil)
  else:
    csleqp.sleqp_log_set_handler(sleqp_python_handler)


class SleqpLogger(logging.Logger):
  def __init__(self, name, level=logging.NOTSET):
    super(SleqpLogger, self).__init__(name, level)

    update_log_handler()


  def setLevel(self, level):
    super(SleqpLogger, self).setLevel(level)

    levels = {
      logging.DEBUG: csleqp.SLEQP_LOG_DEBUG,
      logging.CRITICAL: csleqp.SLEQP_LOG_ERROR,
      logging.ERROR: csleqp.SLEQP_LOG_ERROR,
      logging.WARNING: csleqp.SLEQP_LOG_WARN,
      logging.INFO: csleqp.SLEQP_LOG_INFO,
      logging.NOTSET: csleqp.SLEQP_LOG_DEBUG
    }

    assert level in levels

    csleqp.sleqp_log_set_level(levels[level])


_logging_class = logging.getLoggerClass()
logging.setLoggerClass(SleqpLogger)

sleqp_logger = logging.getLogger('sleqp')

logging.setLoggerClass(_logging_class)
