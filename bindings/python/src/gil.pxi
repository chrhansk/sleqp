cdef bint release_gil = False

cpdef get_release_gil():
  return release_gil

cpdef set_release_gil(bint value):
  global release_gil

  if value == release_gil:
    return

  release_gil = value
  update_log_handler()
  update_callbacks()
