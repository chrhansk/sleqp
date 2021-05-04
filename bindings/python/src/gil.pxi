cdef bint release_gil = False

cpdef get_release_gil():
  """
  Returns whether or not the GIL is released during the solution process

  :rtype: bool
  """
  return release_gil

cpdef set_release_gil(bint value):
  """
  Advising the solver whether or not to release the GIL during
  the solution process (it will be reacquired during function / callback evaluations)

  :param value: whether or not to relase the GIL
  :type value: bool
  """
  global release_gil

  if value == release_gil:
    return

  release_gil = value
  update_log_handler()
  update_func_callbacks()
  update_lsq_func_callbacks()
  update_solver_callbacks()
