cdef class CallbackHandle:
  cdef Solver solver
  cdef object function
  cdef void* function_pointer
  cdef csleqp.SLEQP_SOLVER_EVENT event

  cdef object call_exception

  def __cinit__(self, Solver solver, csleqp.SLEQP_SOLVER_EVENT event, function):
    self.solver = solver
    self.function = function
    self.call_exception = None
    self.event = event

  def unregister(self):
    self.solver.remove_callback(self)


cdef csleqp.SLEQP_RETCODE accepted_iterate(csleqp.SleqpSolver* sol,
                                           csleqp.SleqpIterate* it,
                                           void* callback_data):
  cdef Iterate iterate = Iterate(_create=True)
  cdef CallbackHandle callback_object = None
  iterate._set_iterate(it)

  try:
    callback_object = (<CallbackHandle> callback_data)

    function = callback_object.function
    solver = callback_object.solver

    assert solver.solver == sol

    function(solver, iterate)

  except BaseException as exception:
    callback_object.call_exception = exception
    return csleqp.SLEQP_INTERNAL_ERROR

cdef csleqp.SLEQP_RETCODE accepted_iterate_nogil(csleqp.SleqpSolver* solver,
                                                 csleqp.SleqpIterate* iterate,
                                                 void* callback_data) nogil:
  with gil:
    return accepted_iterate(solver, iterate, callback_data)


cdef void* get_callback_function_pointer(solver_event):
  if solver_event == SolverEvent.AcceptedIterate:
    if get_release_gil():
      return <void*> accepted_iterate_nogil
    else:
      return <void*> accepted_iterate
  else:
    raise Exception("Invalid event: {0}".format(solver_event))
