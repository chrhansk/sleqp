cdef class CallbackHandle:
  cdef Solver solver
  cdef object function
  cdef void* function_pointer
  cdef object event

  cdef object call_exception

  def __cinit__(self, Solver solver, object event, function):
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


cdef csleqp.SLEQP_RETCODE performed_iteration(csleqp.SleqpSolver* sol,
                                              void* callback_data):
  try:
    callback_object = (<CallbackHandle> callback_data)

    function = callback_object.function
    solver = callback_object.solver

    assert solver.solver == sol

    function(solver)

  except BaseException as exception:
    callback_object.call_exception = exception
    return csleqp.SLEQP_INTERNAL_ERROR

cdef csleqp.SLEQP_RETCODE performed_iteration_nogil(csleqp.SleqpSolver* solver,
                                                    void* callback_data) nogil:
  with gil:
    return performed_iteration(solver, callback_data)


cdef  get_callback_function_pointer(solver_event, void** pointer):
  if solver_event == SolverEvent.AcceptedIterate:
    if get_release_gil():
      pointer[0] = <void*> accepted_iterate_nogil
    else:
      pointer[0] = <void*> accepted_iterate
  elif solver_event == SolverEvent.PerformedIteration:
    if get_release_gil():
      pointer[0] = <void*> performed_iteration_nogil
    else:
      pointer[0] = <void*> performed_iteration
  else:
    pointer[0] = NULL
    return csleqp.SLEQP_INTERNAL_ERROR

  return csleqp.SLEQP_OKAY
