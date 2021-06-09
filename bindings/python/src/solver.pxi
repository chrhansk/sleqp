#cython: language_level=3

from libc.stdlib cimport malloc, free

cdef object solvers = weakref.WeakSet()

cdef class Solver:

  cdef dict __dict__

  cdef object __weakref__

  cdef csleqp.SleqpSolver* solver
  cdef object problem
  cdef Params params
  cdef Options options
  cdef csleqp.SleqpSparseVec* residuals

  cdef list callback_handles

  def __cinit__(self,
                object problem,
                Params params,
                Options options,
                np.ndarray primal,
                Scaling scaling=None):

    cdef csleqp.SleqpSparseVec* primal_vec
    cdef _Problem _problem = <_Problem> problem._get_problem()

    self.params = params
    self.options = options
    self.callback_handles = []

    csleqp_call(csleqp.sleqp_sparse_vector_create_empty(&primal_vec,
                                                        problem.num_variables))

    csleqp_call(csleqp.sleqp_sparse_vector_create_empty(&self.residuals, 0))

    array_to_sleqp_sparse_vec(primal, primal_vec)

    csleqp_call(csleqp.sleqp_solver_create(&self.solver,
                                           _problem.cproblem,
                                           params.params,
                                           options.options,
                                           primal_vec,
                                           scaling.scaling if scaling else NULL))

    csleqp_call(csleqp.sleqp_sparse_vector_free(&primal_vec))

    self.problem = problem

    solvers.add(self)

  def __dealloc__(self):
    assert self.solver

    csleqp_call(csleqp.sleqp_solver_release(&self.solver))

    csleqp_call(csleqp.sleqp_sparse_vector_free(&self.residuals))

  cdef _solve(self,
              int max_num_iterations,
              double time_limit):

    cdef csleqp.SleqpSolver* solver = self.solver
    cdef int retcode = csleqp.SLEQP_OKAY

    self.problem.func.call_exception = None

    for obj in self.callback_handles:
      (<CallbackHandle> obj).call_exception = None

    if release_gil:
      with nogil:
        retcode = csleqp.sleqp_solver_solve(self.solver,
                                            max_num_iterations,
                                            time_limit)
    else:
        retcode = csleqp.sleqp_solver_solve(self.solver,
                                            max_num_iterations,
                                            time_limit)

    if retcode != csleqp.SLEQP_OKAY:
      exception = SLEQPError(retcode)
      call_exception = self.problem.func.call_exception
      if call_exception:
        self.problem.func.call_exception = None
        raise exception from call_exception
      else:
        for obj in self.callback_handles:
          call_exception = (<CallbackHandle> obj).call_exception
          if call_exception is not None:
            raise exception from call_exception

        # raise default
        raise exception

  def solve(self,
            max_num_iterations: int = None,
            time_limit: float = None) -> None:

    cdef int max_it = csleqp.SLEQP_NONE
    cdef double time = csleqp.SLEQP_NONE

    if max_num_iterations is not None:
      max_it = max_num_iterations

    if time_limit is not None:
      time = time_limit

    self._solve(max_it, time)

  @property
  def status(self) -> Status:
    return Status(csleqp.sleqp_solver_get_status(self.solver))

  def abort(self):
    """
    Aborts the solver as soon as the current iteration
    finishes.
    """
    csleqp_call(csleqp.sleqp_solver_abort(self.solver))

  @property
  def iterations(self) -> int:
    return csleqp.sleqp_solver_get_iterations(self.solver)

  @property
  def elapsed_seconds(self) -> float:
    return csleqp.sleqp_solver_get_elapsed_seconds(self.solver)

  @property
  def violated_cons(self) -> set:
    num_constraints =  self.problem.num_constraints

    cdef int *violated_cons = <int *> malloc(num_constraints * sizeof(double))
    cdef int num_violated_cons
    cdef csleqp.SleqpIterate* iterate

    try:
      csleqp_call(csleqp.sleqp_solver_get_solution(self.solver, &iterate))

      csleqp_call(csleqp.sleqp_solver_get_violated_constraints(self.solver,
                                                               iterate,
                                                               violated_cons,
                                                               &num_violated_cons))

      violated = set()

      for i in range(num_violated_cons):
        violated.add(violated_cons[i])

      return violated

    finally:
      free(violated_cons)

  @property
  def solution(self) -> Iterate:
    cdef csleqp.SleqpIterate* _iterate
    cdef Iterate iterate = Iterate(_create=True)

    csleqp_call(csleqp.sleqp_solver_get_solution(self.solver, &_iterate))

    iterate._set_iterate(_iterate)

    return iterate


  def add_callback(self, event, function):
    """
    Adds a callback function to the solver.
    """
    cdef CallbackHandle callback_handle = CallbackHandle(self, event, function)

    csleqp_call(get_callback_function_pointer(event,
                                              &callback_handle.function_pointer))

    csleqp_call(csleqp.sleqp_solver_add_callback(self.solver,
                                                 callback_handle.event.value,
                                                 callback_handle.function_pointer,
                                                 <void*> callback_handle))

    self.callback_handles.append(callback_handle)

    return callback_handle

  def remove_callback(self, CallbackHandle callback_handle):

    csleqp_call(csleqp.sleqp_solver_remove_callback(self.solver,
                                                    callback_handle.event.value,
                                                    callback_handle.function_pointer,
                                                    <void*> callback_handle))

    self.callback_handles.remove(callback_handle)

  cdef update_callbacks(self):
    cdef CallbackHandle callback_handle

    for obj in self.callback_handles:
      callback_handle = <CallbackHandle> obj

      csleqp_call(csleqp.sleqp_solver_remove_callback(self.solver,
                                                      callback_handle.event.value,
                                                      callback_handle.function_pointer,
                                                      <void*> callback_handle))

      csleqp_call(get_callback_function_pointer(callback_handle.event,
                                                &callback_handle.function_pointer))

      csleqp_call(csleqp.sleqp_solver_add_callback(self.solver,
                                                   callback_handle.event.value,
                                                   callback_handle.function_pointer,
                                                   <void*> callback_handle))

  @property
  def states(self):

    stat_residuals = None

    csleqp_call(csleqp.sleqp_solver_get_vec_state(self.solver, csleqp.SLEQP_SOLVER_STATE_VEC_SCALED_STAT_RESIDUALS, self.residuals))
    stat_residuals = sleqp_sparse_vec_to_array(self.residuals)

    csleqp_call(csleqp.sleqp_solver_get_vec_state(self.solver, csleqp.SLEQP_SOLVER_STATE_VEC_SCALED_FEAS_RESIDUALS, self.residuals))
    feas_residuals = sleqp_sparse_vec_to_array(self.residuals)

    csleqp_call(csleqp.sleqp_solver_get_vec_state(self.solver, csleqp.SLEQP_SOLVER_STATE_VEC_SCALED_CONS_SLACK_RESIDUALS, self.residuals))
    cons_slack_residuals = sleqp_sparse_vec_to_array(self.residuals)

    csleqp_call(csleqp.sleqp_solver_get_vec_state(self.solver, csleqp.SLEQP_SOLVER_STATE_VEC_SCALED_VAR_SLACK_RESIDUALS, self.residuals))
    vars_slack_residuals = sleqp_sparse_vec_to_array(self.residuals)

    return {
      SolverState.TrustRadius:              self._get_solver_real_state(csleqp.SLEQP_SOLVER_STATE_REAL_TRUST_RADIUS),
      SolverState.LPTrustRadius:            self._get_solver_real_state(csleqp.SLEQP_SOLVER_STATE_REAL_LP_TRUST_RADIUS),
      SolverState.ScaledFuncVal:            self._get_solver_real_state(csleqp.SLEQP_SOLVER_STATE_REAL_SCALED_FUNC_VAL),
      SolverState.ScaledMeritVal:           self._get_solver_real_state(csleqp.SLEQP_SOLVER_STATE_REAL_SCALED_MERIT_VAL),
      SolverState.ScaledFeasRes:            self._get_solver_real_state(csleqp.SLEQP_SOLVER_STATE_REAL_SCALED_FEAS_RES),
      SolverState.ScaledStatRes:            self._get_solver_real_state(csleqp.SLEQP_SOLVER_STATE_REAL_SCALED_STAT_RES),
      SolverState.ScaledSlackRes:           self._get_solver_real_state(csleqp.SLEQP_SOLVER_STATE_REAL_SCALED_SLACK_RES),
      SolverState.PenaltyParameter:         self._get_solver_real_state(csleqp.SLEQP_SOLVER_STATE_REAL_PENALTY_PARAM),
      SolverState.MinRayleigh:              self._get_solver_real_state(csleqp.SLEQP_SOLVER_STATE_REAL_MIN_RAYLEIGH),
      SolverState.MaxRayleigh:              self._get_solver_real_state(csleqp.SLEQP_SOLVER_STATE_REAL_MAX_RAYLEIGH),
      SolverState.LastStepType:             StepType(self._get_solver_int_state(csleqp.SLEQP_SOLVER_STATE_INT_LAST_STEP_TYPE)),
      SolverState.LastStepOnBoundary:       bool(self._get_solver_int_state(csleqp.SLEQP_SOLVER_STATE_INT_LAST_STEP_ON_BDRY)),
      SolverState.ScaledStatResiduals:      stat_residuals,
      SolverState.ScaledFeasResiduals:      feas_residuals,
      SolverState.ScaledConsSlackResiduals: cons_slack_residuals,
      SolverState.ScaledVarSlackResiduals:  vars_slack_residuals
    }

  cpdef double _get_solver_real_state(self, int state):
    cdef double value = 0.
    csleqp_call(csleqp.sleqp_solver_get_real_state(self.solver,
                                                   <csleqp.SLEQP_SOLVER_STATE_REAL> state,
                                                   &value))
    return value

  cpdef int _get_solver_int_state(self, int state):
    cdef int value = 0
    csleqp_call(csleqp.sleqp_solver_get_int_state(self.solver,
                                                  <csleqp.SLEQP_SOLVER_STATE_INT> state,
                                                  &value))
    return value


cdef update_solver_callbacks():
  cdef Solver solver
  for obj in solvers:
    solver = <Solver> obj

    solver.update_callbacks()
