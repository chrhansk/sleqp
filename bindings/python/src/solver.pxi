#cython: language_level=3

from libc.stdlib cimport malloc, free

cdef class Solver:
  cdef csleqp.SleqpSolver* solver
  cdef Problem problem
  cdef Params params
  cdef Options options

  def __cinit__(self,
                Problem problem,
                Params params,
                Options options,
                np.ndarray primal,
                Scaling scaling=None):

    cdef csleqp.SleqpSparseVec* primal_vec

    self.params = params
    self.options = options

    csleqp_call(csleqp.sleqp_sparse_vector_create_empty(&primal_vec,
                                                        problem.num_variables))

    array_to_sleqp_sparse_vec(primal, primal_vec)

    csleqp_call(csleqp.sleqp_solver_create(&self.solver,
                                           problem.problem,
                                           params.params,
                                           options.options,
                                           primal_vec,
                                           scaling.scaling if scaling else NULL))

    csleqp_call(csleqp.sleqp_sparse_vector_free(&primal_vec))

    self.problem = problem

  def __dealloc__(self):
    assert self.solver

    csleqp_call(csleqp.sleqp_solver_release(&self.solver))

  cdef _solve(self,
              int max_num_iterations,
              double time_limit):

    cdef csleqp.SleqpSolver* solver = self.solver
    cdef int retcode = csleqp.SLEQP_OKAY

    self.problem.func.call_exception = None

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
      cdef Iterate iterate = Iterate()

      csleqp_call(csleqp.sleqp_solver_get_solution(self.solver, &_iterate))

      iterate._set_iterate(_iterate)

      return iterate
