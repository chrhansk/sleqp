#cython: language_level=3

cdef class Solver:
  cdef csleqp.SleqpSolver* solver
  cdef csleqp.SleqpSparseVec* x
  cdef Problem problem

  def __cinit__(self,
                Problem problem,
                Params params,
                np.ndarray x,
                Scaling scaling=None):

    csleqp_call(csleqp.sleqp_sparse_vector_create(&self.x,
                                                  problem.num_variables,
                                                  0))

    array_to_sleqp_sparse_vec(x, self.x)

    csleqp_call(csleqp.sleqp_solver_create(&self.solver,
                                           problem.problem,
                                           params.params,
                                           self.x,
                                           scaling.scaling if scaling else NULL))

    array_to_sleqp_sparse_vec(x, self.x)

    self.problem = problem

  cpdef solve(self,
              int max_num_iterations,
              double time_limit):

    try:
      self.problem.func.call_exception = None

      csleqp_call(csleqp.sleqp_solver_solve(self.solver,
                                            max_num_iterations,
                                            time_limit))

    except Exception as exception:
      call_exception = self.problem.func.call_exception
      if call_exception:
        self.problem.func.call_exception = None
        raise exception from call_exception
      else:
        raise exception

  @property
  def status(self):
    return Status(csleqp.sleqp_solver_get_status(self.solver))

  @property
  def iterations(self):
    return csleqp.sleqp_solver_get_iterations(self.solver)

  @property
  def elapsed_seconds(self):
    return csleqp.sleqp_solver_get_elapsed_seconds(self.solver)

  @property
  def primal(self):
    cdef csleqp.SleqpIterate* iterate

    csleqp_call(csleqp.sleqp_solver_get_solution(self.solver, &iterate))

    return sleqp_sparse_vec_to_array(iterate.x)

  @property
  def vars_dual(self):
    cdef csleqp.SleqpIterate* iterate

    csleqp_call(csleqp.sleqp_solver_get_solution(self.solver, &iterate))

    return sleqp_sparse_vec_to_array(iterate.vars_dual)

  @property
  def cons_dual(self):
    cdef csleqp.SleqpIterate* iterate

    csleqp_call(csleqp.sleqp_solver_get_solution(self.solver, &iterate))

    return sleqp_sparse_vec_to_array(iterate.cons_dual)

  def  __dealloc__(self):
    csleqp_call(csleqp.sleqp_solver_free(&self.solver))
