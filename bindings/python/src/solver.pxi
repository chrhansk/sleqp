cdef class Solver:
  cdef csleqp.SleqpSolver* solver

  def __cinit__(self,
                Problem problem,
                Params params,
                np.ndarray x):
    pass

  cpdef solve(self,
              int max_num_iterations,
              double time_limit):
    csleqp.sleqp_solver_solve(self.solver,
                              max_num_iterations,
                              time_limit)

  @property
  def status(self):
      pass

  @property
  def iterations(self):
      return csleqp.sleqp_solver_get_iterations(self.solver)

  @property
  def elapsed_seconds(self):
      return csleqp.sleqp_solver_get_elapsed_seconds(self.solver)

  @property
  def primal_solution(self):
      cdef csleqp.SleqpIterate* iterate

      csleqp.sleqp_solver_get_solution(self.solver, &iterate)

      return sleqp_sparse_vec_to_array(iterate.x)

  @property
  def vars_dual(self):
      cdef csleqp.SleqpIterate* iterate

      csleqp.sleqp_solver_get_solution(self.solver, &iterate)

      return sleqp_sparse_vec_to_array(iterate.vars_dual)

  @property
  def cons_dual(self):
      cdef csleqp.SleqpIterate* iterate

      csleqp.sleqp_solver_get_solution(self.solver, &iterate)

      return sleqp_sparse_vec_to_array(iterate.cons_dual)

  def  __dealloc__(self):
      csleqp.sleqp_solver_free(&self.solver)
