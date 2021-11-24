#cython: language_level=3

cdef class WorkingSet:
  cdef csleqp.SleqpWorkingSet* working_set

  cdef dict __dict__

  cdef _release(self):
    if self.working_set:
      csleqp_call(csleqp.sleqp_working_set_release(&self.working_set))

  cdef _set_working_set(self, csleqp.SleqpWorkingSet* working_set):
    self._release()

    csleqp_call(csleqp.sleqp_working_set_capture(working_set))

    self.working_set = working_set

  @property
  def num_vars(self) -> int:
    cdef csleqp.SleqpProblem* problem

    problem = csleqp.sleqp_working_set_problem(self.working_set)

    return csleqp.sleqp_problem_num_vars(problem)

  @property
  def num_cons(self) -> int:
    cdef csleqp.SleqpProblem* problem

    problem = csleqp.sleqp_working_set_problem(self.working_set)

    return csleqp.sleqp_problem_num_cons(problem)

  @property
  def num_active_vars(self) -> int:
    assert self.working_set

    return csleqp.sleqp_working_set_num_active_vars(self.working_set)

  @property
  def num_active_cons(self) -> int:
    assert self.working_set

    return csleqp.sleqp_working_set_num_active_cons(self.working_set)

  def var_state(self, int i) -> ActiveState:
    assert self.working_set

    state = csleqp.sleqp_working_set_var_state(self.working_set,
                                               i)

    return ActiveState(state)

  def cons_state(self, int i) -> ActiveState:
    assert self.working_set

    state = csleqp.sleqp_working_set_cons_state(self.working_set,
                                                i)

    return ActiveState(state)

  def __str__(self):
    str_val = "Active set, variables: {0}, constraints: {1}\n".format(self.num_vars,
                                                                      self.num_cons)

    str_val += "\n".join(("State of variable {0}: {1}".format(j, str(self.var_state(j)))
                          for j in range(self.num_vars)))

    str_val += "\n"

    str_val += "\n".join(("State of constraint {0}: {1}".format(i, str(self.cons_state(i)))
                          for i in range(self.num_cons)))

    return str_val

  @property
  def size(self) -> int:
    assert self.working_set

    return csleqp.sleqp_working_set_size(self.working_set)

  def __dealloc__(self):
    self._release()
